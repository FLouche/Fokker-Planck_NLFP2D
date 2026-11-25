!*******************************************************************
!*   Resolution of the time-dependent Fokker-Planck Equation       *
!*   with Crank-Nicolson or fully implicit scheme                  *
!*   7-point Fornberg stencil in both v⊥ and v∥                    *
!*                                                                 *
!*   Replaces TimeFP3.f90 / build_matrix_ss-2.f90                 *
!*   by using fd_stencil_2d for higher-order spatial accuracy.     *
!*                                                                 *
!*   Key differences from TimeFP3:                                 *
!*   - The steady-state operator L is built and stored in CSR      *
!*     format (never as a dense matrix), saving O(nbig^2) memory.  *
!*   - Crank-Nicolson: solves  (I - dt/2 * L) f^{n+1}             *
!*                           = (I + dt/2 * L) f^n  + dt * S       *
!*   - Fully implicit:  solves  (I - dt   * L) f^{n+1}            *
!*                           =                 f^n  + dt * S       *
!*   - The LHS matrix is constant (linear problem), so PARDISO     *
!*     performs phases 11+22 once and only phase 33 each step.     *
!*                                                                 *
!*   Version 1.05 - F. Louche                                      *
!*******************************************************************

MODULE mod_timefp_7pt

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: timefp_7pt

CONTAINS

SUBROUTINE timefp_7pt(all00_lin, all10_lin, all01_lin, &
                      all11_lin, all20_lin, all02_lin, &
                      fstart, fout, otime)

  USE mod_fd_stencil_2d               ! provides fd_stencil_2d
  USE pardiso_solver                  ! provides pardiso_handle_t,
                                      !   pardiso_solve_init/step/finalize
  USE shared_grid
  USE shared_plasma
  USE shared_timer
  USE shared_beam
  USE shared_multi, ONLY: tau_E
  USE shared_RF
  USE func_index
  USE mod_conv_diag
  USE time_comps_mod
  USE assemble_FP_lin
  USE coulomb_log_mod
  USE shared_FPterms

  IMPLICIT NONE

  INTEGER, PARAMETER :: dp = KIND(1.0D0)

  !--- Arguments ---------------------------------------------------
  REAL(dp), INTENT(IN),    DIMENSION(nperp,npar) :: &
      all00_lin, all10_lin, all01_lin, &
      all11_lin, all20_lin, all02_lin
  REAL(dp), INTENT(INOUT), DIMENSION(nbig)        :: fstart
  REAL(dp), INTENT(OUT),   DIMENSION(nperp,npar)  :: fout
  REAL(dp), INTENT(IN)                            :: otime

  !--- PARDISO handle ----------------------------------------------
  TYPE(pardiso_handle_t) :: handle_lhs   ! LHS: (I - theta*dt*L)

  !--- CSR storage for the steady-state operator L -----------------
  !    and for the LHS matrix M_lhs = I - theta*dt*L
  INTEGER,  ALLOCATABLE :: ia_L(:), ja_L(:)   ! CSR of L  (1-based)
  REAL(dp), ALLOCATABLE :: aa_L(:)
  INTEGER,  ALLOCATABLE :: ia_lhs(:), ja_lhs(:)
  REAL(dp), ALLOCATABLE :: aa_lhs(:)

  INTEGER :: nnz_L, nnz_max

  !--- Stencil workspace -------------------------------------------
  INTEGER  :: col_idx(49)
  REAL(dp) :: stencil_coeff(49)
  INTEGER  :: n_entries
  REAL(dp) :: rhs_ij
  LOGICAL  :: pattern_ok

  !--- Working vectors ---------------------------------------------
  REAL(dp), ALLOCATABLE :: rhs_vec(:), x_vec(:), Lf(:)

  !--- Scalars -----------------------------------------------------
  REAL(dp) :: theta          ! 0.5 for CN, 1.0 for implicit
  REAL(dp) :: time, dens_tmp, tk, tkperp, tkpar, teff, teff_tmp
  REAL(dp) :: pcoll(nbulk), pRF, psource, plosses, pcoll_self
  REAL(dp) :: ncoll_d(nbulk), nsc_d, nRF_d, nsrc_d, nloss_d   ! dn/dt per operator term
  REAL(dp) :: tau_rf            ! RF tail formation time [s]
  REAL(dp) :: tau_ii, tau_ie    ! effective collisional times [s]
  REAL(dp) :: pcoll_self_perp, pcoll_self_par
  REAL(dp) :: pE_loss, pE_ratio      ! energy-confinement loss power, check ratio

  INTEGER :: ndof, i, j, k, row, ptr, itime, itime_global, iphase, iv, imu, ix
  REAL(dp) :: phase_offset
  INTEGER :: error, ib
  CHARACTER(len=2)   :: ibString
  CHARACTER(len=256) :: dynfname

  EXTERNAL :: time_power_7pt, time_momentum_7pt, time_density_terms_7pt
  EXTERNAL :: time_power_lossE

  real(dp) :: anisotropy
  type(conv_diag_t) :: cdiag        ! Jacobian-weighted convergence diagnostics

  real(dp), dimension(nperp,npar) :: f_init
  REAL(dp), DIMENSION(nperp,npar) :: all00, all10, all01, all11, all20, all02
  REAL(dp), PARAMETER :: gamma0 = 2.390775d-1
  REAL(dp) :: lnab_t, cte0_t, ta_eV
  REAL(dp) :: lnab_arr(nbulk)
  REAL(dp) :: mcoll_perp(nbulk), mcoll_par(nbulk)
  REAL(dp) :: mRF_perp, mRF_par, msrc_perp, msrc_par
  REAL(dp) :: mloss_perp, mloss_par, mSC_perp, mSC_par

  !================================================================
  ! 0.  Setup
  !================================================================
  ndof    = nperp * npar
  nnz_max = ndof * 49          ! upper bound: 7x7 stencil per row

  ! Time-stepping weight
  IF (icn == -1) THEN
    theta = 0.5_dp             ! Crank-Nicolson
  ELSE IF (icn == 1) THEN
    theta = 0.75_dp            ! intermediate
  ELSE
    theta = 1.0_dp             ! fully implicit
  END IF

  ALLOCATE(ia_L(ndof+1), ja_L(nnz_max), aa_L(nnz_max))
  ALLOCATE(ia_lhs(ndof+1), ja_lhs(nnz_max), aa_lhs(nnz_max))
  ALLOCATE(rhs_vec(ndof), x_vec(ndof), Lf(ndof))

  !================================================================
  ! 1.  Build steady-state operator L in CSR
  !     L encodes:  A*f + B*df/dvp + ... + F*d2f/dva2 - taum*f
  !     exactly as in solve_fp_pardiso, but stored in CSR.
  !================================================================

  !--- Pass 1: count non-zeros per row -> ia_L --------------------
  ia_L(1) = 1
  DO i = 1, nperp
    DO j = 1, npar
      row = (i-1)*npar + j
      CALL fd_stencil_2d(i, j, nperp, npar, vperp, dvpar, &
                         all00_lin(i,j), all10_lin(i,j), all01_lin(i,j), &
                         all20_lin(i,j), all11_lin(i,j), all02_lin(i,j), &
                         col_idx, stencil_coeff, n_entries, rhs_ij)
      ia_L(row+1) = ia_L(row) + n_entries
    END DO
  END DO

  nnz_L = ia_L(ndof+1) - 1
  IF (nnz_L > nnz_max) THEN
    WRITE(*,*) 'timefp_7pt: nnz_max exceeded, nnz_L=', nnz_L; STOP
  END IF

  !--- Pass 2: fill ja_L, aa_L ------------------------------------
  ptr = 1
  DO i = 1, nperp
    DO j = 1, npar
      CALL fd_stencil_2d(i, j, nperp, npar, vperp, dvpar, &
                         all00_lin(i,j), all10_lin(i,j), all01_lin(i,j), &
                         all20_lin(i,j), all11_lin(i,j), all02_lin(i,j), &
                         col_idx, stencil_coeff, n_entries, rhs_ij)
      CALL sort_stencil(col_idx, stencil_coeff, n_entries)
      DO k = 1, n_entries
        ja_L(ptr) = col_idx(k)
        aa_L(ptr) = stencil_coeff(k)
        ptr = ptr + 1
      END DO
    END DO
  END DO

  !WRITE(*,'(A,I10,A,F6.2,A)') '  L operator: nnz=', nnz_L, &
   !   '  (', 100.d0*nnz_L/DBLE(ndof)**2, ' %)'

  !================================================================
  ! 2.  Build LHS matrix:  M_lhs = I - theta*dt*L
  !
  !     The sparsity pattern of M_lhs is identical to L plus the
  !     diagonal (identity).  Since the diagonal is already in L
  !     (from the A*f term), the pattern is the same as L.
  !     We construct aa_lhs = -theta*dt * aa_L, then add 1 to
  !     each diagonal entry.
  !================================================================
  ia_lhs = ia_L
  ja_lhs = ja_L

  DO ptr = 1, nnz_L
    aa_lhs(ptr) = -theta * timestep_cur * aa_L(ptr)
  END DO

  ! Add identity: find diagonal entries (col == row) and add 1
  DO row = 1, ndof
    DO ptr = ia_lhs(row), ia_lhs(row+1)-1
      IF (ja_lhs(ptr) == row) THEN
        aa_lhs(ptr) = aa_lhs(ptr) + 1.0_dp
        EXIT
      END IF
    END DO
  END DO

  ! BC rows must be treated as constraints, not time-evolution equations.
  ! Restore them to the original L stencil so the solver enforces
  ! df/dvp=0 (i=1) and f=0 (i=nperp, j=1, j=npar) at every time step.
  DO row = 1, ndof
    CALL index_mat_inv(row, iv, imu)
    IF (iv == 1 .OR. iv == nperp .OR. imu == 1 .OR. imu == npar) THEN
      DO ptr = ia_lhs(row), ia_lhs(row+1)-1
        aa_lhs(ptr) = aa_L(ptr)
      END DO
    END IF
  END DO

  !================================================================
  ! 3.  Open output files (same convention as TimeFP3)
  !================================================================
  IF (otime == 0.d0) THEN
    OPEN(45, file=TRIM(outfile('density_vs_time.txt')),     status='unknown')
    OPEN(46, file=TRIM(outfile('energy_vs_time.txt')),      status='unknown')
    OPEN(47, file=TRIM(outfile('anisotropy_vs_time.txt')), status='unknown')
    IF (nbulk >   1) OPEN(505,file=TRIM(outfile('coulomb_log_vs_time.txt')),     status='unknown')
                     OPEN(507,file=TRIM(outfile('Teff_vs_time.txt')),             status='unknown')
    IF (iplot_pow == -1) THEN
      OPEN(470,file=TRIM(outfile('power_coll_tot_vs_time.txt')), status='unknown')
      DO ib = 1, nbulk
        IF (ib == 1) THEN
          OPEN(471,file=TRIM(outfile('power_coll_e_vs_time.txt')), status='unknown')
        ELSE
          WRITE(ibString,'(i2)') ib-1
          dynfname = 'power_coll_ion'//ibString//'_vs_time.txt'
          OPEN(470+ib, file=TRIM(outfile(dynfname)), status='unknown')
        END IF
      END DO
      IF (irf   == -1) OPEN(480,file=TRIM(outfile('power_RF_vs_time.txt')),        status='unknown')
      IF (irf   == -1) OPEN(515,file=TRIM(outfile('tau_rf_vs_time.txt')),          status='unknown')
      OPEN(516,file=TRIM(outfile('tau_coll_vs_time.txt')),        status='unknown')
      IF (isource==-1) OPEN(490,file=TRIM(outfile('power_NBI_vs_time.txt')),       status='unknown')
      IF (tau_E > 0.d0) OPEN(495,file=TRIM(outfile('power_lossE_vs_time.txt')), status='unknown')
      !DTERMS-OFF! OPEN(520,file=TRIM(outfile('density_terms_vs_time.txt')),   status='unknown')
      !DTERMS-OFF! WRITE(520,'(A)') '# time  total  coll(1:nbulk; 1=e)  self-coll  RF  source  losses   [m^-3 s^-1]'
    END IF
    IF (iplot_mom == -1) THEN
      OPEN(570,file=TRIM(outfile('momentum_coll_tot_vs_time.txt')), status='unknown')
      DO ib = 1, nbulk
        IF (ib == 1) THEN
          OPEN(571,file=TRIM(outfile('momentum_coll_e_vs_time.txt')), status='unknown')
        ELSE
          WRITE(ibString,'(i2)') ib-1
          dynfname = 'momentum_coll_ion'//ibString//'_vs_time.txt'
          OPEN(570+ib, file=TRIM(outfile(dynfname)), status='unknown')
        END IF
      END DO
      IF (irf   == -1) OPEN(580,file=TRIM(outfile('momentum_RF_vs_time.txt')),        status='unknown')
      IF (isource==-1) OPEN(590,file=TRIM(outfile('momentum_NBI_vs_time.txt')),       status='unknown')
    END IF
  ELSE
    OPEN(45, file=TRIM(outfile('density_vs_time.txt')),     status='old', access='append')
    OPEN(46, file=TRIM(outfile('energy_vs_time.txt')),      status='old', access='append')
    OPEN(47, file=TRIM(outfile('anisotropy_vs_time.txt')), status='old', access='append')
    IF (nbulk >   1) OPEN(505,file=TRIM(outfile('coulomb_log_vs_time.txt')),     status='old', access='append')
                     OPEN(507,file=TRIM(outfile('Teff_vs_time.txt')),             status='old', access='append')
    IF (iplot_pow == -1) THEN
      OPEN(470,file=TRIM(outfile('power_coll_tot_vs_time.txt')), status='old', access='append')
      DO ib = 1, nbulk
        IF (ib == 1) THEN
          OPEN(471,file=TRIM(outfile('power_coll_e_vs_time.txt')), status='old', access='append')
        ELSE
          WRITE(ibString,'(i2)') ib-1
          dynfname = 'power_coll_ion'//ibString//'_vs_time.txt'
          OPEN(470+ib, file=TRIM(outfile(dynfname)), status='old', access='append')
        END IF
      END DO
      IF (irf   == -1) OPEN(480,file=TRIM(outfile('power_RF_vs_time.txt')),        status='old', access='append')
      IF (irf   == -1) OPEN(515,file=TRIM(outfile('tau_rf_vs_time.txt')),          status='old', access='append')
      OPEN(516,file=TRIM(outfile('tau_coll_vs_time.txt')),        status='old', access='append')
      IF (isource==-1) OPEN(490,file=TRIM(outfile('power_NBI_vs_time.txt')),       status='old', access='append')
      IF (tau_E > 0.d0) OPEN(495,file=TRIM(outfile('power_lossE_vs_time.txt')), status='old', access='append')
      !DTERMS-OFF! OPEN(520,file=TRIM(outfile('density_terms_vs_time.txt')),   status='old', access='append')
    END IF
    IF (iplot_mom == -1) THEN
      OPEN(570,file=TRIM(outfile('momentum_coll_tot_vs_time.txt')), status='old', access='append')
      DO ib = 1, nbulk
        IF (ib == 1) THEN
          OPEN(571,file=TRIM(outfile('momentum_coll_e_vs_time.txt')), status='old', access='append')
        ELSE
          WRITE(ibString,'(i2)') ib-1
          dynfname = 'momentum_coll_ion'//ibString//'_vs_time.txt'
          OPEN(570+ib, file=TRIM(outfile(dynfname)), status='old', access='append')
        END IF
      END DO
      IF (irf   == -1) OPEN(580,file=TRIM(outfile('momentum_RF_vs_time.txt')),        status='old', access='append')
      IF (isource==-1) OPEN(590,file=TRIM(outfile('momentum_NBI_vs_time.txt')),       status='old', access='append')
    END IF
  END IF

  !================================================================
  ! 4.  Factorise M_lhs once (phases 11 + 22)
  !================================================================
  CALL pardiso_solve_init(handle_lhs, ndof, aa_lhs, ia_lhs, ja_lhs, &
                          a_constant=.FALSE., mtype=11, msglvl=0, error=error, &
                          lu_reuse_L=MERGE(NINT(-LOG10(lu_reuse_tol)), 0, i_lu_reuse == -1))
  IF (error /= 0) THEN
    WRITE(*,*) 'timefp_7pt: pardiso_solve_init failed, error=', error; STOP
  END IF

  !================================================================
  ! 5.  Time loop
  !================================================================
  
  ! Compute initial solution density
  do ix=1,nbig
      call index_mat_inv(ix,iv,imu)
      f_init(iv,imu)=fstart(ix)
  enddo
      
   CALL time_density(f_init, dens_tmp)
   write(*,*) 'Initial density is ',dens_tmp

   call time_energy(f_init, dens_tmp, teff=teff_tmp)
   write(*,*) 'Initial effective temperature is ',teff_tmp


  phase_offset  = 0.0_dp
  itime_global  = 0

  ! Steady-state convergence test (i_ss_check=-1): weights depend only on the
  ! grid, so they are built once here and reused for every step.
  ! nu_ref = 1/tauie (ion-electron collision rate) renders eps dimensionless.
  IF (i_ss_check == -1) &
    CALL conv_diag_init(vteff, 1.0_dp/MAX(tauie, TINY(1.0_dp)), write_hist=.TRUE.)

  phase_loop: DO iphase = 1, 3
    IF (ntimes(iphase) == 0) CYCLE phase_loop
    timestep_cur = timestep(iphase)
    IF (iphase > 1) WRITE(*,'(A,I0,A,ES12.4,A)') &
        '  Phase ', iphase, ': dt = ', timestep_cur, ' s'

  time_loop: DO itime = 1, ntimes(iphase)
    itime_global = itime_global + 1
    time = otime + phase_offset + itime*timestep_cur
   ! WRITE(*,*) 'Time is ', time, ' s'

    !--- Update Coulomb log and rebuild linear operator each step ----
    lnab_arr = 0.0_dp
    IF (.NOT. (isource == -1 .AND. iold == 0 .AND. dens_tmp < 0.05d0 * npart)) THEN
      DO iv = 1, nperp
        DO imu = 1, npar
          f_init(iv,imu) = fstart(index_mat(iv,imu))
        END DO
      END DO
      CALL time_energy(f_init, dens_tmp, teff=teff)
      ta_eV = teff * 1.0d3
      DO ib = 2, nbulk
        CALL coulomb_log_ab(za, aa, ta_eV, npart, &
                            zb(ib-1), ab(ib-1), t(ib), nb(ib), lnab_t)
        lnab_arr(ib) = lnab_t
        cte0_t     = gamma0 * lnab_t * (za/aa)**2
        gammab(ib) = cte0_t * nb(ib) * zb(ib-1)**2
      END DO
      CALL assemble_FP_terms(all00, all10, all01, all20, all11, all02)
    ELSE
      all00 = all00_lin;  all10 = all10_lin;  all01 = all01_lin
      all11 = all11_lin;  all20 = all20_lin;  all02 = all02_lin
    END IF

    !--- Rebuild aa_L values (ja_L pattern unchanged) ---------------
    ! Rows are filled in parallel: row (i,j) starts at ia_L(row) of the fixed
    ! pattern, so no running pointer is needed.  A row whose entry count does
    ! not match the pattern would corrupt the matrix: stop rather than go on.
    pattern_ok = .TRUE.
    !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) &
    !$OMP   PRIVATE(i, j, k, row, col_idx, stencil_coeff, n_entries, rhs_ij) &
    !$OMP   REDUCTION(.AND.:pattern_ok)
    DO i = 1, nperp
      DO j = 1, npar
        row = (i-1)*npar + j
        CALL fd_stencil_2d(i, j, nperp, npar, vperp, dvpar, &
                           all00(i,j), all10(i,j), all01(i,j), &
                           all20(i,j), all11(i,j), all02(i,j), &
                           col_idx, stencil_coeff, n_entries, rhs_ij)
        CALL sort_stencil(col_idx, stencil_coeff, n_entries)
        IF (n_entries == ia_L(row+1) - ia_L(row)) THEN
          DO k = 1, n_entries
            aa_L(ia_L(row) + k - 1) = stencil_coeff(k)
          END DO
        ELSE
          pattern_ok = .FALSE.
        END IF
      END DO
    END DO
    !$OMP END PARALLEL DO
    IF (.NOT. pattern_ok) THEN
      WRITE(*,'(A,I0)') '  *** PATTERN MISMATCH in the matrix refill at step ', itime
      STOP
    END IF

    !--- Rebuild aa_lhs = I - theta*dt*L (with BC restoration) -----
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(row, ptr, iv, imu)
    DO ptr = 1, nnz_L
      aa_lhs(ptr) = -theta * timestep_cur * aa_L(ptr)
    END DO
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(row, ptr, iv, imu)
    DO row = 1, ndof
      DO ptr = ia_lhs(row), ia_lhs(row+1)-1
        IF (ja_lhs(ptr) == row) THEN
          aa_lhs(ptr) = aa_lhs(ptr) + 1.0_dp
          EXIT
        END IF
      END DO
    END DO
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(row, ptr, iv, imu)
    DO row = 1, ndof
      CALL index_mat_inv(row, iv, imu)
      IF (iv == 1 .OR. iv == nperp .OR. imu == 1 .OR. imu == npar) THEN
        DO ptr = ia_lhs(row), ia_lhs(row+1)-1
          aa_lhs(ptr) = aa_L(ptr)
        END DO
      END IF
    END DO
    !$OMP END PARALLEL DO

    !--- Build RHS = (I + (1-theta)*dt*L) * f^n + dt * S ----------
    !
    !    For CN (theta=0.5):     rhs = (I + dt/2 * L)*f^n + dt*S
    !    For implicit (theta=1): rhs =             f^n    + dt*S
    !
    !    Step 1: compute Lf^n = L * fstart  (sparse mat-vec)
    !    Step 2: rhs = fstart + (1-theta)*dt * Lf + dt * source_v

    CALL sparse_matvec_csr(ndof, ia_L, ja_L, aa_L, fstart, Lf)

    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(row, ptr, iv, imu)
    DO row = 1, ndof
      rhs_vec(row) = fstart(row) &
                   + (1.0_dp - theta) * timestep_cur * Lf(row) &
                   + timestep_cur * source_v(row)
    END DO
    !$OMP END PARALLEL DO

    ! BC rows are constraints: RHS must be zero so the solver enforces
    ! the BC exactly (df/dvp=0 or f=0) at every time step.
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(row, ptr, iv, imu)
    DO row = 1, ndof
      CALL index_mat_inv(row, iv, imu)
      IF (iv == 1 .OR. iv == nperp .OR. imu == 1 .OR. imu == npar) THEN
        rhs_vec(row) = 0.0_dp
      END IF
    END DO
    !$OMP END PARALLEL DO

    !--- Solve M_lhs * f^{n+1} = rhs  (phase 33 only) -------------
    CALL pardiso_solve_step(handle_lhs, aa_lhs, ia_lhs, ja_lhs, &
                            rhs_vec, x_vec, &
                            a_changed=.TRUE., error=error)
    IF (error /= 0) THEN
      WRITE(*,*) 'timefp_7pt: pardiso_solve_step failed, error=', error; STOP
    END IF

    !--- Unpack solution into fout ---------------------------------
    DO iv = 1, nperp
      DO imu = 1, npar
        ix = index_mat(iv, imu)
        fout(iv,imu) = x_vec(ix)
      END DO
    END DO

    !--- Optional VDF snapshot (n_snap > 0; off by default) ---------
    !    See write_vdf_snapshot (time_comps_mod).
    CALL write_vdf_snapshot(itime_global, time, fout)

    !--- Diagnostics (identical to TimeFP3) -----------------------
    CALL time_density(fout, dens_tmp)
  !  WRITE(*,*)  'Unnormalised density is ', dens_tmp
    WRITE(45,*) time, dens_tmp
    
    CALL time_energy(fout, dens_tmp, tk, tkperp, tkpar, teff)
    WRITE(46,*) time, tk, tkperp
    anisotropy = merge(100.0_dp*tkperp/(2.0_dp*tk - tkperp), 0.0_dp, tk > 0.0_dp)
    WRITE(47,*) time, anisotropy
    WRITE(507,*) time, teff


    !--- Refresh Teff-dependent coefficients from the just-solved f^{n} ------
    ! The operator used for the time step is built from the start-of-step
    ! (lagged) Teff.  Reporting the collisional/self-collision power with that
    ! lagged operator produces a transient in the power whenever the Teff jumps
    ! a lot in one step -- e.g. at a dt change.  Rebuild gammab from the
    ! post-step teff (and dens_tmp) above so the reported power is consistent
    ! with the current f.
    ! Only the diagnostic coefficients are touched; all** and the solution are
    ! unchanged, and the next step rebuilds everything from fstart regardless.
    ! (This solver runs with isc = 0 only: no self-collision coefficients.)
    IF (.NOT. (isource == -1 .AND. iold == 0 .AND. dens_tmp < 0.05d0 * npart)) THEN
      ta_eV = teff * 1.0d3
      DO ib = 2, nbulk
        CALL coulomb_log_ab(za, aa, ta_eV, npart, &
                            zb(ib-1), ab(ib-1), t(ib), nb(ib), lnab_t)
        lnab_arr(ib) = lnab_t
        gammab(ib)   = gamma0 * lnab_t * (za/aa)**2 * nb(ib) * zb(ib-1)**2
      END DO
    END IF

    CALL time_power_7pt(x_vec, dens_tmp, pcoll, pRF, psource, plosses, &
                        pcoll_self, pcoll_self_perp, pcoll_self_par, &
                        tau_rf, tau_ii, tau_ie)

    IF (iplot_pow == -1) THEN
      WRITE(470,*) time, (SUM(pcoll)+pcoll_self)/1.d6
      DO ib = 1, nbulk
        WRITE(470+ib,*) time, pcoll(ib)/1.d6
      END DO
      IF (irf   == -1) WRITE(480,*) time, pRF/1.d6
      IF (irf   == -1) WRITE(515,*) time, tau_rf
      WRITE(516,*) time, tau_ii, tau_ie
      IF (isource==-1) WRITE(490,*) time, psource/1.d6, plosses/1.d6
      IF (tau_E > 0.d0) THEN
        CALL time_power_lossE(x_vec, dens_tmp, pE_loss, pE_ratio)
        WRITE(495,*) time, pE_loss/1.d6, pE_ratio
      END IF
      !DTERMS-OFF! CALL time_density_terms_7pt(x_vec, ncoll_d, nsc_d, nRF_d, nsrc_d, nloss_d)
      !DTERMS-OFF! WRITE(520,'(ES16.8,*(ES16.7))') time, SUM(ncoll_d)+nsc_d+nRF_d+nsrc_d+nloss_d, &
                                     !DTERMS-OFF! ncoll_d, nsc_d, nRF_d, nsrc_d, nloss_d
    END IF
    IF (nbulk >   1) WRITE(505,*) time, (lnab_arr(ib), ib=2,nbulk)
    IF (iplot_mom == -1) THEN
      CALL time_momentum_7pt(x_vec, dens_tmp, &
                             mcoll_perp, mcoll_par, &
                             mRF_perp,   mRF_par,   &
                             msrc_perp,  msrc_par,  &
                             mloss_perp, mloss_par, &
                             mSC_perp,   mSC_par)
      WRITE(570,*) time, SUM(mcoll_perp)+mSC_perp, SUM(mcoll_par)+mSC_par
      DO ib = 1, nbulk
        WRITE(570+ib,*) time, mcoll_perp(ib), mcoll_par(ib)
      END DO
      IF (irf   == -1) WRITE(580,*) time, mRF_perp, mRF_par
      IF (isource==-1) WRITE(590,'(5ES15.7)') time, msrc_perp, msrc_par, mloss_perp, mloss_par
    END IF

    !--- Steady-state convergence test (single, Jacobian-weighted) --
    ! Rolling window of n_ss_window steps; fout holds f^n.
    if (i_ss_check == -1) then
      call conv_diag_step(itime_global, time, fout, cdiag)
      if (conv_diag_converged(cdiag)) then
        write(*,'(A,F12.5,A)') '  Stopping at t=', time, ' s (conv_diag criteria met).'
        exit phase_loop
      end if
    end if

    !--- Advance solution -----------------------------------------
    fstart = x_vec!*npart/dens_tmp

  END DO time_loop

    phase_offset = phase_offset + ntimes(iphase) * timestep_cur

  END DO phase_loop

  !================================================================
  ! 6.  Finalise
  !================================================================
  IF (i_ss_check == -1) CALL conv_diag_finalize()

  CALL pardiso_solve_finalize(handle_lhs, ia_lhs, ja_lhs, error)
  WRITE(*,*) 'Solve completed.'

  IF (iplot_pow == -1) THEN
    IF (irf     == -1) CLOSE(480)
    IF (irf     == -1) CLOSE(515)
    CLOSE(516)
    IF (isource == -1) CLOSE(490)
    IF (tau_E > 0.d0) CLOSE(495)
    !DTERMS-OFF! CLOSE(520)
    CLOSE(470)
    DO ib = 1, nbulk
      CLOSE(470+ib)
    END DO
  END IF
  IF (nbulk   >   1) CLOSE(505)
  CLOSE(507)
  CLOSE(47); CLOSE(46); CLOSE(45)
  IF (iplot_mom == -1) THEN
    IF (irf     == -1) CLOSE(580)
    IF (isource == -1) CLOSE(590)
    CLOSE(570)
    DO ib = 1, nbulk
      CLOSE(570+ib)
    END DO
  END IF

  !================================================================
  ! 7.  Renormalise (sourceless case)
  !================================================================
  IF (isource == 0) THEN
    WRITE(*,*) 'Renormalizing...'
    fout = fout * npart / dens_tmp
  END IF

  !================================================================
  ! 8.  Write output files (same as TimeFP3)
  !================================================================
  OPEN(40, file=TRIM(outfile('fout.txt')), status='unknown')
  OPEN(42, file=TRIM(outfile('xout.dat')), status='unknown')
  ! Header: time and the total step count, so a restart can continue the
  ! snapshot numbering (main reads "time" alone from older files).
  WRITE(42,*) time, nstep_restart + itime_global
  DO iv = 1, nperp
    DO imu = 1, npar
      WRITE(40,*) vperp(iv), vpar(imu), fout(iv,imu)
      ix = index_mat(iv, imu)
      WRITE(42,*) x_vec(ix)
    END DO
  END DO
  CLOSE(42); CLOSE(40)


  DEALLOCATE(ia_L, ja_L, aa_L, ia_lhs, ja_lhs, aa_lhs, rhs_vec, x_vec, Lf)

CONTAINS

  !----------------------------------------------------------------
  ! Sparse matrix-vector product  y = A_csr * x
  ! CSR format: ia, ja, aa  (1-based)
  !----------------------------------------------------------------
  SUBROUTINE sparse_matvec_csr(n, ia, ja, aa, x, y)
    INTEGER,  INTENT(IN)  :: n, ia(n+1), ja(:)
    REAL(dp), INTENT(IN)  :: aa(:), x(n)
    REAL(dp), INTENT(OUT) :: y(n)
    INTEGER :: row, ptr
    y = 0.0_dp
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(row, ptr)
    DO row = 1, n
      DO ptr = ia(row), ia(row+1)-1
        y(row) = y(row) + aa(ptr) * x(ja(ptr))
      END DO
    END DO
    !$OMP END PARALLEL DO
  END SUBROUTINE sparse_matvec_csr

  !----------------------------------------------------------------
  ! Insertion sort of stencil entries by column index
  ! (required by PARDISO: entries in each row must be sorted)
  !----------------------------------------------------------------
  SUBROUTINE sort_stencil(idx, val, n)
    INTEGER,  INTENT(INOUT) :: idx(n)
    REAL(dp), INTENT(INOUT) :: val(n)
    INTEGER,  INTENT(IN)    :: n
    INTEGER  :: p, q, tmp_i
    REAL(dp) :: tmp_v
    DO p = 2, n
      tmp_i = idx(p);  tmp_v = val(p)
      q = p - 1
      DO WHILE (q >= 1 .AND. idx(q) > tmp_i)
        idx(q+1) = idx(q);  val(q+1) = val(q)
        q = q - 1
      END DO
      idx(q+1) = tmp_i;  val(q+1) = tmp_v
    END DO
  END SUBROUTINE sort_stencil

END SUBROUTINE timefp_7pt

END MODULE mod_timefp_7pt
