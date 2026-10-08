!*******************************************************************
!*   Resolution of the time-dependent Fokker-Planck Equation       *
!*   with Crank-Nicolson or fully implicit scheme                  *
!*   7-point Fornberg stencil in both v⊥ and v∥                    *
!*   Non-linear version: self-collision operator updated each step  *
!*                                                                 *
!*   Differences from timefp_7pt (isc /= -1):                     *
!*   - PARDISO phase 11 (symbolic reordering) done once only.      *
!*   - At each step, main_nlterm(f^n) returns sc**, which are      *
!*     added to all**_lin to form the total operator for that step. *
!*   - aa_L and aa_lhs are rebuilt each step; phases 22+33 are     *
!*     called every step via pardiso_solve_step(a_changed=.TRUE.).  *
!*   - the phi-distance kernel is computed once before the time    *
!*     loop and applied in COMPRESSED form by mod_phi_kernel; it   *
!*     used to be expanded into the dense sum_phi (11.9 GiB at     *
!*     200x200, streamed from RAM on every step).                  *
!*                                                                 *
!*   Version 1.1 - F. Louche                                       *
!*******************************************************************

MODULE mod_timefp_7pt_nl

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: timefp_7pt_nl

CONTAINS

SUBROUTINE timefp_7pt_nl(all00_lin, all10_lin, all01_lin, &
                          all11_lin, all20_lin, all02_lin, &
                          fstart, fout, otime)

  USE mod_fd_stencil_2d
  USE pardiso_solver
  USE shared_grid
  USE shared_plasma
  USE shared_FPterms
  USE shared_timer
  USE shared_beam
  USE shared_RF
  USE func_index
  USE nlterm
  USE mod_grid
  USE mod_conv_diag
  USE time_comps_mod
  USE assemble_FP_lin
  USE coulomb_log_mod
  USE mod_phi_kernel

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
  TYPE(pardiso_handle_t) :: handle_lhs

  !--- CSR storage for total operator L and LHS M = I-theta*dt*L --
  INTEGER,  ALLOCATABLE :: ia_L(:), ja_L(:)
  REAL(dp), ALLOCATABLE :: aa_L(:)
  INTEGER,  ALLOCATABLE :: ia_lhs(:), ja_lhs(:)
  REAL(dp), ALLOCATABLE :: aa_lhs(:)

  INTEGER :: nnz_L, nnz_max

  !--- Stencil workspace -------------------------------------------
  INTEGER  :: col_idx(49)
  REAL(dp) :: stencil_coeff(49)
  INTEGER  :: n_entries
  REAL(dp) :: rhs_ij

  !--- Working vectors ---------------------------------------------
  REAL(dp), ALLOCATABLE :: rhs_vec(:), x_vec(:), Lf(:)

  !--- Total coefficients (linear + self-collision, updated each step) ---
  REAL(dp), DIMENSION(nperp,npar) :: all00, all10, all01, all11, all20, all02

  !--- Scalars and temporaries -------------------------------------
  REAL(dp) :: theta
  REAL(dp) :: time, dens_tmp, tk, tkperp, tkpar,teff,teff_tmp
  REAL(dp) :: Tn          ! density-characteristic temperature (diagnostic)
  REAL(dp) :: pcoll(nbulk), pRF, psource, plosses, pcoll_self
  REAL(dp) :: ncoll_d(nbulk), nsc_d, nRF_d, nsrc_d, nloss_d   ! dn/dt per operator term
  REAL(dp) :: tau_rf            ! RF tail formation time [s]
  REAL(dp) :: tau_ii, tau_ie    ! effective collisional times [s]
  REAL(dp) :: pcoll_self_perp, pcoll_self_par
  REAL(dp) :: t_start, t_end

  INTEGER :: ndof, i, j, k, row, ptr, itime, itime_global, iphase, iv, imu, ix
  REAL(dp) :: phase_offset
  INTEGER :: error, ib
  CHARACTER(len=2)   :: ibString
  CHARACTER(len=256) :: dynfname

  !--- phi-distance kernel cache (phi_kern-CASENAME.dat) ------------
  CHARACTER(len=256) :: kernel_file
  LOGICAL            :: file_exists, need_compute
  INTEGER            :: f_nperp, f_npar
  REAL(dp)           :: f_vperp_min, f_vperp_max, f_vpar_min, f_vpar_max
  REAL(dp)           :: f_vperp_sum      ! SUM(vperp): catches ising/nsing/vbound/p_grid

  REAL(dp), PARAMETER :: gamma0 = 2.390775d-1
  REAL(dp) :: lnab_t, cte0_t, ta_eV, lnaa_t
  REAL(dp) :: lnab_arr(nbulk)
  REAL(dp) :: mcoll_perp(nbulk), mcoll_par(nbulk)
  REAL(dp) :: mRF_perp, mRF_par, msrc_perp, msrc_par
  REAL(dp) :: mloss_perp, mloss_par, mSC_perp, mSC_par

  EXTERNAL :: time_power_7pt, time_momentum_7pt, time_density_terms_7pt

  real(dp) :: anisotropy
  type(conv_diag_t) :: cdiag        ! Jacobian-weighted convergence diagnostics

  REAL(dp), DIMENSION(nperp,npar) :: f_init

  !================================================================
  ! 0.  Setup
  !================================================================
  ndof    = nperp * npar
  nnz_max = ndof * 49

 IF (icn == -1) THEN
      theta = 0.5_dp          ! Crank-Nicolson (may be unstable with NL SC)
  ELSE IF (icn == 1) THEN
      theta = 0.75_dp         ! intermediate
  ELSE
      theta = 1.0_dp          ! fully implicit
  END IF
  
  ALLOCATE(ia_L(ndof+1), ja_L(nnz_max), aa_L(nnz_max))
  ALLOCATE(ia_lhs(ndof+1), ja_lhs(nnz_max), aa_lhs(nnz_max))
  ALLOCATE(rhs_vec(ndof), x_vec(ndof), Lf(ndof))

  !================================================================
  ! 1.  phi-distance kernel (grid-dependent, expensive O(nbig^2)).
  !     Cached per casename in phi_kern-CASENAME.dat, whose header stores
  !     the grid signature (nperp, npar, vperp/vpar min/max).  The kernel
  !     is loaded only when that file exists AND its signature matches the
  !     current namelist; otherwise it is (re)computed and saved.  Fully
  !     self-managing:
  !       - file missing               -> compute + save
  !       - file present, grid matches -> load
  !       - file present, grid differs -> compute + save
  !================================================================
  ! The cache file is named phi_kern.dat, NOT sum_phi.dat.  The old file holds
  ! the dense matrix and carries the same grid signature as this one, so a
  ! stale sum_phi.dat would pass the signature check and then be read as a
  ! kernel -- silently, since the first npar*nperp^2 doubles of the dense
  ! matrix are perfectly readable numbers.  Renaming makes old caches simply
  ! not found.  Any leftover sum_phi-*.dat can be deleted; they are ~11.9 GiB
  ! each at 200x200 and are no longer used.
  ALLOCATE(phi_kern(0:npar-1, nperp, nperp))
  kernel_file  = TRIM(outfile('phi_kern.dat'))
  need_compute = .TRUE.
  INQUIRE(file=kernel_file, exist=file_exists)
  IF (file_exists) THEN
    OPEN(55, file=kernel_file, status='old', form='unformatted', access='stream', iostat=error)
    IF (error == 0) THEN
      ! SUM(vperp) closes a hole in the signature: nperp and the box limits do
      ! not change with ising, nsing, vbound or p_grid, so a kernel for another
      ! v_perp grid used to load as 'matching'.  An older cache has a kernel
      ! value where the sum now sits, fails the test and is simply recomputed.
      READ(55, iostat=error) f_nperp, f_npar, f_vperp_min, f_vperp_max, f_vpar_min, f_vpar_max, &
                             f_vperp_sum
      IF (error == 0 .AND. f_nperp == nperp .AND. f_npar == npar .AND. &
          f_vperp_min == vperp_min .AND. f_vperp_max == vperp_max .AND. &
          f_vpar_min  == vpar_min  .AND. f_vpar_max  == vpar_max  .AND. &
          f_vperp_sum == SUM(vperp)) THEN
        READ(55, iostat=error) phi_kern
        IF (error == 0) THEN
          need_compute = .FALSE.
          WRITE(*,'(A)') '  phi-distance kernel loaded from '//TRIM(kernel_file)//' (grid matches).'
        ELSE
          WRITE(*,'(A)') '  '//TRIM(kernel_file)//' is truncated -- recomputing kernel.'
        END IF
      ELSE
        WRITE(*,'(A)') '  grid signature in '//TRIM(kernel_file)//' differs -- recomputing kernel.'
      END IF
      CLOSE(55)
    END IF
  ELSE
    WRITE(*,'(A)') '  '//TRIM(kernel_file)//' not found -- computing kernel.'
  END IF

  IF (need_compute) THEN
    WRITE(*,*) 'Computing phi-distance kernel for non-linear self-collisions...'
    CALL cpu_time(t_start)
    CALL distance_v_gauss_legendre
    CALL cpu_time(t_end)
    WRITE(*,'(A,F10.3,A)') '  Done. CPU time = ', t_end - t_start, ' s'
    OPEN(55, file=kernel_file, status='replace', form='unformatted', access='stream')
    WRITE(55) nperp, npar, vperp_min, vperp_max, vpar_min, vpar_max, SUM(vperp)
    WRITE(55) phi_kern
    CLOSE(55)
    WRITE(*,'(A)') '  phi-distance kernel saved to '//TRIM(kernel_file)
  END IF

  ! Circulant transform of the kernel, used by the FFT matvec.  Cheap
  ! (nperp^2 transforms of length phi_m) and done once, whether the kernel was
  ! just computed or loaded from cache.
  CALL cpu_time(t_start)
  CALL phi_kernel_transform
  CALL cpu_time(t_end)
  WRITE(*,'(A,I0,A,F8.3,A)') '  phi kernel transformed (circulant length ', &
       phi_m, '), CPU time = ', t_end - t_start, ' s'
  WRITE(*,'(A,F10.1,A,F10.1,A)') '  kernel storage: ', &
       DBLE(npar)*DBLE(nperp)**2*8.d0/1024.d0**2, ' MiB + ', &
       DBLE(nperp)**2*DBLE(phi_nf)*8.d0/1024.d0**2, &
       ' MiB transform  (dense sum_phi would be '// &
       'nbig^2*8 B)'

  !================================================================
  ! 2.  Build sparsity pattern of L using linear coefficients.
  !     The pattern is fixed: sc** enters through the same stencil
  !     positions as all**_lin, so ia_L/ja_L never change.
  !================================================================

  !--- Pass 1: count non-zeros per row ----------------------------
  ! Use sentinel 1.0 for E_ij when all11_lin=0 so the pattern always
  ! includes mixed-derivative off-diagonal entries.  sc11 is added in
  ! step 5b and may be non-zero there even when all11_lin==0; without
  ! the sentinel, step 5c would generate more entries than ja_L holds.
  ia_L(1) = 1
  DO i = 1, nperp
    DO j = 1, npar
      row = (i-1)*npar + j
      CALL fd_stencil_2d(i, j, nperp, npar, vperp, dvpar, &
                         all00_lin(i,j), all10_lin(i,j), all01_lin(i,j), &
                         all20_lin(i,j), &
                         MERGE(all11_lin(i,j), 1.0_dp, all11_lin(i,j) /= 0.0_dp), &
                         all02_lin(i,j), &
                         col_idx, stencil_coeff, n_entries, rhs_ij)
      ia_L(row+1) = ia_L(row) + n_entries
    END DO
  END DO

  nnz_L = ia_L(ndof+1) - 1
  IF (nnz_L > nnz_max) THEN
    WRITE(*,*) 'timefp_7pt_nl: nnz_max exceeded, nnz_L=', nnz_L; STOP
  END IF

  !--- Pass 2: fill ja_L and aa_L ---------------------------------
  ptr = 1
  DO i = 1, nperp
    DO j = 1, npar
      CALL fd_stencil_2d(i, j, nperp, npar, vperp, dvpar, &
                         all00_lin(i,j), all10_lin(i,j), all01_lin(i,j), &
                         all20_lin(i,j), &
                         MERGE(all11_lin(i,j), 1.0_dp, all11_lin(i,j) /= 0.0_dp), &
                         all02_lin(i,j), &
                         col_idx, stencil_coeff, n_entries, rhs_ij)
      CALL sort_stencil(col_idx, stencil_coeff, n_entries)
      DO k = 1, n_entries
        ja_L(ptr) = col_idx(k)
        aa_L(ptr) = stencil_coeff(k)
        ptr = ptr + 1
      END DO
    END DO
  END DO

  !WRITE(*,'(A,I10,A,F6.2,A)') '  L operator (7-pt NL): nnz=', nnz_L, &
   !   '  (', 100.d0*nnz_L/DBLE(ndof)**2, ' %)'

  !--- Copy sparsity pattern to LHS arrays (values filled per step)
  ia_lhs = ia_L
  ja_lhs = ja_L

  !================================================================
  ! 3.  PARDISO phase 11 only (symbolic reordering, once).
  !     Pass a_constant=.FALSE. so phase 22 is deferred to each step.
  !     The aa_lhs values here are just to give PARDISO a valid array;
  !     the actual values are updated before every solve.
  !================================================================
  aa_lhs = -theta * timestep_cur * aa_L
  DO row = 1, ndof
    DO ptr = ia_lhs(row), ia_lhs(row+1)-1
      IF (ja_lhs(ptr) == row) THEN
        aa_lhs(ptr) = aa_lhs(ptr) + 1.0_dp
        EXIT
      END IF
    END DO
  END DO

  CALL pardiso_solve_init(handle_lhs, ndof, aa_lhs, ia_lhs, ja_lhs, &
                          a_constant=.FALSE., mtype=11, msglvl=0, error=error, &
                          lu_reuse_L=MERGE(NINT(-LOG10(lu_reuse_tol)), 0, i_lu_reuse == -1))
  IF (error /= 0) THEN
    WRITE(*,*) 'timefp_7pt_nl: pardiso_solve_init failed, error=', error; STOP
  END IF

  !================================================================
  ! 4.  Open output files
  !================================================================
  IF (otime == 0.d0) THEN
    OPEN(45, file=TRIM(outfile('density_vs_time.txt')),     status='unknown')
    OPEN(46, file=TRIM(outfile('energy_vs_time.txt')),      status='unknown')
    OPEN(47, file=TRIM(outfile('anisotropy_vs_time.txt')), status='unknown')
    IF (nbulk > 1) OPEN(505, file=TRIM(outfile('coulomb_log_vs_time.txt')), status='unknown')
    IF (isc /= 0)  OPEN(506, file=TRIM(outfile('coulomb_log_self_vs_time.txt')), status='unknown')
                   OPEN(507, file=TRIM(outfile('Teff_vs_time.txt')),              status='unknown')
                   OPEN(514, file=TRIM(outfile('Tn_vs_time.txt')),                status='unknown')
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
      IF (irf    == -1) OPEN(480,file=TRIM(outfile('power_RF_vs_time.txt')),        status='unknown')
      IF (irf    == -1) OPEN(515,file=TRIM(outfile('tau_rf_vs_time.txt')),          status='unknown')
      OPEN(516,file=TRIM(outfile('tau_coll_vs_time.txt')),        status='unknown')
      IF (isource== -1) OPEN(490,file=TRIM(outfile('power_NBI_vs_time.txt')),       status='unknown')
      IF (isc /= 0)     OPEN(500,file=TRIM(outfile('power_coll_self_vs_time.txt')), status='unknown')
      OPEN(520,file=TRIM(outfile('density_terms_vs_time.txt')),   status='unknown')
      WRITE(520,'(A)') '# time  total  coll(1:nbulk; 1=e)  self-coll  RF  source  losses   [m^-3 s^-1]'
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
      IF (irf    == -1) OPEN(580,file=TRIM(outfile('momentum_RF_vs_time.txt')),        status='unknown')
      IF (isource== -1) OPEN(590,file=TRIM(outfile('momentum_NBI_vs_time.txt')),       status='unknown')
      IF (isc /= 0)     OPEN(600,file=TRIM(outfile('momentum_coll_self_vs_time.txt')), status='unknown')
    END IF
  ELSE
    OPEN(45, file=TRIM(outfile('density_vs_time.txt')),     status='old', access='append')
    OPEN(46, file=TRIM(outfile('energy_vs_time.txt')),      status='old', access='append')
    OPEN(47, file=TRIM(outfile('anisotropy_vs_time.txt')), status='old', access='append')
    IF (nbulk > 1) OPEN(505, file=TRIM(outfile('coulomb_log_vs_time.txt')), status='old', access='append')
    IF (isc /= 0)  OPEN(506, file=TRIM(outfile('coulomb_log_self_vs_time.txt')), status='unknown', position='append')
                   OPEN(507, file=TRIM(outfile('Teff_vs_time.txt')),              status='old',     access='append')
                   OPEN(514, file=TRIM(outfile('Tn_vs_time.txt')),                status='old',     access='append')
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
      IF (irf    == -1) OPEN(480,file=TRIM(outfile('power_RF_vs_time.txt')),        status='old', access='append')
      IF (irf    == -1) OPEN(515,file=TRIM(outfile('tau_rf_vs_time.txt')),          status='old', access='append')
      OPEN(516,file=TRIM(outfile('tau_coll_vs_time.txt')),        status='old', access='append')
      IF (isource== -1) OPEN(490,file=TRIM(outfile('power_NBI_vs_time.txt')),       status='old', access='append')
      IF (isc /= 0)     OPEN(500,file=TRIM(outfile('power_coll_self_vs_time.txt')), status='old', access='append')
      OPEN(520,file=TRIM(outfile('density_terms_vs_time.txt')),   status='old', access='append')
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
      IF (irf    == -1) OPEN(580,file=TRIM(outfile('momentum_RF_vs_time.txt')),        status='old', access='append')
      IF (isource== -1) OPEN(590,file=TRIM(outfile('momentum_NBI_vs_time.txt')),       status='old', access='append')
      IF (isc /= 0)     OPEN(600,file=TRIM(outfile('momentum_coll_self_vs_time.txt')), status='old', access='append')
    END IF
  END IF

  !================================================================
  ! 5.  Time loop  (all phases in one uninterrupted loop)
  !================================================================
  DO ix = 1, nbig
    CALL index_mat_inv(ix, iv, imu)
    f_init(iv,imu) = fstart(ix)
  END DO
  CALL time_density(f_init, dens_tmp)
  WRITE(*,*) 'Initial density is ', dens_tmp

  CALL time_energy(f_init, dens_tmp, teff=teff_tmp)
  WRITE(*,*) 'Initial effective temperature is ', teff_tmp

  teff        = teff_tmp
  lnaa_t      = 0.0_dp
  phase_offset = 0.0_dp
  itime_global = 0

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

    !--- 5a. Self-collision coefficients from f^n ------------------
    ! When starting from zero (iold=0, isource=-1), skip SC while beam
    ! density is still negligible.  The Rosenbluth potential psi is then
    ! dominated by floating-point noise; regularise_axis_3 detects nearly
    ! every near-axis row as "bad" (ratio test fires on noise/noise) and
    ! overwrites phi with a spurious polynomial, producing garbage sc**.
    ! dens_tmp holds the density of fstart (initialised to 0 before the
    ! loop for istart=0, updated at the end of every step thereafter).
    IF (isource == -1 .AND. iold == 0 .AND. dens_tmp < 0.05d0 * npart) THEN
      sc00 = 0.0_dp;  sc10 = 0.0_dp;  sc01 = 0.0_dp
      sc20 = 0.0_dp;  sc11 = 0.0_dp;  sc02 = 0.0_dp
    ELSE
      CALL main_nlterm(fstart, teff,sc00, sc10, sc01, sc20, sc11, sc02)
    END IF

    !--- 5b. Update Coulomb log and recompute linear coefficients ----
    lnab_arr = 0.0_dp
    IF (.NOT. (isource == -1 .AND. iold == 0 .AND. dens_tmp < 0.05d0 * npart)) THEN
      DO iv = 1, nperp
        DO imu = 1, npar
          f_init(iv,imu) = fstart(index_mat(iv,imu))
        END DO
      END DO
      CALL time_energy(f_init, dens_tmp, teff=teff)
      ta_eV = teff * 1.0d3
      IF (isc /= 0) CALL coulomb_log_ab(za, aa, ta_eV, npart, za, aa, ta_eV, npart, lnaa_t)
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
    !--- 5b'. Add self-collision coefficients -----------------------
    all00 = all00 + sc00;  all10 = all10 + sc10;  all01 = all01 + sc01
    all11 = all11 + sc11;  all20 = all20 + sc20;  all02 = all02 + sc02

    !--- 5c. Rebuild aa_L values only (ja_L pattern unchanged) -----
    ptr = 1
    DO i = 1, nperp
      DO j = 1, npar
        CALL fd_stencil_2d(i, j, nperp, npar, vperp, dvpar, &
                           all00(i,j), all10(i,j), all01(i,j), &
                           all20(i,j), all11(i,j), all02(i,j), &
                           col_idx, stencil_coeff, n_entries, rhs_ij)
        CALL sort_stencil(col_idx, stencil_coeff, n_entries)
        DO k = 1, n_entries
          aa_L(ptr) = stencil_coeff(k)
          ptr = ptr + 1
        END DO
      END DO
    END DO
    IF (itime == 1) THEN
      IF (ptr-1 /= nnz_L) THEN
        WRITE(*,'(A,I0,A,I0)') '  *** PATTERN MISMATCH: ptr-1=', ptr-1, ' nnz_L=', nnz_L
      ELSE
        WRITE(*,'(A,I0)')      '  Pattern OK: nnz_L=', nnz_L
      END IF
   !   WRITE(*,'(A,2ES14.5)')   '  aa_L min/max:', MINVAL(aa_L(1:nnz_L)), MAXVAL(aa_L(1:nnz_L))
    END IF

    !--- 5d. Rebuild aa_lhs = I - theta*dt*L ----------------------
    DO ptr = 1, nnz_L
      aa_lhs(ptr) = -theta * timestep_cur * aa_L(ptr)
    END DO
    DO row = 1, ndof
      DO ptr = ia_lhs(row), ia_lhs(row+1)-1
        IF (ja_lhs(ptr) == row) THEN
          aa_lhs(ptr) = aa_lhs(ptr) + 1.0_dp
          EXIT
        END IF
      END DO
    END DO
    ! BC rows must remain constraints: restore them to the L stencil.
    DO row = 1, ndof
      CALL index_mat_inv(row, iv, imu)
      IF (iv == 1 .OR. iv == nperp .OR. imu == 1 .OR. imu == npar) THEN
        DO ptr = ia_lhs(row), ia_lhs(row+1)-1
          aa_lhs(ptr) = aa_L(ptr)
        END DO
      END IF
    END DO
    !IF (itime == 1) THEN
    !  WRITE(*,'(A,2ES14.5)') '  aa_lhs min/max:', MINVAL(aa_lhs(1:nnz_L)), MAXVAL(aa_lhs(1:nnz_L))
    !END IF

    !--- 5e. Build RHS: (I + (1-theta)*dt*L)*f^n + dt*S -----------
    CALL sparse_matvec_csr(ndof, ia_L, ja_L, aa_L, fstart, Lf)
    DO row = 1, ndof
      rhs_vec(row) = fstart(row) &
                   + (1.0_dp - theta) * timestep_cur * Lf(row) &
                   + timestep_cur * source_v(row)
    END DO
    ! BC rows are constraints: zero RHS so the BC is enforced exactly.
    DO row = 1, ndof
      CALL index_mat_inv(row, iv, imu)
      IF (iv == 1 .OR. iv == nperp .OR. imu == 1 .OR. imu == npar) THEN
        rhs_vec(row) = 0.0_dp
      END IF
    END DO
    !IF (itime == 1) THEN
    !  WRITE(*,'(A,2ES14.5)') '  ||fstart||, ||Lf||:', SQRT(SUM(fstart**2)), SQRT(SUM(Lf**2))
    !  WRITE(*,'(A,2ES14.5)') '  ||rhs||, ||rhs-f||:', SQRT(SUM(rhs_vec**2)), SQRT(SUM((rhs_vec-fstart)**2))
    !  WRITE(*,'(A,2ES14.5)') '  Lf min/max:', MINVAL(Lf), MAXVAL(Lf)
    !END IF

    !--- 5f. Numerical factorisation + solve (phases 22 + 33) ------
    CALL pardiso_solve_step(handle_lhs, aa_lhs, ia_lhs, ja_lhs, &
                            rhs_vec, x_vec, a_changed=.TRUE., error=error)
    IF (error /= 0) THEN
      WRITE(*,*) 'timefp_7pt_nl: pardiso_solve_step failed, error=', error; STOP
    END IF

    !--- Unpack solution into fout ---------------------------------
    DO iv = 1, nperp
      DO imu = 1, npar
        ix = index_mat(iv, imu)
        fout(iv,imu) = x_vec(ix)
      END DO
    END DO

    !--- Optional VDF snapshot (n_snap > 0; off by default) ---------
    !    Used to follow moments of f through the transient (e.g. the
    !    drift-frame anisotropy that sets the sign of the self-collision
    !    power split).  See write_vdf_snapshot (time_comps_mod).
    CALL write_vdf_snapshot(itime_global, time, fout)

    !--- Diagnostics -----------------------------------------------
    CALL time_density(fout, dens_tmp)
  !  WRITE(*,*)  'Unnormalised density is ', dens_tmp
    WRITE(45,*) time, dens_tmp

    CALL time_energy(fout, dens_tmp, tk, tkperp, tkpar,teff)
    WRITE(46,*) time, tk, tkperp
    anisotropy = merge(100.0_dp*tkperp/(2.0_dp*tk - tkperp), 0.0_dp, tk > 0.0_dp)
    WRITE(47,*) time, anisotropy
    WRITE(507,*) time, teff
    CALL time_Tn(fout, dens_tmp, Tn)
    WRITE(514,*) time, Tn

    CALL time_power_7pt(x_vec, dens_tmp, pcoll, pRF, psource, plosses, &
                        pcoll_self, pcoll_self_perp, pcoll_self_par, &
                        tau_rf, tau_ii, tau_ie)

    IF (iplot_pow == -1) THEN
      WRITE(470,*) time, (SUM(pcoll)+pcoll_self)/1.d6
      DO ib = 1, nbulk
        WRITE(470+ib,*) time, pcoll(ib)/1.d6
      END DO
      IF (irf    == -1) WRITE(480,*) time, pRF/1.d6
      IF (irf    == -1) WRITE(515,*) time, tau_rf
      WRITE(516,*) time, tau_ii, tau_ie
      IF (isource== -1) WRITE(490,*) time, psource/1.d6, plosses/1.d6
      IF (isc /= 0)     WRITE(500,*) time, pcoll_self/1.d6, &
                                         pcoll_self_perp/1.d6, pcoll_self_par/1.d6
      CALL time_density_terms_7pt(x_vec, ncoll_d, nsc_d, nRF_d, nsrc_d, nloss_d)
      WRITE(520,'(ES16.8,*(ES16.7))') time, SUM(ncoll_d)+nsc_d+nRF_d+nsrc_d+nloss_d, &
                                     ncoll_d, nsc_d, nRF_d, nsrc_d, nloss_d
    END IF
    IF (nbulk > 1)  WRITE(505,*) time, (lnab_arr(ib), ib=2,nbulk)
    IF (isc /= 0)   WRITE(506,*) time, lnaa_t

    IF (iplot_mom == -1) THEN
      CALL time_momentum_7pt(x_vec, dens_tmp, mcoll_perp, mcoll_par, &
                             mRF_perp, mRF_par, msrc_perp, msrc_par, &
                             mloss_perp, mloss_par, mSC_perp, mSC_par)
      WRITE(570,*) time, SUM(mcoll_perp)+mSC_perp, SUM(mcoll_par)+mSC_par
      DO ib = 1, nbulk
        WRITE(570+ib,*) time, mcoll_perp(ib), mcoll_par(ib)
      END DO
      IF (irf    == -1) WRITE(580,*) time, mRF_perp, mRF_par
      IF (isource== -1) WRITE(590,'(5ES15.7)') time, msrc_perp, msrc_par, mloss_perp, mloss_par
      IF (isc /= 0)     WRITE(600,*) time, mSC_perp, mSC_par
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
    fstart = x_vec

  END DO time_loop

    phase_offset = phase_offset + ntimes(iphase) * timestep_cur

  END DO phase_loop

  !================================================================
  ! 6.  Finalise
  !================================================================
  IF (i_ss_check == -1) CALL conv_diag_finalize()

  CALL pardiso_solve_finalize(handle_lhs, ia_lhs, ja_lhs, error)
  WRITE(*,*) 'Solve completed.'

  DEALLOCATE(phi_kern)
  IF (ALLOCATED(phi_khat)) DEALLOCATE(phi_khat)

  IF (iplot_pow == -1) THEN
    IF (irf     == -1) CLOSE(480)
    IF (irf     == -1) CLOSE(515)
    CLOSE(516)
    IF (isource == -1) CLOSE(490)
    IF (isc /= 0)      CLOSE(500)
    CLOSE(520)
    CLOSE(470)
    DO ib = 1, nbulk; CLOSE(470+ib); END DO
  END IF
  IF (nbulk > 1) CLOSE(505)
  IF (isc /= 0)  CLOSE(506)
                 CLOSE(507); CLOSE(514)
  CLOSE(47); CLOSE(46); CLOSE(45)
  IF (iplot_mom == -1) THEN
    IF (irf     == -1) CLOSE(580)
    IF (isource == -1) CLOSE(590)
    IF (isc /= 0)      CLOSE(600)
    CLOSE(570)
    DO ib = 1, nbulk; CLOSE(570+ib); END DO
  END IF

  !================================================================
  ! 7.  Renormalise (sourceless case)
  !================================================================
  IF (isource == 0) THEN
    WRITE(*,*) 'Renormalizing...'
    fout = fout * npart / dens_tmp
  END IF

  !================================================================
  ! 8.  Write output files
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

  SUBROUTINE sparse_matvec_csr(n, ia, ja, aa, x, y)
    INTEGER,  INTENT(IN)  :: n, ia(n+1), ja(:)
    REAL(dp), INTENT(IN)  :: aa(:), x(n)
    REAL(dp), INTENT(OUT) :: y(n)
    INTEGER :: row, ptr
    y = 0.0_dp
    DO row = 1, n
      DO ptr = ia(row), ia(row+1)-1
        y(row) = y(row) + aa(ptr) * x(ja(ptr))
      END DO
    END DO
  END SUBROUTINE sparse_matvec_csr

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

END SUBROUTINE timefp_7pt_nl

END MODULE mod_timefp_7pt_nl
