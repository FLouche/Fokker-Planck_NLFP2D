!*******************************************************************
!*   Resolution of the time-dependent Fokker-Planck Equation       *
!*   with Crank-Nicolson or fully implicit scheme                  *
!*   3-point (2nd-order) stencil in v⊥ and v∥                     *
!*                                                                 *
!*   Updated replacement for TimeFP3 / build_matrix_td:           *
!*   - Operator L assembled directly in CSR format (no dense       *
!*     matrix), saving O(nbig^2) memory.                           *
!*   - Crank-Nicolson: (I - dt/2*L) f^{n+1} = (I + dt/2*L) f^n   *
!*                                           + dt * S              *
!*   - Fully implicit: (I - dt*L)   f^{n+1} =            f^n      *
!*                                           + dt * S              *
!*   - LHS is constant: PARDISO phases 11+22 done once,           *
!*     only phase 33 repeated each step.                           *
!*   - Uses fd_stencil_2d_3 for 3-point 2nd-order FD weights.     *
!*                                                                 *
!*   Version 1.0 - F. Louche                                       *
!*******************************************************************

MODULE mod_timefp3_upd

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: timefp_upd

CONTAINS

SUBROUTINE timefp_upd(all00_lin, all10_lin, all01_lin, &
                      all11_lin, all20_lin, all02_lin, &
                      fstart, fout, otime)

  USE mod_fd_stencil_2d_3              ! provides fd_stencil_2d_3
  USE pardiso_solver                   ! provides pardiso_handle_t,
                                       !   pardiso_solve_init/step/finalize
  USE shared_grid
  USE shared_plasma
  USE shared_timer
  USE shared_beam
  USE shared_RF
  USE func_index
  USE mod_ss_check

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

  !--- CSR storage for operator L and LHS matrix M = I-theta*dt*L --
  INTEGER,  ALLOCATABLE :: ia_L(:), ja_L(:)
  REAL(dp), ALLOCATABLE :: aa_L(:)
  INTEGER,  ALLOCATABLE :: ia_lhs(:), ja_lhs(:)
  REAL(dp), ALLOCATABLE :: aa_lhs(:)

  INTEGER :: nnz_L, nnz_max

  !--- Stencil workspace -------------------------------------------
  INTEGER  :: col_idx(9)
  REAL(dp) :: stencil_coeff(9)
  INTEGER  :: n_entries
  REAL(dp) :: rhs_ij

  !--- Working vectors ---------------------------------------------
  REAL(dp), ALLOCATABLE :: rhs_vec(:), x_vec(:), Lf(:)

  !--- Scalars -----------------------------------------------------
  REAL(dp) :: theta
  REAL(dp) :: time, dens_tmp, tk, tkperp, tkpar, teff
  REAL(dp) :: pcoll(nbulk), pRF, psource, plosses, pcoll_self

  INTEGER :: ndof, i, j, k, row, ptr, itime, iv, imu, ix
  INTEGER :: error, ib
  CHARACTER(len=2)   :: ibString
  CHARACTER(len=256) :: dynfname

  EXTERNAL :: time_density, time_energy, time_power_7pt

  logical  :: ss_converged
  real(dp) :: p_net_ss, p_drive_ss, anisotropy

  REAL(dp), DIMENSION(nperp,npar) :: f_init

  !================================================================
  ! 0.  Setup
  !================================================================
  ndof    = nperp * npar
  nnz_max = ndof * 9          ! upper bound: 3x3 stencil per row

  IF (icn == -1) THEN
    theta = 0.5_dp             ! Crank-Nicolson
  ELSE
    theta = 1.0_dp             ! fully implicit
  END IF

  ALLOCATE(ia_L(ndof+1), ja_L(nnz_max), aa_L(nnz_max))
  ALLOCATE(ia_lhs(ndof+1), ja_lhs(nnz_max), aa_lhs(nnz_max))
  ALLOCATE(rhs_vec(ndof), x_vec(ndof), Lf(ndof))

  !================================================================
  ! 1.  Build steady-state operator L in CSR
  !     Two-pass assembly (count then fill), same pattern as
  !     timefp_7pt but calling fd_stencil_2d_3 instead of
  !     fd_stencil_2d.
  !================================================================

  !--- Pass 1: count non-zeros per row -> ia_L --------------------
  ia_L(1) = 1
  DO i = 1, nperp
    DO j = 1, npar
      row = (i-1)*npar + j
      CALL fd_stencil_2d_3(i, j, nperp, npar, vperp, dvpar, &
                            all00_lin(i,j), all10_lin(i,j), all01_lin(i,j), &
                            all20_lin(i,j), all11_lin(i,j), all02_lin(i,j), &
                            col_idx, stencil_coeff, n_entries, rhs_ij)
      ia_L(row+1) = ia_L(row) + n_entries
    END DO
  END DO

  nnz_L = ia_L(ndof+1) - 1
  IF (nnz_L > nnz_max) THEN
    WRITE(*,*) 'timefp_upd: nnz_max exceeded, nnz_L=', nnz_L; STOP
  END IF

  !--- Pass 2: fill ja_L, aa_L ------------------------------------
  ! Note: fd_stencil_2d_3 returns entries in ascending column order
  ! for any npar>=4 (natural loop order on 3x3 stencil), so no
  ! sort is needed. sort_stencil is called defensively.
  ptr = 1
  DO i = 1, nperp
    DO j = 1, npar
      CALL fd_stencil_2d_3(i, j, nperp, npar, vperp, dvpar, &
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

  WRITE(*,'(A,I10,A,F6.2,A)') '  L operator (3-pt): nnz=', nnz_L, &
      '  (', 100.d0*nnz_L/DBLE(ndof)**2, ' %)'

  !================================================================
  ! 2.  Build LHS:  M_lhs = I - theta*dt*L
  !     Sparsity pattern identical to L (diagonal already in L
  !     via the A*f term); just scale values and shift diagonal.
  !================================================================
  ia_lhs = ia_L
  ja_lhs = ja_L

  DO ptr = 1, nnz_L
    aa_lhs(ptr) = -theta * timestep * aa_L(ptr)
  END DO

  DO row = 1, ndof
    DO ptr = ia_lhs(row), ia_lhs(row+1)-1
      IF (ja_lhs(ptr) == row) THEN
        aa_lhs(ptr) = aa_lhs(ptr) + 1.0_dp
        EXIT
      END IF
    END DO
  END DO

  !================================================================
  ! 3.  Open output files (same convention as timefp / timefp_7pt)
  !================================================================
  IF (otime == 0.d0) THEN
    OPEN(45, file=TRIM(outfile('density_vs_time.txt')),        status='unknown')
    OPEN(46, file=TRIM(outfile('energy_vs_time.txt')),         status='unknown')
    OPEN(47, file=TRIM(outfile('anisotropy_vs_time.txt')),    status='unknown')
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
    IF (isource== -1) OPEN(490,file=TRIM(outfile('power_NBI_vs_time.txt')),       status='unknown')
    IF (isc    /=  0) OPEN(500,file=TRIM(outfile('power_coll_self_vs_time.txt')), status='unknown')
  ELSE
    OPEN(45, file=TRIM(outfile('density_vs_time.txt')),        status='old', access='append')
    OPEN(46, file=TRIM(outfile('energy_vs_time.txt')),         status='old', access='append')
    OPEN(47, file=TRIM(outfile('anisotropy_vs_time.txt')),    status='old', access='append')
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
    IF (isource== -1) OPEN(490,file=TRIM(outfile('power_NBI_vs_time.txt')),       status='old', access='append')
    IF (isc    /=  0) OPEN(500,file=TRIM(outfile('power_coll_self_vs_time.txt')), status='old', access='append')
  END IF

  !================================================================
  ! 4.  Factorise M_lhs once (phases 11 + 22)
  !================================================================
  CALL pardiso_solve_init(handle_lhs, ndof, aa_lhs, ia_lhs, ja_lhs, &
                          a_constant=.TRUE., mtype=11, msglvl=0, error=error)
  IF (error /= 0) THEN
    WRITE(*,*) 'timefp_upd: pardiso_solve_init failed, error=', error; STOP
  END IF

  !================================================================
  ! 5.  Time loop
  !================================================================

  ! Print initial density
  DO ix = 1, nbig
    CALL index_mat_inv(ix, iv, imu)
    f_init(iv,imu) = fstart(ix)
  END DO
  CALL time_density(f_init, dens_tmp)
  WRITE(*,*) 'Initial density is ', dens_tmp

  time_loop: DO itime = 1, ntimes

    time = otime + itime*timestep
    WRITE(*,*) 'Time is ', time, ' s'

    !--- Build RHS = (I + (1-theta)*dt*L)*f^n + dt*S ---------------
    ! CN (theta=0.5): rhs = f^n + dt/2 * L*f^n + dt*S
    ! Implicit:       rhs = f^n                 + dt*S
    CALL sparse_matvec_csr(ndof, ia_L, ja_L, aa_L, fstart, Lf)

    DO row = 1, ndof
      rhs_vec(row) = fstart(row) &
                   + (1.0_dp - theta) * timestep * Lf(row) &
                   + timestep * source_v(row)
    END DO

    !--- Solve M_lhs * f^{n+1} = rhs  (phase 33 only) -------------
    CALL pardiso_solve_step(handle_lhs, aa_lhs, ia_lhs, ja_lhs, &
                            rhs_vec, x_vec, &
                            a_changed=.FALSE., error=error)
    IF (error /= 0) THEN
      WRITE(*,*) 'timefp_upd: pardiso_solve_step failed, error=', error; STOP
    END IF

    !--- Unpack solution into fout ---------------------------------
    DO iv = 1, nperp
      DO imu = 1, npar
        ix = index_mat(iv, imu)
        fout(iv,imu) = x_vec(ix)
      END DO
    END DO

    !--- Diagnostics -----------------------------------------------
    CALL time_density(fout, dens_tmp)
    WRITE(*,*)  'Unnormalised density is ', dens_tmp
    WRITE(45,*) time, dens_tmp

    CALL time_energy(fout, dens_tmp, tk, tkperp, tkpar, teff)
    WRITE(46,*) time, tk, tkperp
    anisotropy = merge(100.0_dp*(tkperp/tk - 2.0_dp/3.0_dp)/(2.0_dp/3.0_dp), 0.0_dp, tk > 0.0_dp)
    WRITE(47,*) time, anisotropy


    CALL time_power_7pt(x_vec, dens_tmp, pcoll, pRF, psource, plosses, pcoll_self)

    WRITE(470,*) time, (SUM(pcoll)+pcoll_self)/1.d6
    DO ib = 1, nbulk
      WRITE(470+ib,*) time, pcoll(ib)/1.d6
    END DO
    IF (irf    == -1) WRITE(480,*) time, pRF/1.d6
    IF (isource== -1) WRITE(490,*) time, psource/1.d6, plosses/1.d6
    IF (isc    /=  0) WRITE(500,*) time, pcoll_self/1.d6

    !--- Steady-state convergence check (optional) ----------------
    if (i_ss_check == -1) then
      p_net_ss   = sum(pcoll(1:nbulk)) + pcoll_self + pRF + psource + plosses
      p_drive_ss = max(abs(pRF), abs(pcoll_self), abs(psource))
      do ib = 1, nbulk
        p_drive_ss = max(p_drive_ss, abs(pcoll(ib)))
      end do
      p_drive_ss = max(p_drive_ss, 1.0_dp)
      call ss_check(itime, tk, tkperp, pRF, p_net_ss, p_drive_ss, ss_converged)
      if (ss_converged) then
        write(*,'(A,F12.5,A)') '  Stopping at t=', time, ' s (steady state reached).'
        exit time_loop
      end if
    end if

    !--- Advance solution -----------------------------------------
    fstart = x_vec

  END DO time_loop

  !================================================================
  ! 6.  Finalise PARDISO and close files
  !================================================================
  CALL pardiso_solve_finalize(handle_lhs, ia_lhs, ja_lhs, error)
  WRITE(*,*) 'Solve completed.'

  IF (irf     == -1) CLOSE(480)
  IF (isource == -1) CLOSE(490)
  IF (isc     /=  0) CLOSE(500)
  CLOSE(470); CLOSE(47); CLOSE(46); CLOSE(45)

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
  WRITE(42,*) time
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
  ! Sparse matrix-vector product  y = A_csr * x  (CSR, 1-based)
  !----------------------------------------------------------------
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

  !----------------------------------------------------------------
  ! Insertion sort of stencil entries by column index
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

END SUBROUTINE timefp_upd

END MODULE mod_timefp3_upd
