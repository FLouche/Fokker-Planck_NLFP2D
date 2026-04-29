SUBROUTINE solve_fp_pardiso(A, B, C, D, E, F,          &
                             steady_no_beam,             &
                             f_sol)
  !-----------------------------------------------------------------
  ! Assembles and solves the 2D Fokker-Planck equation:
  !
  !   A*f + B*df/dvp + C*df/dva + D*d2f/dvp2
  !       + E*d2f/dvp_dva + F*d2f/dva2 = 0
  !
  ! Boundary conditions:
  !   i=1      (vp=0)      : df/dvp = 0  (Neumann)
  !   i=nperp  (vp=vpmax)  : f = 0       (Dirichlet)
  !   j=1,npar (vpa=+-max) : f = 0       (Dirichlet)
  !
  ! Solvability constraint (activated when steady_no_beam=.TRUE.):
  !   f(i_ref, j_ref) = 1  (fixes the normalisation of the null-space
  !   solution; the caller should rescale f_sol afterward to enforce
  !   the physical normalisation, e.g. integral f dv = n).
  !   i_ref, j_ref are compile-time PARAMETERs defined below.
  !
  ! Input:
  !   steady_no_beam : .TRUE.  -> solvability constraint applied
  !                    .FALSE. -> standard solve, no constraint
  !-----------------------------------------------------------------
 
    use shared_grid
    use shared_beam
    
    USE mod_fd_stencil_2d, ONLY: fd_stencil_2d
  
    IMPLICIT NONE

  INTEGER, PARAMETER :: dp = KIND(1.0D0)
  
  !--- compile-time reference point for solvability constraint -----
  ! v∥ index: defined in grid.f90 --> j_ref= jmid

  REAL(dp), INTENT(IN)  :: A(nperp,npar), B(nperp,npar), C(nperp,npar)
  REAL(dp), INTENT(IN)  :: D(nperp,npar), E(nperp,npar), F(nperp,npar)
  LOGICAL,  INTENT(IN)  :: steady_no_beam
  REAL(dp), INTENT(OUT) :: f_sol(nperp,npar)

  !--- sparse matrix (CSR, 1-based) --------------------------------
  INTEGER  :: ndof, nnz, nnz_max, ptr, row_ref
  INTEGER,  ALLOCATABLE :: ia(:), ja(:)
  REAL(dp), ALLOCATABLE :: aa(:), rhs_vec(:), x_vec(:)

  !--- stencil workspace -------------------------------------------
  INTEGER  :: col_idx(49)
  REAL(dp) :: stencil_coeff(49)
  INTEGER  :: n_entries
  REAL(dp) :: rhs_ij

  !--- loop indices ------------------------------------------------
  INTEGER  :: i, j, row, k

  !--- PARDISO variables -------------------------------------------
  INTEGER(8) :: pt(64)
  INTEGER    :: iparm(64)
  INTEGER    :: mtype, phase, nrhs, msglvl, error, maxfct, mnum, idum
  REAL(dp)   :: ddum

  !================================================================
  ! 0. Sizes
  !================================================================
  ndof    = nperp * npar
  nnz_max = ndof * 49
  row_ref = (imid-1)*npar + jmid 
  write(*,*) 'row_ref = ',row_ref! global row of the constraint

  ALLOCATE(ia(ndof+1), ja(nnz_max), aa(nnz_max))
  ALLOCATE(rhs_vec(ndof), x_vec(ndof))

  if(steady_no_beam) then
      rhs_vec = 0.0_dp
  else
      DO i = 1, nperp
        DO j = 1, npar
             rhs_vec((i-1)*npar + j) = -source(i,j)
        END DO
      END DO
  endif
  
  x_vec   = 0.0_dp

  !================================================================
  ! 1. First pass: count non-zeros per row -> build ia(:)
  !================================================================
  ia(1) = 1
  DO i = 1, nperp
    DO j = 1, npar
      row = (i-1)*npar + j

      !--- solvability constraint row: single diagonal entry -------
      IF (steady_no_beam .AND. row == row_ref) THEN
        ia(row+1) = ia(row) + 1
        CYCLE
      END IF

      CALL fd_stencil_2d(i, j, nperp, npar, vperp, dvpar, &
                   A(i,j), B(i,j), C(i,j),           &
                   D(i,j), E(i,j), F(i,j),           &
                   col_idx, stencil_coeff, n_entries, rhs_ij)
      
      ia(row+1) = ia(row) + n_entries
    END DO
  END DO

  nnz = ia(ndof+1) - 1
  IF (nnz > nnz_max) THEN
    WRITE(*,*) 'solve_fp_pardiso: nnz_max exceeded, nnz=', nnz
    STOP
  END IF

  !================================================================
  ! 2. Second pass: fill ja(:), aa(:), rhs_vec(:)
  !================================================================
  ptr = 1
  DO i = 1, nperp
    DO j = 1, npar
      row = (i-1)*npar + j

      !--- solvability constraint: f(i_ref,j_ref) = 1 -------------
      IF (steady_no_beam .AND. row == row_ref) THEN
          write(*,*) 'Solvability constraint-> we impose f = 1 somewhere'
        ja(ptr)      = row_ref
        aa(ptr)      = 1.0_dp
        rhs_vec(row) = 1.0_dp    ! fix f at this point = 1
        ptr          = ptr + 1
        CYCLE
      END IF

      CALL fd_stencil_2d(i, j, nperp, npar, vperp, dvpar, &
                   A(i,j), B(i,j), C(i,j),           &
                   D(i,j), E(i,j), F(i,j),           &
                   col_idx, stencil_coeff, n_entries, rhs_ij)

      CALL sort_stencil(col_idx, stencil_coeff, n_entries)

      DO k = 1, n_entries
        ja(ptr) = col_idx(k)
        aa(ptr) = stencil_coeff(k)
        ptr     = ptr + 1
      END DO

      rhs_vec(row) = rhs_vec(row) + rhs_ij
    END DO
  END DO

  !================================================================
  ! 3. Initialise PARDISO
  !================================================================
  mtype  = 11
  nrhs   = 1
  msglvl = 0
  maxfct = 1
  mnum   = 1
  pt     = 0
  iparm  = 0

  iparm(1)  = 1
  iparm(2)  = 3
  iparm(3)  = 1
  iparm(4)  = 0
  iparm(5)  = 0
  iparm(6)  = 0
  iparm(8)  = 10
  iparm(10) = 13
  iparm(11) = 1
  iparm(13) = 1
  iparm(18) = -1
  iparm(19) = -1
  iparm(35) = 0

  !================================================================
  ! 4. Phase 11: reordering and symbolic factorization
  !================================================================
  phase = 11
  CALL pardiso(pt, maxfct, mnum, mtype, phase,          &
               ndof, aa, ia, ja, idum, nrhs,             &
               iparm, msglvl, ddum, ddum, error)
  IF (error /= 0) THEN
    WRITE(*,'(A,I4)') 'PARDISO phase 11 error: ', error
    CALL pardiso_error_msg(error);  STOP
  END IF

  !================================================================
  ! 5. Phase 22: numerical factorization
  !================================================================
  phase = 22
  CALL pardiso(pt, maxfct, mnum, mtype, phase,          &
               ndof, aa, ia, ja, idum, nrhs,             &
               iparm, msglvl, ddum, ddum, error)
  IF (error /= 0) THEN
    WRITE(*,'(A,I4)') 'PARDISO phase 22 error: ', error
    CALL pardiso_error_msg(error);  STOP
  END IF

  !================================================================
  ! 6. Phase 33: back-substitution
  !================================================================
  phase = 33
  CALL pardiso(pt, maxfct, mnum, mtype, phase,          &
               ndof, aa, ia, ja, idum, nrhs,             &
               iparm, msglvl, rhs_vec, x_vec, error)
  IF (error /= 0) THEN
    WRITE(*,'(A,I4)') 'PARDISO phase 33 error: ', error
    CALL pardiso_error_msg(error);  STOP
  END IF

  IF (msglvl == 1) THEN
    WRITE(*,'(A,I12)') '  Non-zeros in factors:     ', iparm(18)
    WRITE(*,'(A,I12)') '  Mflops for factorization: ', iparm(19)
  END IF

  !================================================================
  ! 7. Phase -1: release PARDISO memory
  !================================================================
  phase = -1
  CALL pardiso(pt, maxfct, mnum, mtype, phase,          &
               ndof, ddum, ia, ja, idum, nrhs,           &
               iparm, msglvl, ddum, ddum, error)

  !================================================================
  ! 8. Reshape solution -> 2D, then rescale if constraint was used
  !================================================================
  DO i = 1, nperp
    DO j = 1, npar
      f_sol(i,j) = x_vec((i-1)*npar + j)
    END DO
  END DO

  !--- Rescale so that the physical normalisation is preserved -----
  ! The constraint fixed f(i_ref,j_ref)=1 arbitrarily. The caller
  ! should rescale to enforce e.g. 2*pi * integral f vp dvp dvpa = n.
  ! A reminder is printed so the caller does not forget.
  IF (steady_no_beam) THEN
    WRITE(*,'(A)') '  solve_fp_pardiso: solvability constraint applied.'
    WRITE(*,'(A)') '  f_sol is normalised so that f(i_ref,j_ref)=1.'
    WRITE(*,'(A)') '  Caller must rescale to physical normalisation.'
  END IF

  DEALLOCATE(ia, ja, aa, rhs_vec, x_vec)

CONTAINS

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

  SUBROUTINE pardiso_error_msg(err)
    INTEGER, INTENT(IN) :: err
    SELECT CASE (err)
    CASE(-1);  WRITE(*,*) '  Input inconsistent'
    CASE(-2);  WRITE(*,*) '  Not enough memory'
    CASE(-3);  WRITE(*,*) '  Reordering problem'
    CASE(-4);  WRITE(*,*) '  Zero pivot'
    CASE(-5);  WRITE(*,*) '  Unclassified internal error'
    CASE(-6);  WRITE(*,*) '  Reordering failed'
    CASE(-7);  WRITE(*,*) '  Diagonal matrix is singular'
    CASE(-8);  WRITE(*,*) '  32-bit integer overflow'
    CASE(-9);  WRITE(*,*) '  Not enough memory for OOC'
    CASE(-10); WRITE(*,*) '  Error opening OOC files'
    CASE(-11); WRITE(*,*) '  Read/write error with OOC files'
    CASE DEFAULT; WRITE(*,'(A,I4)') '  Unknown error: ', err
    END SELECT
  END SUBROUTINE pardiso_error_msg

END SUBROUTINE solve_fp_pardiso