!*******************************************************************
!*   Eigenvalue analysis of the steady-state FP operator bigm_ss   *
!*                                                                 *
!*   SUBSPACE ITERATION WITH DEFLATION                              *
!*                                                                 *
!*   Finds the k eigenvalues of B closest to a user-supplied shift *
!*   sigma, by orthogonal iteration on (sigma*I - B)^(-1).          *
!*                                                                 *
!*   Method:                                                        *
!*     1. Start with k random orthonormal vectors V = [v1 ... vk]  *
!*     2. Repeat until convergence:                                 *
!*         W = (sigma*I - B)^(-1) * V   (k PARDISO solves)          *
!*         QR decompose W = Q * R to orthonormalise                 *
!*         Compute Rayleigh quotient H = V^T (sigma*I - B)^(-1) V  *
!*         Diagonalise H to get Ritz values                         *
!*         V <- Q                                                   *
!*     3. Recover B-eigenvalues from Ritz values:                   *
!*         lambda_j = sigma - 1 / mu_j                              *
!*                                                                 *
!*   This is the standard block generalisation of power iteration. *
!*   It finds multiple eigenvalues simultaneously and handles      *
!*   close pairs / complex conjugate pairs naturally.              *
!*                                                                 *
!*   Uses only PARDISO and LAPACK (dgeqrf, dorgqr, dgeev).         *
!*                                                                 *
!*   Call this routine AFTER build_ss has filled bigm_ss.          *
!*                                                                 *
!*   Arguments:                                                    *
!*     bigm       - dense steady-state matrix (n x n)              *
!*     n          - dimension (= nbig)                             *
!*     sigma      - shift parameter                                 *
!*     k_req      - number of eigenvalues to find (typ. 5-20)       *
!*     max_iter   - maximum iterations (e.g. 300)                  *
!*     tol        - convergence tolerance (e.g. 1d-8)               *
!*     write_vec  - if .TRUE., writes eigenvectors to file         *
!*                                                                 *
!*******************************************************************

MODULE mod_eigen_bigm

  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY : IEEE_IS_FINITE

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: eigen_analysis

CONTAINS

SUBROUTINE eigen_analysis(bigm, n, sigma, k_req, max_iter, tol, write_vec)

  USE pardiso_solver

  IMPLICIT NONE

  INTEGER, PARAMETER :: dp = KIND(1.0D0)

  !--- Arguments ---------------------------------------------------
  INTEGER,  INTENT(IN) :: n, k_req, max_iter
  REAL(dp), INTENT(IN) :: bigm(n,n)
  REAL(dp), INTENT(IN) :: sigma, tol
  LOGICAL,  INTENT(IN) :: write_vec

  !--- Local --------------------------------------------------------
  TYPE(pardiso_handle_t) :: handle
  INTEGER  :: iter, i, j, nnz, error_p, info
  INTEGER  :: k              ! active block size

  REAL(dp), ALLOCATABLE :: V(:,:), W(:,:)        ! n x k iterates
  REAL(dp), ALLOCATABLE :: H(:,:)                ! k x k Rayleigh matrix
  REAL(dp), ALLOCATABLE :: mu_r(:), mu_i(:)      ! Ritz values of H
  REAL(dp), ALLOCATABLE :: lam_r(:), lam_i(:)    ! B eigenvalues
  REAL(dp), ALLOCATABLE :: lam_prev(:)
  REAL(dp), ALLOCATABLE :: VL(:,:), VR(:,:)      ! eigenvectors from dgeev
  REAL(dp), ALLOCATABLE :: tau(:), work(:)
  INTEGER               :: lwork
  REAL(dp) :: conv, rnorm

  !--- CSR storage --------------------------------------------------
  INTEGER,  ALLOCATABLE :: ia(:), ja(:)
  REAL(dp), ALLOCATABLE :: aa(:)

  !=================================================================
  ! 0. Header
  !=================================================================
  k = MIN(k_req, n/2)
  k = MAX(k, 1)

  WRITE(*,'(A)') ' ======================================================'
  WRITE(*,'(A)') '  Subspace iteration for eigenvalues of bigm_ss'
  WRITE(*,'(A)') ' ======================================================'
  WRITE(*,'(A,I8)')     '   matrix size    n       = ', n
  WRITE(*,'(A,ES12.4)') '   shift          sigma   = ', sigma
  WRITE(*,'(A,I4)')     '   block size     k       = ', k
  WRITE(*,'(A,I4)')     '   max iterations          = ', max_iter
  WRITE(*,'(A,ES10.2)') '   tolerance              = ', tol
  WRITE(*,*)

  !=================================================================
  ! 1. Build M = sigma*I - B in CSR
  !=================================================================
  CALL dense_to_csr_shift(bigm, n, sigma, ia, ja, aa, nnz)

  !=================================================================
  ! 2. Factorise M once
  !=================================================================
  CALL pardiso_solve_init(handle, n, aa, ia, ja, &
                          a_constant=.TRUE., mtype=11, msglvl=0, error=error_p)
  IF (error_p /= 0) THEN
    WRITE(*,'(A,I4)') ' PARDISO init failed, error = ', error_p
    IF (error_p == -4) WRITE(*,'(A)') ' (zero pivot: sigma coincides with an eigenvalue)'
    STOP
  END IF

  !=================================================================
  ! 3. Set up iterates
  !=================================================================
  ALLOCATE(V(n,k), W(n,k), H(k,k))
  ALLOCATE(mu_r(k), mu_i(k), lam_r(k), lam_i(k), lam_prev(k))
  ALLOCATE(VL(1,k), VR(k,k), tau(k))

  ! Deterministic pseudo-random initial vectors
  DO j = 1, k
    DO i = 1, n
      V(i,j) = SIN(DBLE(i)*0.17D0 + DBLE(j)*2.7D0) &
             + 0.3D0*COS(DBLE(i)*0.041D0 + DBLE(j)*1.3D0) &
             + 0.1D0*DBLE(j)
    END DO
  END DO

  ! Orthonormalise via QR
  lwork = 4*k*n
  ALLOCATE(work(lwork))
  CALL dgeqrf(n, k, V, n, tau, work, lwork, info)
  CALL dorgqr(n, k, k, V, n, tau, work, lwork, info)
  DEALLOCATE(work)

  lam_prev = 0.0D0

  !=================================================================
  ! 4. Subspace iteration loop
  !=================================================================
  WRITE(*,'(A)') ' Iteration  |  max conv  |  leading Ritz values (B eigenvalues)'
  WRITE(*,'(A)') ' ---------  +  --------  +  -------------------------------------'

  iter_loop: DO iter = 1, max_iter

    !-- W <- (sigma*I - B)^(-1) * V   : k PARDISO solves
    DO j = 1, k
      CALL pardiso_solve_step(handle, aa, ia, ja, V(:,j), W(:,j), &
                              a_changed=.FALSE., error=error_p)
      IF (error_p /= 0) THEN
        WRITE(*,*) ' PARDISO solve failed at column j=', j
        EXIT iter_loop
      END IF
    END DO

    !-- H <- V^T * W  (Rayleigh quotient of M^-1 in the current subspace)
    H = MATMUL(TRANSPOSE(V), W)

    !-- Diagonalise H to get Ritz values mu = mu_r + i*mu_i
    lwork = 8*k
    ALLOCATE(work(lwork))
    CALL dgeev('N', 'V', k, H, k, mu_r, mu_i, VL, 1, VR, k, work, lwork, info)
    DEALLOCATE(work)
    IF (info /= 0) THEN
      WRITE(*,*) ' dgeev failed, info=', info
      EXIT iter_loop
    END IF

    !-- Convert mu (of M^-1) to eigenvalue of B: lambda = sigma - 1/mu
    !   Complex: lambda = sigma - conj(mu) / |mu|^2
    DO j = 1, k
      IF (mu_i(j) == 0.0D0) THEN
        IF (ABS(mu_r(j)) < 1.0D-30) THEN
          lam_r(j) = 0.0D0;  lam_i(j) = 0.0D0   ! undefined
        ELSE
          lam_r(j) = sigma - 1.0D0 / mu_r(j)
          lam_i(j) = 0.0D0
        END IF
      ELSE
        rnorm = mu_r(j)**2 + mu_i(j)**2
        lam_r(j) = sigma - mu_r(j) / rnorm
        lam_i(j) = mu_i(j) / rnorm
      END IF
    END DO

    !-- Sort by distance to sigma (closest first)
    CALL sort_by_distance(lam_r, lam_i, k, sigma)

    !-- Convergence: max change in real parts of the first min(k,5) eigenvalues
    conv = 0.0D0
    DO j = 1, MIN(k,5)
      conv = MAX(conv, ABS(lam_r(j) - lam_prev(j)))
    END DO
    lam_prev = lam_r

    !-- Display progress
    IF (MOD(iter,5) == 0 .OR. iter <= 3) THEN
      WRITE(*,'(I6,6X,ES10.2,3X,4(ES13.5,:,", "))') &
          iter, conv, (lam_r(j), j=1,MIN(4,k))
      IF (k > 4) WRITE(*,'(27X,4(ES13.5,:,", "))') (lam_r(j), j=5,MIN(8,k))
    END IF

    !-- Check convergence
    IF (iter > 10 .AND. conv < tol) THEN
      WRITE(*,'(A,I0,A)') ' Converged in ', iter, ' iterations.'
      EXIT iter_loop
    END IF

    !-- Orthonormalise W and copy to V for next iteration
    lwork = 4*k*n
    ALLOCATE(work(lwork))
    CALL dgeqrf(n, k, W, n, tau, work, lwork, info)
    CALL dorgqr(n, k, k, W, n, tau, work, lwork, info)
    DEALLOCATE(work)
    V = W

  END DO iter_loop

  !=================================================================
  ! 5. Final report
  !=================================================================
  WRITE(*,*)
  WRITE(*,'(A)') ' ================  EIGENVALUE RESULTS  ================'
  WRITE(*,'(A)') '   #      Re(lambda)         Im(lambda)       |lambda-sigma|    note'
  WRITE(*,'(A)') ' ----   --------------    --------------    --------------   ------'
  DO j = 1, k
    rnorm = SQRT((lam_r(j) - sigma)**2 + lam_i(j)**2)
    IF (lam_r(j) > 1.0D-10) THEN
      WRITE(*,'(I4,3X,ES14.6,4X,ES14.6,4X,ES14.6,3X,A)') &
           j, lam_r(j), lam_i(j), rnorm, 'UNSTABLE (Re>0)'
    ELSE IF (ABS(lam_r(j)) <= 1.0D-10) THEN
      WRITE(*,'(I4,3X,ES14.6,4X,ES14.6,4X,ES14.6,3X,A)') &
           j, lam_r(j), lam_i(j), rnorm, 'null (Maxwellian)'
    ELSE
      WRITE(*,'(I4,3X,ES14.6,4X,ES14.6,4X,ES14.6,3X,A)') &
           j, lam_r(j), lam_i(j), rnorm, 'stable'
    END IF
  END DO
  WRITE(*,'(A)') ' ======================================================='

  IF (ANY(lam_r > 1.0D-10)) THEN
    WRITE(*,*)
    WRITE(*,'(A)') ' *** Positive real parts found - unstable modes ***'
    DO j = 1, k
      IF (lam_r(j) > 1.0D-10) THEN
        WRITE(*,'(A,I2,A,F12.2,A)') '   mode #', j, &
            ': e-folding time = ', 1.0D0/lam_r(j), ' s'
      END IF
    END DO
  END IF

  !=================================================================
  ! 6. Write eigenvectors to file
  !=================================================================
  IF (write_vec) THEN
    OPEN(77, file='eigenvectors.dat', status='unknown')
    WRITE(77,'(A)')         '# Eigenvectors of bigm_ss (subspace iteration)'
    WRITE(77,'(A,ES14.6)')  '# sigma = ', sigma
    WRITE(77,'(A,I4)')      '# k     = ', k
    WRITE(77,'(A,I8)')      '# n     = ', n
    WRITE(77,'(A)')         '# Columns: row_index, then v_j(row) for j=1..k'

    ! Eigenvectors in original space: W = V * VR
    W = MATMUL(V, VR)

    ! Eigenvalue header lines
    WRITE(77,'(A)', ADVANCE='NO') '# lambda_r:'
    DO j = 1, k
      WRITE(77,'(ES16.8)',ADVANCE='NO') lam_r(j)
    END DO
    WRITE(77,*)
    WRITE(77,'(A)', ADVANCE='NO') '# lambda_i:'
    DO j = 1, k
      WRITE(77,'(ES16.8)',ADVANCE='NO') lam_i(j)
    END DO
    WRITE(77,*)

    ! Data: one line per row, giving row index and then eigenvector components
    DO i = 1, n
      WRITE(77,'(I8)', ADVANCE='NO') i
      DO j = 1, k
        WRITE(77,'(ES16.8)', ADVANCE='NO') W(i,j)
      END DO
      WRITE(77,*)
    END DO
    CLOSE(77)
    WRITE(*,'(A,I4,A)') ' Wrote ', k, ' eigenvectors to eigenvectors.dat'
  END IF

  !=================================================================
  ! 7. Cleanup
  !=================================================================
  CALL pardiso_solve_finalize(handle, ia, ja, error_p)
  DEALLOCATE(V, W, H, mu_r, mu_i, lam_r, lam_i, lam_prev, VL, VR, tau)
  DEALLOCATE(ia, ja, aa)

CONTAINS

  !----------------------------------------------------------------
  ! Build M = sigma*I - B in CSR (1-based)
  !----------------------------------------------------------------
  SUBROUTINE dense_to_csr_shift(A, nmat, sig, ia_out, ja_out, aa_out, nnz_out)
    REAL(dp), INTENT(IN)  :: A(nmat,nmat), sig
    INTEGER,  INTENT(IN)  :: nmat
    INTEGER,  ALLOCATABLE, INTENT(OUT) :: ia_out(:), ja_out(:)
    REAL(dp), ALLOCATABLE, INTENT(OUT) :: aa_out(:)
    INTEGER,  INTENT(OUT) :: nnz_out

    INTEGER  :: ii, jj, pt
    REAL(dp) :: val

    nnz_out = 0
    DO ii = 1, nmat
      DO jj = 1, nmat
        val = -A(ii,jj)
        IF (ii == jj) val = val + sig
        IF (val /= 0.0D0) nnz_out = nnz_out + 1
      END DO
    END DO

    ALLOCATE(ia_out(nmat+1), ja_out(nnz_out), aa_out(nnz_out))

    pt = 1
    ia_out(1) = 1
    DO ii = 1, nmat
      DO jj = 1, nmat
        val = -A(ii,jj)
        IF (ii == jj) val = val + sig
        IF (val /= 0.0D0) THEN
          ja_out(pt) = jj
          aa_out(pt) = val
          pt = pt + 1
        END IF
      END DO
      ia_out(ii+1) = pt
    END DO
  END SUBROUTINE dense_to_csr_shift

  !----------------------------------------------------------------
  ! Sort eigenvalues by distance to sigma (closest first).
  ! Simple insertion/selection sort; k is small (<30 in practice).
  !----------------------------------------------------------------
  SUBROUTINE sort_by_distance(lr, li, kk, sig)
    INTEGER,  INTENT(IN)    :: kk
    REAL(dp), INTENT(INOUT) :: lr(kk), li(kk)
    REAL(dp), INTENT(IN)    :: sig

    INTEGER  :: ii, jj
    REAL(dp) :: d_i, d_j, tr, ti

    DO ii = 1, kk-1
      DO jj = ii+1, kk
        d_i = (lr(ii)-sig)**2 + li(ii)**2
        d_j = (lr(jj)-sig)**2 + li(jj)**2
        IF (d_j < d_i) THEN
          tr = lr(ii); lr(ii) = lr(jj); lr(jj) = tr
          ti = li(ii); li(ii) = li(jj); li(jj) = ti
        END IF
      END DO
    END DO
  END SUBROUTINE sort_by_distance

END SUBROUTINE eigen_analysis

END MODULE mod_eigen_bigm
