!==============================================================================!
! MODULE pardiso_solver
!
! A structured interface to Intel MKL PARDISO that exploits phased calling:
!
!   Phase 11  – Reordering & symbolic factorisation  (sparsity pattern fixed)
!   Phase 22  – Numerical factorisation               (values of A)
!   Phase 33  – Back-substitution / solve             (rhs b)
!   Phase -1  – Release internal PARDISO memory
!
! Three public solve entry points:
!
!   pardiso_solve_steady   – one-shot: phases 11 + 22 + 33 then release
!   pardiso_solve_init     – time-loop setup: phase 11 (+ optionally 22)
!   pardiso_solve_step     – one time step: phase 22 (if A changed) + 33
!   pardiso_solve_finalize – release PARDISO memory after the time loop
!
! LU reuse (pardiso_solve_init argument lu_reuse_L > 0): when A has changed,
! pardiso_solve_step first tries phase 23 with iparm(4) = 10*L + 1, i.e. CGS
! preconditioned with the LU factors of the LAST factorisation, to a relative
! residual of 10^-L.  The result is accepted only if its true residual, checked
! here, is within max(10^-L, 10 x that of the latest direct solve); otherwise A
! is refactorised (phase 22 + 33).  The time-step operator changes slowly, so
! one factorisation typically serves hundreds of steps with ~3 CGS iterations
! each: 2.7-3.4x faster full runs, results within 10^-6 of refactorising every
! step (Simulations/REPORTS_ALL/Parallelisation).  Stats: pardiso_solve_finalize.
!
! Conventions (CSR, 1-based, upper-triangular for symmetric problems):
!   ia(n+1)  – row pointers  (1-based)
!   ja(nnz)  – column indices (1-based)
!   a(nnz)   – non-zero values
!   b(n)     – right-hand side
!   x(n)     – solution vector
!
! Matrix type (mtype):
!    1 = real structurally symmetric
!    2 = real symmetric positive definite
!   -2 = real symmetric indefinite
!   11 = real non-symmetric  (general)
!  (see PARDISO documentation for the full list)
!==============================================================================!
MODULE pardiso_solver

  USE iso_fortran_env, ONLY: dp => real64
  IMPLICIT NONE
  PRIVATE

  !-- Public API ---------------------------------------------------------------
  PUBLIC :: pardiso_handle_t
  PUBLIC :: pardiso_solve_steady
  PUBLIC :: pardiso_solve_init
  PUBLIC :: pardiso_solve_step
  PUBLIC :: pardiso_solve_finalize

  !-- Opaque handle that carries all PARDISO internal state --------------------
  TYPE :: pardiso_handle_t
    PRIVATE
    INTEGER(8)         :: pt(64)     = 0   ! PARDISO internal pointer array
    INTEGER            :: iparm(64)  = 0   ! integer parameters
    REAL(dp)           :: dparm(64)  = 0.0_dp ! double parameters (DSS path)
    INTEGER            :: mtype      = 11  ! matrix type (default: non-symmetric)
    INTEGER            :: n          = 0   ! matrix dimension
    INTEGER            :: nrhs       = 1   ! number of right-hand sides
    INTEGER            :: maxfct     = 1
    INTEGER            :: mnum       = 1
    INTEGER            :: msglvl     = 0   ! 0=silent, 1=verbose
    LOGICAL            :: symbolic_done  = .FALSE.
    LOGICAL            :: numeric_done   = .FALSE.
    ! LU reuse (0 = off) and its statistics
    INTEGER            :: lu_reuse_L     = 0
    INTEGER            :: n_steps = 0, n_fact = 0, n_cgs_ok = 0, n_cgs_fail = 0
    INTEGER(8)         :: n_cgs_iter = 0
    REAL(dp)           :: res_reuse_max = 0.0_dp, res_direct_max = 0.0_dp
    REAL(dp)           :: res_direct_last = 0.0_dp   ! residual of the latest direct solve
  END TYPE pardiso_handle_t

  !-- PARDISO external interface (MKL signature) --------------------------------
  INTERFACE
    SUBROUTINE pardiso(pt, maxfct, mnum, mtype, phase, n, &
                       a, ia, ja, perm, nrhs, iparm, msglvl, &
                       b, x, error)
      USE iso_fortran_env, ONLY: dp => real64
      IMPLICIT NONE
      INTEGER(8), INTENT(INOUT) :: pt(64)
      INTEGER,    INTENT(IN)    :: maxfct, mnum, mtype, phase, n, nrhs, msglvl
      REAL(dp),   INTENT(INOUT) :: a(*)      ! PARDISO may scale/modify a in place
      INTEGER,    INTENT(INOUT) :: ia(*), ja(*) ! may be reordered internally
      INTEGER,    INTENT(INOUT) :: perm(*)
      INTEGER,    INTENT(INOUT) :: iparm(64)
      REAL(dp),   INTENT(INOUT) :: b(*), x(*)
      INTEGER,    INTENT(OUT)   :: error
    END SUBROUTINE pardiso
  END INTERFACE

CONTAINS

  !=============================================================================
  ! pardiso_handle_init  (private helper)
  !   Initialise iparm defaults and set matrix type.
  !=============================================================================
  SUBROUTINE handle_init(h, n, mtype, msglvl)
    TYPE(pardiso_handle_t), INTENT(INOUT) :: h
    INTEGER, INTENT(IN) :: n, mtype, msglvl

    INTEGER :: error!, dummy_perm(1)
  !  REAL(dp) :: dummy(1)

    h%n      = n
    h%mtype  = mtype
    h%msglvl = msglvl
    h%pt     = 0

    ! Let PARDISO set default iparm values (phase = -1 with pt=0 triggers init)
    h%iparm     = 0
    h%iparm(1)  = 1   ! Do not use default values; we set them ourselves below.

    ! Call pardisoinit to get defaults for this matrix type
    ! (Some MKL versions expose pardisoinit; fall back to manual defaults.)
    CALL pardisoinit_safe(h, error)
    IF (error /= 0) THEN
      WRITE(*,'(A,I0)') 'pardiso_solver: pardisoinit error = ', error
      STOP
    END IF

    ! Common overrides
    h%iparm(1)  = 1      ! non-default settings
    h%iparm(2)  = 3      ! METIS reordering
    h%iparm(8)  = 10     ! max iterative refinement steps
    h%iparm(10) = 13     ! pivot perturbation (10^{-13})
    h%iparm(11) = 1      ! scaling (useful for non-symmetric)
    h%iparm(13) = 1      ! improved accuracy for non-symmetric
    h%iparm(27) = 1      ! check the matrix
    h%iparm(35) = 0      ! 1-based indexing
    h%nrhs      = 1

  END SUBROUTINE handle_init

  !-- Thin wrapper around pardisoinit; if unavailable, sets iparm manually -----
  SUBROUTINE pardisoinit_safe(h, error)
    TYPE(pardiso_handle_t), INTENT(INOUT) :: h
    INTEGER, INTENT(OUT) :: error
    error = 0
    ! pardisoinit(pt, mtype, iparm)
    ! If your MKL version has it, uncomment the next line and remove the manual
    ! iparm block below.
    ! CALL pardisoinit(h%pt, h%mtype, h%iparm)

    ! Manual iparm initialisation (safe fallback)
    h%iparm     = 0
    h%iparm(1)  = 1
    h%iparm(2)  = 3      ! METIS reordering
    h%iparm(4)  = 0      ! direct algorithm
    h%iparm(5)  = 0      ! user permutation ignored
    h%iparm(6)  = 0      ! write solution to x
    h%iparm(8)  = 10     ! iterative refinement steps
    h%iparm(10) = 13     ! pivot perturbation exponent
    h%iparm(11) = 1      ! enable scaling
    h%iparm(13) = 1      ! improved accuracy
    h%iparm(18) = -1     ! report number of non-zeros in factors
    h%iparm(19) = -1     ! report Mflops for factorisation
    h%iparm(27) = 1      ! matrix checker
    h%iparm(35) = 0      ! 1-based CSR indexing
  END SUBROUTINE pardisoinit_safe

  !-- ||b - A x||_inf / ||b||_inf for a 1-based CSR matrix ---------------------
  FUNCTION rel_residual(n, a, ia, ja, b, x) RESULT(res)
    INTEGER,  INTENT(IN) :: n, ia(*), ja(*)
    REAL(dp), INTENT(IN) :: a(*), b(*), x(*)
    REAL(dp) :: res, r, rmax, bmax
    INTEGER  :: i, p
    rmax = 0.0_dp; bmax = 0.0_dp
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(i, p, r) REDUCTION(MAX:rmax, bmax)
    DO i = 1, n
      r = b(i)
      DO p = ia(i), ia(i+1) - 1
        r = r - a(p) * x(ja(p))
      END DO
      rmax = MAX(rmax, ABS(r))
      bmax = MAX(bmax, ABS(b(i)))
    END DO
    !$OMP END PARALLEL DO
    res = rmax / MAX(bmax, TINY(1.0_dp))
  END FUNCTION rel_residual

  !=============================================================================
  ! PHASE 11 – symbolic factorisation (reordering)
  !   Only depends on ia, ja (sparsity pattern).  Called once.
  !=============================================================================
  SUBROUTINE phase_symbolic(h, a, ia, ja, error)
    TYPE(pardiso_handle_t), INTENT(INOUT) :: h
    REAL(dp), INTENT(INOUT) :: a(*)
    INTEGER,  INTENT(INOUT) :: ia(*), ja(*)
    INTEGER,  INTENT(OUT) :: error

    INTEGER  :: perm(h%n)
    REAL(dp) :: dummy(1)

    CALL pardiso(h%pt, h%maxfct, h%mnum, h%mtype, &
                 11, &                   ! phase 11: analysis + reordering
                 h%n, a, ia, ja, perm, h%nrhs, h%iparm, h%msglvl, &
                 dummy, dummy, error)
    IF (error /= 0) THEN
      WRITE(*,'(A,I0)') 'pardiso_solver: phase 11 error = ', error
      RETURN
    END IF
    h%symbolic_done = .TRUE.
  END SUBROUTINE phase_symbolic

  !=============================================================================
  ! PHASE 22 – numerical factorisation
  !   Depends on the values in a.  Re-called whenever A changes.
  !=============================================================================
  SUBROUTINE phase_numeric(h, a, ia, ja, error)
    TYPE(pardiso_handle_t), INTENT(INOUT) :: h
    REAL(dp), INTENT(INOUT) :: a(*)
    INTEGER,  INTENT(INOUT) :: ia(*), ja(*)
    INTEGER,  INTENT(OUT) :: error

    INTEGER  :: perm(h%n)
    REAL(dp) :: dummy(1)

    IF (.NOT. h%symbolic_done) THEN
      WRITE(*,'(A)') 'pardiso_solver: phase 22 called before phase 11'
      error = -100
      RETURN
    END IF

    CALL pardiso(h%pt, h%maxfct, h%mnum, h%mtype, &
                 22, &                   ! phase 22: numerical factorisation
                 h%n, a, ia, ja, perm, h%nrhs, h%iparm, h%msglvl, &
                 dummy, dummy, error)
    IF (error /= 0) THEN
      WRITE(*,'(A,I0)') 'pardiso_solver: phase 22 error = ', error
      RETURN
    END IF
    h%numeric_done = .TRUE.
  END SUBROUTINE phase_numeric

  !=============================================================================
  ! PHASE 33 – back-substitution / solve
  !   Depends only on the factored matrix (stored internally by PARDISO) and b.
  !=============================================================================
  SUBROUTINE phase_solve(h, a, ia, ja, b, x, error)
    TYPE(pardiso_handle_t), INTENT(INOUT) :: h
    REAL(dp), INTENT(INOUT) :: a(*)
    INTEGER,  INTENT(INOUT) :: ia(*), ja(*)
    REAL(dp), INTENT(INOUT) :: b(*)   ! INOUT: PARDISO may use b as workspace
    REAL(dp), INTENT(INOUT) :: x(*)
    INTEGER,  INTENT(OUT)   :: error

    INTEGER :: perm(h%n)
    INTEGER :: nrefine_requested

    IF (.NOT. h%numeric_done) THEN
      WRITE(*,'(A)') 'pardiso_solver: phase 33 called before phase 22'
      error = -101
      RETURN
    END IF

    ! iparm(8) is dual-use: INPUT = max refinement steps requested;
    ! OUTPUT = actual steps taken (PARDISO overwrites it).
    ! Save and restore so every call uses the same requested value,
    ! otherwise later time steps silently get 0 refinement steps.
    nrefine_requested = h%iparm(8)

    CALL pardiso(h%pt, h%maxfct, h%mnum, h%mtype, &
                 33, &                   ! phase 33: solve + iterative refinement
                 h%n, a, ia, ja, perm, h%nrhs, h%iparm, h%msglvl, &
                 b, x, error)

    h%iparm(8) = nrefine_requested   ! restore for next call

    IF (error /= 0) THEN
      WRITE(*,'(A,I0)') 'pardiso_solver: phase 33 error = ', error
    END IF
  END SUBROUTINE phase_solve

  !=============================================================================
  ! PUBLIC ROUTINES
  !=============================================================================

  !-----------------------------------------------------------------------------
  ! pardiso_solve_steady
  !
  !   One-shot solve for the steady-state case:
  !     Phase 11 + 22 + 33 + release
  !   The handle h is initialised inside this routine.
  !-----------------------------------------------------------------------------
  SUBROUTINE pardiso_solve_steady(n, a, ia, ja, b, x, &
                                  mtype, msglvl, error)
    INTEGER,  INTENT(IN)    :: n
    REAL(dp), INTENT(INOUT) :: a(*)         ! PARDISO may modify a (scaling)
    INTEGER,  INTENT(INOUT) :: ia(n+1)      ! CSR row pointers
    INTEGER,  INTENT(INOUT) :: ja(*)        ! CSR column indices
    REAL(dp), INTENT(INOUT) :: b(n)         ! rhs (PARDISO may use as workspace)
    REAL(dp), INTENT(OUT)   :: x(n)         ! solution
    INTEGER,  INTENT(IN), OPTIONAL :: mtype    ! matrix type (default 11)
    INTEGER,  INTENT(IN), OPTIONAL :: msglvl   ! verbosity (default 0)
    INTEGER,  INTENT(OUT), OPTIONAL :: error

    TYPE(pardiso_handle_t) :: h
    INTEGER :: mt, ml, ierr, perm(n)
    REAL(dp) :: dummy(1)

    mt = 11 ; IF (PRESENT(mtype))  mt = mtype
    ml = 0  ; IF (PRESENT(msglvl)) ml = msglvl

    CALL handle_init(h, n, mt, ml)

    !-- Phase 11: symbolic ---------------------------------------------------
    CALL phase_symbolic(h, a, ia, ja, ierr)
    IF (ierr /= 0) GOTO 99

    !-- Phase 22: numeric ----------------------------------------------------
    CALL phase_numeric(h, a, ia, ja, ierr)
    IF (ierr /= 0) GOTO 99

    !-- Phase 33: solve -------------------------------------------------------
    CALL phase_solve(h, a, ia, ja, b, x, ierr)

99  CONTINUE
    !-- Release memory --------------------------------------------------------
    CALL pardiso(h%pt, h%maxfct, h%mnum, h%mtype, &
                 -1, h%n, dummy, ia, ja, perm, h%nrhs, h%iparm, h%msglvl, &
                 dummy, dummy, ierr)

    IF (PRESENT(error)) error = ierr
  END SUBROUTINE pardiso_solve_steady

  !-----------------------------------------------------------------------------
  ! pardiso_solve_init
  !
  !   Initialise a handle for a time-dependent problem.
  !   Always performs phase 11 (symbolic).
  !   If a_constant=.TRUE., also performs phase 22 (numeric) right away,
  !   so that pardiso_solve_step only needs to do phase 33.
  !   If a_constant=.FALSE., phase 22 is deferred to pardiso_solve_step.
  !-----------------------------------------------------------------------------
  SUBROUTINE pardiso_solve_init(h, n, a, ia, ja, a_constant, &
                                mtype, msglvl, error, lu_reuse_L)
    TYPE(pardiso_handle_t), INTENT(OUT) :: h
    INTEGER,  INTENT(IN)    :: n
    REAL(dp), INTENT(INOUT) :: a(*)
    INTEGER,  INTENT(INOUT) :: ia(n+1), ja(*)
    LOGICAL,  INTENT(IN)  :: a_constant   ! .TRUE. => A never changes
    INTEGER,  INTENT(IN),  OPTIONAL :: mtype
    INTEGER,  INTENT(IN),  OPTIONAL :: msglvl
    INTEGER,  INTENT(OUT), OPTIONAL :: error
    INTEGER,  INTENT(IN),  OPTIONAL :: lu_reuse_L  ! >0: reuse LU, CGS tol 10^-L

    INTEGER :: mt, ml, ierr

    mt = 11 ; IF (PRESENT(mtype))  mt = mtype
    ml = 0  ; IF (PRESENT(msglvl)) ml = msglvl

    CALL handle_init(h, n, mt, ml)
    IF (PRESENT(lu_reuse_L)) h%lu_reuse_L = MAX(0, lu_reuse_L)

    !-- Phase 11 (sparsity pattern, always once) --------------------------------
    CALL phase_symbolic(h, a, ia, ja, ierr)
    IF (ierr /= 0) GOTO 99

    !-- Phase 22 (values) – only if A is constant --------------------------------
    IF (a_constant) THEN
      CALL phase_numeric(h, a, ia, ja, ierr)
    END IF

99  IF (PRESENT(error)) error = ierr
  END SUBROUTINE pardiso_solve_init

  !-----------------------------------------------------------------------------
  ! pardiso_solve_step
  !
  !   Solve one time step:
  !   - If a_changed=.TRUE., re-runs phase 22 before phase 33.
  !   - If a_changed=.FALSE., skips phase 22 (reuses existing factorisation).
  !
  !   Call after pardiso_solve_init.
  !-----------------------------------------------------------------------------
  SUBROUTINE pardiso_solve_step(h, a, ia, ja, b, x, a_changed, error)
    TYPE(pardiso_handle_t), INTENT(INOUT) :: h
    REAL(dp), INTENT(INOUT) :: a(*)
    INTEGER,  INTENT(INOUT) :: ia(*), ja(*)
    REAL(dp), INTENT(INOUT) :: b(*)         ! rhs (PARDISO may use as workspace)
    REAL(dp), INTENT(INOUT) :: x(*)
    LOGICAL,  INTENT(IN)    :: a_changed    ! .TRUE. => redo phase 22
    INTEGER,  INTENT(OUT), OPTIONAL :: error

    INTEGER  :: ierr, nref, perm(h%n)
    REAL(dp) :: res

    h%n_steps = h%n_steps + 1

    !-- LU reuse: CGS preconditioned with the last LU factors (phase 23) -------
    !   PARDISO's own verdict (iparm(20) > 0) is not enough: it reports success
    !   even for tolerances it cannot reach (10^-60 "converges" in 13
    !   iterations).  So the step is accepted only if the TRUE relative residual
    !   ||b - A x||_inf / ||b||_inf, computed here (one sparse mat-vec), is
    !   <= MAX(10^-L, 10 * residual of the latest direct solve); otherwise A is
    !   refactorised and solved directly below.  The second bound matters: on
    !   some problems the direct solve itself only reaches ~5e-10 (NBI fill-up),
    !   and a fixed 10^-10 then rejected almost every step although the CGS
    !   results matched the direct ones to 1e-13.  With iparm(6) = 0 PARDISO
    !   does not overwrite b, so the direct solve still sees the caller's rhs.
    IF (a_changed .AND. h%lu_reuse_L > 0 .AND. h%numeric_done) THEN
      nref = h%iparm(8)
      h%iparm(4) = 10*h%lu_reuse_L + 1
      CALL pardiso(h%pt, h%maxfct, h%mnum, h%mtype, 23, h%n, a, ia, ja, perm, &
                   h%nrhs, h%iparm, h%msglvl, b, x, ierr)
      h%iparm(4) = 0
      h%iparm(8) = nref                 ! iparm(8) is also an output; restore it
      IF (ierr == 0 .AND. h%iparm(20) > 0) THEN
        res = rel_residual(h%n, a, ia, ja, b, x)
        IF (res <= MAX(10.0_dp**(-h%lu_reuse_L), 10.0_dp*h%res_direct_last)) THEN
          h%n_cgs_ok   = h%n_cgs_ok + 1
          h%n_cgs_iter = h%n_cgs_iter + h%iparm(20)
          h%res_reuse_max = MAX(h%res_reuse_max, res)
          GOTO 99
        END IF
      END IF
      h%n_cgs_fail = h%n_cgs_fail + 1
    END IF

    !-- Optional phase 22: only when matrix values have changed ----------------
    IF (a_changed) THEN
      CALL phase_numeric(h, a, ia, ja, ierr)
      IF (ierr /= 0) GOTO 99
      h%n_fact = h%n_fact + 1
    END IF

    !-- Phase 33: solve with current b ----------------------------------------
    CALL phase_solve(h, a, ia, ja, b, x, ierr)
    IF (h%lu_reuse_L > 0 .AND. ierr == 0) THEN
      h%res_direct_last = rel_residual(h%n, a, ia, ja, b, x)
      h%res_direct_max  = MAX(h%res_direct_max, h%res_direct_last)
    END IF

99  IF (PRESENT(error)) error = ierr
  END SUBROUTINE pardiso_solve_step

  !-----------------------------------------------------------------------------
  ! pardiso_solve_finalize
  !
  !   Release all PARDISO internal memory.  Call once after the time loop ends.
  !-----------------------------------------------------------------------------
  SUBROUTINE pardiso_solve_finalize(h, ia, ja, error)
    TYPE(pardiso_handle_t), INTENT(INOUT) :: h
    INTEGER,  INTENT(INOUT) :: ia(*), ja(*)
    INTEGER,  INTENT(OUT), OPTIONAL :: error

    INTEGER  :: perm(1), ierr
    REAL(dp) :: dummy(1)

    CALL pardiso(h%pt, h%maxfct, h%mnum, h%mtype, &
                 -1, &                  ! phase -1: release memory
                 h%n, dummy, ia, ja, perm, h%nrhs, h%iparm, h%msglvl, &
                 dummy, dummy, ierr)

    h%symbolic_done = .FALSE.
    h%numeric_done  = .FALSE.

    IF (h%n_steps > 0) THEN
      IF (h%lu_reuse_L > 0) THEN
        WRITE(*,'(A,I0,A,I0,A,I0,A,F5.1,A,I0,A)') ' PARDISO: ', h%n_steps, &
          ' steps, ', h%n_fact, ' factorisations; LU reused on ', h%n_cgs_ok, &
          ' steps (', DBLE(h%n_cgs_iter)/MAX(1,h%n_cgs_ok), &
          ' CGS iterations on average, ', h%n_cgs_fail, ' rejected)'
        WRITE(*,'(A,ES9.2,A,ES9.2,A,ES9.2)') '          max relative residual: LU reuse ', &
          h%res_reuse_max, ', direct ', h%res_direct_max, ';  tolerance ', &
          10.0_dp**(-h%lu_reuse_L)
      ELSE
        WRITE(*,'(A,I0,A,I0,A)') ' PARDISO: ', h%n_steps, ' steps, ', &
          h%n_fact, ' factorisations (LU reuse off)'
      END IF
    END IF

    IF (PRESENT(error)) error = ierr
  END SUBROUTINE pardiso_solve_finalize

END MODULE pardiso_solver
