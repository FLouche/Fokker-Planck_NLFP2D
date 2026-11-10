! ---------------------------------------------------------------------------
! DEAD CODE, kept for reference only. Not in the Visual Studio project.
! Moved here on 2026-10-27 from 'pardiso_solver (2).f90': subroutine pardiso_solve_steady, a contained routine
! of MODULE pardiso_solver. Its only caller was sparse_solve_ss (also in Dead/).
! ---------------------------------------------------------------------------

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

