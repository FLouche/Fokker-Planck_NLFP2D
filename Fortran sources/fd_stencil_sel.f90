! TEMPORARY: 3-point stencil comparison study — remove after 2026-09-30
!
! Runtime dispatcher between the production 7-point Fornberg stencil
! (mod_fd_stencil_2d, UNCHANGED) and the recovered 3-point stencil
! (mod_fd_stencil_2d_3), selected by the namelist flag `stencil`:
!
!     stencil = 7   -> fd_stencil_2d      (default; production behaviour)
!     stencil = 3   -> fd_stencil_2d_3    (study arm only)
!
! Every caller that previously called fd_stencil_2d directly now calls
! fd_stencil_sel with an identical argument list, so removing this study
! is a mechanical rename back to fd_stencil_2d. See STENCIL_STUDY_REMOVAL.md.
!
MODULE mod_fd_stencil_sel

  USE shared_grid,        ONLY: stencil
  USE mod_fd_stencil_2d,  ONLY: fd_stencil_2d
  USE mod_fd_stencil_2d_3, ONLY: fd_stencil_2d_3

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: fd_stencil_sel, stencil_nnz_per_row

  INTEGER, PARAMETER :: dp = KIND(1.0D0)

CONTAINS

SUBROUTINE fd_stencil_sel(i, j, nperp, npar, vperp, dvpar, &
                          A_ij, B_ij, C_ij, D_ij, E_ij, F_ij, &
                          col_idx, coeff, n_entries, rhs)
  INTEGER,  INTENT(IN)  :: i, j, nperp, npar
  REAL(dp), INTENT(IN)  :: vperp(nperp), dvpar
  REAL(dp), INTENT(IN)  :: A_ij, B_ij, C_ij, D_ij, E_ij, F_ij
  INTEGER,  INTENT(OUT) :: col_idx(49)
  REAL(dp), INTENT(OUT) :: coeff(49)
  INTEGER,  INTENT(OUT) :: n_entries
  REAL(dp), INTENT(OUT) :: rhs

  IF (stencil == 3) THEN
    CALL fd_stencil_2d_3(i, j, nperp, npar, vperp, dvpar, &
                         A_ij, B_ij, C_ij, D_ij, E_ij, F_ij, &
                         col_idx, coeff, n_entries, rhs)
  ELSE
    CALL fd_stencil_2d(i, j, nperp, npar, vperp, dvpar, &
                       A_ij, B_ij, C_ij, D_ij, E_ij, F_ij, &
                       col_idx, coeff, n_entries, rhs)
  END IF

END SUBROUTINE fd_stencil_sel

!-------------------------------------------------------------------
! Upper bound on non-zeros per row for the active stencil.
! Used to size ja/aa so the 3-point arm's memory footprint is measured
! honestly rather than being padded to the 7-point allocation.
!-------------------------------------------------------------------
PURE FUNCTION stencil_nnz_per_row() RESULT(n)
  INTEGER :: n
  IF (stencil == 3) THEN
    n = 9
  ELSE
    n = 49
  END IF
END FUNCTION stencil_nnz_per_row

END MODULE mod_fd_stencil_sel
