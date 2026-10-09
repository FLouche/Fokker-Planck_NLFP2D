module mod_apply_operator
    
    contains
    
      !================================================================
  ! SUBROUTINE apply_operator
  !
  ! Computes  Lf = L(a00,a10,...) * f  without forming the matrix.
  !
  ! For each row (i,j), fd_stencil_2d returns the stencil entries
  ! {col_idx(k), coeff(k)}.  The matrix-vector product is:
  !
  !   Lf(row) = sum_k  coeff(k) * f(col_idx(k))
  !
  ! taum is set by the caller before entering this routine.
  ! The stencil call uses dvpar from shared_grid (uniform v∥).
  ! Replace dvpar -> vpar if using the non-uniform vpar version.
  !================================================================
  SUBROUTINE apply_operator(a20, a02, a11, a10, a01, a00, f_in, Lf_out)
  
  USE mod_fd_stencil_2d     ! provides fd_stencil_2d
  USE shared_grid
  
    IMPLICIT NONE

  INTEGER, PARAMETER :: dp = KIND(1.0D0)

    REAL(dp), INTENT(IN)  :: a20(nperp,npar), a02(nperp,npar)
    REAL(dp), INTENT(IN)  :: a11(nperp,npar), a10(nperp,npar)
    REAL(dp), INTENT(IN)  :: a01(nperp,npar), a00(nperp,npar)
    REAL(dp), INTENT(IN)  :: f_in(nbig)
    REAL(dp), INTENT(OUT) :: Lf_out(nbig)

    INTEGER  :: col_idx_loc(49)
    REAL(dp) :: coeff_loc(49)
    INTEGER  :: n_ent, row_loc, k_loc
    REAL(dp) :: rhs_loc
    INTEGER  :: i_loc, j_loc

    Lf_out = 0.0_dp

    ! Rows are independent: each writes only Lf_out(row_loc).  Called ~8 times
    ! per step by the power / density-term diagnostics, which made it the
    ! largest serial cost once the linear solve was cheap.
    !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) &
    !$OMP   PRIVATE(i_loc, j_loc, row_loc, k_loc, n_ent, rhs_loc, col_idx_loc, coeff_loc)
    DO i_loc = 1, nperp
      DO j_loc = 1, npar

        row_loc = (i_loc-1)*npar + j_loc

        CALL fd_stencil_2d(i_loc, j_loc, nperp, npar, vperp, dvpar, &
                           a00(i_loc,j_loc), a10(i_loc,j_loc), a01(i_loc,j_loc), &
                           a20(i_loc,j_loc), a11(i_loc,j_loc), a02(i_loc,j_loc), &
                           col_idx_loc, coeff_loc, n_ent, rhs_loc)

        DO k_loc = 1, n_ent
          Lf_out(row_loc) = Lf_out(row_loc) &
                          + coeff_loc(k_loc) * f_in(col_idx_loc(k_loc))
        END DO

      END DO
    END DO
    !$OMP END PARALLEL DO

  END SUBROUTINE apply_operator
  
    end module mod_apply_operator
    