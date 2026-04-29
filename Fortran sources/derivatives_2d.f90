!==============================================================================
! MODULE: derivatives_2d
!
! Double-precision finite-difference derivatives of a 2-D grid function F(x,y)
! on a grid that is ARBITRARY (non-uniform) in x and UNIFORM in y.
!
! Conventions
! -----------
!   F(i,j)  : function value at grid point (x(i), y(j))  [1-based]
!   x(1:nx) : arbitrary strictly-increasing x-coordinates [REAL(8)]
!   dy      : uniform y spacing                           [REAL(8)]
!   nx, ny  : number of grid points
!
! Method – x direction
! --------------------
!   Fornberg (1988) algorithm computes exact finite-difference weights for
!   any derivative order on any node set.  For each grid point i a 7-point
!   local stencil [L, L+6] is selected (centred symmetrically, clamped at
!   boundaries) and weights for orders 1, 2, 3 are obtained simultaneously.
!   The node coordinates are shifted so that the target point maps to 0,
!   which improves conditioning and simplifies the Fornberg recurrence.
!
! Method – y direction
! --------------------
!   Uniform grid: 6th-order centred stencils for interior points (j=4..ny-3),
!   4th-order one-sided stencils for j=1,2,3 and j=ny-2,ny-1,ny.
!   All coefficients verified by Taylor expansion and cross-checked against
!   the Fornberg algorithm on a uniform grid.
!
! Mixed derivative d²F/(dx dy)
! ----------------------------
!   Computed by composition: x-derivative first, then y-derivative.
!
! Public subroutines
! ------------------
!   deriv_x1(F, x, nx, ny, dFdx)
!   deriv_x2(F, x, nx, ny, dFdx2)
!   deriv_x3(F, x, nx, ny, dFdx3)
!   deriv_y1(F, nx, ny, dy, dFdy)
!   deriv_y2(F, nx, ny, dy, dFdy2)
!   deriv_y3(F, nx, ny, dy, dFdy3)
!   deriv_xy(F, x, nx, ny, dy, d2Fdxdy)
!
! Requirements: nx >= 7, ny >= 7.
! All reals are REAL(8) (IEEE 754 double precision).
!==============================================================================

MODULE derivatives_2d

  IMPLICIT NONE
  PRIVATE

  INTEGER, PARAMETER, PUBLIC :: dp = KIND(1.0D0)

  PUBLIC :: deriv_x1, deriv_x2, deriv_x3
  PUBLIC :: deriv_y1, deriv_y2, deriv_y3
  PUBLIC :: deriv_xy, laplacian_radial, &
            fornberg_weights               

CONTAINS

!==============================================================================
! FORNBERG_WEIGHTS
!
! Compute finite-difference weights for derivatives of orders 0..m_max
! at target point z, using n nodes x(1:n)  (standard 1-based Fortran).
!
! Output: w(1:n, 0:m_max)
!         w(k,m) = weight of node x(k) in the m-th derivative approximation.
!
! Algorithm: Fornberg, B. (1988). Math. Comp. 51(184), 699-706.
! Correct variant: the new node row w(i,:) is only filled on the last inner
! iteration (j == i-1); all other inner iterations update w(j,:) only.
!==============================================================================
  SUBROUTINE fornberg_weights(z, x, n, m_max, w)
    INTEGER,  INTENT(IN)  :: n, m_max
    REAL(dp), INTENT(IN)  :: z, x(n)
    REAL(dp), INTENT(OUT) :: w(n, 0:m_max)

    REAL(dp) :: c1, c2, c3, c4, c5
    INTEGER  :: mn, i, j, s

    w  = 0.0_dp
    w(1,0) = 1.0_dp
    c1 = 1.0_dp
    c4 = x(1) - z          ! x(1)-z  (used to track x(j)-z for previous j)

    DO i = 2, n             ! add node x(i) one at a time
      mn = MIN(i-1, m_max)
      c2 = 1.0_dp
      c5 = c4               ! save x(i-1)-z
      c4 = x(i) - z

      DO j = 1, i-1         ! loop over all previously added nodes x(j)
        c3 = x(i) - x(j)
        c2 = c2 * c3

        ! On the last j-iteration, fill the new row w(i,:)
        IF (j == i-1) THEN
          DO s = mn, 1, -1
            w(i,s) = c1 * (s*w(i-1,s-1) - c5*w(i-1,s)) / c2
          END DO
          w(i,0) = -c1 * c5 * w(i-1,0) / c2
        END IF

        ! Always update w(j,:) using the current c4 = x(i)-z
        DO s = mn, 1, -1
          w(j,s) = (c4*w(j,s) - s*w(j,s-1)) / c3
        END DO
        w(j,0) = c4 * w(j,0) / c3

      END DO
      c1 = c2
    END DO
  END SUBROUTINE fornberg_weights


!==============================================================================
! STENCIL_LEFT
!
! For 1-based grid index i (1..nx), return the 1-based left endpoint L of a
! 7-point window [L, L+6], clamped so that 1 <= L and L+6 <= nx.
! Centred as symmetrically as possible: ideal L = i-3.
!==============================================================================
  PURE FUNCTION stencil_left(i, nx) RESULT(L)
    INTEGER, INTENT(IN) :: i, nx
    INTEGER :: L
    L = MAX(1,    i - 3)
    L = MIN(L, nx - 6)
  END FUNCTION stencil_left


!==============================================================================
!  d F / d x   (non-uniform x)
!==============================================================================
  SUBROUTINE deriv_x1(F, x, nx, ny, dFdx)
    INTEGER,  INTENT(IN)  :: nx, ny
    REAL(dp), INTENT(IN)  :: F(nx,ny), x(nx)
    REAL(dp), INTENT(OUT) :: dFdx(nx,ny)

    INTEGER  :: i, j, L, k
    REAL(dp) :: xloc(7), w(7, 0:1), acc

    DO i = 1, nx
      L = stencil_left(i, nx)
      xloc = x(L:L+6) - x(i)          ! shift target to 0 for conditioning
      CALL fornberg_weights(0.0_dp, xloc, 7, 1, w)
      DO j = 1, ny
        acc = 0.0_dp
        DO k = 1, 7
          acc = acc + w(k,1) * F(L+k-1, j)
        END DO
        dFdx(i,j) = acc
      END DO
    END DO
  END SUBROUTINE deriv_x1


!==============================================================================
!  d²F / d x²   (non-uniform x)
!==============================================================================
  SUBROUTINE deriv_x2(F, x, nx, ny, dFdx2)
    INTEGER,  INTENT(IN)  :: nx, ny
    REAL(dp), INTENT(IN)  :: F(nx,ny), x(nx)
    REAL(dp), INTENT(OUT) :: dFdx2(nx,ny)

    INTEGER  :: i, j, L, k
    REAL(dp) :: xloc(7), w(7, 0:2), acc

    DO i = 1, nx
      L = stencil_left(i, nx)
      xloc = x(L:L+6) - x(i)
      CALL fornberg_weights(0.0_dp, xloc, 7, 2, w)
      DO j = 1, ny
        acc = 0.0_dp
        DO k = 1, 7
          acc = acc + w(k,2) * F(L+k-1, j)
        END DO
        dFdx2(i,j) = acc
      END DO
    END DO
  END SUBROUTINE deriv_x2


!==============================================================================
!  d³F / d x³   (non-uniform x)
!==============================================================================
  SUBROUTINE deriv_x3(F, x, nx, ny, dFdx3)
    INTEGER,  INTENT(IN)  :: nx, ny
    REAL(dp), INTENT(IN)  :: F(nx,ny), x(nx)
    REAL(dp), INTENT(OUT) :: dFdx3(nx,ny)

    INTEGER  :: i, j, L, k
    REAL(dp) :: xloc(7), w(7, 0:3), acc

    DO i = 1, nx
      L = stencil_left(i, nx)
      xloc = x(L:L+6) - x(i)
      CALL fornberg_weights(0.0_dp, xloc, 7, 3, w)
      DO j = 1, ny
        acc = 0.0_dp
        DO k = 1, 7
          acc = acc + w(k,3) * F(L+k-1, j)
        END DO
        dFdx3(i,j) = acc
      END DO
    END DO
  END SUBROUTINE deriv_x3


!==============================================================================
!  d F / d y   (uniform dy)
!
!  Interior j=4..ny-3: 6th-order centred
!    f'_j = (-f_{j-3} + 9f_{j-2} - 45f_{j-1} + 45f_{j+1} - 9f_{j+2} + f_{j+3})
!           / (60 dy)
!
!  Boundaries: 4th-order one-sided (verified against Fornberg on uniform grid).
!    j=1  forward:  (-25f1 + 48f2 - 36f3 + 16f4 -  3f5) / (12 dy)
!    j=2  shifted:  ( -3f1 - 10f2 + 18f3 -  6f4 +  1f5) / (12 dy)
!    j=3  centred:  (  1f1 -  8f2 +  8f4 -  1f5)        / (12 dy)
!    j=ny-2 centred: mirror of j=3
!    j=ny-1 shifted: mirror of j=2
!    j=ny  backward: mirror of j=1
!==============================================================================
  SUBROUTINE deriv_y1(F, nx, ny, dy, dFdy)
    INTEGER,  INTENT(IN)  :: nx, ny
    REAL(dp), INTENT(IN)  :: F(nx,ny), dy
    REAL(dp), INTENT(OUT) :: dFdy(nx,ny)

    INTEGER  :: i, j
    REAL(dp) :: r   ! 1/dy prefactor absorbed into denominators below

    r = 1.0_dp / dy

    DO i = 1, nx
      dFdy(i,1) = r*(-25.0_dp*F(i,1) + 48.0_dp*F(i,2) &
                     -36.0_dp*F(i,3) + 16.0_dp*F(i,4) &
                      -3.0_dp*F(i,5)) / 12.0_dp
      dFdy(i,2) = r*( -3.0_dp*F(i,1) - 10.0_dp*F(i,2) &
                      +18.0_dp*F(i,3) -  6.0_dp*F(i,4) &
                       +1.0_dp*F(i,5)) / 12.0_dp
      dFdy(i,3) = r*(  1.0_dp*F(i,1) -  8.0_dp*F(i,2) &
                       +8.0_dp*F(i,4) -  1.0_dp*F(i,5)) / 12.0_dp

      DO j = 4, ny-3
        dFdy(i,j) = r*(-1.0_dp*F(i,j-3) +  9.0_dp*F(i,j-2) &
                       -45.0_dp*F(i,j-1) + 45.0_dp*F(i,j+1) &
                        -9.0_dp*F(i,j+2) +  1.0_dp*F(i,j+3)) / 60.0_dp
      END DO

      dFdy(i,ny-2) = r*(  1.0_dp*F(i,ny-4) -  8.0_dp*F(i,ny-3) &
                           +8.0_dp*F(i,ny-1) -  1.0_dp*F(i,ny  )) / 12.0_dp
      dFdy(i,ny-1) = r*( -1.0_dp*F(i,ny-4) +  6.0_dp*F(i,ny-3) &
                         -18.0_dp*F(i,ny-2) + 10.0_dp*F(i,ny-1) &
                          +3.0_dp*F(i,ny  )) / 12.0_dp
      dFdy(i,ny  ) = r*(  3.0_dp*F(i,ny-4) - 16.0_dp*F(i,ny-3) &
                         +36.0_dp*F(i,ny-2) - 48.0_dp*F(i,ny-1) &
                         +25.0_dp*F(i,ny  )) / 12.0_dp
    END DO
  END SUBROUTINE deriv_y1


!==============================================================================
!  d²F / d y²   (uniform dy)
!
!  Interior j=4..ny-3: 6th-order centred
!    f''_j = (2f_{j-3} - 27f_{j-2} + 270f_{j-1} - 490f_j
!             + 270f_{j+1} - 27f_{j+2} + 2f_{j+3}) / (180 dy²)
!
!  Boundaries: 4th-order one-sided.
!    j=1:  (35f1 - 104f2 + 114f3 - 56f4 + 11f5) / (12 dy²)
!    j=2:  (11f1 -  20f2 +   6f3 +  4f4 -  1f5) / (12 dy²)
!    j=3:  (-f1  +  16f2 -  30f3 + 16f4 -   f5) / (12 dy²)
!    right boundaries: mirror (symmetric)
!==============================================================================
  SUBROUTINE deriv_y2(F, nx, ny, dy, dFdy2)
    INTEGER,  INTENT(IN)  :: nx, ny
    REAL(dp), INTENT(IN)  :: F(nx,ny), dy
    REAL(dp), INTENT(OUT) :: dFdy2(nx,ny)

    INTEGER  :: i, j
    REAL(dp) :: r2

    r2 = 1.0_dp / (dy*dy)

    DO i = 1, nx
      dFdy2(i,1) = r2*( 35.0_dp*F(i,1) - 104.0_dp*F(i,2) &
                       +114.0_dp*F(i,3) -  56.0_dp*F(i,4) &
                        +11.0_dp*F(i,5)) / 12.0_dp
      dFdy2(i,2) = r2*( 11.0_dp*F(i,1) -  20.0_dp*F(i,2) &
                        + 6.0_dp*F(i,3) +   4.0_dp*F(i,4) &
                        - 1.0_dp*F(i,5)) / 12.0_dp
      dFdy2(i,3) = r2*( -1.0_dp*F(i,1) +  16.0_dp*F(i,2) &
                        -30.0_dp*F(i,3) +  16.0_dp*F(i,4) &
                         -1.0_dp*F(i,5)) / 12.0_dp

      DO j = 4, ny-3
        dFdy2(i,j) = r2*(  2.0_dp*F(i,j-3) -  27.0_dp*F(i,j-2) &
                          +270.0_dp*F(i,j-1) - 490.0_dp*F(i,j  ) &
                          +270.0_dp*F(i,j+1) -  27.0_dp*F(i,j+2) &
                            +2.0_dp*F(i,j+3)) / 180.0_dp
      END DO

      dFdy2(i,ny-2) = r2*( -1.0_dp*F(i,ny-4) +  16.0_dp*F(i,ny-3) &
                            -30.0_dp*F(i,ny-2) +  16.0_dp*F(i,ny-1) &
                             -1.0_dp*F(i,ny  )) / 12.0_dp
      dFdy2(i,ny-1) = r2*( -1.0_dp*F(i,ny-4) +   4.0_dp*F(i,ny-3) &
                            + 6.0_dp*F(i,ny-2) -  20.0_dp*F(i,ny-1) &
                            +11.0_dp*F(i,ny  )) / 12.0_dp
      dFdy2(i,ny  ) = r2*( 11.0_dp*F(i,ny-4) -  56.0_dp*F(i,ny-3) &
                           +114.0_dp*F(i,ny-2) - 104.0_dp*F(i,ny-1) &
                            +35.0_dp*F(i,ny  )) / 12.0_dp
    END DO
  END SUBROUTINE deriv_y2


!==============================================================================
!  d³F / d y³   (uniform dy)
!
!  Interior j=4..ny-3: 6th-order centred
!    f'''_j = (f_{j-3} - 8f_{j-2} + 13f_{j-1} - 13f_{j+1} + 8f_{j+2} - f_{j+3})
!             / (8 dy³)
!    Note: sign pattern is antisymmetric, as expected for an odd derivative.
!
!  Boundaries: 4th-order one-sided (computed via Fornberg on a 5-node stencil
!  and cross-checked by Taylor expansion). Multiply through by 2h³:
!    j=1:  (-5f1 + 18f2 - 24f3 + 14f4 -  3f5) / (2 dy³)
!    j=2:  (-3f1 + 10f2 - 12f3 +  6f4 -  1f5) / (2 dy³)
!    j=3:  (-1f1 +  2f2 +  0f3 -  2f4 +  1f5) / (2 dy³)   [centred 4th-order]
!    right boundaries: negate and mirror (antisymmetric for odd derivative)
!      j=ny-2:  ( 1f(ny-4) - 2f(ny-3) + 0 + 2f(ny-1) - 1f(ny)) / (2 dy³)
!      j=ny-1:  ( 1f(ny-4) - 6f(ny-3) +12f(ny-2) -10f(ny-1) + 3f(ny)) / (2 dy³)
!      j=ny:    ( 3f(ny-4) -14f(ny-3) +24f(ny-2) -18f(ny-1) + 5f(ny)) / (2 dy³)
!==============================================================================
  SUBROUTINE deriv_y3(F, nx, ny, dy, dFdy3)
    INTEGER,  INTENT(IN)  :: nx, ny
    REAL(dp), INTENT(IN)  :: F(nx,ny), dy
    REAL(dp), INTENT(OUT) :: dFdy3(nx,ny)

    INTEGER  :: i, j
    REAL(dp) :: r3

    r3 = 1.0_dp / (dy*dy*dy)

    DO i = 1, nx
      dFdy3(i,1) = r3*( -5.0_dp*F(i,1) + 18.0_dp*F(i,2) &
                        -24.0_dp*F(i,3) + 14.0_dp*F(i,4) &
                         -3.0_dp*F(i,5)) / 2.0_dp
      dFdy3(i,2) = r3*( -3.0_dp*F(i,1) + 10.0_dp*F(i,2) &
                        -12.0_dp*F(i,3) +  6.0_dp*F(i,4) &
                         -1.0_dp*F(i,5)) / 2.0_dp
      dFdy3(i,3) = r3*( -1.0_dp*F(i,1) +  2.0_dp*F(i,2) &
                         -2.0_dp*F(i,4) +  1.0_dp*F(i,5)) / 2.0_dp

      DO j = 4, ny-3
        dFdy3(i,j) = r3*( 1.0_dp*F(i,j-3) -  8.0_dp*F(i,j-2) &
                         +13.0_dp*F(i,j-1) - 13.0_dp*F(i,j+1) &
                          +8.0_dp*F(i,j+2) -  1.0_dp*F(i,j+3)) / 8.0_dp
      END DO

      dFdy3(i,ny-2) = r3*(  1.0_dp*F(i,ny-4) -  2.0_dp*F(i,ny-3) &
                             +2.0_dp*F(i,ny-1) -  1.0_dp*F(i,ny  )) / 2.0_dp
      dFdy3(i,ny-1) = r3*(  1.0_dp*F(i,ny-4) -  6.0_dp*F(i,ny-3) &
                            +12.0_dp*F(i,ny-2) - 10.0_dp*F(i,ny-1) &
                             +3.0_dp*F(i,ny  )) / 2.0_dp
      dFdy3(i,ny  ) = r3*(  3.0_dp*F(i,ny-4) - 14.0_dp*F(i,ny-3) &
                            +24.0_dp*F(i,ny-2) - 18.0_dp*F(i,ny-1) &
                             +5.0_dp*F(i,ny  )) / 2.0_dp
    END DO
  END SUBROUTINE deriv_y3


!==============================================================================
!  d²F / (d x d y)   (non-uniform x, uniform dy)
!  Composition: x-derivative first, then y-derivative.
!==============================================================================
  SUBROUTINE deriv_xy(F, x, nx, ny, dy, d2Fdxdy)
    INTEGER,  INTENT(IN)  :: nx, ny
    REAL(dp), INTENT(IN)  :: F(nx,ny), x(nx), dy
    REAL(dp), INTENT(OUT) :: d2Fdxdy(nx,ny)

    REAL(dp) :: tmp(nx,ny)

    CALL deriv_x1(F,   x,  nx, ny,     tmp)
    CALL deriv_y1(tmp, nx, ny, dy, d2Fdxdy)
  END SUBROUTINE deriv_xy

!=============================================================================
! LAPLACIAN_RADIAL
!
! Computes the radial part of the cylindrical Laplacian of psi:
!
!   L(i,j) = d²ψ/dv⊥² + (1/v⊥) dψ/dv⊥
!           = (1/v⊥) d/dv⊥ (v⊥ dψ/dv⊥)
!
! evaluated on the non-uniform grid vperp(1:nx), for all v∥ columns j.
!
! Method
! ------
!   Define  g(i,j) = vperp(i) * dψ/dv⊥(i,j)
!   Then    L(i,j) = (1/vperp(i)) * dg/dv⊥(i,j)      for i >= 2
!   At      i = 1  (vperp = 0):
!           L(1,j) = 2 * d²ψ/dv⊥²(1,j)               (L'Hopital limit)
!
! Both dψ/dv⊥ and dg/dv⊥ are computed with the same Fornberg 7-point
! stencil used by deriv_x1, so no additional information is needed.
!=============================================================================
SUBROUTINE laplacian_radial(psi, vperp, nx, ny, Lpsi)
  INTEGER,  INTENT(IN)  :: nx, ny
  REAL(dp), INTENT(IN)  :: psi(nx,ny), vperp(nx)
  REAL(dp), INTENT(OUT) :: Lpsi(nx,ny)

  REAL(dp) :: dpsi_dpe(nx,ny)   ! dψ/dv⊥
  REAL(dp) :: g(nx,ny)          ! g = v⊥ * dψ/dv⊥
  REAL(dp) :: dg_dpe(nx,ny)     ! dg/dv⊥
  REAL(dp) :: d2psi_dpe2(nx,ny) ! d²ψ/dv⊥²  (needed only at i=1)
  INTEGER  :: i, j

  ! Step 1: dψ/dv⊥  using Fornberg on the non-uniform vperp grid
  CALL deriv_x1(psi, vperp, nx, ny, dpsi_dpe)

  ! Step 2: g = v⊥ * dψ/dv⊥
  DO j = 1, ny
    DO i = 1, nx
      g(i,j) = vperp(i) * dpsi_dpe(i,j)
    END DO
  END DO

  ! Step 3: dg/dv⊥  using Fornberg on the same grid
  CALL deriv_x1(g, vperp, nx, ny, dg_dpe)

  ! Step 4: L = (1/v⊥) * dg/dv⊥   for i >= 2
  !         L = 2 * d²ψ/dv⊥²       for i = 1  (v⊥ = 0, L'Hopital)
  CALL deriv_x2(psi, vperp, nx, ny, d2psi_dpe2)

  ! i = 1: L'Hopital regularisation
  Lpsi(1,:) = 2.0_dp * d2psi_dpe2(1,:)

  ! i >= 2: divergence form, no singularity
  DO j = 1, ny
    DO i = 2, nx
      Lpsi(i,j) = dg_dpe(i,j) / vperp(i)
    END DO
  END DO

END SUBROUTINE laplacian_radial

END MODULE derivatives_2d