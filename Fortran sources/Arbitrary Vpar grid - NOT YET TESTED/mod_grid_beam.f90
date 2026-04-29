!**********************************************************
!
! MODULE mod_grid_beam
!
! Provides:
!   1. make_grid_vperp : combined quadratic + beam-refined v⊥ grid
!   2. make_grid_vpar  : combined uniform  + beam-refined v∥ grid
!   3. integrate_2d    : 2D numerical integration over (v⊥, v∥)
!                        using the trapezoidal rule with full
!                        non-uniform spacing in both directions,
!                        including the cylindrical v⊥ factor.
!
! Grid design principle (density function method):
!   The node positions are defined as the inverse cumulative integral
!   of a density function rho(v):
!
!     v⊥: rho(v) = 1/(v + v0)  +  A * exp(-(v-v_beam)^2 / 2*sr^2)
!         Base term gives quadratic spacing near v=0.
!
!     v∥: rho(v) = 1/v0_pa     +  A * exp(-(v-v_beam)^2 / 2*sr^2)
!         Base term gives uniform spacing.
!
!   In both cases: rho_beam amplitude A is chosen so that n_refine
!   extra nodes are concentrated in the beam region.
!   All spacing ratios h(i+1)/h(i) are <= 1.4 by construction
!   (safe for 7-point Fornberg stencil).
!
!   by Fabrice Louche
!
!**********************************************************

MODULE mod_grid_beam

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: make_grid_vperp, make_grid_vpar, integrate_2d

  INTEGER,  PARAMETER :: dp = KIND(1.0D0)
  REAL(dp), PARAMETER :: pi = 3.141592653589793_dp
  INTEGER,  PARAMETER :: n_fine = 100000   ! resolution of density function

CONTAINS

  !================================================================
  ! SUBROUTINE make_grid_vperp
  !
  ! Builds the v⊥ grid on [0, vperp_max] with:
  !   - quadratic-like spacing near v⊥ = 0  (1/(v+v0) base density)
  !   - local refinement of n_refine extra nodes around vpe_s
  !
  ! Input:
  !   vperp_max          : maximum v⊥ [m/s]
  !   nperp              : number of grid points
  !   vpe_s              : beam centre in v⊥ [m/s]
  !   sigma_vpe          : beam width in v⊥ [m/s]
  !   n_refine           : number of extra nodes in beam region
  !   sigma_refine_factor: refinement width = factor * sigma_vpe
  !                        (0.8 recommended: smooth transitions)
  !
  ! Output:
  !   vperp_grid(nperp)  : node positions [m/s], vperp_grid(1) = 0
  !================================================================
  SUBROUTINE make_grid_vperp(vperp_max, nperp, vpe_s, sigma_vpe, &
                              n_refine, sigma_refine_factor, vperp_grid)

    REAL(dp), INTENT(IN)  :: vperp_max, vpe_s, sigma_vpe, sigma_refine_factor
    INTEGER,  INTENT(IN)  :: nperp, n_refine
    REAL(dp), INTENT(OUT) :: vperp_grid(nperp)

    REAL(dp) :: v_fine(n_fine), rho(n_fine), cum(n_fine)
    REAL(dp) :: rho_base, A_beam, C_base, sigma_r, v0
    REAL(dp) :: dv, cum_total, cum_target
    INTEGER  :: k, inode

    ! Fine mesh
    dv = vperp_max / REAL(n_fine-1, dp)
    DO k = 1, n_fine
      v_fine(k) = REAL(k-1, dp) * dv
    END DO

    ! Base density: 1/(v + v0) -> quadratic near v=0
    v0 = vperp_max / REAL(nperp - n_refine, dp)

    ! Beam refinement width
    sigma_r = sigma_refine_factor * sigma_vpe

    ! Compute base integral C_base via trapezoidal rule
    C_base = 0.0_dp
    DO k = 2, n_fine
      C_base = C_base + 0.5_dp * (1.0_dp/(v_fine(k)+v0) + 1.0_dp/(v_fine(k-1)+v0)) * dv
    END DO

    ! Beam amplitude: n_refine extra nodes in beam region
    A_beam = REAL(n_refine, dp) * C_base &
           / (REAL(nperp - n_refine, dp) * sigma_r * SQRT(2.0_dp*pi))

    ! Build density and cumulative integral
    cum(1) = 0.0_dp
    DO k = 2, n_fine
      rho(k-1) = 1.0_dp/(v_fine(k-1)+v0) &
               + A_beam * EXP(-0.5_dp*((v_fine(k-1)-vpe_s)/sigma_r)**2)
      rho(k)   = 1.0_dp/(v_fine(k)  +v0) &
               + A_beam * EXP(-0.5_dp*((v_fine(k)  -vpe_s)/sigma_r)**2)
      cum(k)   = cum(k-1) + 0.5_dp*(rho(k)+rho(k-1)) * dv
    END DO
    cum_total = cum(n_fine)

    ! Invert: uniform spacing in cumulative space -> node positions
    vperp_grid(1)     = 0.0_dp
    vperp_grid(nperp) = vperp_max
    DO inode = 2, nperp-1
      cum_target = REAL(inode-1, dp) / REAL(nperp-1, dp) * cum_total
      ! Binary search in cum(:)
      vperp_grid(inode) = interp_inverse(cum, v_fine, n_fine, cum_target)
    END DO

  END SUBROUTINE make_grid_vperp


  !================================================================
  ! SUBROUTINE make_grid_vpar
  !
  ! Builds the v∥ grid on [-vpar_max, +vpar_max] with:
  !   - uniform base spacing
  !   - local refinement of n_refine extra nodes around vpa_s
  !
  ! Input:
  !   vpar_max           : maximum |v∥| [m/s]
  !   npar               : number of grid points
  !   vpa_s              : beam centre in v∥ [m/s]
  !   sigma_vpa          : beam width in v∥ [m/s]
  !   n_refine           : number of extra nodes in beam region
  !   sigma_refine_factor: refinement width = factor * sigma_vpa
  !
  ! Output:
  !   vpar_grid(npar)    : node positions [m/s],
  !                        vpar_grid(1) = -vpar_max, vpar_grid(npar) = +vpar_max
  !================================================================
  SUBROUTINE make_grid_vpar(vpar_max, npar, vpa_s, sigma_vpa, &
                             n_refine, sigma_refine_factor, vpar_grid)

    REAL(dp), INTENT(IN)  :: vpar_max, vpa_s, sigma_vpa, sigma_refine_factor
    INTEGER,  INTENT(IN)  :: npar, n_refine
    REAL(dp), INTENT(OUT) :: vpar_grid(npar)

    REAL(dp) :: v_fine(n_fine), rho(n_fine), cum(n_fine)
    REAL(dp) :: A_beam, C_base, sigma_r, v0_pa
    REAL(dp) :: dv, cum_total, cum_target
    REAL(dp) :: v_lo, v_hi
    INTEGER  :: k, inode

    v_lo = -vpar_max;  v_hi = vpar_max
    dv   = (v_hi - v_lo) / REAL(n_fine-1, dp)
    DO k = 1, n_fine
      v_fine(k) = v_lo + REAL(k-1, dp) * dv
    END DO

    ! Uniform base density
    v0_pa  = (v_hi - v_lo) / REAL(npar - n_refine, dp)
    sigma_r = sigma_refine_factor * sigma_vpa

    ! Base integral
    C_base = (v_hi - v_lo) / v0_pa   ! exact for uniform

    ! Beam amplitude
    A_beam = REAL(n_refine, dp) * C_base &
           / (REAL(npar - n_refine, dp) * sigma_r * SQRT(2.0_dp*pi))

    ! Build density and cumulative integral
    cum(1) = 0.0_dp
    DO k = 2, n_fine
      rho(k-1) = 1.0_dp/v0_pa &
               + A_beam * EXP(-0.5_dp*((v_fine(k-1)-vpa_s)/sigma_r)**2)
      rho(k)   = 1.0_dp/v0_pa &
               + A_beam * EXP(-0.5_dp*((v_fine(k)  -vpa_s)/sigma_r)**2)
      cum(k)   = cum(k-1) + 0.5_dp*(rho(k)+rho(k-1)) * dv
    END DO
    cum_total = cum(n_fine)

    ! Invert
    vpar_grid(1)    = v_lo
    vpar_grid(npar) = v_hi
    DO inode = 2, npar-1
      cum_target = REAL(inode-1, dp) / REAL(npar-1, dp) * cum_total
      vpar_grid(inode) = interp_inverse(cum, v_fine, n_fine, cum_target)
    END DO

  END SUBROUTINE make_grid_vpar


  !================================================================
  ! SUBROUTINE integrate_2d
  !
  ! Computes:
  !   I = 2*pi * integral_{v⊥=0}^{vperp_max}
  !             integral_{v∥=-vpar_max}^{+vpar_max}
  !             f(v⊥, v∥) * v⊥  dv⊥ dv∥
  !
  ! using the 2D trapezoidal rule on the non-uniform grid.
  ! Boundaries (Dirichlet f=0) are excluded automatically.
  !
  ! This gives: n = 2*pi * integral f * vperp dvperp dvpar
  !
  ! Input:
  !   f(nperp,npar)      : distribution function on the grid
  !   vperp(nperp)       : v⊥ grid [m/s], vperp(1)=0
  !   vpar(npar)         : v∥ grid [m/s]
  !   nperp, npar        : grid dimensions
  !
  ! Output:
  !   integral           : value of I
  !================================================================
  SUBROUTINE integrate_2d(f, vperp, vpar, nperp, npar, integral)

    INTEGER,  INTENT(IN)  :: nperp, npar
    REAL(dp), INTENT(IN)  :: f(nperp,npar), vperp(nperp), vpar(npar)
    REAL(dp), INTENT(OUT) :: integral

    REAL(dp) :: dvp, dva, fij, vp
    REAL(dp) :: w_vp(nperp), w_va(npar)
    INTEGER  :: iv, ip

    ! Trapezoidal weights in v⊥ (non-uniform)
    ! w(1)     = (vperp(2) - vperp(1)) / 2
    ! w(i)     = (vperp(i+1) - vperp(i-1)) / 2   for i=2..nperp-1
    ! w(nperp) = (vperp(nperp) - vperp(nperp-1)) / 2
    w_vp(1)     = 0.5_dp * (vperp(2) - vperp(1))
    w_vp(nperp) = 0.5_dp * (vperp(nperp) - vperp(nperp-1))
    DO iv = 2, nperp-1
      w_vp(iv) = 0.5_dp * (vperp(iv+1) - vperp(iv-1))
    END DO

    ! Trapezoidal weights in v∥ (non-uniform)
    w_va(1)    = 0.5_dp * (vpar(2) - vpar(1))
    w_va(npar) = 0.5_dp * (vpar(npar) - vpar(npar-1))
    DO ip = 2, npar-1
      w_va(ip) = 0.5_dp * (vpar(ip+1) - vpar(ip-1))
    END DO

    ! Double sum: 2*pi * sum_{iv,ip} f(iv,ip) * vperp(iv) * w_vp(iv) * w_va(ip)
    integral = 0.0_dp
    DO iv = 1, nperp
      vp = vperp(iv)
      DO ip = 1, npar
        integral = integral + f(iv,ip) * vp * w_vp(iv) * w_va(ip)
      END DO
    END DO
    integral = 2.0_dp * pi * integral

  END SUBROUTINE integrate_2d


  !================================================================
  ! FUNCTION interp_inverse  (private)
  !
  ! Given a monotone increasing array cum(1:n) and corresponding
  ! positions x(1:n), returns the x value where cum = target,
  ! via binary search + linear interpolation.
  !================================================================
  FUNCTION interp_inverse(cum, x, n, target) RESULT(v)
    INTEGER,  INTENT(IN) :: n
    REAL(dp), INTENT(IN) :: cum(n), x(n), target
    REAL(dp) :: v

    INTEGER  :: lo, hi, mid
    REAL(dp) :: frac

    lo = 1;  hi = n
    DO WHILE (hi - lo > 1)
      mid = (lo + hi) / 2
      IF (cum(mid) <= target) THEN
        lo = mid
      ELSE
        hi = mid
      END IF
    END DO

    IF (cum(hi) == cum(lo)) THEN
      v = x(lo)
    ELSE
      frac = (target - cum(lo)) / (cum(hi) - cum(lo))
      v    = x(lo) + frac * (x(hi) - x(lo))
    END IF

  END FUNCTION interp_inverse

END MODULE mod_grid_beam
