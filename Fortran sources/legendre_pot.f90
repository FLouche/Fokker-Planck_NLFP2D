!***********************************************************************
!*  mod_legendre_pot -- Rosenbluth potentials of one species by        *
!*  Legendre expansion, evaluated anywhere (any species' grid), and    *
!*  the Fokker-Planck coefficients they give for any pair.             *
!*                                                                     *
!*  Stage 2 of docs/Multispecies_design.pdf.                           *
!*                                                                     *
!*  SOURCE SIDE (legpot_build, on the grid of the source species b,    *
!*  loaded in shared_grid):                                            *
!*    f_l(v) = (2l+1)/2 int_{-1}^{1} f(v sqrt(1-xi^2), v xi) P_l(xi) dxi *
!*  on a uniform radial grid, xi by Gauss-Legendre, f by 4x4 Lagrange  *
!*  interpolation (f even in v_perp across the axis, f = 0 outside the *
!*  box).  Between radial nodes f_l is taken linear, and the four      *
!*  cumulative moments                                                 *
!*    I_k(v) = int_0^v f_l v'^k dv'   (k = l+2, l+4)                   *
!*    J_k(v) = int_v^inf f_l v'^k dv' (k = 1-l, 3-l)                   *
!*  are integrated exactly, so they are known at any v.                *
!*                                                                     *
!*  TARGET SIDE (legpot_eval, pair_coefs): with x = v/v0, K=4pi/(2l+1) *
!*    h_l = K [x^{-l-1} I_{l+2} + x^l J_{1-l}]                         *
!*    g_l = K [x^{-l-1} I_{l+4}/(2l+3) - x^{1-l} I_{l+2}/(2l-1)        *
!*             + x^{l+2} J_{1-l}/(2l+3) - x^l J_{3-l}/(2l-1)]          *
!*  and their v-derivatives analytically (the boundary terms of the    *
!*  cumulative integrals cancel).  The code's potentials are           *
!*    psi = -g/8pi,   phi = laplacian(psi) = -h/4pi,                   *
!*  and the Cartesian (v_perp, v_par) derivatives follow from the      *
!*  chain rule in (v, xi = v_par/v).                                   *
!*                                                                     *
!*  COEFFICIENTS for target a, source b, mu = m_a/m_b:                 *
!*    C(f_a,f_b) = coef [psi_ij f_a,ij + (1-mu) phi_j f_a,j            *
!*                       - mu laplacian(phi) f_a],  lap(phi) = f_b,     *
!*    coef = -4 pi gamma0 lnLambda_ab (Z_a/A_a)^2 Z_b^2,               *
!*  i.e. main_nlterm's expansion with the friction terms times mu.     *
!*                                                                     *
!*  20/11/2026                                                          *
!***********************************************************************

module mod_legendre_pot

  implicit none
  private

  public :: legpot_t, legpot_build, legpot_free, legpot_eval, pair_coefs, &
            pair_coef_factor, legpot_density

  double precision, parameter :: pi = 3.141592653589793d0
  double precision, parameter :: gamma0 = 2.390775d-1   ! as in the solvers

  type legpot_t
    integer :: lmax = -1, nr = 0, nxi = 0
    double precision :: v0 = 0.d0          ! velocity scale: x = v/v0, x_nr = 1
    double precision :: dx = 0.d0          ! radial step in x
    double precision, allocatable :: fl(:,:)            ! (0:nr,0:lmax) f_l at the nodes
    double precision, allocatable :: i2(:,:), i4(:,:)   ! I_{l+2}, I_{l+4} at the nodes
    double precision, allocatable :: j1(:,:), j3(:,:)   ! J_{1-l}, J_{3-l} at the nodes (m>=1)
  end type legpot_t

  ! Everything legpot_eval returns at one point (physical units).
  type, public :: legval_t
    double precision :: psi, psi_pe, psi_pa, psi_pepe, psi_pepa, psi_papa
    double precision :: psi_pe_on_vpe2     ! psi_pe / v_perp^2
    double precision :: phi, phi_pe, phi_pa
    double precision :: fb                 ! f_b itself (= laplacian(phi))
  end type legval_t

contains

!-----------------------------------------------------------------------
! Coefficient factor of a pair: -4 pi gamma0 lnL (Z_a/A_a)^2 Z_b^2.
!-----------------------------------------------------------------------
  double precision function pair_coef_factor(za_t, aa_t, zb_s, lnl)
    double precision, intent(in) :: za_t, aa_t, zb_s, lnl
    pair_coef_factor = -4.d0 * pi * gamma0 * lnl * (za_t / aa_t)**2 * zb_s**2
  end function pair_coef_factor

!-----------------------------------------------------------------------
! int_a^b x^k dx
!-----------------------------------------------------------------------
  pure double precision function pw(k, a, b)
    integer, intent(in) :: k
    double precision, intent(in) :: a, b
    if (k == -1) then
      pw = log(b / a)
    else
      pw = (b**(k+1) - a**(k+1)) / dble(k+1)
    end if
  end function pw

!-----------------------------------------------------------------------
! int_a^b f(x) x^k dx with f linear through (xa,fa), (xb,fb).
!-----------------------------------------------------------------------
  pure double precision function cellint(k, fa, fb, xa, xb, a, b)
    integer, intent(in) :: k
    double precision, intent(in) :: fa, fb, xa, xb, a, b
    double precision :: alpha, beta
    if (b <= a) then
      cellint = 0.d0
      return
    end if
    beta  = (fb - fa) / (xb - xa)
    alpha = fa - beta * xa
    cellint = alpha * pw(k, a, b) + beta * pw(k+1, a, b)
  end function cellint

!-----------------------------------------------------------------------
! Gauss-Legendre nodes and weights on [-1,1].
!-----------------------------------------------------------------------
  subroutine gauss_legendre(n, x, w)
    integer, intent(in) :: n
    double precision, intent(out) :: x(n), w(n)
    integer :: i, j, it
    double precision :: z, z1, p1, p2, p3, pp
    do i = 1, (n + 1) / 2
      z = cos(pi * (dble(i) - 0.25d0) / (dble(n) + 0.5d0))
      do it = 1, 100
        p1 = 1.d0;  p2 = 0.d0
        do j = 1, n
          p3 = p2;  p2 = p1
          p1 = ((2.d0*j - 1.d0) * z * p2 - (j - 1.d0) * p3) / dble(j)
        end do
        pp = dble(n) * (z * p1 - p2) / (z*z - 1.d0)
        z1 = z
        z  = z1 - p1 / pp
        if (abs(z - z1) < 1.d-15) exit
      end do
      x(i) = -z;  x(n+1-i) = z
      w(i) = 2.d0 / ((1.d0 - z*z) * pp*pp);  w(n+1-i) = w(i)
    end do
  end subroutine gauss_legendre

!-----------------------------------------------------------------------
! f at (vpe, vpa) on the grid in shared_grid: 4x4 Lagrange interpolation.
! f is even in v_perp (nodes mirrored across the axis), zero outside the
! box (the walls are Dirichlet f = 0).
!-----------------------------------------------------------------------
  double precision function interp_f(f, vpe, vpa)
    use shared_grid
    double precision, intent(in) :: f(nperp,npar), vpe, vpa
    integer :: i0, j0, lo, hi, mid, ii, jj, k, m
    integer :: iidx(4), jidx(4)
    double precision :: xn(4), wi(4), wj(4), t

    interp_f = 0.d0
    if (vpa < vpar(1) .or. vpa > vpar(npar) .or. vpe > vperp(nperp)) return

    ! v_par: uniform; j0 such that vpar(j0) <= vpa < vpar(j0+1)
    t  = (vpa - vpar(1)) / dvpar
    j0 = min(max(int(t) + 1, 1), npar - 1)
    jj = min(max(j0 - 1, 1), npar - 3)
    do k = 1, 4
      jidx(k) = jj + k - 1
      xn(k)   = vpar(jidx(k))
    end do
    call lagw(xn, vpa, wj)

    ! v_perp: i0 such that vp(i0) <= vpe < vp(i0+1), vp(i) = vperp(i) for
    ! i >= 1 and -vperp(1-i) for i <= 0 (mirror).
    if (vpe < vperp(1)) then
      i0 = 0
    else
      lo = 1;  hi = nperp
      do while (hi - lo > 1)
        mid = (lo + hi) / 2
        if (vperp(mid) <= vpe) then
          lo = mid
        else
          hi = mid
        end if
      end do
      i0 = lo
    end if
    ii = min(i0 - 1, nperp - 3)
    do k = 1, 4
      m = ii + k - 1
      if (m >= 1) then
        iidx(k) = m;      xn(k) = vperp(m)
      else
        iidx(k) = 1 - m;  xn(k) = -vperp(1 - m)
      end if
    end do
    call lagw(xn, vpe, wi)

    do k = 1, 4
      do m = 1, 4
        interp_f = interp_f + wi(k) * wj(m) * f(iidx(k), jidx(m))
      end do
    end do

  contains
    subroutine lagw(xn, x, w)
      double precision, intent(in)  :: xn(4), x
      double precision, intent(out) :: w(4)
      integer :: a, b
      do a = 1, 4
        w(a) = 1.d0
        do b = 1, 4
          if (b /= a) w(a) = w(a) * (x - xn(b)) / (xn(a) - xn(b))
        end do
      end do
    end subroutine lagw
  end function interp_f

!-----------------------------------------------------------------------
! Build the Legendre potentials of f, given on the grid now in shared_grid.
!   lmax   highest harmonic
!   nr_in  radial intervals (0: auto, step dvpar/4)
!   nxi_in Gauss points in xi (0: auto, max(64, 4 lmax + 8))
!-----------------------------------------------------------------------
  subroutine legpot_build(p, f, lmax, nr_in, nxi_in)
    use shared_grid
    type(legpot_t), intent(inout) :: p
    double precision, intent(in) :: f(nperp,npar)
    integer, intent(in) :: lmax, nr_in, nxi_in

    integer :: nr, nxi, m, k, l
    double precision :: vmax, v, fv, xa, xb
    double precision, allocatable :: xi(:), wxi(:), pl(:,:)

    call legpot_free(p)

    vmax = sqrt(vperp(nperp)**2 + max(vpar(1)**2, vpar(npar)**2))
    if (nr_in > 0) then
      nr = nr_in
    else
      nr = min(max(ceiling(vmax / (0.25d0 * dvpar)), 64), 40000)
    end if
    if (nxi_in > 0) then
      nxi = nxi_in
    else
      nxi = max(64, 4*lmax + 8)
    end if

    p%lmax = lmax;  p%nr = nr;  p%nxi = nxi
    p%v0 = vmax;    p%dx = 1.d0 / dble(nr)

    allocate(xi(nxi), wxi(nxi), pl(nxi,0:lmax))
    call gauss_legendre(nxi, xi, wxi)
    do k = 1, nxi
      pl(k,0) = 1.d0
      if (lmax >= 1) pl(k,1) = xi(k)
      do l = 1, lmax - 1
        pl(k,l+1) = ((2*l+1) * xi(k) * pl(k,l) - l * pl(k,l-1)) / dble(l+1)
      end do
    end do

    allocate(p%fl(0:nr,0:lmax), p%i2(0:nr,0:lmax), p%i4(0:nr,0:lmax), &
             p%j1(0:nr,0:lmax), p%j3(0:nr,0:lmax))
    p%fl = 0.d0
    do m = 0, nr
      v = dble(m) * p%dx * vmax
      do k = 1, nxi
        fv = interp_f(f, v * sqrt(max(0.d0, 1.d0 - xi(k)**2)), v * xi(k))
        do l = 0, lmax
          p%fl(m,l) = p%fl(m,l) + wxi(k) * fv * pl(k,l)
        end do
      end do
    end do
    do l = 0, lmax
      p%fl(:,l) = p%fl(:,l) * 0.5d0 * dble(2*l+1)
    end do
    p%fl(nr,:) = 0.d0          ! the box corner is a wall point

    ! Cumulative moments, exact for piecewise-linear f_l.
    p%i2 = 0.d0;  p%i4 = 0.d0;  p%j1 = 0.d0;  p%j3 = 0.d0
    do l = 0, lmax
      do m = 0, nr - 1
        xa = m * p%dx;  xb = (m+1) * p%dx
        p%i2(m+1,l) = p%i2(m,l) + cellint(l+2, p%fl(m,l), p%fl(m+1,l), xa, xb, xa, xb)
        p%i4(m+1,l) = p%i4(m,l) + cellint(l+4, p%fl(m,l), p%fl(m+1,l), xa, xb, xa, xb)
      end do
      do m = nr - 1, 1, -1
        xa = m * p%dx;  xb = (m+1) * p%dx
        p%j1(m,l) = p%j1(m+1,l) + cellint(1-l, p%fl(m,l), p%fl(m+1,l), xa, xb, xa, xb)
        p%j3(m,l) = p%j3(m+1,l) + cellint(3-l, p%fl(m,l), p%fl(m+1,l), xa, xb, xa, xb)
      end do
    end do

    deallocate(xi, wxi, pl)
  end subroutine legpot_build

  subroutine legpot_free(p)
    type(legpot_t), intent(inout) :: p
    if (allocated(p%fl)) deallocate(p%fl, p%i2, p%i4, p%j1, p%j3)
    p%lmax = -1;  p%nr = 0
  end subroutine legpot_free

!-----------------------------------------------------------------------
! Density of the expansion, 4 pi int f_0 v^2 dv (to compare with the 2D
! quadrature of f on the source grid).
!-----------------------------------------------------------------------
  double precision function legpot_density(p)
    type(legpot_t), intent(in) :: p
    legpot_density = 4.d0 * pi * p%v0**3 * p%i2(p%nr, 0)
  end function legpot_density

!-----------------------------------------------------------------------
! Potentials and their derivatives at (vpe, vpa), any point (v > 0).
!-----------------------------------------------------------------------
  subroutine legpot_eval(p, vpe, vpa, r)
    type(legpot_t), intent(in) :: p
    double precision, intent(in) :: vpe, vpa
    type(legval_t), intent(out) :: r

    integer :: l, m
    double precision :: v, x, xa, xb, s, c, a, b, v0
    double precision :: i2, i4, j1, j3, flx, fa, fb_, kl, dl
    double precision :: hx, hxp, gx, gxp, gxpp
    double precision :: psil, psilp, psilpp, phil, philp
    double precision :: pl, plm1, plp, plpm1, plpp, plppm1, tmp, newp, newpp
    double precision :: fv, fvv, fxi, fvxi, fxixi, gv, gxi, fbsum, psis, phis
    ! running powers (x^l, x^-l, xa^l, xb^-l) and fixed small powers: the loop
    ! over l needs no x**k (the former ~15 power calls per degree and point)
    double precision :: pl_x, ql_x, al_x, bl_x, xm1, xm2, xm3, x2, x3, x4, x5, x6
    double precision :: xa3, xa4, xa5, xa6, xb2, xb3, xb4, xb5, xbinv, alpha, beta
    double precision :: c4, c3, c2, c1

    v0 = p%v0
    v  = sqrt(vpe**2 + vpa**2)
    x  = v / v0
    s  = vpe / v
    c  = vpa / v
    a  = -c * s / v          ! d xi / d v_perp
    b  =  s * s / v          ! d xi / d v_par

    if (x >= 1.d0) then
      m = p%nr
    else
      m = min(int(x / p%dx), p%nr - 1)
    end if
    xa = m * p%dx;  xb = (m + 1) * p%dx

    xm1 = 1.d0 / x;  xm2 = xm1 * xm1;  xm3 = xm2 * xm1
    x2 = x * x;  x3 = x2 * x;  x4 = x3 * x;  x5 = x4 * x;  x6 = x5 * x
    xa3 = xa**3;  xa4 = xa3 * xa;  xa5 = xa4 * xa;  xa6 = xa5 * xa
    xb2 = xb * xb;  xb3 = xb2 * xb;  xb4 = xb3 * xb;  xb5 = xb4 * xb;  xbinv = 1.d0 / xb
    pl_x = 1.d0;  ql_x = 1.d0;  al_x = 1.d0;  bl_x = 1.d0
    c4 = v0**4 / (8.d0*pi);  c3 = v0**3 / (8.d0*pi);  c2 = v0**2 / (8.d0*pi)
    c1 = v0**2 / (4.d0*pi)

    fv = 0.d0; fvv = 0.d0; fxi = 0.d0; fvxi = 0.d0; fxixi = 0.d0
    gv = 0.d0; gxi = 0.d0; fbsum = 0.d0; psis = 0.d0; phis = 0.d0

    ! Legendre recurrences: P_l, P_l', P_l''
    plm1 = 0.d0;  pl = 1.d0
    plpm1 = 0.d0; plp = 0.d0
    plppm1 = 0.d0; plpp = 0.d0

    do l = 0, p%lmax
      if (m >= p%nr) then
        i2 = p%i2(p%nr,l);  i4 = p%i4(p%nr,l);  j1 = 0.d0;  j3 = 0.d0;  flx = 0.d0
      else
        fa = p%fl(m,l);  fb_ = p%fl(m+1,l)
        ! f_l = alpha + beta x on the cell; partial integrals of the cell
        ! (= cellint, with the powers from the running products)
        beta  = (fb_ - fa) / (xb - xa)
        alpha = fa - beta * xa
        i2 = p%i2(m,l) + alpha * (pl_x * x3 - al_x * xa3) / dble(l+3) &
                       + beta  * (pl_x * x4 - al_x * xa4) / dble(l+4)
        i4 = p%i4(m,l) + alpha * (pl_x * x5 - al_x * xa5) / dble(l+5) &
                       + beta  * (pl_x * x6 - al_x * xa6) / dble(l+6)
        j1 = p%j1(m+1,l) + alpha * pint(2-l, bl_x * xb2, ql_x * x2, xb, x) &
                         + beta  * pint(3-l, bl_x * xb3, ql_x * x3, xb, x)
        j3 = p%j3(m+1,l) + alpha * pint(4-l, bl_x * xb4, ql_x * x4, xb, x) &
                         + beta  * pint(5-l, bl_x * xb5, ql_x * x5, xb, x)
        flx = fa + (fb_ - fa) * (x - xa) / (xb - xa)
      end if

      kl = 4.d0 * pi / dble(2*l+1)
      dl = dble(l)
      ! x^(-l-k) = ql_x * x^-k,  x^(l+k) = pl_x * x^k
      hx  = kl * (ql_x * xm1 * i2 + pl_x * j1)
      hxp = kl * (-(dl+1.d0) * ql_x * xm2 * i2 + dl * pl_x * xm1 * j1)
      gx  = kl * ( ql_x * xm1 * i4 / (2*l+3) - ql_x * x * i2 / (2*l-1) &
                 + pl_x * x2 * j1 / (2*l+3) - pl_x * j3 / (2*l-1) )
      gxp = kl * ( -(dl+1.d0) * ql_x * xm2 * i4 / (2*l+3) + (dl-1.d0) * ql_x * i2 / (2*l-1) &
                 + (dl+2.d0) * pl_x * x * j1 / (2*l+3) - dl * pl_x * xm1 * j3 / (2*l-1) )
      gxpp = kl * ( (dl+1.d0)*(dl+2.d0) * ql_x * xm3 * i4 / (2*l+3) &
                  - dl*(dl-1.d0) * ql_x * xm1 * i2 / (2*l-1) &
                  + (dl+1.d0)*(dl+2.d0) * pl_x * j1 / (2*l+3) &
                  - dl*(dl-1.d0) * pl_x * xm2 * j3 / (2*l-1) )

      ! physical units: h ~ v0^2, g ~ v0^4, one 1/v0 per derivative
      psil   = -gx   * c4
      psilp  = -gxp  * c3
      psilpp = -gxpp * c2
      phil   = -hx   * c1
      philp  = -hxp  * c1 / v0

      pl_x = pl_x * x;  ql_x = ql_x * xm1;  al_x = al_x * xa;  bl_x = bl_x * xbinv

      psis  = psis  + psil   * pl
      fv    = fv    + psilp  * pl
      fvv   = fvv   + psilpp * pl
      fxi   = fxi   + psil   * plp
      fvxi  = fvxi  + psilp  * plp
      fxixi = fxixi + psil   * plpp
      phis  = phis  + phil   * pl
      gv    = gv    + philp  * pl
      gxi   = gxi   + phil   * plp
      fbsum = fbsum + flx    * pl

      ! advance to l+1:  P_{l+1} = ((2l+1) xi P_l - l P_{l-1})/(l+1),
      ! P'_{l+1} = P'_{l-1} + (2l+1) P_l,  P''_{l+1} = P''_{l-1} + (2l+1) P'_l
      tmp = ((2*l+1) * c * pl - dl * plm1) / (dl + 1.d0)
      if (l == 0) then
        plpp = 0.d0;  plppm1 = 0.d0          ! P''_1 = 0
        plpm1 = plp;  plp = 1.d0             ! P'_1 = 1
      else
        newp  = plpm1  + (2*l+1) * pl
        newpp = plppm1 + (2*l+1) * plp
        plpm1 = plp;    plp  = newp
        plppm1 = plpp;  plpp = newpp
      end if
      plm1 = pl;  pl = tmp
    end do

    r%psi    = psis
    r%phi    = phis
    r%fb     = fbsum
    r%psi_pe = fv * s + fxi * a
    r%psi_pa = fv * c + fxi * b
    r%psi_pe_on_vpe2 = (fv - fxi * c / v) / (v * vpe)
    r%psi_pepe = fvv * s*s + 2.d0 * fvxi * a * s + fxixi * a*a + fv * c*c / v &
               + fxi * (c / v**2) * (3.d0*s*s - 1.d0)
    r%psi_papa = fvv * c*c + 2.d0 * fvxi * b * c + fxixi * b*b + fv * s*s / v &
               - fxi * 3.d0 * s*s * c / v**2
    r%psi_pepa = fvv * s*c + fvxi * (b*s + a*c) + fxixi * a*b - fv * s*c / v &
               + fxi * (s / v**2) * (3.d0*c*c - 1.d0)
    r%phi_pe = gv * s + gxi * a
    r%phi_pa = gv * c + gxi * b

  contains

    ! int_lo^hi t^(e-1) dt from the powers hi^e, lo^e (log when e = 0)
    pure double precision function pint(e, hip, lop, hi, lo)
      integer, intent(in) :: e
      double precision, intent(in) :: hip, lop, hi, lo
      if (e == 0) then
        pint = log(hi / lo)
      else
        pint = (hip - lop) / dble(e)
      end if
    end function pint
  end subroutine legpot_eval

!-----------------------------------------------------------------------
! FP coefficients on the grid now in shared_grid (target species a) from
! the potentials p of source species b:  mu = m_a/m_b,
! coef = pair_coef_factor(Z_a, A_a, Z_b, lnLambda_ab).
!   c20 d2f/dvpe2 + c11 d2f/dvpe dvpa + c02 d2f/dvpa2
!   + c10 df/dvpe + c01 df/dvpa + c00 f
!-----------------------------------------------------------------------
  subroutine pair_coefs(p, mu, coef, c20, c02, c11, c10, c01, c00)
    use shared_grid
    type(legpot_t), intent(in) :: p
    double precision, intent(in) :: mu, coef
    double precision, dimension(nperp,npar), intent(out) :: c20, c02, c11, c10, c01, c00
    integer :: i, j
    type(legval_t) :: r

    !$OMP PARALLEL DO PRIVATE(i, j, r) SCHEDULE(static)
    do j = 1, npar
      do i = 1, nperp
        call legpot_eval(p, vperp(i), vpar(j), r)
        c20(i,j) = coef * r%psi_pepe
        c02(i,j) = coef * r%psi_papa
        c11(i,j) = coef * 2.d0 * r%psi_pepa
        c10(i,j) = coef * (r%psi_pe_on_vpe2 + (1.d0 - mu) * r%phi_pe)
        c01(i,j) = coef * (1.d0 - mu) * r%phi_pa
        c00(i,j) = -coef * mu * r%fb
      end do
    end do
    !$OMP END PARALLEL DO
  end subroutine pair_coefs

end module mod_legendre_pot
