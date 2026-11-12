!====================================================================
!
!   coul_fv.f90
!
!   FPColl_2D
!
!   Coulomb collision operator in conservative (flux) form, discretised
!   by finite volumes with Chang-Cooper weights.  Selected by the namelist
!   flag i_cc = -1; the default i_cc = 0 keeps the 7-point Fornberg path.
!
!   The operator of cblin (Maxwellian background b) is exactly
!
!     L_b f = (1/vp) d/dvp [vp G_p] + d/dva G_a ,
!     G     = A f + D.grad f ,  A = Psi v ,
!     D_pp  = Theta + va^2 Phi,  D_aa = Theta + vp^2 Phi,  D_pa = -vp va Phi
!
!   and since Psi = (m_a/T_b) Theta, the flux is G = M_b D.grad(f/M_b) with
!   M_b = exp(-kappa_b v^2), kappa_b = maonmb(b)/(2 vt(b)^2) = m_a/(2 T_b).
!
!   Discretisation (cells around the nodes, faces half way between them):
!   - diagonal part of each face flux: Chang-Cooper / Scharfetter-Gummel,
!       G = (D/h) [B(-s) f_+ - B(s) f_-],  B(x) = x/(e^x-1),
!     s = ln M_- - ln M_+ the jump of ln M across the face;
!   - cross part: M d(f/M)/dv at the two nodes of the face, averaged, by
!     Scharfetter-Gummel two-point weights in the transverse direction
!     (zero on M; bounded at high cell Peclet number, see sg_weights).
!   So the discrete flux vanishes exactly on the nodal M_b (cross term
!   included), and the scheme is conservative: the net fluxes telescope.
!   The cross term is not sign-definite: positivity is not guaranteed.
!
!   Cells: v_perp faces at 0 (the axis: zero flux) and (vp_i+vp_i+1)/2,
!   the last cell ends at vp(nperp); v_par faces half way, the first and
!   last cells are half cells.  Walls (i_cc_wall): -1 zero flux (closed
!   box, default), 0 Dirichlet f = 0 at i = nperp, j = 1, j = npar.
!
!   The per-species matrices do not change in time (grid and T_b fixed):
!   they are built once, and each step only rescales them by gammab(b).
!
!   Maxwellian self-collisions (isc = 1, 2, 3) are the same operator with
!   the minority as background: mass ratio 1, thermal speed vteff_t, factor
!   gammaa (self_coll_max).  For isc = 2, 3 vteff_t follows Teff or Tn, so
!   cfv_set_sc recomputes that matrix whenever vteff_t changes.
!
!   The quasilinear RF operator (cfv_set_rf) is a rank-one diffusion with
!   no drag, discretised at the cell corners by a symmetric scheme: see
!   there.  It does not change in time and is built once.
!
!====================================================================
module mod_coul_fv

  use shared_grid
  use shared_plasma, only: nbulk, vt, maonmb, gammab

  implicit none
  private
  public :: cfv_build, cfv_fill, cfv_apply_species, cfv_density, bc_row, &
            coll_apply_species, coll_apply_total, cfv_set_sc, sc_apply, &
            cfv_normalise_source, cfv_set_nl, cfv_nl_off, cfv_set_rf, rf_apply
  public :: cfv_nnz, cfv_ia, cfv_ja

  integer, parameter :: dp = kind(1.d0)

  logical,  save :: built = .false.
  integer,  save :: cfv_nnz = 0
  integer,  allocatable, save :: cfv_ia(:), cfv_ja(:), idiag(:)
  real(dp), allocatable, save :: val(:,:)      ! (nnz, nbulk): L_b without gammab
  real(dp), allocatable, save :: vol(:)        ! cell volume 2 pi R_i W_j
  real(dp), allocatable, save :: eface(:)      ! v_perp faces e(0:nperp)
  ! Maxwellian self-collisions: L_sc without gammaa, for vth_sc
  real(dp), allocatable, save :: val_sc(:)
  real(dp), save :: vth_sc = -1.0_dp, gam_sc = 0.0_dp
  ! Nonlinear self-collisions (isc = -1): nodal D and A from main_nlterm
  ! (factor included), and the matrix built from them each step
  logical,  save :: nl_on = .false.
  real(dp), allocatable, save :: val_nl(:)
  real(dp), allocatable, save :: n_dpp(:,:), n_dpa(:,:), n_daa(:,:), n_ap(:,:), n_aa(:,:)
  ! Quasilinear RF (irf = -1): matrix built once by cfv_set_rf
  real(dp), allocatable, save :: val_rf(:)

contains

  !------------------------------------------------------------------
  ! Rows that are constraints (f = 0 or df/dvp = 0) rather than time
  ! evolution equations.  Flux form: only Dirichlet walls, if chosen.
  !------------------------------------------------------------------
  logical function bc_row(iv, imu)
    integer, intent(in) :: iv, imu
    if (i_cc == -1) then
      bc_row = (i_cc_wall == 0) .and. (iv == nperp .or. imu == 1 .or. imu == npar)
    else
      bc_row = (iv == 1 .or. iv == nperp .or. imu == 1 .or. imu == npar)
    end if
  end function bc_row

  !------------------------------------------------------------------
  ! Theta and Phi of cblin at an arbitrary point (v_perp, v_par).
  ! Series below u = v/(sqrt2 vt) = 0.1, where (erf - 3G)/(2 v^3) loses
  ! digits to cancellation (they agree to ~1e-13 at u = 0.1).
  !------------------------------------------------------------------
  subroutine theta_phi(vpe, vpa, vth, theta, phi)
    real(dp), intent(in)  :: vpe, vpa, vth
    real(dp), intent(out) :: theta, phi
    real(dp), parameter :: pi = 3.141592653589793238462643_dp
    real(dp) :: v, u, u2, c, e, de, g, s2v
    v   = sqrt(vpe**2 + vpa**2)
    s2v = sqrt(2.0_dp) * vth
    u   = v / s2v
    c   = 2.0_dp / sqrt(pi)
    if (u < 0.1_dp) then
      u2    = u*u
      theta = c*(1.0_dp/3 - u2/5 + u2**2/14 - u2**3/54 + u2**4/264 - u2**5/1560) / s2v
      phi   = c*(4.0_dp/15 - 4*u2/35 + 2*u2**2/63 - 2*u2**3/297 + u2**4/858) &
              / (2.0_dp * s2v**3)
    else
      e     = erf(u)
      de    = c * exp(-u*u)
      g     = (e - u*de) / (2.0_dp*u*u)            ! Chandrasekhar function
      theta = g / v
      phi   = (e - 3.0_dp*g) / (2.0_dp*v**3)
    end if
  end subroutine theta_phi

  !------------------------------------------------------------------
  ! ln B(x), B(x) = x/(e^x - 1) > 0, without overflow
  !------------------------------------------------------------------
  real(dp) function log_bern(x)
    real(dp), intent(in) :: x
    if (abs(x) < 1.0e-3_dp) then
      log_bern = log(1.0_dp - x/2 + x*x/12 - x**4/720)
    else if (x > 0.0_dp) then
      log_bern = log(x) - x - log(1.0_dp - exp(-x))
    else
      log_bern = log(-x) - log(1.0_dp - exp(x))
    end if
  end function log_bern

  real(dp) function lnM(i, j, kap)
    integer,  intent(in) :: i, j
    real(dp), intent(in) :: kap
    lnM = -kap * (vperp(i)**2 + vpar(j)**2)
  end function lnM

  !------------------------------------------------------------------
  ! Flux through the v_perp face between (i,j) and (i+1,j), for a
  ! Maxwellian background of thermal speed vth and mass ratio mrat:
  ! G = sum c(a,l) f(i+a, j+l), a = 0,1, l = -1..1
  !------------------------------------------------------------------
  subroutine perp_face(i, j, vth, mrat, c)
    integer,  intent(in)  :: i, j
    real(dp), intent(in)  :: vth, mrat
    real(dp), intent(out) :: c(0:1,-1:1)
    real(dp) :: kap, h, th, ph, dpp, dpa, s, lb, lm0, w(-1:1), ws(-1:1)
    integer  :: a
    c   = 0.0_dp
    kap = mrat / (2.0_dp * vth**2)
    h   = vperp(i+1) - vperp(i)
    call theta_phi(eface(i), vpar(j), vth, th, ph)
    dpp = th + vpar(j)**2 * ph
    dpa = -eface(i) * vpar(j) * ph
    lm0 = lnM(i, j, kap)
    s   = lm0 - lnM(i+1, j, kap)
    lb  = log_bern(s)
    ! diagonal (Chang-Cooper): (D/h) M~ (g+ - g-), M~ = B(s) M_-
    c(1,0) = c(1,0) + dpp/h * exp(lb + lm0 - lnM(i+1, j, kap))
    c(0,0) = c(0,0) - dpp/h * exp(lb)
    ! cross: D_pa M d(f/M)/dva at the two nodes of the face, averaged;
    ! Scharfetter-Gummel weights along va (see sg_weights)
    call dva_weights(j, w)
    do a = 0, 1
      call sg_weights(w, [lnM(i+a, max(j-1,1), kap), lnM(i+a, j, kap), &
                          lnM(i+a, min(j+1,npar), kap)], ws)
      c(a,:) = c(a,:) + 0.5_dp * dpa * ws
    end do
  end subroutine perp_face

  !------------------------------------------------------------------
  ! Two-point difference w (weights -1/H at node lo, +1/H at node hi)
  ! turned into the Scharfetter-Gummel approximation of M d(f/M)/dv:
  !   [B(-s) f_hi - B(s) f_lo] / H,   s = ln M_lo - ln M_hi,
  ! zero on M exactly, -> the central difference as s -> 0.  The former
  ! form M~ (g_hi - g_lo)/H, g = f/M, has weights exp(ln M~ - ln M_node)
  ! ~ e^|s| at high cell Peclet number; these grow only like |s|.
  !------------------------------------------------------------------
  subroutine sg_weights(w, lm, ws)
    real(dp), intent(in)  :: w(-1:1), lm(-1:1)
    real(dp), intent(out) :: ws(-1:1)
    integer  :: lo, hi, k
    real(dp) :: s
    lo = 0;  hi = 0
    do k = -1, 1
      if (w(k) < 0.0_dp) lo = k
      if (w(k) > 0.0_dp) hi = k
    end do
    s = lm(lo) - lm(hi)
    ws = 0.0_dp
    ws(hi) = w(hi) * exp(log_bern(-s))
    ws(lo) = w(lo) * exp(log_bern(s))
  end subroutine sg_weights

  !------------------------------------------------------------------
  ! Flux through the v_par face between (i,j) and (i,j+1):
  ! G = sum c(k,m) f(i+k, j+m), k = -1..1, m = 0,1
  !------------------------------------------------------------------
  subroutine par_face(i, j, vth, mrat, c)
    integer,  intent(in)  :: i, j
    real(dp), intent(in)  :: vth, mrat
    real(dp), intent(out) :: c(-1:1,0:1)
    real(dp) :: kap, vf, th, ph, daa, dpa, s, lb, lm0, w(-1:1), ws(-1:1)
    integer  :: m
    c   = 0.0_dp
    kap = mrat / (2.0_dp * vth**2)
    vf  = 0.5_dp * (vpar(j) + vpar(j+1))
    call theta_phi(vperp(i), vf, vth, th, ph)
    daa = th + vperp(i)**2 * ph
    dpa = -vperp(i) * vf * ph
    lm0 = lnM(i, j, kap)
    s   = lm0 - lnM(i, j+1, kap)
    lb  = log_bern(s)
    c(0,1) = c(0,1) + daa/dvpar * exp(lb + lm0 - lnM(i, j+1, kap))
    c(0,0) = c(0,0) - daa/dvpar * exp(lb)
    ! cross: D_pa M d(f/M)/dvp at the two nodes of the face, averaged
    ! (at i = 1 the mirror node -vp(1) is node 1 itself: same M)
    call dvp_weights(i, w)
    do m = 0, 1
      call sg_weights(w, [lnM(max(i-1,1), j+m, kap), lnM(i, j+m, kap), &
                          lnM(min(i+1,nperp), j+m, kap)], ws)
      c(:,m) = c(:,m) + 0.5_dp * dpa * ws
    end do
  end subroutine par_face

  !------------------------------------------------------------------
  ! Chang-Cooper flux A f + D df/dv across one face of width h, for local
  ! coefficients: G = cm f_- + cp f_+, with s = A h / D (exact for constant
  ! A and D).  D <= 0 (should not happen) falls back to upwind drag.
  !------------------------------------------------------------------
  subroutine cc_pair(d, a, h, cm, cp)
    real(dp), intent(in)  :: d, a, h
    real(dp), intent(out) :: cm, cp
    real(dp) :: s
    if (d > 0.0_dp) then
      s  = a * h / d
      cm = -d/h * exp(log_bern(s))
      cp =  d/h * exp(log_bern(-s))
    else if (a >= 0.0_dp) then
      cm = 0.0_dp;  cp = a
    else
      cm = a;       cp = 0.0_dp
    end if
  end subroutine cc_pair

  ! Nonlinear self-collisions: v_perp face between (i,j) and (i+1,j),
  ! coefficients averaged from the two nodes; same layout as perp_face
  subroutine perp_face_nl(i, j, c)
    integer,  intent(in)  :: i, j
    real(dp), intent(out) :: c(0:1,-1:1)
    real(dp) :: dpp, dpa, ap, w(-1:1)
    integer  :: a, l
    c   = 0.0_dp
    dpp = 0.5_dp * (n_dpp(i,j) + n_dpp(i+1,j))
    dpa = 0.5_dp * (n_dpa(i,j) + n_dpa(i+1,j))
    ap  = 0.5_dp * (n_ap(i,j)  + n_ap(i+1,j))
    call cc_pair(dpp, ap, vperp(i+1) - vperp(i), c(0,0), c(1,0))
    call dva_weights(j, w)
    do a = 0, 1
      do l = -1, 1
        c(a,l) = c(a,l) + 0.5_dp * dpa * w(l)
      end do
    end do
  end subroutine perp_face_nl

  ! Nonlinear self-collisions: v_par face between (i,j) and (i,j+1)
  subroutine par_face_nl(i, j, c)
    integer,  intent(in)  :: i, j
    real(dp), intent(out) :: c(-1:1,0:1)
    real(dp) :: daa, dpa, apa, w(-1:1)
    integer  :: k, m
    c   = 0.0_dp
    daa = 0.5_dp * (n_daa(i,j) + n_daa(i,j+1))
    dpa = 0.5_dp * (n_dpa(i,j) + n_dpa(i,j+1))
    apa = 0.5_dp * (n_aa(i,j)  + n_aa(i,j+1))
    call cc_pair(daa, apa, dvpar, c(0,0), c(0,1))
    call dvp_weights(i, w)
    do m = 0, 1
      do k = -1, 1
        c(k,m) = c(k,m) + 0.5_dp * dpa * w(k)
      end do
    end do
  end subroutine par_face_nl

  ! d/dvpar at node column j: central inside, one-sided at the walls
  subroutine dva_weights(j, w)
    integer,  intent(in)  :: j
    real(dp), intent(out) :: w(-1:1)
    w = 0.0_dp
    if (j == 1) then
      w(0) = -1.0_dp/dvpar;  w(1) = 1.0_dp/dvpar
    else if (j == npar) then
      w(-1) = -1.0_dp/dvpar; w(0) = 1.0_dp/dvpar
    else
      w(-1) = -0.5_dp/dvpar; w(1) = 0.5_dp/dvpar
    end if
  end subroutine dva_weights

  ! d/dvperp at node row i.  i = 1: difference across the axis between the
  ! mirror node -vp(1) and vp(2) (g is even in vp); D_pa ~ vp makes it
  ! vanish when vp(1) = 0.
  subroutine dvp_weights(i, w)
    integer,  intent(in)  :: i
    real(dp), intent(out) :: w(-1:1)
    w = 0.0_dp
    if (i == 1) then
      w(0) = -1.0_dp/(vperp(2) + vperp(1));  w(1) = 1.0_dp/(vperp(2) + vperp(1))
    else if (i == nperp) then
      w(-1) = -1.0_dp/(vperp(i) - vperp(i-1)); w(0) = 1.0_dp/(vperp(i) - vperp(i-1))
    else
      w(-1) = -1.0_dp/(vperp(i+1) - vperp(i-1)); w(1) = 1.0_dp/(vperp(i+1) - vperp(i-1))
    end if
  end subroutine dvp_weights

  !------------------------------------------------------------------
  ! Matrix values, on the fixed pattern, of the operator for a Maxwellian
  ! background (vth, mrat), without its gamma factor; or, with nl, of the
  ! nonlinear self-collision operator from the nodal n_* arrays (vth and
  ! mrat unused).  Rows are independent.
  !------------------------------------------------------------------
  subroutine op_values(vth, mrat, v, nl)
    real(dp), intent(in)  :: vth, mrat
    real(dp), intent(out) :: v(:)
    logical,  intent(in), optional :: nl
    integer  :: i, j, row, p, ip, jp, a, m
    real(dp) :: acc(-1:1,-1:1), cp(0:1,-1:1), ca(-1:1,0:1), ri, wj
    logical  :: use_nl
    use_nl = .false.
    if (present(nl)) use_nl = nl
    !$OMP PARALLEL DO SCHEDULE(STATIC) &
    !$OMP   PRIVATE(i, j, row, p, ip, jp, a, m, acc, cp, ca, ri, wj)
    do i = 1, nperp
      ri = 0.5_dp * (eface(i)**2 - eface(i-1)**2)       ! int vp dvp over the cell
      do j = 1, npar
        row = (i-1)*npar + j
        if (bc_row(i, j)) then                     ! Dirichlet: set in cfv_fill
          v(cfv_ia(row):cfv_ia(row+1)-1) = 0.0_dp
          cycle
        end if
        wj = dvpar
        if (j == 1 .or. j == npar) wj = 0.5_dp * dvpar
        acc = 0.0_dp
        if (i < nperp) then                         ! outer v_perp face
          if (use_nl) then
            call perp_face_nl(i, j, cp)
          else
            call perp_face(i, j, vth, mrat, cp)
          end if
          do a = 0, 1
            acc(a,:) = acc(a,:) + eface(i)/ri * cp(a,:)
          end do
        end if
        if (i > 1) then                             ! inner v_perp face
          if (use_nl) then
            call perp_face_nl(i-1, j, cp)
          else
            call perp_face(i-1, j, vth, mrat, cp)
          end if
          do a = 0, 1
            acc(a-1,:) = acc(a-1,:) - eface(i-1)/ri * cp(a,:)
          end do
        end if
        if (j < npar) then                          ! upper v_par face
          if (use_nl) then
            call par_face_nl(i, j, ca)
          else
            call par_face(i, j, vth, mrat, ca)
          end if
          do m = 0, 1
            acc(:,m) = acc(:,m) + ca(:,m)/wj
          end do
        end if
        if (j > 1) then                             ! lower v_par face
          if (use_nl) then
            call par_face_nl(i, j-1, ca)
          else
            call par_face(i, j-1, vth, mrat, ca)
          end if
          do m = 0, 1
            acc(:,m-1) = acc(:,m-1) - ca(:,m)/wj
          end do
        end if
        do p = cfv_ia(row), cfv_ia(row+1)-1
          ip = (cfv_ja(p)-1)/npar + 1
          jp = cfv_ja(p) - (ip-1)*npar
          v(p) = acc(ip-i, jp-j)
        end do
      end do
    end do
    !$OMP END PARALLEL DO
  end subroutine op_values

  !------------------------------------------------------------------
  ! Check that the Maxwellian of (vth, mrat) is a null vector of v, row by
  ! row (rows reaching a Dirichlet node excluded)
  !------------------------------------------------------------------
  subroutine check_null(v, vth, mrat, label)
    real(dp), intent(in) :: v(:), vth, mrat
    character(*), intent(in) :: label
    real(dp), allocatable :: fm(:)
    real(dp) :: kap, res, scl, rmax
    integer  :: i, j, row, p
    allocate(fm(nbig))
    kap = mrat / (2.0_dp * vth**2)
    do i = 1, nperp
      do j = 1, npar
        fm((i-1)*npar + j) = exp(lnM(i, j, kap))
      end do
    end do
    rmax = 0.0_dp
    do row = 1, nbig
      i = (row-1)/npar + 1
      j = row - (i-1)*npar
      if (bc_row(i, j)) cycle
      if (i_cc_wall == 0 .and. (i >= nperp-1 .or. j <= 2 .or. j >= npar-1)) cycle
      res = 0.0_dp;  scl = 0.0_dp
      do p = cfv_ia(row), cfv_ia(row+1)-1
        res = res + v(p) * fm(cfv_ja(p))
        scl = scl + abs(v(p) * fm(cfv_ja(p)))
      end do
      if (scl > 0.0_dp) rmax = max(rmax, abs(res)/scl)
    end do
    write(*,'(A,A,A,ES9.2)') '  Coulomb flux form: ', label, &
          ', max |L M| / sum |terms| on its Maxwellian = ', rmax
    deallocate(fm)
  end subroutine check_null

  !------------------------------------------------------------------
  ! Build the fixed CSR pattern (3x3 per row, sorted) and the per-species
  ! values; check that each species' Maxwellian is an exact null vector.
  !------------------------------------------------------------------
  subroutine cfv_build()
    integer  :: i, j, b, row, p, ip, jp
    real(dp) :: ri, wj
    character(len=16) :: lab

    if (built) return
    allocate(eface(0:nperp))
    eface(0) = 0.0_dp
    do i = 1, nperp-1
      eface(i) = 0.5_dp * (vperp(i) + vperp(i+1))
    end do
    eface(nperp) = vperp(nperp)

    allocate(vol(nbig), cfv_ia(nbig+1), idiag(nbig))
    ! pattern: every neighbour in the 3x3 block inside the grid; Dirichlet
    ! rows keep their diagonal only and the other rows drop wall columns
    cfv_ia(1) = 1
    do i = 1, nperp
      do j = 1, npar
        row = (i-1)*npar + j
        p = 0
        do ip = max(1,i-1), min(nperp,i+1)
          do jp = max(1,j-1), min(npar,j+1)
            if (keep(i, j, ip, jp)) p = p + 1
          end do
        end do
        cfv_ia(row+1) = cfv_ia(row) + p
      end do
    end do
    cfv_nnz = cfv_ia(nbig+1) - 1
    allocate(cfv_ja(cfv_nnz), val(cfv_nnz, nbulk))

    do i = 1, nperp
      ri = 0.5_dp * (eface(i)**2 - eface(i-1)**2)
      do j = 1, npar
        wj = dvpar
        if (j == 1 .or. j == npar) wj = 0.5_dp * dvpar
        row = (i-1)*npar + j
        vol(row) = 2.0_dp * 3.141592653589793238462643_dp * ri * wj
        p = cfv_ia(row)
        do ip = max(1,i-1), min(nperp,i+1)
          do jp = max(1,j-1), min(npar,j+1)
            if (.not. keep(i, j, ip, jp)) cycle
            cfv_ja(p) = (ip-1)*npar + jp
            if (ip == i .and. jp == j) idiag(row) = p
            p = p + 1
          end do
        end do
      end do
    end do
    built = .true.

    do b = 1, nbulk
      call op_values(vt(b), maonmb(b), val(:,b))
      write(lab,'(A,I0)') 'species ', b
      call check_null(val(:,b), vt(b), maonmb(b), trim(lab))
    end do
  end subroutine cfv_build

  ! is column (ip,jp) part of row (i,j)?
  logical function keep(i, j, ip, jp)
    integer, intent(in) :: i, j, ip, jp
    if (bc_row(i, j)) then
      keep = (ip == i .and. jp == j)
    else
      keep = .not. bc_row(ip, jp)
    end if
  end function keep

  !------------------------------------------------------------------
  ! Maxwellian self-collisions: thermal speed vth (vteff_t) and factor gam
  ! (gammaa).  The matrix is recomputed only when vth changes; the first
  ! one is checked against its Maxwellian.
  !------------------------------------------------------------------
  subroutine cfv_set_sc(vth, gam)
    real(dp), intent(in) :: vth, gam
    logical :: first
    call cfv_build()
    first = .not. allocated(val_sc)
    if (first) allocate(val_sc(cfv_nnz))
    if (vth /= vth_sc) then
      call op_values(vth, 1.0_dp, val_sc)
      vth_sc = vth
      if (first) call check_null(val_sc, vth, 1.0_dp, 'self-collisions')
    end if
    gam_sc = gam
  end subroutine cfv_set_sc

  !------------------------------------------------------------------
  ! Nonlinear self-collisions: nodal D (pp, pa, aa) and drag A (p, a) from
  ! main_nlterm, factor included.  The operator is conservative whatever
  ! the coefficients; its steady state is Maxwellian only as far as the
  ! numerical potentials make A = (m/T) D.v.
  !------------------------------------------------------------------
  subroutine cfv_set_nl(dpp, dpa, daa, ap, apa)
    real(dp), intent(in), dimension(nperp,npar) :: dpp, dpa, daa, ap, apa
    logical, save :: warned = .false.
    integer :: nneg
    call cfv_build()
    if (.not. allocated(val_nl)) then
      allocate(val_nl(cfv_nnz))
      allocate(n_dpp(nperp,npar), n_dpa(nperp,npar), n_daa(nperp,npar), &
               n_ap(nperp,npar), n_aa(nperp,npar))
    end if
    n_dpp = dpp;  n_dpa = dpa;  n_daa = daa;  n_ap = ap;  n_aa = apa
    call op_values(1.0_dp, 1.0_dp, val_nl, nl=.true.)
    nl_on = .true.
    nneg = count(dpp <= 0.0_dp) + count(daa <= 0.0_dp)
    if (nneg > 0 .and. .not. warned) then
      write(*,'(A,I0,A)') '  Nonlinear SC, flux form: ', nneg, &
            ' nodal D_pp or D_aa <= 0 (faces there use upwind drag); reported once'
      warned = .true.
    end if
  end subroutine cfv_set_nl

  !------------------------------------------------------------------
  ! Quasilinear RF operator (irf = -1) in flux form.  qlrfterm's rf**
  ! arrays are, expanded, the exact cylindrical divergence
  !   L f = (1/vp) d/dvp [vp G_p] + d/dva G_a ,  G = K T b (b . grad f),
  !   b = (vph - va, vp),  K = rfcte delta(va - vres),  T = |Theta|^2,
  ! a rank-one diffusion with no drag.  Face-centred differences of the
  ! cross term (as for the collisions) would not keep this operator
  ! dissipative, so it is discretised at the cell corners c = (e_i, j+1/2),
  ! i = 1..nperp-1, j = 1..npar-1 (symmetric scheme, Guenter et al. 2005):
  ! - corner gradient from its four nodes,
  !     G_p f = [(f(i+1,j) - f(i,j)) + (f(i+1,j+1) - f(i,j+1))] / (2 hp),
  !     G_a f = [rp (f(i,j+1) - f(i,j)) + rm (f(i+1,j+1) - f(i+1,j))] / om,
  !   rp = (e_i^2 - vp_i^2)/2, rm = (vp_i+1^2 - e_i^2)/2 the two halves of
  !   int vp dvp, om = e_i hp dva (= (rp + rm) dva);
  ! - corner flux G_c = K_c T_c b_c (b_c . G f), delivered to the four
  !   quarter faces that meet at c (v_perp faces: e_i dva/2 each; v_par
  !   faces: rp and rm).
  ! These weights give, exactly,
  !   sum_cells V g (L f) = - sum_c om_c K_c T_c (b.G g)(b.G f),
  ! so V L is symmetric and negative semidefinite: conservative and
  ! dissipative for any grid and any anisotropy.  T is evaluated on the
  ! v_perp faces, K is the mean of the Gaussian delta over [va_j, va_j+1]
  ! (so sum_c dva K_c = rfcte int delta exactly).  The walls carry no flux
  ! (corners are interior); the strip [0, vp_1] of the first cell's v_par
  ! faces is not covered (empty when vp_1 = 0; G_a ~ vp there anyway).
  !------------------------------------------------------------------
  subroutine cfv_set_rf()
    use shared_rf,    only: rfcte, vph, vres, delta_rf
    use mod_qlrfterm, only: rf_theta2
    real(dp), allocatable :: beta(:,:,:,:), wc(:,:), tf(:), kc(:)
    real(dp) :: hp, e, rp, rm, om, bp, ba, ri, wj, vr, acc(-1:1,-1:1)
    integer  :: i, j, ic, jc, a, l, a2, l2, row, p, ip, jp, nneg

    call cfv_build()
    if (.not. allocated(val_rf)) allocate(val_rf(cfv_nnz))
    allocate(beta(0:1,0:1,nperp-1,npar-1), wc(nperp-1,npar-1), tf(nperp-1), kc(npar-1))

    do ic = 1, nperp-1
      tf(ic) = rf_theta2(eface(ic))                ! serial: CBESJ
    end do
    do jc = 1, npar-1
      kc(jc) = rfcte * (erf((vpar(jc+1) - vres)/delta_rf) - erf((vpar(jc) - vres)/delta_rf)) &
               / (2.0_dp * dvpar)
    end do

    do jc = 1, npar-1
      do ic = 1, nperp-1
        hp = vperp(ic+1) - vperp(ic)
        e  = eface(ic)
        rp = 0.5_dp * (e**2 - vperp(ic)**2)
        rm = 0.5_dp * (vperp(ic+1)**2 - e**2)
        om = e * hp * dvpar
        bp = vph - 0.5_dp * (vpar(jc) + vpar(jc+1))
        ba = e
        beta(0,0,ic,jc) = -bp/(2.0_dp*hp) - ba*rp/om
        beta(1,0,ic,jc) =  bp/(2.0_dp*hp) - ba*rm/om
        beta(0,1,ic,jc) = -bp/(2.0_dp*hp) + ba*rp/om
        beta(1,1,ic,jc) =  bp/(2.0_dp*hp) + ba*rm/om
        wc(ic,jc) = om * kc(jc) * tf(ic)
      end do
    end do
    nneg = count(wc < 0.0_dp)
    if (nneg > 0) write(*,'(A,I0,A)') '  RF flux form: WARNING, ', nneg, &
          ' corners with K |Theta|^2 < 0 (operator not dissipative there)'

    !$OMP PARALLEL DO SCHEDULE(STATIC) &
    !$OMP   PRIVATE(i, j, ic, jc, a, l, a2, l2, row, p, ip, jp, ri, wj, vr, acc)
    do i = 1, nperp
      ri = 0.5_dp * (eface(i)**2 - eface(i-1)**2)
      do j = 1, npar
        row = (i-1)*npar + j
        if (bc_row(i, j)) then
          val_rf(cfv_ia(row):cfv_ia(row+1)-1) = 0.0_dp
          cycle
        end if
        wj = dvpar
        if (j == 1 .or. j == npar) wj = 0.5_dp * dvpar
        vr = ri * wj                                ! cell volume / 2 pi
        acc = 0.0_dp
        do ic = max(1, i-1), min(nperp-1, i)
          do jc = max(1, j-1), min(npar-1, j)
            a = i - ic;  l = j - jc                 ! this node's place in corner c
            do a2 = 0, 1
              do l2 = 0, 1
                acc(ic+a2-i, jc+l2-j) = acc(ic+a2-i, jc+l2-j) &
                    - wc(ic,jc) * beta(a,l,ic,jc) * beta(a2,l2,ic,jc) / vr
              end do
            end do
          end do
        end do
        do p = cfv_ia(row), cfv_ia(row+1)-1
          ip = (cfv_ja(p)-1)/npar + 1
          jp = cfv_ja(p) - (ip-1)*npar
          val_rf(p) = acc(ip-i, jp-j)
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    deallocate(beta, wc, tf, kc)
    call check_rf()
  end subroutine cfv_set_rf

  !------------------------------------------------------------------
  ! Start-of-run check of the RF matrix: V L symmetric, and (closed box)
  ! its columns sum to zero, i.e. sum V L f = 0 for any f
  !------------------------------------------------------------------
  subroutine check_rf()
    real(dp), allocatable :: colsum(:), colabs(:)
    real(dp) :: asym, scl, vrow, x, y
    integer  :: row, col, p, q, i, j, ic, jcol
    allocate(colsum(nbig), colabs(nbig))
    colsum = 0.0_dp;  colabs = 0.0_dp
    asym = 0.0_dp;  scl = 0.0_dp
    do row = 1, nbig
      i = (row-1)/npar + 1
      j = row - (i-1)*npar
      if (bc_row(i, j)) cycle
      vrow = vol(row)
      do p = cfv_ia(row), cfv_ia(row+1)-1
        col = cfv_ja(p)
        x = vrow * val_rf(p)
        colsum(col) = colsum(col) + x
        colabs(col) = colabs(col) + abs(x)
        scl = max(scl, abs(x))
        ! transposed entry (col, row)
        ic = (col-1)/npar + 1
        jcol = col - (ic-1)*npar
        if (bc_row(ic, jcol)) cycle
        y = 0.0_dp
        do q = cfv_ia(col), cfv_ia(col+1)-1
          if (cfv_ja(q) == row) then
            y = vol(col) * val_rf(q);  exit
          end if
        end do
        asym = max(asym, abs(x - y))
      end do
    end do
    if (scl > 0.0_dp) then
      write(*,'(A,ES9.2)') '  RF flux form: max |V L - (V L)^T| / max |V L| = ', asym/scl
      if (i_cc_wall == -1) write(*,'(A,ES9.2)') &
          '  RF flux form: max |sum_rows V L| / sum |terms| per column = ', &
          maxval(abs(colsum) / max(colabs, tiny(1.0_dp)))
    end if
    deallocate(colsum, colabs)
  end subroutine check_rf

  !------------------------------------------------------------------
  ! RF operator applied to f, through the scheme in use: what the
  ! diagnostics call instead of apply_operator(rf**)
  !------------------------------------------------------------------
  subroutine rf_apply(f, Lf)
    use mod_apply_operator, only: apply_operator
    use shared_FPterms
    real(dp), intent(in)  :: f(nbig)
    real(dp), intent(out) :: Lf(nbig)
    real(dp), allocatable :: rf00(:,:)
    if (i_cc == -1) then
      call csr_apply(val_rf, 1.0_dp, f, Lf)
    else
      allocate(rf00(nperp,npar))
      rf00 = 0.0_dp
      call apply_operator(rf20, rf02, rf11, rf10, rf01, rf00, f, Lf)
      deallocate(rf00)
    end if
  end subroutine rf_apply

  ! no nonlinear self-collision term this step (NBI fill-up from f = 0)
  subroutine cfv_nl_off()
    nl_on = .false.
  end subroutine cfv_nl_off

  !------------------------------------------------------------------
  ! aa = sum_b gammab(b) L_b (+ gammaa L_sc, + L_nl) - taum on the diagonal;
  ! Dirichlet rows = 1
  !------------------------------------------------------------------
  subroutine cfv_fill(aa, taum)
    real(dp), intent(out) :: aa(:)
    real(dp), intent(in)  :: taum
    integer :: p, b, row, i, j
    logical :: sc, rf
    sc = allocated(val_sc)
    rf = allocated(val_rf)
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(p, b)
    do p = 1, cfv_nnz
      aa(p) = 0.0_dp
      do b = 1, nbulk
        aa(p) = aa(p) + gammab(b) * val(p,b)
      end do
      if (sc) aa(p) = aa(p) + gam_sc * val_sc(p)
      if (nl_on) aa(p) = aa(p) + val_nl(p)
      if (rf) aa(p) = aa(p) + val_rf(p)
    end do
    !$OMP END PARALLEL DO
    do row = 1, nbig
      i = (row-1)/npar + 1
      j = row - (i-1)*npar
      if (bc_row(i, j)) then
        aa(idiag(row)) = 1.0_dp
      else
        aa(idiag(row)) = aa(idiag(row)) - taum
      end if
    end do
  end subroutine cfv_fill

  ! Lf = gam * V f (0 on constraint rows: their values are zero)
  subroutine csr_apply(v, gam, f, Lf)
    real(dp), intent(in)  :: v(:), gam, f(nbig)
    real(dp), intent(out) :: Lf(nbig)
    integer :: row, p
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(row, p)
    do row = 1, nbig
      Lf(row) = 0.0_dp
      do p = cfv_ia(row), cfv_ia(row+1)-1
        Lf(row) = Lf(row) + v(p) * f(cfv_ja(p))
      end do
      Lf(row) = gam * Lf(row)
    end do
    !$OMP END PARALLEL DO
  end subroutine csr_apply

  !------------------------------------------------------------------
  ! Lf = gammab(b) L_b f (0 on constraint rows)
  !------------------------------------------------------------------
  subroutine cfv_apply_species(b, f, Lf)
    integer,  intent(in)  :: b
    real(dp), intent(in)  :: f(nbig)
    real(dp), intent(out) :: Lf(nbig)
    call csr_apply(val(:,b), gammab(b), f, Lf)
  end subroutine cfv_apply_species

  !------------------------------------------------------------------
  ! NBI source: beam_source normalises it analytically (exact integral
  ! npart/taus), but the rate the scheme injects is its sum over the cells,
  ! off by the quadrature error (+4.6e-5 at 151^2, JET beam 1).  Rescale it
  ! so that sum V S = npart/taus exactly: then the closed-box steady state
  ! is n = npart on the cells.
  !------------------------------------------------------------------
  subroutine cfv_normalise_source()
    use shared_beam,   only: source, source_v, taus
    use shared_plasma, only: npart
    real(dp) :: rate, fac
    call cfv_build()
    rate = sum(vol * source_v)
    if (rate <= 0.0_dp) return
    fac = (npart / taus) / rate
    source_v = source_v * fac
    source   = source   * fac
    write(*,'(A,ES10.3)') '  NBI source rescaled to npart/taus on the finite-volume cells, factor - 1 = ', &
          fac - 1.0_dp
  end subroutine cfv_normalise_source

  ! density on the finite-volume cells: 2 pi sum R_i W_j f_ij
  real(dp) function cfv_density(f)
    real(dp), intent(in) :: f(nbig)
    cfv_density = sum(vol * f)
  end function cfv_density

  !------------------------------------------------------------------
  ! Collision operator of background species ib applied to f, through the
  ! scheme in use: what the diagnostics call instead of apply_operator.
  !------------------------------------------------------------------
  subroutine coll_apply_species(ib, f, Lf)
    use mod_apply_operator, only: apply_operator
    use shared_FPterms
    integer,  intent(in)  :: ib
    real(dp), intent(in)  :: f(nbig)
    real(dp), intent(out) :: Lf(nbig)
    if (i_cc == -1) then
      call cfv_build()
      call cfv_apply_species(ib, f, Lf)
    else
      call apply_operator(gammab(ib)*colin20_sp(:,:,ib), &
                          gammab(ib)*colin02_sp(:,:,ib), &
                          gammab(ib)*colin11_sp(:,:,ib), &
                          gammab(ib)*colin10_sp(:,:,ib), &
                          gammab(ib)*colin01_sp(:,:,ib), &
                          gammab(ib)*colin00_sp(:,:,ib), &
                          f, Lf)
    end if
  end subroutine coll_apply_species

  ! all background species together (the colin** sum)
  subroutine coll_apply_total(f, Lf)
    use mod_apply_operator, only: apply_operator
    use shared_FPterms
    real(dp), intent(in)  :: f(nbig)
    real(dp), intent(out) :: Lf(nbig)
    real(dp), allocatable :: Lb(:)
    integer :: ib
    if (i_cc == -1) then
      allocate(Lb(nbig))
      Lf = 0.0_dp
      do ib = 1, nbulk
        call coll_apply_species(ib, f, Lb)
        Lf = Lf + Lb
      end do
      deallocate(Lb)
    else
      call apply_operator(colin20, colin02, colin11, colin10, colin01, colin00, f, Lf)
    end if
  end subroutine coll_apply_total

  !------------------------------------------------------------------
  ! Self-collision operator applied to f, through the scheme in use (the
  ! one set by the last cfv_set_sc / self_coll_max call)
  !------------------------------------------------------------------
  subroutine sc_apply(f, Lf)
    use mod_apply_operator, only: apply_operator
    use shared_FPterms
    real(dp), intent(in)  :: f(nbig)
    real(dp), intent(out) :: Lf(nbig)
    if (i_cc == -1) then
      if (nl_on) then                              ! isc = -1
        call csr_apply(val_nl, 1.0_dp, f, Lf)
      else if (allocated(val_sc)) then             ! isc = 1, 2, 3
        call csr_apply(val_sc, gam_sc, f, Lf)
      else
        Lf = 0.0_dp
      end if
    else
      call apply_operator(sc20, sc02, sc11, sc10, sc01, sc00, f, Lf)
    end if
  end subroutine sc_apply

end module mod_coul_fv
