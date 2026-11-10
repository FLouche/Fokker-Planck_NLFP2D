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
!   - cross part: central differences of g = f/M, averaged onto the face,
!     scaled by the same face value of M (B(s) M_-).
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
!====================================================================
module mod_coul_fv

  use shared_grid
  use shared_plasma, only: nbulk, vt, maonmb, gammab

  implicit none
  private
  public :: cfv_build, cfv_fill, cfv_apply_species, cfv_density, bc_row, &
            coll_apply_species, coll_apply_total, cfv_set_sc, sc_apply
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
    real(dp) :: kap, h, th, ph, dpp, dpa, s, lb, lm0, w(-1:1)
    integer  :: a, l
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
    ! cross: D_pa M~ dg/dva, dg/dva averaged over the two nodes of the face
    call dva_weights(j, w)
    do a = 0, 1
      do l = -1, 1
        if (w(l) == 0.0_dp) cycle
        c(a,l) = c(a,l) + 0.5_dp * dpa * w(l) * exp(lb + lm0 - lnM(i+a, j+l, kap))
      end do
    end do
  end subroutine perp_face

  !------------------------------------------------------------------
  ! Flux through the v_par face between (i,j) and (i,j+1):
  ! G = sum c(k,m) f(i+k, j+m), k = -1..1, m = 0,1
  !------------------------------------------------------------------
  subroutine par_face(i, j, vth, mrat, c)
    integer,  intent(in)  :: i, j
    real(dp), intent(in)  :: vth, mrat
    real(dp), intent(out) :: c(-1:1,0:1)
    real(dp) :: kap, vf, th, ph, daa, dpa, s, lb, lm0, w(-1:1)
    integer  :: k, m
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
    call dvp_weights(i, w)
    do m = 0, 1
      do k = -1, 1
        if (w(k) == 0.0_dp) cycle
        c(k,m) = c(k,m) + 0.5_dp * dpa * w(k) * exp(lb + lm0 - lnM(i+k, j+m, kap))
      end do
    end do
  end subroutine par_face

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
  ! background (vth, mrat), without its gamma factor.  Rows are independent.
  !------------------------------------------------------------------
  subroutine op_values(vth, mrat, v)
    real(dp), intent(in)  :: vth, mrat
    real(dp), intent(out) :: v(:)
    integer  :: i, j, row, p, ip, jp, a, m
    real(dp) :: acc(-1:1,-1:1), cp(0:1,-1:1), ca(-1:1,0:1), ri, wj
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
          call perp_face(i, j, vth, mrat, cp)
          do a = 0, 1
            acc(a,:) = acc(a,:) + eface(i)/ri * cp(a,:)
          end do
        end if
        if (i > 1) then                             ! inner v_perp face
          call perp_face(i-1, j, vth, mrat, cp)
          do a = 0, 1
            acc(a-1,:) = acc(a-1,:) - eface(i-1)/ri * cp(a,:)
          end do
        end if
        if (j < npar) then                          ! upper v_par face
          call par_face(i, j, vth, mrat, ca)
          do m = 0, 1
            acc(:,m) = acc(:,m) + ca(:,m)/wj
          end do
        end if
        if (j > 1) then                             ! lower v_par face
          call par_face(i, j-1, vth, mrat, ca)
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
  ! aa = sum_b gammab(b) L_b (+ gammaa L_sc) - taum on the diagonal;
  ! Dirichlet rows = 1
  !------------------------------------------------------------------
  subroutine cfv_fill(aa, taum)
    real(dp), intent(out) :: aa(:)
    real(dp), intent(in)  :: taum
    integer :: p, b, row, i, j
    logical :: sc
    sc = allocated(val_sc)
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(p, b)
    do p = 1, cfv_nnz
      aa(p) = 0.0_dp
      do b = 1, nbulk
        aa(p) = aa(p) + gammab(b) * val(p,b)
      end do
      if (sc) aa(p) = aa(p) + gam_sc * val_sc(p)
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
      if (allocated(val_sc)) then
        call csr_apply(val_sc, gam_sc, f, Lf)
      else
        Lf = 0.0_dp
      end if
    else
      call apply_operator(sc20, sc02, sc11, sc10, sc01, sc00, f, Lf)
    end if
  end subroutine sc_apply

end module mod_coul_fv
