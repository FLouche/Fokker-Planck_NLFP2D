!***********************************************************************
!*  mod_ms_checks -- stage-2 checks of the Legendre potentials and of   *
!*  the any-pair coefficients (docs/Multispecies_design.pdf).           *
!*  Namelist i_ms_check = -1 (with i_multi = -1).                       *
!*                                                                     *
!*  ms_check_cblin  (after species_setup, every pair (a <- b)):         *
!*     the initial f of every species is a Maxwellian, so the pair      *
!*     coefficients from the Legendre potentials of f_b must equal      *
!*     cblin(vth_b, m_a/m_b) x gamma0 (Z_a/A_a)^2 Z_b^2 n_b on the grid *
!*     of a (lnLambda = 1 on both sides).                               *
!*                                                                     *
!*  ms_check_kernel (end of timefp_7pt_nl, phi kernel still alive):     *
!*     psi of the final (RF-anisotropic) f from the Legendre expansion  *
!*     for increasing lmax against compute_psi (FFT kernel), and the    *
!*     self-collision coefficients against main_nlterm.                 *
!*                                                                     *
!*  Results on stdout and in ms_check.txt.                              *
!*                                                                     *
!*  20/11/2026                                                          *
!***********************************************************************

module mod_ms_checks

  use mod_legendre_pot
  implicit none
  private
  public :: ms_check_cblin, ms_check_kernel

  integer, parameter :: uchk = 496

contains

!-----------------------------------------------------------------------
! Relative errors of six coefficient arrays over a mask:
!   diffusion (20, 02, 11) normalised by max(|ref20|,|ref02|) on the mask,
!   friction  (10, 01)     by max(|ref10|,|ref01|),
!   c00                    by max|ref00|.
!-----------------------------------------------------------------------
  subroutine coef_errors(n1, n2, mask, a20,a02,a11,a10,a01,a00, &
                         r20,r02,r11,r10,r01,r00, err)
    integer, intent(in) :: n1, n2
    logical, intent(in) :: mask(n1,n2)
    double precision, dimension(n1,n2), intent(in) :: a20,a02,a11,a10,a01,a00
    double precision, dimension(n1,n2), intent(in) :: r20,r02,r11,r10,r01,r00
    double precision, intent(out) :: err(6)
    double precision :: sd, sf, s0
    sd = max(maxval(abs(r20), mask), maxval(abs(r02), mask), tiny(1.d0))
    sf = max(maxval(abs(r10), mask), maxval(abs(r01), mask), tiny(1.d0))
    s0 = max(maxval(abs(r00), mask), tiny(1.d0))
    err(1) = maxval(abs(a20 - r20), mask) / sd
    err(2) = maxval(abs(a02 - r02), mask) / sd
    err(3) = maxval(abs(a11 - r11), mask) / sd
    err(4) = maxval(abs(a10 - r10), mask) / sf
    err(5) = maxval(abs(a01 - r01), mask) / sf
    err(6) = maxval(abs(a00 - r00), mask) / s0
  end subroutine coef_errors

!-----------------------------------------------------------------------
  subroutine ms_check_cblin
    use shared_grid
    use shared_multi, only: lmax_leg, nr_leg, nxi_leg
    use mod_species
    external cblin

    integer :: s, t, i, j
    double precision :: mu, coef, gref, nleg, vth_t
    double precision :: c20,c02,c11,c10,c01,c00, err(6), errb(6)
    double precision, allocatable, dimension(:,:) :: a20,a02,a11,a10,a01,a00
    double precision, allocatable, dimension(:,:) :: r20,r02,r11,r10,r01,r00
    logical, allocatable :: mall(:,:), mbulk(:,:)
    type(legpot_t) :: pot
    character(len=4) :: tg_t, tg_s

    open(uchk, file='ms_check.txt', status='unknown')
    write(*,'(/,A)') ' STAGE-2 CHECK 1: pair coefficients from the Legendre potentials of'
    write(*,'(A)')   '   the initial Maxwellians against cblin (lnLambda = 1).  Errors: max|diff|'
    write(*,'(A)')   '   over the grid of a, normalised by the diffusion (20,02,11), friction'
    write(*,'(A)')   '   (10,01) and c00 scales; "bulk" = |v| <= 4 vth_a.'
    write(uchk,'(A)') '# Stage-2 check 1: Legendre pair coefficients vs cblin (Maxwellian sources)'
    write(uchk,'(A)') '# a  b  mu  n_leg/n_b  then err(20,02,11,10,01,00) all grid, then bulk |v|<=4vth_a'

    do s = 1, nspec
      call select_species(s)
      call legpot_build(pot, spec(s)%f, 0, nr_leg, nxi_leg)
      nleg = legpot_density(pot) / spec(s)%dens

      do t = 1, nspec
        call select_species(t)
        allocate(a20(nperp,npar), a02(nperp,npar), a11(nperp,npar), a10(nperp,npar), &
                 a01(nperp,npar), a00(nperp,npar), r20(nperp,npar), r02(nperp,npar), &
                 r11(nperp,npar), r10(nperp,npar), r01(nperp,npar), r00(nperp,npar), &
                 mall(nperp,npar), mbulk(nperp,npar))
        mu   = spec(t)%aa / spec(s)%aa
        coef = pair_coef_factor(spec(t)%za, spec(t)%aa, spec(s)%za, 1.d0)
        call pair_coefs(pot, mu, coef, a20, a02, a11, a10, a01, a00)

        gref  = 2.390775d-1 * (spec(t)%za / spec(t)%aa)**2 * spec(s)%za**2 * spec(s)%dens
        vth_t = spec(t)%vth0
        do j = 1, npar
          do i = 1, nperp
            call cblin(i, j, spec(s)%vth0, mu, c20, c02, c11, c10, c01, c00)
            r20(i,j) = gref*c20;  r02(i,j) = gref*c02;  r11(i,j) = gref*c11
            r10(i,j) = gref*c10;  r01(i,j) = gref*c01;  r00(i,j) = gref*c00
            mall(i,j)  = .true.
            mbulk(i,j) = vperp(i)**2 + vpar(j)**2 <= (4.d0 * vth_t)**2
          end do
        end do
        call coef_errors(nperp, npar, mall,  a20,a02,a11,a10,a01,a00, r20,r02,r11,r10,r01,r00, err)
        call coef_errors(nperp, npar, mbulk, a20,a02,a11,a10,a01,a00, r20,r02,r11,r10,r01,r00, errb)

        tg_t = spec(t)%tag;  if (t == 1) tg_t = 'a'
        tg_s = spec(s)%tag;  if (s == 1) tg_s = 'a'
        if (s == 1 .and. t == 1) write(*,'(A)') &
          '    a <- b        mu    n_leg/n_b |  all:  20       02       11       10       01       00'// &
          '     | bulk: 20       02       11       10       01       00'
        write(*,'(2X,A4,A,A4,ES10.3,F10.6,A,6ES9.2,A,6ES9.2)') tg_t, ' <- ', tg_s, mu, nleg, ' |', err, ' |', errb
        write(uchk,'(2A6,ES12.4,F12.8,12ES11.3)') tg_t, tg_s, mu, nleg, err, errb

        deallocate(a20,a02,a11,a10,a01,a00,r20,r02,r11,r10,r01,r00,mall,mbulk)
      end do
    end do
    call legpot_free(pot)
    call select_species(1)
    close(uchk)
    write(*,*)
  end subroutine ms_check_cblin

!-----------------------------------------------------------------------
! On the heated species' grid (selected), at the end of timefp_7pt_nl.
!-----------------------------------------------------------------------
  subroutine ms_check_kernel(xvec, teff)
    use shared_grid
    use shared_plasma
    use shared_multi, only: lmax_leg, nr_leg, nxi_leg
    use func_index
    use nlterm
    use coulomb_log_mod

    double precision, intent(in) :: xvec(nbig), teff

    integer :: i, j, il, l, nl
    integer, parameter :: llist(9) = (/ 0, 1, 2, 4, 8, 16, 24, 32, 48 /)
    double precision :: lnaa, coef, err(6), errt(6), epsi, epsit, pmax, fmax, ta_ev
    double precision, allocatable, dimension(:,:) :: f, psik, psil, &
         s20,s02,s11,s10,s01,s00, a20,a02,a11,a10,a01,a00
    logical, allocatable :: mall(:,:), mtail(:,:)
    type(legpot_t) :: pot
    type(legval_t) :: r

    allocate(f(nperp,npar), psik(nperp,npar), psil(nperp,npar), &
             s20(nperp,npar), s02(nperp,npar), s11(nperp,npar), s10(nperp,npar), &
             s01(nperp,npar), s00(nperp,npar), a20(nperp,npar), a02(nperp,npar), &
             a11(nperp,npar), a10(nperp,npar), a01(nperp,npar), a00(nperp,npar), &
             mall(nperp,npar), mtail(nperp,npar))
    do i = 1, nperp
      do j = 1, npar
        f(i,j) = xvec(index_mat(i,j))
      end do
    end do
    fmax = maxval(f)
    mall  = .true.
    mtail = f > 1.d-6 * fmax

    call compute_psi(f, psik)
    pmax = maxval(abs(psik))

    open(uchk, file='ms_check.txt', status='old', position='append')
    write(*,'(/,A)') ' STAGE-2 CHECK 2: psi of the final f, Legendre expansion vs FFT kernel'
    write(*,'(A)')   '   (relative to max|psi|; "f>1e-6" = where f > 1e-6 max f)'
    write(*,'(A)')   '    lmax    all grid      f>1e-6'
    write(uchk,'(A)') '# Stage-2 check 2: psi Legendre(lmax) vs kernel: lmax err_all err_f>1e-6'

    do il = 1, size(llist)
      l = llist(il)
      if (l > lmax_leg) exit
      call legpot_build(pot, f, l, nr_leg, nxi_leg)
      !$OMP PARALLEL DO PRIVATE(i, j, r)
      do j = 1, npar
        do i = 1, nperp
          call legpot_eval(pot, vperp(i), vpar(j), r)
          psil(i,j) = r%psi
        end do
      end do
      !$OMP END PARALLEL DO
      epsi  = maxval(abs(psil - psik)) / pmax
      epsit = maxval(abs(psil - psik), mtail) / pmax
      write(*,'(I8,2ES12.3)') l, epsi, epsit
      write(uchk,'(I6,2ES12.4)') l, epsi, epsit
      nl = l
    end do

    ! Coefficients: Legendre (lmax_leg) vs main_nlterm (kernel + finite differences)
    if (nl /= lmax_leg) then
      call legpot_build(pot, f, lmax_leg, nr_leg, nxi_leg)
      do j = 1, npar
        do i = 1, nperp
          call legpot_eval(pot, vperp(i), vpar(j), r)
          psil(i,j) = r%psi
        end do
      end do
    end if
    call main_nlterm(xvec, teff, s00, s10, s01, s20, s11, s02)
    ta_ev = teff * 1.d3
    call coulomb_log_ab(za, aa, ta_ev, npart, za, aa, ta_ev, npart, lnaa)
    coef = pair_coef_factor(za, aa, za, lnaa)
    call pair_coefs(pot, 1.d0, coef, a20, a02, a11, a10, a01, a00)
    ! main_nlterm's sc11 multiplies d2f/dvpe dvpa like c11
    call coef_errors(nperp, npar, mall,  a20,a02,a11,a10,a01,a00, s20,s02,s11,s10,s01,s00, err)
    call coef_errors(nperp, npar, mtail, a20,a02,a11,a10,a01,a00, s20,s02,s11,s10,s01,s00, errt)
    write(*,'(A,I0,A)') '   self-collision coefficients, Legendre (lmax=', lmax_leg, &
                        ') vs main_nlterm:'
    write(*,'(A)') '               20       02       11       10       01       00'
    write(*,'(A,6ES9.2)') '   all     ', err
    write(*,'(A,6ES9.2)') '   f>1e-6  ', errt
    write(uchk,'(A)') '# self coefficients vs main_nlterm: err(20,02,11,10,01,00) all, then f>1e-6'
    write(uchk,'(12ES11.3)') err, errt
    ! c00 = -coef mu laplacian(phi) = -coef f exactly: where does each side depart?
    s20 = -coef * f
    block
      integer :: loc(2)
      double precision :: sc0
      integer :: loc2
      integer, parameter :: ilist(4) = (/ 2, 3, 5, 9 /)
      sc0 = maxval(abs(s20))
      loc = maxloc(abs(s00 - s20))
      write(*,'(A,ES9.2,A,2I5,A,2ES11.3,A)') '   c00 vs -coef*f: main_nlterm ', &
           maxval(abs(s00 - s20)) / sc0, ' (max at i,j =', loc, ', v =', vperp(loc(1)), vpar(loc(2)), ')'
      loc = maxloc(abs(a00 - s20))
      write(*,'(A,ES9.2,A,2I5,A,2ES11.3,A)') '                   Legendre    ', &
           maxval(abs(a00 - s20)) / sc0, ' (max at i,j =', loc, ', v =', vperp(loc(1)), vpar(loc(2)), ')'
      write(*,'(A,4ES10.2)') '   main_nlterm c00 error off the axis, i >= 2, 3, 5, 9: ', &
           (maxval(abs(s00(ilist(loc2):,:) - s20(ilist(loc2):,:))) / sc0, loc2 = 1, 4)
      write(uchk,'(A,2ES11.3)') '# c00 vs -coef*f (grid f): main_nlterm, Legendre', &
           maxval(abs(s00 - s20)) / sc0, maxval(abs(a00 - s20)) / sc0
    end block

    ! Anisotropic part alone.  psi and the coefficients are dominated by the
    ! isotropic part of f, so the errors above say little about the l > 0
    ! terms.  Both methods are linear in f: apply each to f and to its
    ! isotropic part f_0(v) sampled on the grid, and compare the differences,
    ! in which the common discretisation error of the l = 0 part cancels.
    block
      type(legpot_t) :: pot0
      double precision, allocatable, dimension(:,:) :: fiso, psiki, psi0, &
           t20,t02,t11,t10,t01,t00, b20,b02,b11,b10,b01,b00
      double precision, allocatable :: xiso(:)
      double precision :: dk, ek(7), sz(7)
      integer :: k
      allocate(fiso(nperp,npar), psiki(nperp,npar), psi0(nperp,npar), xiso(nbig), &
               t20(nperp,npar), t02(nperp,npar), t11(nperp,npar), t10(nperp,npar), &
               t01(nperp,npar), t00(nperp,npar), b20(nperp,npar), b02(nperp,npar), &
               b11(nperp,npar), b10(nperp,npar), b01(nperp,npar), b00(nperp,npar))
      call legpot_build(pot0, f, 0, nr_leg, nxi_leg)
      do j = 1, npar
        do i = 1, nperp
          call legpot_eval(pot0, vperp(i), vpar(j), r)
          fiso(i,j) = r%fb
          xiso(index_mat(i,j)) = fiso(i,j)
        end do
      end do
      ! f_0 sampled on the grid is cut by the box (the shells of the RF tail
      ! leave it at |v_par| > vpar_max), so it is no longer isotropic: both
      ! methods must be applied to this same grid function, with all l.
      call legpot_build(pot0, fiso, lmax_leg, nr_leg, nxi_leg)
      do j = 1, npar
        do i = 1, nperp
          call legpot_eval(pot0, vperp(i), vpar(j), r)
          psi0(i,j) = r%psi
        end do
      end do
      call compute_psi(fiso, psiki)
      ! s** of f were overwritten for the c00 test above: recompute
      call main_nlterm(xvec, teff, s00, s10, s01, s20, s11, s02)
      call main_nlterm(xiso, teff, t00, t10, t01, t20, t11, t02)
      call pair_coefs(pot0, 1.d0, coef, b20, b02, b11, b10, b01, b00)

      dk  = maxval(abs(psik - psiki), mtail)
      ek(1) = maxval(abs((psil - psi0) - (psik - psiki)), mtail) / max(dk, tiny(1.d0))
      sz(1) = dk / pmax
      call aniso(a20, b20, s20, t20, ek(2), sz(2))
      call aniso(a02, b02, s02, t02, ek(3), sz(3))
      call aniso(a11, b11, s11, t11, ek(4), sz(4))
      call aniso(a10, b10, s10, t10, ek(5), sz(5))
      ! c00: against the exact -coef (f - f_iso) (main_nlterm's c00 is noisy on
      ! the axis); c01 is zero by construction for mu = 1 and is left out.
      t20 = -coef * f;  t02 = -coef * fiso
      call aniso(a00, b00, t20, t02, ek(6), sz(6))
      write(*,'(A)') '   anisotropic part (f - f_0), where f > 1e-6 max f:'
      write(*,'(A)') '                       psi      20       02       11       10     00(exact)'
      write(*,'(A,6ES9.2)') '   size / max|full|  ', sz(1:6)
      write(*,'(A,6ES9.2)') '   Legendre error    ', ek(1:6)
      write(*,'(A)') '   (error relative to the kernel''s anisotropic part)'
      write(uchk,'(A)') '# anisotropic part: size/max|full| (psi,20,02,11,10,00), then error/size'
      write(uchk,'(12ES11.3)') sz(1:6), ek(1:6)
      call legpot_free(pot0)
      deallocate(fiso, psiki, psi0, xiso, t20,t02,t11,t10,t01,t00, b20,b02,b11,b10,b01,b00)
    end block

    close(uchk)
    write(*,*)

    call legpot_free(pot)
    deallocate(f, psik, psil, s20,s02,s11,s10,s01,s00, a20,a02,a11,a10,a01,a00, mall, mtail)

  contains

    subroutine aniso(aL, a0, sK, sK0, e, sz_)
      double precision, dimension(nperp,npar), intent(in) :: aL, a0, sK, sK0
      double precision, intent(out) :: e, sz_
      double precision :: d
      d   = maxval(abs(sK - sK0), mtail)
      e   = maxval(abs((aL - a0) - (sK - sK0)), mtail) / max(d, tiny(1.d0))
      sz_ = d / max(maxval(abs(sK), mtail), tiny(1.d0))
    end subroutine aniso
  end subroutine ms_check_kernel

end module mod_ms_checks
