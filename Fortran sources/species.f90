!***********************************************************************
!*  mod_species -- species containers for the multispecies extension   *
!*                                                                     *
!*  Every species (the heated one, s = 1, and each background species, *
!*  s = 1 + ib, ib = 1 electrons, 2.. ions) has its own velocity grid    *
!*  and distribution function.  The rest of the code works on ONE grid, *
!*  held in shared_grid; select_species(s) loads species s into it, so  *
!*  fd_stencil_2d, the moments, the diagnostics, ... run unchanged on   *
!*  whichever species is selected.  The grid arrays are copied exactly, *
!*  so switching back reproduces the heated species' grid to the bit.   *
!*                                                                     *
!*  Stage 1 (docs/Multispecies_design.pdf): the containers are built    *
!*  and the initial Maxwellians written.  With i_coupled = 0 the        *
!*  backgrounds stay frozen (the heated species sees them through       *
!*  cblin); with i_coupled = -1 every species evolves (TimeFP_ms.f90).  *
!*                                                                     *
!*  20/11/2026                                                          *
!***********************************************************************

module mod_species

  use shared_grid
  implicit none
  private

  public :: species_t, spec, nspec, species_setup, select_species, &
            store_grid, load_grid, species_current, bulk_grid

  integer, parameter :: dp = kind(1.0d0)

  ! Electron to proton mass ratio, as in consts (vt(1), maonmb(1)).
  double precision, parameter, public :: mp_on_me = 1.8362d3

  type species_t
    character(len=8) :: tag = ''          ! output tag ('' heated, 'e', 'i1', ...)
    integer          :: kind = 0          ! 0 heated, 1 electrons, 2 background ion
    integer          :: ib = 0            ! bulk index (0 for the heated species)
    double precision :: aa = 0.d0         ! mass in proton masses
    double precision :: za = 0.d0         ! charge number (magnitude)
    double precision :: dens = 0.d0       ! density [m^-3]
    double precision :: t0 = 0.d0         ! initial temperature [eV]
    double precision :: vth0 = 0.d0       ! sqrt(T0/m) [m/s]
    ! grid (mirror of shared_grid)
    integer          :: nperp = 0, npar = 0, nbig = 0, ising = 0, nsing = 0
    integer          :: jmid = 0, imid = 0
    double precision :: vperp_min = 0.d0, vperp_max = 0.d0
    double precision :: vpar_min = 0.d0, vpar_max = 0.d0
    double precision :: vbound = 0.d0, p_grid = 2.d0
    double precision :: dvperp = 0.d0, dvpar = 0.d0, dvleft = 0.d0, dvright = 0.d0
    double precision :: dv2 = 0.d0, dmu2 = 0.d0
    double precision, allocatable :: vperp(:), vpar(:), jacob(:,:)
    double precision, allocatable :: f(:,:)        ! distribution [m^-6 s^3]
  end type species_t

  integer :: nspec = 0
  integer :: species_current = 0
  type(species_t), allocatable, target :: spec(:)

contains

!-----------------------------------------------------------------------
! Copy the grid now in shared_grid into species s.
!-----------------------------------------------------------------------
  subroutine store_grid(s)
    integer, intent(in) :: s
    associate (p => spec(s))
      p%nperp = nperp;  p%npar = npar;  p%nbig = nbig
      p%ising = ising;  p%nsing = nsing;  p%jmid = jmid;  p%imid = imid
      p%vperp_min = vperp_min;  p%vperp_max = vperp_max
      p%vpar_min  = vpar_min;   p%vpar_max  = vpar_max
      p%vbound = vbound;  p%p_grid = p_grid
      p%dvperp = dvperp;  p%dvpar = dvpar;  p%dvleft = dvleft;  p%dvright = dvright
      p%dv2 = dv2;  p%dmu2 = dmu2
      if (allocated(p%vperp)) deallocate(p%vperp, p%vpar, p%jacob)
      allocate(p%vperp(nperp), p%vpar(npar), p%jacob(nperp,npar))
      p%vperp = vperp;  p%vpar = vpar;  p%jacob = jacob
    end associate
  end subroutine store_grid

!-----------------------------------------------------------------------
! Load the grid of species s into shared_grid.
!-----------------------------------------------------------------------
  subroutine load_grid(s)
    integer, intent(in) :: s
    associate (p => spec(s))
      if (allocated(vperp)) then
        if (size(vperp) /= p%nperp .or. size(vpar) /= p%npar) &
          deallocate(vperp, vpar, jacob)
      end if
      if (.not. allocated(vperp)) allocate(vperp(p%nperp), vpar(p%npar), jacob(p%nperp,p%npar))
      nperp = p%nperp;  npar = p%npar;  nbig = p%nbig
      ising = p%ising;  nsing = p%nsing;  jmid = p%jmid;  imid = p%imid
      vperp_min = p%vperp_min;  vperp_max = p%vperp_max
      vpar_min  = p%vpar_min;   vpar_max  = p%vpar_max
      vbound = p%vbound;  p_grid = p%p_grid
      dvperp = p%dvperp;  dvpar = p%dvpar;  dvleft = p%dvleft;  dvright = p%dvright
      dv2 = p%dv2;  dmu2 = p%dmu2
      vperp = p%vperp;  vpar = p%vpar;  jacob = p%jacob
    end associate
  end subroutine load_grid

!-----------------------------------------------------------------------
! Make species s the current one: its grid in shared_grid, its tag on the
! output file names.
!-----------------------------------------------------------------------
  subroutine select_species(s)
    use shared_timer, only: species_tag
    use mod_fd_stencil_2d, only: fd_wcache_sync
    integer, intent(in) :: s
    call load_grid(s)
    ! Serial rebuild of the stencil-weight cache for this grid (thread safety).
    call fd_wcache_sync(nperp, npar, vperp, dvpar)
    species_tag = spec(s)%tag
    species_current = s
  end subroutine select_species

!-----------------------------------------------------------------------
! Grid of a bulk species s: the heated species' v_perp grid type (power
! law, nsing, vbound) rescaled so that v_perp <= vhalf, and a uniform
! v_par grid on [vcenter - vhalf, vcenter + vhalf].  Built in shared_grid
! and stored in spec(s); the caller selects a species afterwards.  Used
! at setup (vcenter = 0) and by the regridding of the coupled loop.
!-----------------------------------------------------------------------
  subroutine bulk_grid(s, np_b, nq_b, vhalf, vcenter)
    use mod_grid, only: make_grid
    integer,          intent(in) :: s, np_b, nq_b
    double precision, intent(in) :: vhalf, vcenter
    double precision :: scale
    scale = vhalf / spec(1)%vperp_max
    nperp = np_b;  npar = nq_b
    ising = spec(1)%ising;  nsing = spec(1)%nsing;  p_grid = spec(1)%p_grid
    vperp_min = spec(1)%vperp_min * scale
    vperp_max = vhalf
    vbound    = spec(1)%vbound * scale
    vpar_min  = vcenter - vhalf
    vpar_max  = vcenter + vhalf
    if (allocated(vperp)) deallocate(vperp, vpar, jacob)
    allocate(vperp(nperp), vpar(npar), jacob(nperp,npar))
    call make_grid
    call store_grid(s)
  end subroutine bulk_grid

!-----------------------------------------------------------------------
! Build every species: the heated one from the grid already made by main
! (make_grid) and the Stix Maxwellian (consts); each background from its
! own box, |v| <= vbox_b * vth0, with an initial Maxwellian at t(ib) and
! density nb(ib).  Called once, after consts.  Leaves the heated species
! selected, its grid restored exactly.
!-----------------------------------------------------------------------
  subroutine species_setup

    use shared_plasma
    use shared_multi
    use shared_timer, only: fstix, outfile
    use mod_ncint

    integer          :: s, ib, i, j, np_b, nq_b
    double precision :: vbox, dens_num, cinit, arg, vth
    double precision, allocatable :: fint(:,:)
    ! the heated species' grid size (default for the background grids)
    integer          :: nperp_a, npar_a

    double precision, parameter :: twopi = 6.283185307179586d0

    nspec = 1 + nbulk
    if (allocated(spec)) deallocate(spec)
    allocate(spec(nspec))

    !--- s = 1: heated species, grid and Stix Maxwellian as built by main/consts
    spec(1)%tag  = ''
    spec(1)%kind = 0
    spec(1)%ib   = 0
    spec(1)%aa   = aa
    spec(1)%za   = za
    spec(1)%dens = npart
    spec(1)%vth0 = vteff
    spec(1)%t0   = aa * (vteff / 9.79d3)**2
    call store_grid(1)
    allocate(spec(1)%f(nperp,npar))
    spec(1)%f = fstix

    nperp_a = nperp;  npar_a = npar

    !--- s = 1+ib: background species
    do ib = 1, nbulk
      s = 1 + ib
      spec(s)%ib = ib
      spec(s)%t0 = t(ib)
      spec(s)%dens = nb(ib)
      if (ib == 1) then
        spec(s)%tag  = 'e'
        spec(s)%kind = 1
        spec(s)%aa   = 1.d0 / mp_on_me
        spec(s)%za   = 1.d0
        spec(s)%vth0 = vt(1)              ! 9.79d3*sqrt(1836.2*Te), as in consts
      else
        write(spec(s)%tag, '(A,I0)') 'i', ib - 1
        spec(s)%kind = 2
        spec(s)%aa   = ab(ib-1)
        spec(s)%za   = zb(ib-1)
        spec(s)%vth0 = vt(ib)             ! 9.79d3*sqrt(T/A)
      end if

      ! Grid: the heated species' v_perp grid type, rescaled to the box.
      vbox = vbox_b(ib);  if (vbox <= 0.d0) vbox = vbox_default
      np_b = nperp_b(ib); if (np_b <= 0) np_b = nperp_a
      nq_b = npar_b(ib);  if (nq_b <= 0) nq_b = npar_a
      call bulk_grid(s, np_b, nq_b, vbox * spec(s)%vth0, 0.d0)

      ! Initial Maxwellian at T0, normalised to the density on this grid
      ! (as fstix in consts).
      allocate(spec(s)%f(nperp,npar), fint(nperp,npar))
      vth   = spec(s)%vth0
      cinit = spec(s)%dens / (twopi * vth**2)**1.5d0
      do i = 1, nperp
        do j = 1, npar
          arg = (vperp(i)**2 + vpar(j)**2) / (2.d0 * vth**2)
          spec(s)%f(i,j) = cinit * exp(-arg)
          fint(i,j) = spec(s)%f(i,j) * jacob(i,j)
        end do
      end do
      call ncint_2d(fint, dens_num)
      if (dens_num > 0.d0) spec(s)%f = spec(s)%f * spec(s)%dens / dens_num
      deallocate(fint)
    end do

    !--- Summary and initial distributions
    write(*,'(/,A)') ' MULTISPECIES (i_multi = -1): species and grids'
    write(*,'(A)')   '   s  tag kind     A          Z    n [m^-3]    T0 [eV]   vth0 [m/s]  '// &
                     'vmax/vth0 nperp x npar  vperp_max   vpar_min    vpar_max [m/s]'
    do s = 1, nspec
      associate (p => spec(s))
        write(*,'(I4,2X,A4,I4,ES12.4,F6.1,4ES12.4,I6,A,I4,3ES12.4)') s, p%tag, p%kind, p%aa, p%za, &
             p%dens, p%t0, p%vth0, p%vperp_max/p%vth0, p%nperp, ' x', p%npar, &
             p%vperp_max, p%vpar_min, p%vpar_max
      end associate
    end do
    write(*,'(A)') '   Initial distributions written (finit-*.txt).'
    if (i_coupled /= 0 .and. i_bulk_frozen == -1) then
      write(*,'(A)')   '   Coupled loop (i_coupled = -1) with i_bulk_frozen = -1: only the heated'
      write(*,'(A,/)') '   species evolves; the bulk species keep their initial Maxwellians.'
    else if (i_coupled /= 0) then
      write(*,'(A,/)') '   Coupled loop (i_coupled = -1): every species evolves on its own grid.'
    else
      write(*,'(A)')   '   i_coupled = 0: the background distributions stay frozen; the heated'
      write(*,'(A,/)') '   species sees them through cblin (v3.2 physics).'
    end if

    do s = 1, nspec
      call select_species(s)
      open(40, file=trim(outfile('finit.txt')), status='unknown')
      do i = 1, nperp
        do j = 1, npar
          write(40,*) vperp(i), vpar(j), spec(s)%f(i,j)
        end do
      end do
      close(40)
    end do

    call select_species(1)

  end subroutine species_setup

end module mod_species
