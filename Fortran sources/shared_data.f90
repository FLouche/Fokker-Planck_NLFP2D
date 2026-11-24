!
!   shared_data.f90
!
!   FPColl_2D
    ! vperp/vpar version
!   
!   Created by Fabrice Louche on 02/07/25.
!   Copyright 2025 LPP-ERM/KMS. All rights reserved.
!   
module shared_grid

implicit none

save

integer :: nperp, npar, nbig
double precision :: vperp_min,vperp_max,vpar_min,vpar_max
double precision :: vbound,dvleft,dvright,dvperp,dvpar
double precision :: dv2, dmu2

double precision, allocatable, dimension(:) :: vperp,vpar
double precision, allocatable, dimension(:,:) :: jacob

! phi-distance kernel for the non-linear self-collision term.
!
! This used to be the dense matrix sum_phi(nbig,nbig), which is
! block-Toeplitz with Toeplitz blocks and therefore holds only
! npar*nperp^2 distinct values in nperp^2*npar^2 slots (11.9 GiB against
! 61 MiB at 200x200).  Only the compressed kernel is stored now, plus its
! circulant transform used by the FFT matvec.  See mod_phi_kernel.
double precision, allocatable, dimension(:,:,:) :: phi_kern   ! (0:npar-1, nperp, nperp)
double precision, allocatable, dimension(:,:,:) :: phi_khat   ! (nperp, nperp, 0:phi_nf-1)
integer :: phi_m = 0      ! circulant length (power of two, >= 2*npar-2)
integer :: phi_nf = 0     ! independent frequencies, phi_m/2 + 1

integer :: ising,nsing
double precision :: p_grid = 2.d0  ! ising=+1: vperp(i) = vperp_min + vperp_max*((i-1/2)/(nperp-1))**p_grid
                                   ! p_grid=2 is the original quadratic grid; a larger p_grid puts
                                   ! more points in the thermal bulk at the cost of a coarser outer
                                   ! (tail) region, the spacing still vanishing at the axis for p_grid>1.
integer :: jmid,imid
integer :: i_upwind = 0   ! v_perp convection-diffusion scheme in fd_stencil_2d:
                          !  0: central Fornberg drag + Fornberg diffusion (default)
                          ! -1: Peclet-hybrid first-order upwind of the drag (B) term
                          !     in cells with |B|dv/D > 2; central diffusion
                          !  1: Patankar (1980) power-law scheme, applied only in
                          !     cells with |B|dv/D > 2 (Fornberg kept in the bulk)

integer :: n_snap = 0     ! VDF snapshots: write f every n_snap steps to
                          ! vdf_snap_<step>.txt (0 = off, the default).  Meant
                          ! for time-resolved moment diagnostics -- e.g. the
                          ! drift-frame anisotropy behind the sign of the
                          ! self-collision power split -- not for production
                          ! runs: a 151x151 snapshot is ~0.7 MB.

integer :: i_lu_reuse = -1          ! Time-step linear solves (pardiso_solver):
                                    ! -1: reuse the LU factors of the last
                                    !     factorisation as the preconditioner of a
                                    !     CGS iteration (PARDISO phase 23); refactorise
                                    !     only when CGS fails (default).
                                    !  0: refactorise every step (the former path).
double precision :: lu_reuse_tol = 1.d-10  ! CGS stopping tolerance (relative
                                    ! residual) for i_lu_reuse=-1; rounded to a
                                    ! power of ten, 10^-L, as PARDISO requires.

integer :: nstep_restart = 0  ! Time steps already done before this run: the
                              ! step count stored on the first line of xout.dat
                              ! by the run being restarted (iold=-1), 0 for a
                              ! fresh run or an older xout.dat that holds only
                              ! the time.  Snapshots are numbered
                              ! nstep_restart + itime_global, so a restart
                              ! continues the previous run's numbering instead
                              ! of overwriting its files.

end module shared_grid

!***************************************

module shared_plasma

implicit none

save

double precision, dimension(10) :: nb,maonmb,vt,gammab,t
double precision :: ne,npart
double precision :: gammaa
double precision, dimension(9) :: ab,zb,xb
!double precision :: vcr
double precision :: xpart,aa,za
double precision ::  vteff,spit,tauie
double precision :: core_frac   ! isc=3 Tn log-slope: fit cells with f > core_frac*max(f)
	   
integer nbulk
!integer isource
	   	
    end module shared_plasma

! ***************************************
    
module shared_beam

implicit none

save

double PRECISION :: beam_ekin,beam_angle, beam_angle_deg
double PRECISION ::  beam_v,beam_vperp,beam_vpar
double precision :: beam_dvperp,beam_dvpar
double precision :: taus, taum
double precision, allocatable, dimension(:,:) :: source
double precision, allocatable, dimension(:) :: source_v

integer isource


    end module shared_beam

!! ***************************************

    
module shared_timer
!
implicit none
!
save
!
! ntimes(k) and timestep(k): up to 3 sequential phases.
! Phase k runs ntimes(k) steps of timestep(k) seconds.
! ntimes(2) and ntimes(3) default to 0 (phase skipped).
! ntimes(1)=0 → steady-state run (phases 2 and 3 ignored).
integer, dimension(3) :: ntimes   = [0, 0, 0]
double precision, dimension(3) :: timestep = [0.0d0, 0.0d0, 0.0d0]
! Per-phase scalars used internally by the time solvers (set by main before each solver call).
integer :: ntimes_cur   = 0
double precision :: timestep_cur = 0.0d0
integer :: iold,icn,isc
integer :: i_ss_check  = 0       ! 0: disabled; -1: auto-stop when SS reached, using the
                                 !    Jacobian-weighted rolling-window test in mod_conv_diag
integer :: n_ss_window = 50      ! SS test evaluated every n_ss_window steps; rates are
                                 ! measured over that window (mod_conv_diag)
! Steady-state tolerances (mod_conv_diag).  All are rates PER SECOND measured over the
! window, so like the criterion they replace they are independent of the time step.
! eps and eps_tail are compared after division by the reference rate nu_ref, so these are
! dimensionless "fraction of the physical rate" thresholds.
double precision :: ss_tol_eps    = 1.0d-3  ! tolerance on eps/nu_ref      (bulk L2 rate)
double precision :: ss_tol_tail   = 1.0d-2  ! tolerance on eps_tail/nu_ref (tail-weighted)
double precision :: ss_tol_moment = 1.0d-3  ! tolerance on the density/flow/energy drifts
integer :: i_conv_shape = 0      ! convergence criterion (mod_conv_diag):
                                 !  0: amplitude -- eps, eps_tail and the density/flow/energy
                                 !     drifts.  Correct when f reaches a true steady state.
                                 ! -1: SHAPE -- f is normalised to unit norm before
                                 !     differencing, so a uniformly draining solution reads as
                                 !     converged once its shape stops changing.  Use for
                                 !     sourceless runs, where particles absorbed at the
                                 !     Dirichlet boundaries put a floor under the amplitude
                                 !     rate that eps can never fall below.  The moment tests
                                 !     then use the INTENSIVE moments (u_par, E/n); the
                                 !     density drift is skipped, being pure amplitude.
integer :: istart = 1            ! TD initial condition: 0=zero(beam only) 1=Stix 2=SS no-SC 3=SS Maxw-SC
integer :: iplot_pow = -1        ! -1: write power vs time files; 0: skip
integer :: iplot_mom = 0         ! -1: write momentum vs time files; 0: skip
integer :: notxt     = 0         ! .txt output: 0=write all; 1=suppress all; 2=keep only the minimal test/diagnostic set
double precision, dimension(:,:), allocatable :: fstix
character(len=64) :: casename = ''
character(len=8)  :: species_tag = ''  ! multispecies: tag of the species whose grid is
                                       ! loaded (mod_species, select_species).  Empty for
                                       ! the heated species, so its files keep the v3.2
                                       ! names; 'e', 'i1', 'i2', ... for the backgrounds,
                                       ! which outfile appends to the case name.
!
contains

  ! Returns 'stem-casename.ext' when casename is non-empty,
  ! or 'stem.ext' unchanged when casename = ''.
  function outfile(name) result(fname)
    character(len=*), intent(in) :: name
    character(len=256) :: fname
    integer :: idot
    logical :: keep
    character(len=64) :: stem
    character(len=80) :: suffix
    ! .txt output control via notxt: 0 = write all; 1 = suppress all; 2 = keep
    ! only the minimal test/diagnostic set below.  Suppressed .txt files are
    ! routed to the OS null device.  (.dat restart/kernel files are unaffected.)
    if (notxt /= 0) then
      idot = index(name, '.', back=.true.)
      if (idot > 0 .and. name(idot:len_trim(name)) == '.txt') then
        keep = .false.
        if (notxt == 2) then
          stem = name(1:idot-1)
          keep =      trim(stem) == 'Ekin_perp'          &
                 .or. trim(stem) == 'Ekin_perp_at_vpar0' &
                 .or. trim(stem) == 'fout_at_vpar0'      &
                 .or. trim(stem) == 'Teff_vs_time'       &
                 .or. trim(stem) == 'anisotropy_vs_time'
        end if
        if (.not. keep) then
          fname = 'NUL'
          return
        end if
      end if
    end if
    ! Suffix = casename, plus '_<species_tag>' for a background species
    ! (or the tag alone when casename is empty).
    suffix = casename
    if (len_trim(species_tag) > 0) then
      if (len_trim(suffix) > 0) then
        suffix = trim(suffix) // '_' // trim(species_tag)
      else
        suffix = species_tag
      end if
    end if
    if (len_trim(suffix) == 0) then
      fname = trim(name)
    else
      idot = index(name, '.', back=.true.)
      if (idot > 0) then
        fname = name(1:idot-1) // '-' // trim(suffix) // name(idot:len_trim(name))
      else
        fname = trim(name) // '-' // trim(suffix)
      end if
    end if
  end function outfile

    end module shared_timer

!! ***************************************


module shared_FPterms

implicit none

save

!! ----- COLLISIONS

!  --> Total
double precision , dimension(:,:), allocatable :: colin00,colin10,colin01
double precision , dimension(:,:), allocatable :: colin11,colin20,colin02

!  --> Per species
double precision , dimension(:,:,:), allocatable :: colin10_sp,colin01_sp,colin00_sp
double precision , dimension(:,:,:), allocatable :: colin11_sp,colin20_sp,colin02_sp

! --> Self-collisions

double precision , dimension(:,:), allocatable :: sc10,sc01,sc00
double precision , dimension(:,:), allocatable :: sc11,sc20,sc02

! --> RF terms

double precision , dimension(:,:), allocatable :: rf10,rf01
double precision , dimension(:,:), allocatable :: rf11,rf20,rf02


    end module shared_FPterms
    
!! ***************************************

    
module shared_multi
!
! Multispecies extension (docs/Multispecies_design.pdf).  Defaults reproduce v3.2.
!
implicit none

save

integer :: i_multi = 0              ! 0: single species, the backgrounds are fixed
                                    !    Maxwellians (v3.2).
                                    ! -1: multispecies.  Every background species gets
                                    !    its own grid and distribution (mod_species).
                                    !    Stage 1: built and written, still frozen.
double precision :: tau_E = 0.d0    ! Energy confinement time [s] of the loss term
                                    !   L_E f = (1/2 tau_E) div(v f),
                                    ! applied to every species, and only when there is a
                                    ! source (irf=-1 or isource=-1).  0 = no loss term.

! Background-species grids, indexed like the bulk species (1 = electrons, 2.. = ions).
! The box is |v| <= vbox_b * sqrt(T_b/m_b) at the initial temperature, in both v_perp
! and v_par; the v_perp grid type (ising, p_grid, nsing) is the heated species' one.
integer,          dimension(10) :: nperp_b = 0    ! 0 -> nperp of the heated species
integer,          dimension(10) :: npar_b  = 0    ! 0 -> npar  of the heated species
double precision, dimension(10) :: vbox_b  = 0.d0 ! 0 -> vbox_default
double precision, parameter     :: vbox_default = 8.d0  ! f(wall)/f(0) = exp(-32) = 1.3e-14

! Legendre potentials (mod_legendre_pot), stage 2.
integer :: lmax_leg = 16     ! highest harmonic of the heated species' potentials
integer :: lmax_bulk = 4     ! coupled loop: highest Legendre degree of the bulk species'
                             ! potentials.  L = 0 (isotropic) breaks momentum conservation:
                             ! a drifting species then collides with a non-drifting copy
                             ! of itself (e-e drift damped at ~0.3 nu_ee).  L >= 1 needed.
integer :: nr_leg   = 0      ! radial intervals of the expansion (0: auto, step dvpar/4)
integer :: nxi_leg  = 0      ! Gauss points in cos(theta) (0: auto, max(64, 4 lmax_leg + 8))
integer :: i_ms_check = 0    ! -1: stage-2 checks (Legendre pair coefficients vs cblin at
                             !     start; vs the phi kernel at the end of timefp_7pt_nl)

! Stage 3: coupled evolution of all species (mod_timefp_ms).
integer :: i_coupled = 0     !  0: background species frozen (stages 1-2, v3.2 physics)
                             ! -1: every species evolves on its own grid; the collision
                             !     operator of each pair is built from the Legendre
                             !     potentials of the background species, lagged by one
                             !     step (heated species' self-collisions: main_nlterm).
                             !     Needs i_multi = -1, isc = -1, iold = 0.
double precision :: ta_init = 0.d0  ! initial temperature [eV] of the heated species in
                             ! the coupled run (Maxwellian); 0 = the namelist/Stix start.
double precision :: ua_init = 0.d0  ! initial parallel drift [m/s] of the heated species (with ta_init)
integer :: n_picard = 0      ! coupled loop: Picard iterations per step.  0: potentials
                             ! lagged (f^n); k > 0: potentials, coefficients and
                             ! corrections re-evaluated k times at
                             ! theta f^(n+1) + (1-theta) f^n (energy exchange
                             ! antisymmetric at convergence).
integer :: i_ms_cons = 2     ! conservation of the coupled step (implicit moment projection
                             ! onto the responses to f, div(v f), df/dv_par; TimeFP_ms.f90):
                             !  0: none (raw discrete pair operators)
                             !  1: particle number of every species
                             !  2: 1 + total energy and total parallel momentum (self-pair
                             !     rates removed, inter-species defects shared n- and
                             !     rho-weighted): all conserved to round-off
integer :: i_bulk_frozen = 0 ! coupled loop, stage-4 verification:
                             !  0: every species evolves
                             ! -1: the bulk species keep their initial Maxwellians (no
                             !     step, no L_E, no projection): the coupled loop then
                             !     solves the v3.2 problem with its own operators
                             !     (Legendre a <- b coefficients instead of cblin)

! f*f term of the nonlinear self-collisions (main_nlterm), c00 = -coef laplacian(phi).
integer :: i_sc00 = 0        !  0: finite-difference laplacian(phi) after the axis
                             !     regularisation (v3.2)
                             ! -1: exact laplacian(phi) = f, i.e. c00 = -coef f

end module shared_multi

! ***************************************

module shared_rf

implicit none

save

integer :: nharm

double precision :: b0,frek,kpar,delta_rf
double precision :: omega,vph,omc,vres,rfcte

complex*16 :: eplus,emin,kperp
complex*16 :: ceplus,cemin,eta

integer irf


end module shared_rf

! ***************************************

!
!module energy_td
!
!implicit none
!
!save
!
!double precision :: densi, convold, correct
!double precision :: ener,enpa,enpe,enav
!double precision :: temp
!double precision :: coll,pcolnl
!double precision, dimension(10) :: colpd
!double precision :: pbal
!
!double precision odensi,oenpa,oenpe,oenav,oener, &
!            otemp,opassrf,otime,oconv1,oconv2
!			
!logical ex,ex2
!
!end module energy_td
!
!! ***************************************
!