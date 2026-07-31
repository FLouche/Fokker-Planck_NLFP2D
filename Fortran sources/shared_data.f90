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

double precision, allocatable, dimension(:,:) :: sum_phi

integer :: ising,nsing
integer :: jmid,imid
integer :: i_upwind = 0   ! v_perp convection-diffusion scheme in fd_stencil_2d:
                          !  0: central Fornberg drag + Fornberg diffusion (default)
                          ! -1: Peclet-hybrid first-order upwind of the drag (B) term
                          !     in cells with |B|dv/D > 2; central diffusion
                          !  1: Patankar (1980) power-law scheme, applied only in
                          !     cells with |B|dv/D > 2 (Fornberg kept in the bulk)

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
integer :: idiag     = 0         ! 0: power balance only; -1: all balances (density+momentum+power)
integer :: notxt     = 0         ! .txt output: 0=write all; 1=suppress all; 2=keep only the minimal test/diagnostic set
double precision, dimension(:,:), allocatable :: fstix
character(len=64) :: casename = ''
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
    if (len_trim(casename) == 0) then
      fname = trim(name)
    else
      idot = index(name, '.', back=.true.)
      if (idot > 0) then
        fname = name(1:idot-1) // '-' // trim(casename) // name(idot:len_trim(name))
      else
        fname = trim(name) // '-' // trim(casename)
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