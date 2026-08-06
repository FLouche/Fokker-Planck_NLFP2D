program FP_Coll_2D

! ====================================================
!  Non-linear Fokker-Planck code
!    - 2D in vperp/vpar
!    - steady-state solution OR 
!    - time-dependant solver (Cranck-Nicholson method)
!    - collisions with Maxwellian background
!    - self-collisions accounted for:
!           _ either Maxwellian background (time-independant)
!           _ or general background (time-dependant)
!    - finite differences
!    - homogeneous grid
!    - beam 
!
! ====================================================
!
    ! Version 2.5 - 30th July 2026 (FL)
    
    ! Convergence in time assessed with L2 norm:  we added a Jacobian-weighted convergence 
    !      diagnostics module (mod_conv_diag) 
    
    ! Version 2.4 - 28th July 2026 (FL)
    
    !  Patankar power-law convection–diffusion scheme, activated by i_upwind = 1
    
    ! Version 2.3 - 3rd July 2026 (FL)
    
    ! Parameter i_upwind introduced to alleviate sign oscillations
    ! of fout at vperp max
    
    ! Version 2.2 - 23 June 2026 (FL)
!
!   New definiton of isc=3 


!     Version 2.1  - 03/06/2026 (FL)
!
!    Various corrections + new definition of the anisotropy factor
!
!     Version 2 - 29/05/2026 (FL)
!
!  Self-collision term (non-linear) consistently considers the variations
!   of the Coulomb logarithm
!
!   Version 1.7 - 29 May 2026 (FL):

!   new module to build linear FP terms -> for varying Coulomb log vs time
!       the treatment of varying Coulomb log is accounted for
!         in for the options isc=0, isc=1, and isc=2
!   isc=2: Maxwellian SC background at varying temperature (starts at Tstix)
    


!    
!    Fabrice Louche
!
!     MAIN PROGRAM
!     ------------
!
! Initialisation
! --------------

use shared_grid
use shared_plasma
use shared_beam
use shared_timer
use shared_rf
use mod_qlrfterm

use shared_FPterms

use mod_grid
use mod_beam

use assemble_FP_lin
use coulomb_log_mod
!
use mod_linear
use mod_timefp_7pt
use mod_timefp_7pt_nl
!
use mod_anal
!
use func_index
!
use mod_ncint
use time_comps_mod        ! explicit interface for time_energy(..., teff=)

implicit none

external cblin, consts, self_coll_max


! Elements of the Fokker-Planck linearized collision operator

double precision :: c20,c02,c10,c01,c00,c11

! Factors in front of the terms of the equation (total)

double precision , dimension(:,:), allocatable :: all00,all10,all01
double precision , dimension(:,:), allocatable :: all11,all20,all02

double precision, dimension (:), allocatable :: fin

double precision, dimension (:,:), allocatable :: fout


double precision time1
double precision :: teff_sc, lnaa_sc
double precision, parameter :: gamma0_sc = 2.390775d-1

! do loops indexes

integer :: ib,iv,imu,ix

! Other variables

!double PRECISION :: dens_tmp

! Math constants

double precision :: pi,twopi

common/mathcons/pi,twopi

data pi/3.141592653589793238462643d0/

double precision start_time,end_time

!external derf


!====================================================================
!
! Data's : Input Variables
! ------
!
! Definition of some variables:
!  icn: choice of time-differencing scheme:
!           +1: Crank-Nicholson
!            0: midpoint leap-frog (explicit)
!           -1: implicit scheme

!  beam_ekin: beam kinetic energy in keV
!
!  isc: treatment of self-collisions:
!           -1: Consistent self-collisions (non-linear operator)
!            0: no self-collisions
!           +1: Maxwellian background at fixed Tstix (Stix solution without RF)
!           +2: Maxwellian background at varying Teff (energy-weighted; starts at Tstix)
!           +3: Maxwellian SC background at the density-characteristic (cold-bulk)
!          temperature Tn, from a log-slope fit of ln(f) vs v^2 over the thermal
!          core (starts at Tstix). core_frac (namelist, default 3.8d-3) sets the
!          core threshold f > core_frac*max(f); 3.8d-3 calibrates isc=3 to isc=-1.
!  ising: homogeneity of the grid in vperp:
!            0: homogeneous grid
!           -1: inhomogeneous grid made of two domains (vperp<vbound and vperp>vbound) with different meshings
!                      for vperp <= vbound: nsing points (increase the density for small vperp)
!           +1: Quadratic spacing for higher resolution near vperp=0 (or vperp_min)
!               
!
! Convergence for time-dependent simulation:
!
!1. With i_ss_check=-1 a line like
!     [conv] t= 2.5E-02 s  eps= 3.1E-01 /s  epsT= 8.4E-01 /s  eps/nu= 3.1E-03 ...
!   appears on stdout once every n_ss_window steps (mod_conv_diag).
!2. When all five criteria are met (eps/nu < ss_tol_eps, epsT/nu < ss_tol_tail, and the
!   three moment drifts < ss_tol_moment) the run stops early and prints a CONVERGED line.
!3. After early exit, all output files (fout.txt, xout.dat, *_vs_time.txt) should be complete and the end-of-run plots
!should still be produced.
!
!The test is evaluated every n_ss_window steps; each rate is divided by the physical
!duration that window spans, so the tolerances are independent of timestep. If they are too tight
!(run never converges) or too loose (stops too early), adjust them together with n_ss_window; a
!wider window is more immune to short-term fluctuations and proportionally cheaper, but coarsens
!the time resolution of the diagnostic and delays the first verdict (which comes at step 2*n_ss_window).
!
! Initial solution for time-dependent solver:
!
!  - istart = 0 -> empty solution, f = 0 everywhere; only possible when isource = -1. The code should include a test at the start to guarantee that isource = -1
!  - istart = 1 ->  Stix's solution as initial solution 
!  - istart = 2 -> initial solution is the steady-state solution of the linear time independent code, computed without the self-collisions
!  - istart = 3 -> initial solution is the steady-state solution of the linear time independent code, computed with a Maxwellian background for the self-collisions

namelist /INPUT/ casename, &
                 nperp,npar,vperp_min,vperp_max,vpar_min,vpar_max,&
                ising,nsing,vbound,&
               nbulk,t,aa,ab,za,zb,ne,xpart,xb, &
                isource,beam_ekin,beam_angle_deg, & 
                beam_dvperp, beam_dvpar, taus, &
                irf,eplus,emin,kperp, &
                kpar,frek,delta_RF,b0,nharm, &
                icn, ntimes, timestep, iold, istart, isc, core_frac, &
                i_ss_check, n_ss_window, &
                ss_tol_eps, ss_tol_tail, ss_tol_moment, i_conv_shape, &
                iplot_pow, iplot_mom, idiag, notxt, &
                i_upwind


!write(*,*) 'Read namelist'

! Default for the isc=3 Tn log-slope core fraction (overridable via namelist).
! 3.8d-3 calibrates isc=3 to the rigorous isc=-1 reference (JET RF case5:
! Teff ~43.8 keV, grid-independent); see SC_models_grid_convergence report.
!core_frac = 3.8d-3

read(5,INPUT)

! isc=-1, isc=2 and isc=3 require a time-dependent run
if ((isc == -1 .or. isc == 2 .or. isc == 3) .and. ntimes(1) == 0) then
    write(*,'(A,I0,A)') 'ERROR: isc=', isc, ' requires a time-dependent run (ntimes > 0).'
    write(*,*) 'Steady-state solver cannot be used with a time-varying self-collision operator.'
    stop
endif

! Coherence checks for istart (only relevant for a fresh TD run)
if (ntimes(1) /= 0 .and. iold /= -1) then
    if (istart == 0 .and. isource /= -1) then
        write(*,*) 'ERROR: istart=0 (zero initial condition) requires isource=-1 (beam source).'
        write(*,*) 'Without a source, starting from f=0 gives a trivial zero solution.'
        stop
    endif
    if (istart < 0 .or. istart > 3) then
        write(*,'(A,I0,A)') 'ERROR: istart=', istart, ' is not valid. Use 0, 1, 2 or 3.'
        stop
    endif
endif

twopi=2.d0*pi

call cpu_time(start_time)

!====================================================================
!
! Building grid
! -------------

allocate(vperp(nperp),vpar(npar),jacob(nperp,npar))

call make_grid

call consts

!====================================================================
!
! Coefficients in the linear part of the collision operator 
! ---------------------------------------------------------

allocate(colin00_sp(nperp,npar,nbulk))
allocate(colin10_sp,colin01_sp,colin20_sp,colin02_sp,colin11_sp,mold=colin00_sp)
!
allocate(colin00(nperp,npar))
allocate(colin10,colin01,colin20,colin02,colin11,mold=colin00)

!

vperp_loop: do iv = 1,nperp
    
    vpar_loop: do imu = 1,npar
    
        background_loop: do ib=1,nbulk 
            
            call cblin(iv,imu,vt(ib),maonmb(ib),c20,c02,c11,c10,c01,c00)
            
            colin20_sp(iv,imu,ib) = c20
            colin02_sp(iv,imu,ib) = c02
            colin11_sp(iv,imu,ib) = c11
            colin10_sp(iv,imu,ib) = c10
            colin01_sp(iv,imu,ib) = c01
            colin00_sp(iv,imu,ib) = c00
                        
        enddo background_loop        
                
            
    enddo vpar_loop
    
enddo vperp_loop


!!====================================================================
!!
!!  Beam source
!
allocate(source(nperp,npar),source_v(nbig))

if (isource == -1) then

    call beam_source
    
    open(40,file=TRIM(outfile('beam.txt')), status='unknown')


do iv=1,nperp
        do imu=1,npar
        write(40,*) vperp(iv), vpar(imu), source(iv,imu)
        enddo
enddo

close(40)
    
else
        
    source = 0.d0
    source_v = 0.d0
        
endif
!
!!====================================================================
!
! RF TERM
! -------

if (irf == -1) then
    
    allocate(rf01,rf10,rf20,rf11,rf02,mold=colin00)
    
    call qlrfterm
    
endif


!!====================================================================
!!
! We allocate the factors in front of each derivative 

allocate(all00(nperp,npar))
allocate(all10,all01,all11,all20,all02,mold=all00)

call assemble_FP_terms(all00,all10,all01,all20,all11,all02)

!
!!====================================================================
!

! --> If self-collisions are taken into account, we consider two options:
!     1) we assume a Maxwellian background with Stix's effective temperature
!         --> only for steady-state solution OR isc = +1
!     2) we compute the time-dependant self-collision operator assuming
!        slowing-down on the pitch-angle averaged vdf
!         --> only in the TD module if isc = -1 (self-collisions neglected otherwise)
    
    if(isc /= 0) then
        allocate(sc00(nperp,npar))
        allocate(sc20,sc02,sc11,sc10,sc01,mold=sc00)
        if (isc == 1 .OR. isc == 2 .OR. isc == 3) then
           ! Initial SC term: Maxwellian background at Tstix (same starting point for isc=1,2,3)
           teff_sc = aa * (vteff / 9.79d3)**2     ! convert vteff [m/s] to T [eV]
           call coulomb_log_ab(za, aa, teff_sc, npart, za, aa, teff_sc, npart, lnaa_sc)
           gammaa = gamma0_sc * lnaa_sc * (za/aa)**2 * npart * za**2
           call self_coll_max(vteff, gammaa, sc20, sc02, sc11, sc10, sc01, sc00)

             all20 = all20+sc20
             all02 = all02+sc02
             all10 = all10+sc10
             all01 = all01+sc01
             all11 = all11+sc11
             all00 = all00+sc00
             if (isc == 1) write(*,*) 'Self-collisions: Maxwellian background at fixed Tstix'
             if (isc == 2) write(*,*) 'Self-collisions: Maxwellian background at varying Teff (initial: Tstix)'
             if (isc == 3) write(*,*) 'Self-collisions: Maxwellian background at varying Tn (initial: Tstix)'
        !else
        !    allocate(sum_phi(nbig,nbig))
        !    write(*,*) 'Starting distance evaluation'
        !    call cpu_time(start_time)
        !     call distance_v_gauss_legendre
        !    !call distance_v_elliptic_E
        !     write(*,*) 'Distance evaluation one'
        !     call cpu_time(end_time)
        !     write(*,*) 'Evaluation time: ',end_time-start_time,'seconds'
         endif
    endif    

allocate(fout(nperp,npar))
!
!!====================================================================
!! Steady-state case...
!! --------------------
!    
steady_state: if(ntimes(1) == 0) then
!
	write(*,*) ' '
    write(*,*) 'STEADY-STATE SOLUTION OF THE LINEAR FP EQUATION'
	write(*,*) '***********************************************'
	write(*,*) ' '

    call FP_steady_state(all20,all02,all11,all10,all01,all00,fout)
    
    
!
else steady_state

!====================================================================
!            
!  TIME-DEPENDANT SOLVER
!  ---------------------
    

!  --> Initial distribution construction 
!      ---------------------------------
    
    allocate(fin(nbig))
    fin = 0.d0
    
    if(iold == -1) then ! restart from previously stored solution

        open(40,file=TRIM(outfile('xout.dat')),status='old')
        read(40,*) time1
        do ix=1,nbig
            read(40,*) fin(ix)
        enddo

    else

        time1=0.d0

        select case (istart)

        case (0)
            ! f = 0 everywhere — beam-driven case only (isource=-1 guaranteed above)
            write(*,*) 'Initial condition: f = 0 (beam-driven).'
            fin = 0.d0

        case (1)
            ! Stix Maxwellian (sourceless) or zero (beam-driven) — original behaviour
            write(*,*) 'Initial condition: Stix Maxwellian (sourceless) or zero (beam).'
            do iv = 1,nperp
                do imu = 1,npar
                    ix = index_mat(iv,imu)
                    if (isource == 0) fin(ix) = fstix(iv,imu)
                enddo
            enddo

        case (2)
            ! Steady-state of linear code without self-collisions
            write(*,*) 'Initial condition: computing SS solution (no SC)...'
            if (isc == 1 .OR. isc == 2) then
                ! all** currently includes Maxwellian SC; strip it back to Coulomb+RF only
                block
                    double precision, dimension(nperp,npar) :: &
                        c00_ns,c10_ns,c01_ns,c11_ns,c20_ns,c02_ns
                    if (irf == -1) then
                        c20_ns=colin20+rf20; c02_ns=colin02+rf02; c11_ns=colin11+rf11
                        c10_ns=colin10+rf10; c01_ns=colin01+rf01; c00_ns=colin00
                    else
                        c20_ns=colin20; c02_ns=colin02; c11_ns=colin11
                        c10_ns=colin10; c01_ns=colin01; c00_ns=colin00
                    end if
                    call FP_steady_state(c20_ns,c02_ns,c11_ns,c10_ns,c01_ns,c00_ns,fout)
                end block
            else
                ! isc=0 or isc=-1: all** already excludes SC
                call FP_steady_state(all20,all02,all11,all10,all01,all00,fout)
            end if
            do iv = 1,nperp
                do imu = 1,npar
                    ix = index_mat(iv,imu)
                    fin(ix) = fout(iv,imu)
                enddo
            enddo
            write(*,*) 'Initial condition: SS (no SC) done.'

        case (3)
            ! Steady-state of linear code with Maxwellian SC background
            write(*,*) 'Initial condition: computing SS solution (Maxwellian SC)...'
            block
                double precision, dimension(nperp,npar) :: &
                    c00_sc,c10_sc,c01_sc,c11_sc,c20_sc,c02_sc
                double precision, dimension(nperp,npar) :: &
                    sc00t,sc10t,sc01t,sc11t,sc20t,sc02t
                double precision :: teff_sc_t, lnaa_sc_t, gammaa_sc_t
                if (isc == 1 .OR. isc == 2) then
                    ! all** already includes Maxwellian SC — use as-is
                    c20_sc=all20; c02_sc=all02; c11_sc=all11
                    c10_sc=all10; c01_sc=all01; c00_sc=all00
                else
                    ! Compute Maxwellian SC and add temporarily
                    teff_sc_t = aa * (vteff / 9.79d3)**2
                    call coulomb_log_ab(za, aa, teff_sc_t, npart, &
                                        za, aa, teff_sc_t, npart, lnaa_sc_t)
                    gammaa_sc_t = gamma0_sc * lnaa_sc_t * (za/aa)**2 * npart * za**2
                    call self_coll_max(vteff, gammaa_sc_t, &
                                       sc20t, sc02t, sc11t, sc10t, sc01t, sc00t)
                    c20_sc=all20+sc20t; c02_sc=all02+sc02t; c11_sc=all11+sc11t
                    c10_sc=all10+sc10t; c01_sc=all01+sc01t; c00_sc=all00+sc00t
                end if
                call FP_steady_state(c20_sc,c02_sc,c11_sc,c10_sc,c01_sc,c00_sc,fout)
            end block
            do iv = 1,nperp
                do imu = 1,npar
                    ix = index_mat(iv,imu)
                    fin(ix) = fout(iv,imu)
                enddo
            enddo
            write(*,*) 'Initial condition: SS (Maxwellian SC) done.'

        end select

    endif
            
    ! All phases run inside the solver in one uninterrupted loop.
    ! ntimes_cur/timestep_cur are no longer used here; the solvers
    ! iterate over ntimes(1..3)/timestep(1..3) directly.
    if (isc == -1) then
        call timefp_7pt_nl(all00,all10,all01,all11,all20,all02,fin,fout,time1)
    else
        call timefp_7pt(all00,all10,all01,all11,all20,all02,fin,fout,time1)
    end if
    time1 = time1 + sum(ntimes * timestep)
    
 !  
endif steady_state
!        
!
!!====================================================================
!
!! We compute the various moments of the vdf
!
call analysis(fout)

!====================================================================
!
! Domain adequacy check
! ---------------------
!
! The Dirichlet walls impose f=0, while the true solution there is
!   f(wall)/f(peak) = exp(-(V/vth)^2/2).
! Whatever that ratio is, it is the floor on the solution error, and it does
! NOT improve with grid refinement.  The check must use the FINAL (heated)
! temperature, not the cold background one: an RF or beam case can be perfectly
! well sized for its initial Maxwellian and badly under-sized for the
! distribution it actually evolves into.
!
! Target ratio 1e-12 (i.e. |v|max ~ 7.4 vth) was measured to give a shape error
! ~2e-8, against ~3e-5 at 5 vth.

domain_check: block

    double precision :: dens_chk, teff_chk, vth_chk, v_req
    double precision :: r_perp, r_par
    double precision, parameter :: f_target = 1.d-12

    call time_density(fout, dens_chk)
    call time_energy(fout, dens_chk, teff=teff_chk)

    if (teff_chk > 0.d0) then

        vth_chk = 9.79d3*dsqrt(teff_chk*1.d3/aa)
        v_req   = vth_chk*dsqrt(2.d0*dlog(1.d0/f_target))

        r_perp = vperp_max/vth_chk
        r_par  = min(dabs(vpar_min), dabs(vpar_max))/vth_chk

        write(*,*) ' '
        write(*,*) 'DOMAIN ADEQUACY CHECK (Dirichlet wall)'
        write(*,*) '--------------------------------------'
        write(*,'(A,F9.3,A,ES10.3,A)') '  final Teff  = ', teff_chk, &
             ' keV   ->  vth = ', vth_chk, ' m/s'
        write(*,'(A,ES9.2,A,ES10.3,A)') '  for f(wall)/f(peak) < ', f_target, &
             '  recommend |v|max > ', v_req, ' m/s'
        write(*,'(A,ES10.3,A,F6.2,A,ES9.2,A)') '  vperp_max   = ', vperp_max, &
             '  (', r_perp, ' vth)  f(wall)/f(peak) = ', &
             dexp(-0.5d0*r_perp**2), '  '//trim(merge('OK       ', &
             'TOO SMALL', vperp_max >= v_req))
        write(*,'(A,ES10.3,A,F6.2,A,ES9.2,A)') '  |vpar|max   = ', &
             min(dabs(vpar_min),dabs(vpar_max)), &
             '  (', r_par, ' vth)  f(wall)/f(peak) = ', &
             dexp(-0.5d0*r_par**2), '  '//trim(merge('OK       ', &
             'TOO SMALL', min(dabs(vpar_min),dabs(vpar_max)) >= v_req))

        if (vperp_max < v_req .or. &
            min(dabs(vpar_min),dabs(vpar_max)) < v_req) then
            write(*,*) ' '
            write(*,'(A,ES10.3,A)') '  ==> the domain truncates the tail. '// &
                 'Re-run with |v|max >= ', v_req, ' m/s;'
            write(*,*) '      refining the grid will NOT reduce this error.'
        end if

    end if

end block domain_check

if(isc /= 0) deallocate(sc20,sc02,sc11,sc10,sc01,sc00)

! TEST: plot the solution at vpar = 0

open(40,file=TRIM(outfile('fout_at_vpar0.txt')),status='unknown')
!open(41,file=TRIM(outfile('fstix_at_vpar0.txt')),status='unknown')

do iv = 1,nperp
    write(40,*) vperp(iv),fout(iv,jmid)
 !   write(41,*) vperp(iv),fstix(iv,jmid)
enddo
!close(41)
close(40)

open(40,file=TRIM(outfile('fout_at_vperp0.txt')),status='unknown')
!open(41,file=TRIM(outfile('fstix_at_vperp0.txt')),status='unknown')
!open(42,file=TRIM(outfile('fout_at_vperpmax.txt')),status='unknown')
!open(43,file=TRIM(outfile('fstix_at_vperpmax.txt')),status='unknown')

do iv = 1,npar
    write(40,*) vpar(iv),fout(1,iv)
 !   write(41,*) vpar(iv),fstix(1,iv)
 !   write(42,*) vpar(iv),fout(nperp,iv)
 !   write(43,*) vpar(iv),fstix(nperp,iv)
enddo
!close(43)
!close(42)
!close(41)
close(40)

!open(40,file=TRIM(outfile('fout_at_vparmax.txt')),status='unknown')
!open(41,file=TRIM(outfile('fstix_at_vparmax.txt')),status='unknown')	
!do iv = 1,nperp
!    write(40,*) vperp(iv),fout(iv,1)
!    write(41,*) vperp(iv),fstix(iv,1)
!enddo
!close(41)
!close(40)

call cpu_time(end_time)
write(*,*) ' '
    write(*,*) 'Simulation duration: ',end_time-start_time,'seconds'

end program FP_Coll_2D
