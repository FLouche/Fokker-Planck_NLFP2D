program FP_Coll_2D

! ====================================================
! First attempt at a Non-linear Fokker-Planck code
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
!   Version 1.65 - 29 May 2026 (FL):
    
!   new module to build linear FP terms -> for varying Coulomb log vs time
!       the treatment of varying Coulomb log is only accounted for 
!         in for the optionms isc=0 and isc=1 (no SC or MAxwellian SC)

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
!           +1: Stix's Maxwellian solution without RF used as collisional background (constant temperature)
!  ising: homogeneity of the grid in vperp:
!            0: homogeneous grid
!           -1: inhomogeneous grid made of two domains (vperp<vbound and vperp>vbound) with different meshings
!                      for vperp <= vbound: nsing points (increase the density for small vperp)
!           +1: Quadratic spacing for higher resolution near vperp=0 (or vperp_min)
!                 ==> FD scheme needs to be adapted ===> DO NOT USE !!!!
!
! Convergence for time-dependent simulation:
!
!1. Every n_ss_window steps, a line like [SS] step=50  dE/E= 1.23E-02  dn/n= 4.56E-03  dP/Pd= 7.89E-03 [tol= 1.00E-03]
!should appear on stdout.
!2. When all three criteria drop below ss_tol, the run stops early and prints [SS] CONVERGED at step NNN ....
!3. After early exit, all output files (fout.txt, xout.dat, *_vs_time.txt) should be complete and the end-of-run plots
!should still be produced.
!
!If ss_tol is too tight (run never converges) or too loose (stops too early), adjust it together with n_ss_window. A
!wider window is more immune to short-term fluctuations.
!
! Initial solution for time-dependent solver:
!
!  - istart = 0 -> empty solution, f = 0 everywhere; only possible when isource = -1. The code should include a test at the start to guarantee that isource = -1
!  - istart = 1 ->  Stix's solution as initial solution 
!  - istart = 2 -> initial solution is the steady-state solution of the linear time independent code, computed without the self-collisions
!  - istart = 3 -> initial solution is the steady-state solution of the linear time independent code, computed with a Maxwellian background for the self-collisions

namelist /INPUT/ casename, &
                 new_grid, &
                 nperp,npar,vperp_min,vperp_max,vpar_min,vpar_max,&
                ising,nsing,vbound,&
               nbulk,t,aa,ab,za,zb,ne,xpart,xb, &
                isource,beam_ekin,beam_angle_deg, & 
                beam_dvperp, beam_dvpar, taus, &
                irf,eplus,emin,kperp, &
                kpar,frek,delta_RF,b0,nharm, &
                icn, ntimes, timestep, iold, istart, isc, &
                i_ss_check, n_ss_window, ss_tol
                

!write(*,*) 'Read namelist'

read(5,INPUT)

! Coherence check: restarting from a previous solution requires the same grid
! (only relevant for NLSC runs where the sum_phi kernel is cached on disk)
if (new_grid == -1 .and. iold == -1 .and. isc == -1) then
    write(*,*) 'ERROR: new_grid=-1 (new grid) is incompatible with iold=-1 (restart).'
    write(*,*) 'A restart uses the solution from a previous run, which requires the same grid.'
    write(*,*) 'Set new_grid=0 to reuse the existing grid, or iold=0 to start fresh.'
    stop
endif

! Coherence checks for istart (only relevant for a fresh TD run)
if (ntimes /= 0 .and. iold /= -1) then
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
        if (isc == 1) then
           ! Compute self-collision Coulomb log using Stix effective temperature
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
             write(*,*) 'Self-collisions accounted (Maxwellian approximation)'
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
steady_state: if(ntimes == 0) then
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
            if (isc == 1) then
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
                if (isc == 1) then
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
            
    if(isc == -1) then
        call timefp_7pt_nl(all00,all10,all01,all11,all20,all02,fin,fout,time1)
    else
        call timefp_7pt(all00,all10,all01,all11,all20,all02,fin,fout,time1)
    endif
    
 !  
endif steady_state
!        
!
!!====================================================================
!
!! We compute the various moments of the vdf
!
call analysis(fout)

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
