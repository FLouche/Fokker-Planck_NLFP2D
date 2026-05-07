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
!   Version 1.4 - 4 May 2026

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
!
use mod_linear
!use mod_timefp3!3
use mod_timefp3_upd
use mod_timefp_7pt
use mod_timefp3_nl
use mod_timefp_7pt_nl
!
use mod_anal
use mod_dislin_plots
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

! do loops indexes

integer :: ib,iv,imu,ix

! Other variables

!double PRECISION :: dens_tmp

! Math constants

double precision :: pi,twopi

common/mathcons/pi,twopi

data pi/3.141592653589793238462643d0/

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
!           -1: Stix's Maxwellian background
!            0: no self-collisions
!           +1: pitch-angle averaged solution used as collisional background
!  ising: homogeneity of the grid in vperp:
!            0: homogeneous grid
!           -1: inhomogeneous grid made of two domains (vperp<vbound and vperp>vbound) with different meshings
!                      for vperp <= vbound: nsing points (increase the density for small vperp)
!           +1: Quadratic spacing for higher resolution near vperp=0 (or vperp_min)
!                 ==> FD scheme needs to be adapted ===> DO NOT USE !!!!

namelist /INPUT/ casename, &
                 new_grid, &
                 nperp,npar,vperp_min,vperp_max,vpar_min,vpar_max,&
                ising,nsing,vbound,&
               nbulk,t,aa,ab,za,zb,ne,xpart,xb, &
                isource,beam_ekin,beam_angle_deg, & 
                beam_dvperp, beam_dvpar, taus, &
                irf,eplus,emin,kperp, &
                kpar,frek,delta_RF,b0,nharm, &
                icn, ntimes, timestep, iold, isc, ifd7, iplot_traces

!write(*,*) 'Read namelist'

read(5,INPUT)

! Coherence check: restarting from a previous solution requires the same grid
if (new_grid == -1 .and. iold == -1) then
    write(*,*) 'ERROR: new_grid=-1 (new grid) is incompatible with iold=-1 (restart).'
    write(*,*) 'A restart uses the solution from a previous run, which requires the same grid.'
    write(*,*) 'Set new_grid=0 to reuse the existing grid, or iold=0 to start fresh.'
    stop
endif

twopi=2.d0*pi

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
            
            call cblin(iv,imu,vt(ib),gammab(ib),maonmb(ib),c20,c02,c11,c10,c01,c00)
            
            colin20_sp(iv,imu,ib) = c20
            colin02_sp(iv,imu,ib) = c02
            colin11_sp(iv,imu,ib) = c11
            colin10_sp(iv,imu,ib) = c10
            colin01_sp(iv,imu,ib) = c01
            colin00_sp(iv,imu,ib) = c00
                        
        enddo background_loop        
                
            
    enddo vpar_loop
    
enddo vperp_loop


colin20 = sum(colin20_sp, DIM = 3)
colin02 = sum(colin02_sp, DIM = 3)
colin11 = sum(colin11_sp, DIM = 3)
colin10 = sum(colin10_sp, DIM = 3)
colin01 = sum(colin01_sp, DIM = 3)
colin00 = sum(colin00_sp, DIM = 3)

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
! We collect the factors in front of each derivative 

allocate(all00(nperp,npar))
allocate(all10,all01,all11,all20,all02,mold=all00)

if (irf == -1) then
    all00 = colin00
    all10 = colin10+rf10
    all01 = colin01+rf01
    all20 = colin20+rf20
    all02 = colin02+rf02
    all11 = colin11+rf11

else
    all00 = colin00
    all10 = colin10
    all01 = colin01
    all20 = colin20
    all02 = colin02
    all11 = colin11
endif

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
           call self_coll_max(vteff,sc20,sc02,sc11,sc10,sc01,sc00)
                      
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

    if (ifd7 == -1) then
        call FP_steady_state(all20,all02,all11,all10,all01,all00,fout)
    else
        call linear(all20,all02,all11,all10,all01,all00,fout)
    endif
    
    
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
    
    if(iold == -1) then !we start from previousy stored solution
    
        open(40,file=TRIM(outfile('xout.dat')),status='old')
        read(40,*) time1
        do ix=1,nbig
            read(40,*) fin(ix)
        enddo
        
    else
        
        time1=0.d0
    
! Stix's Maxwellian solution for sourceless case or zero function for driven case
    
    !open(40,file='fstart.txt',status='unknown')
    !
    do iv = 1,nperp
        do imu = 1,npar
            
            ix = index_mat(iv,imu)
            
            if (isource == 0) then

                fin(ix) = fstix(iv,imu)
                
            endif
            
        enddo
    enddo
    
    endif
    
    !close(40)
            
    if(ifd7 == -1) then
        if(isc == -1) then
            call timefp_7pt_nl(all00,all10,all01,all11,all20,all02,fin,fout,time1)
        else
            call timefp_7pt(all00,all10,all01,all11,all20,all02,fin,fout,time1)
        endif
    else
        if(isc == -1) then
            call timefp_nl(all00,all10,all01,all11,all20,all02,fin,fout,time1)
        else
            call timefp_upd(all00,all10,all01,all11,all20,all02,fin,fout,time1)
        endif
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
open(41,file=TRIM(outfile('fstix_at_vpar0.txt')),status='unknown')

do iv = 1,nperp
    write(40,*) vperp(iv),fout(iv,jmid)
    write(41,*) vperp(iv),fstix(iv,jmid)
enddo
close(41)
close(40)

open(40,file=TRIM(outfile('fout_at_vperp0.txt')),status='unknown')
open(41,file=TRIM(outfile('fstix_at_vperp0.txt')),status='unknown')
open(42,file=TRIM(outfile('fout_at_vperpmax.txt')),status='unknown')
open(43,file=TRIM(outfile('fstix_at_vperpmax.txt')),status='unknown')

do iv = 1,npar
    write(40,*) vpar(iv),fout(1,iv)
    write(41,*) vpar(iv),fstix(1,iv)
    write(42,*) vpar(iv),fout(nperp,iv)
    write(43,*) vpar(iv),fstix(nperp,iv)
enddo
close(43)
close(42)
close(41)
close(40)

open(40,file=TRIM(outfile('fout_at_vparmax.txt')),status='unknown')
open(41,file=TRIM(outfile('fstix_at_vparmax.txt')),status='unknown')	
do iv = 1,nperp
    write(40,*) vperp(iv),fout(iv,1)
    write(41,*) vperp(iv),fstix(iv,1)
enddo
close(41)
close(40)


call plot_endof_run()
if (ntimes /= 0 .and. iplot_traces == -1) call plot_time_traces()

end program FP_Coll_2D
