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
integer :: ifd7

end module shared_grid

!***************************************

module shared_plasma

implicit none

save

double precision, dimension(10) :: nb,maonmb,vt,gammab,t
double precision :: ne,npart
double precision :: gammaa
double precision, dimension(9) :: ab,zb,xb
double precision :: vcr
double precision :: xpart,aa,za
double precision ::  vteff,spit,tauie
	   
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
integer :: ntimes,iold,icn,isc
double precision :: timestep
double precision, dimension(:,:), allocatable :: fstix
!
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