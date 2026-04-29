!
!   beam_source.f90
!   bechsi2011
!   
!   Created by Fabrice Louche on 23/11/11.
!   Modified by FL on 18/08/25 :
!       - vperp, vpar system
    
!   Copyright 2011-2025 LPP-ERM/KMS. All rights reserved.
!
!*******************************************************

! Computes the source term for NBI
! --------------------------------
module mod_beam

contains


subroutine beam_source

! Common information

use shared_plasma
use shared_grid
use shared_beam
use mod_ncint

use delta_dirac

use func_index


implicit none 

intrinsic dsqrt

! Working variables

integer :: i,j,ix1

double precision :: x,S0, fac
double precision, dimension(nperp) :: svel
double precision, dimension(npar) :: spa

double PRECISION :: pi,twopi
double precision :: dens_beam

double precision pmass, kev_in_J

common/mathcons/pi,twopi
data pmass/1.6726d-27/ !proton mass in kg
data kev_in_J/1.60218d-16/ !convert keV to Joule
    

! =====================================================

! Beam velocity is computed from beam injection energy
    
    
    beam_v = DSQRT(2.d0*beam_ekin*kev_in_J/(pmass*aa))
    
    beam_vperp = beam_v*dsin(beam_angle)
    beam_vpar = beam_v*dcos(beam_angle)

    write(*,*) 'Beam injection perp. velocity = ',beam_vperp,' m/s'
    write(*,*) 'Beam injection par. velocity = ',beam_vpar,' m/s'



!=================================================

! Approximation of delta functions
! --------------------------------

S0 = npart/taus

do i=1,nperp
	x = vperp(i)-beam_vperp
	svel(i) = delta_d(x,beam_dvperp)
enddo


do i=1,npar
	x = vpar(i)-beam_vpar
	spa(i) = delta_d(x,beam_dvpar)
enddo


!=================================================

! Source term
! -----------

! Normalising factor

fac = twopi*beam_vperp

do i=1, nperp
	do j=1,npar
	
		source(i,j)=S0*svel(i)*spa(j)/fac
        ix1 = index_mat(i,j)
        source_v(ix1)=source(i,j) ! useful for matrix formulations
        enddo
    enddo

call ncint_2D(source*jacob,dens_beam)

write(*,*) 'Beam density ', dens_beam, '/m3/s'
write(*,*) 'Total density ',npart,'/m3'

   
end subroutine beam_source


!*******************************************************

end module mod_beam