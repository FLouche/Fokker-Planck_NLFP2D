!
!   qlrfterm.f90
!   bechsi2011
!   
!   Created by Fabrice Louche on 20/10/11.
!   Copyright 2011 LPP-ERM/KMS. All rights reserved.
!

! Evaluation of the coefficients of the quasilinear RF term
!  at the limit of pure perpendicular diffusion
!
!*********************************************************

module mod_qlrfterm

contains

subroutine qlrfterm

! Initialisation
! --------------

use shared_plasma
use shared_grid
use shared_rf

use shared_FPterms

use shared_timer, only: notxt   ! RF_dirac.txt written only when notxt==0

use delta_dirac

use Complex_Bessel

implicit none

integer :: i,j,n
integer :: nz,ierr

double precision :: x,comfac1,comfac2
double precision, allocatable, dimension(:) :: dirac
complex*16, dimension(nharm+3) :: cy
double precision, allocatable, dimension(:) :: dn,ddndvperp

complex*16 :: arg,jnm1,jnp1,cjnm1,cjnp1,djnm1,djnp1
complex*16, dimension(0:nharm+2) :: jbes
complex*16, dimension(0:nharm+1) :: djbes

!=========================================================
!
! Preliminary computations
! ------------------------
!
!  Gaussian approximation of Delta function
!

allocate(dirac(npar))

if (notxt == 0) open (40,file='RF_dirac.txt',status='unknown')

do i = 1,npar

   x = vpar(i)-vres
   dirac(i) = delta_d(x,delta_RF)
   if (notxt == 0) write(40,*) vpar(i),dirac(i)

enddo

if (notxt == 0) close(40)

! Bessel functions
! ----------------

allocate(dn(nperp),ddndvperp(nperp))

do i = 1,nperp

   arg = eta*vperp(i)
   
   call cbesj(arg,0.d0,1,nharm+3,cy,nz,ierr)

   if(ierr /= 0) then
	write(*,*) 'Error in CBESJ #',ierr
   endif	
   ! We start at N=0 => cyr(n) = Real(J(n-1))
   ! The last value of the series is cyr(nharm+3) = Real(J(nharm+2))
   
    do n = 0,nharm+2
	    jbes(n) = cy(n+1)
    enddo

    jnm1=jbes(nharm-1)
	jnp1=jbes(nharm+1)

   ! Complex conjugates
   
   cjnm1=dconjg(jnm1)
   cjnp1=dconjg(jnp1)
   
   ! 1st derivative of Bessel functions
   ! J'(n) = (J(n-1)-j(n+1))/2

	do n = 0,nharm+1
		if(n == 0) then
			djbes(n)=-jbes(1)
		else
			djbes(n)=(jbes(n-1)-jbes(n+1))/2.d0
		endif
	enddo
	
	djnm1=djbes(nharm-1)
	djnp1=djbes(nharm+1)
	   
   ! Function (Theta)^2 appearing in ql term
   ! Theta = E+ J(n-1)+E- J(n+1)
   
   dn(i)=(eplus*jnm1+emin*jnp1)*(ceplus*cjnm1+cemin*cjnp1)
   
   ! Derivative of (Theta)^2 wrt Vperp
   
   ddndvperp(i) = 2.d0*dreal(eta*(eplus*djnm1+emin*djnp1)*(ceplus*cjnm1+cemin*cjnp1))
   
   
enddo

!=========================================================
!
! Coefficients of the FP equation
! -------------------------------

do j=1,npar
    
       comfac1 = vph-vpar(j)
       comfac2 = vpar(j)-vres
   
   do i=1,nperp
          
!      Coefficient of df/dvperp
!      ------------------------

rf10(i,j) = rfcte*dirac(j)*(comfac1*(comfac1*(dn(i)/vperp(i)+ddndvperp(i)) &
    -2.d0*vperp(i)/delta_RF**2*comfac2*dn(i))-vperp(i)*dn(i))
			  
!      Coefficient of d2f/dvperp2
!      --------------------------

rf20(i,j) = rfcte*dirac(j)*comfac1**2*dn(i)

!      Coefficient of df/dvpar
!      ------------------------

rf01(i,j) = rfcte*dirac(j)*(comfac1*(2.d0*dn(i)+vperp(i)*ddndvperp(i)) &
    -2.d0*comfac2*vperp(i)**2*dn(i)/delta_RF**2)

!      Coefficient of d2f/dvpar2
!      --------------------------

rf02(i,j) = rfcte*dirac(j)*vperp(i)**2*dn(i)

!      Coefficient of d2f/dvperpdvpar
!      ------------------------------

rf11(i,j) = rfcte*dirac(j)*2.d0*comfac1*vperp(i)*dn(i)



	enddo
enddo

!=========================================================

deallocate(dn,ddndvperp)

!*********************************************************

end subroutine qlrfterm

end module mod_qlrfterm