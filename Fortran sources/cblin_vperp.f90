  !*******************************************************************
!   Computation of the Fokker-Planck collision operator in
!     (vperp,vpar) coordinates for a population of Maxwellian 
!        background species
!
!     Created by Fabrice Louche on 18/11/25.
!   Copyright 2025 LPP-ERM/KMS. All rights reserved.
!*******************************************************************

subroutine cblin(iperp,ipar,vth,gamma,maonmb,c20,c02,c11,c10,c01,c00)

! Initialisation
! --------------
!use shared_plasma
use shared_grid

implicit none

double precision, intent(in) :: vth, gamma,maonmb
double precision, intent(out) :: c20,c02,c10,c01,c00,c11
integer, intent(in) :: iperp,ipar

double precision pi,twopi,sq2
double precision arg,arg2,arg3,arg5,coef1,v1,v2,v3
double precision func1,derfarg,chandra,chandrap,dGonv

double precision Theta, Phi, Psi, dPsi, dTheta

common/mathcons/pi,twopi

!  DERF is the error function (double precision)

!external derf

!=======================================================================

! The Fokker-Planck diffusion/friction terms for a Maxwellian background
!  are only function of the velocity magnitude
!
! Ref. Karney Comp. Phys. Rep. 4(3-4), (Aug. 1986)
! arXiv: physics/0501066v1

! We use the formalism of Van Eester, Plasma Phys. Control. Fusion 36 (1994)

sq2 = dsqrt(2.d0)


v1=dSQRT(vperp(iperp)**2+vpar(ipar)**2)
v2=vperp(iperp)**2+vpar(ipar)**2
v3=v1**3
    
! Main mathematical functions

arg=v1/sq2/vth
arg2=v2/2.d0/vth**2
arg3=arg**3
arg5=arg**5

coef1=2.d0/dsqrt(pi)

func1 = derf(arg) ! Erf[u]
derfarg = coef1*dexp(-arg2) ! Erf'[u]
    
! Chandrasekhar function is defined as G=[Erf-uERf']/(2u^2)

!if (arg < 1d-1) then
!    chandra = (2*arg/3.d0-2*arg2/5.d0+arg5/7.d0)/DSQRT(pi)
!else
    
chandra = (func1-arg*derfarg)/(2*arg**2)

!endif

! First derivative of G wrt. its argument

chandrap = derfarg-2*chandra/arg

! First derivative of G/v wrt. v

dGonv = (arg*chandrap-chandra)/v2

!========================

!      Theta
!      -----

Theta = gamma*chandra/v1

!      Phi
!      ---

Phi = gamma/2.d0/v3*(func1-3.d0*chandra)

!      Psi
!      ---

Psi = gamma/v1/vth**2*maonmb*chandra

!     dPsi/dv
!     -------------------

dPsi = gamma/vth**2*maonmb*dGonv

!     dTheta/dv
!     -------------------

dTheta = gamma*dGonv

!=======================================================================
! 
! Assembling the coefficients of the derivative of f
!
!=====


c20 = Theta+vpar(ipar)**2*Phi

c02 = Theta+vperp(iperp)**2*Phi

c11 = -2.d0*+vpar(ipar)*vperp(iperp)*Phi

c10 = vperp(iperp)*(-2.d0*Phi+dTheta/v1+Psi+(Theta+v2*Phi)/vperp(iperp)**2)

c01 = vpar(ipar)*(-2.d0*Phi+dTheta/v1+Psi)

c00 = 3.d0*Psi+v1*dPsi


!===================================================================


end subroutine cblin

!*******************************************************************

