!
!   ncint.f
!   nlfp2011
!   
!   Created by Fabrice Louche on 20/01/11.
!   Copyright 2011 LPP-ERM/KMS. All rights reserved.
!   
!
!***********************************************
!*          NEWTON-COTES FORMULAS              *
!* on a grid of arbitrary distributed points   *
!*                             *
!***********************************************

module mod_ncint

contains

subroutine ncint_2D(f,s)

! Initialisation
! -------------- 

use shared_grid
       
implicit none

double precision, intent(in) :: f(nperp,npar)
double precision, intent(out):: s

double precision dvp_local
integer i,j

! ===========================================================================

! F is the values of the function to be integrated on the complete grid
! s  : output value of the integral

s=0.d0

! Compute 2π ∫∫ f * vperp * dvperp * dvpar
 
DO i = 1, nperp-1
    dvp_local = vperp(i+1) - vperp(i)   ! works for any grid
    DO j = 2, npar-1                                ! skip Dirichlet boundaries
        s = s + f(i,j) * dvp_local * dvpar
    END DO
END DO

!
!============================================================================
!
end subroutine ncint_2D
	   


!
!============================================================================

end module mod_ncint