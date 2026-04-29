!
!   self_coll_max.f
!   FP2D_QLRF 
!   
!   Created by Fabrice Louche on 17/12/25.
!   Copyright 2025 LPP-ERM/KMS. All rights reserved.
!   
!
!***********************************************
!*  Computation of the self-collision FP term  *
!*    assuming Maxwellian background           *
!***********************************************

    subroutine self_coll_max(vth,sc20,sc02,sc11,sc10,sc01,sc00)

! Initialisation
! -------------- 

use shared_grid
use shared_plasma
       
implicit none

double precision, intent(in):: vth
double precision, intent(out), dimension(nperp,npar) :: sc20,sc02,sc11,sc10,sc01,sc00
double precision :: c20,c02,c11,c10,c01,c00

integer i,j

external cblin

! ===========================================================================

do i=1,nperp
    do j=1,npar
        call cblin(i,j,vth,gammaa,1.d0,c20,c02,c11,c10,c01,c00)        

        sc20(i,j) = c20
        sc02(i,j) = c02
        sc10(i,j) = c10
        sc01(i,j) = c01
        sc00(i,j) = c00
        sc11(i,j) = c11
        
        
    enddo
enddo


! ===========================================================================


end subroutine self_coll_max