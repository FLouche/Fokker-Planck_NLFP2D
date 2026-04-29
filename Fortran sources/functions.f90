!
!   functions.f90
!
!   FPColl_2D
!   
!   Created by Fabrice Louche on 02/07/25.
!   Copyright 2025 LPP-ERM/KMS. All rights reserved.

!   
module func_index

contains 

double precision function index_mat(i,j)
    
use shared_grid

implicit none
    
integer, intent(in):: i,j

index_mat = (i-1)*npar+j

return 

end function index_mat

!***************************************


subroutine index_mat_inv(ix,ipe,ipa)

use shared_grid

implicit none

integer, intent(in):: ix
integer, intent(out) :: ipe,ipa

ipa = modulo(ix,npar)

if (ipa == 0) ipa=npar

ipe = (ix-ipa)/npar+1

end subroutine index_mat_inv

!***************************************


end module func_index

!***************************************



    module delta_dirac
    
    contains
    
    double PRECISION function delta_d(x,deltax)
    
    
    implicit none
    
    double PRECISION, intent(in):: x, deltax
    double PRECISION pi
    
    data pi/3.141592653589793238462643d0/

    
    delta_d = dexp(-(x/deltax)**2)/dsqrt(pi)/deltax
    
    end function delta_d
    
    end module delta_dirac

    
    