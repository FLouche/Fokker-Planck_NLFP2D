!
!   sparse_solve.f
!   
!   Created by Fabrice Louche on 19/01/11.
!   Modifications:
    ! - 30/03/2026: new pardiso_solver module
!
!   Copyright 2011-2026 LPP-ERM/KMS. All rights reserved.
!   

! ***********************************************************************************
! 
!  The system is solved taking into account the sparsity of the matrix.
!
!  The PARDISO (Parallel Direct Sparse Solver) interface is used. It's provided 
!   by Intel(c) Math Kernel Libray (MKL). 
!   See "Intel(R) Math Kernel Library Reference Manual", ch. 8 p. 2495
!    or http://www.pardiso.org 
!
! ***********************************************************************************
!
module mod_sparse_solve

contains

subroutine sparse_solve_ss(flag,bigm,b,x)

! Case of a steady-state one-shot computation

! Initialisation
! --------------
use shared_grid
!use MKL_PARDISO
use pardiso_solver

implicit none

!include 'mkl_pardiso.f90'

	
integer, intent(in) :: flag
integer :: nz

double precision, intent(in), dimension(nbig,nbig) :: bigm
double precision, intent(inout), dimension(nbig) :: b,x

integer i,j,index

!	 
!  PARDISO input parameters
!  ------------------------  
! Other variables
integer n, nrhs, error, flag1
integer, dimension(:), allocatable :: ia,ja
double precision, dimension(:), allocatable :: a ! b is actually bigv

! ***********************************************************************************
!
! FIRST COMPUTATION(S)
! =================

!first_comp: if (flag == -1) then

! Set up PARDISO control parameter

n = nbig
nrhs = 1

! Settings
! --------       
      
! --> computes the number of non-zeros elements of Bigm
!      
nz=0         		
do i=1,nbig
	do j=1,nbig
	 if(bigm(i,j)/=0.d0)  then
		       nz=nz+1
	 endif
	enddo
enddo

write(*,*) 'Number of non-zeros in Bigm',nz, ' = ',100.d0*nz/nbig**2, ' %'
	   
! --> allocates memory space to working arrays
!
allocate(a(1:nz))
!allocate(x(1:nbig))
allocate(ia(1:(nbig+1)))
allocate(ja(1:nz))

! --> non-zeros elements of sparse Bigm matrix are stored
!
index=0
do i=1,nbig
	flag1 = 0
	do j=1,nbig
		if(bigm(i,j) /= 0.d0) then
			index=index+1
			if(flag1 == 0) then
				flag1 = 1
				ia(i)= index
			endif
		    a(index)=bigm(i,j)
			ja(index)=j
		endif
	enddo
enddo 
	   
ia(nbig+1) = nz+1

CALL pardiso_solve_steady(n, a, ia, ja, b, x, &
                            mtype=11, msglvl=0, error=error)
  

! ***********************************************************************************

       end subroutine sparse_solve_ss
	   
! ***********************************************************************************

end module mod_sparse_solve