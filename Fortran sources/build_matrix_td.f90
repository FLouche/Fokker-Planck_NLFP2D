!**********************************************************
!
! This routine builds the matrix containing the terms
!  of the FP equation written on a 2D grid using a
!   finite-differences differentiation scheme and an 
!   implicit time-difference scheme
!
!   Preliminary version includes:
!
!    - homogeneous grid only (ising = 0)
!    - only Coulomb collisions (no df/dvdmu term)
!    - beam source
!    - Various time differencing algorithms are available:
!           - Cranck-Nicholson method
!           - implicit  method (03/09)
!
!   by Fabrice Louche - version 1.0 (20/08/2025)
!                     - version 1.1 (03/09/2025):
!                       + implicit scheme
!
!     (c) LPP ERM-KMS
!**********************************************************

module mod_build_td

contains

subroutine build_cn(fin,bigm_ss,bigm_cn,bigv)

! Crank-Nicholson method

! Initialisation
! --------------

use shared_grid
use shared_plasma
use shared_beam 
use shared_timer

use func_index

!use mkl_blas

implicit none

external sparse_matrix_vect_mult

double precision,  dimension (nbig,nbig), intent(in) :: bigm_ss

double precision,  dimension (nbig,nbig), intent(out) :: bigm_cn

double precision, dimension(nbig):: fin
double precision, dimension (nbig),intent(out) :: bigv
double PRECISION, dimension(:,:), allocatable :: idmat,mat_wrk
double PRECISION, dimension(:), allocatable :: f_wrk

integer j

! Preliminary: we build the identity matrix

allocate(idmat(nbig,nbig))

idmat = 0.d0
    
forall(j= 1:nbig) idmat(j,j) = 1.d0

!====================================================================
!
! The left-hand side is computed
!
! We first express the source term in vector form
!
    
    
    allocate(mat_wrk(nbig,nbig),f_wrk(nbig))
    
    mat_wrk = timestep*bigm_ss/2.d0+idmat
                
    call sparse_matrix_vect_mult(mat_wrk,fin,f_wrk)
    
    bigv = f_wrk+timestep*source_v    
            
! New rhs
    
    bigm_cn = idmat-timestep*bigm_ss/2.d0
    
    deallocate(mat_wrk,f_wrk,idmat)
        
                    
!====================================================================


end subroutine build_cn

!====================================================================


subroutine build_imp(fin,bigm_ss,bigm_imp,bigv)

! Implicit scheme

! Initialisation
! --------------

use shared_grid
use shared_plasma
use shared_beam 
use shared_timer

use func_index

!use mkl_blas

implicit none

!external dgemv

double precision,  dimension (nbig,nbig), intent(in) :: bigm_ss

double precision,  dimension (nbig,nbig), intent(out) :: bigm_imp

double precision, dimension(nbig):: fin
double precision, dimension (nbig),intent(out) :: bigv
double PRECISION, dimension(:,:), allocatable :: idmat

integer j
!integer iv, ip, ix1, j

! Preliminary: we build the identity matrix

allocate(idmat(nbig,nbig))

idmat = 0.d0
    
forall(j= 1:nbig) idmat(j,j) = 1.d0

!====================================================================
!
! The left-hand side is computed
!
! We first express the source term in vector form
!
    
    !do iv=1,nv
    !    do ip=1,npa
    !        ix1 = index_mat(iv,ip)
    !        bigv(ix1)=source(iv,ip)
    !    enddo
    !enddo

    bigm_imp = idmat-timestep*bigm_ss 
    
    bigv = fin+timestep*source_v
    
            
!====================================================================


end subroutine build_imp

!====================================================================

!====================================================================

!**********************************************************


end module mod_build_td
