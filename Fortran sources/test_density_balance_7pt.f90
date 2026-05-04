! TEST: DOES THE GLOBAL EQUATION RESPECT PARTICLE DENSITY CONSERVATION
    
    !   FPColl_2D
!   
!   Created by Fabrice Louche on 04/09/25.
!   Modified on 19/11/2025: cylindrical coordinates

!   Copyright 2025 LPP-ERM/KMS. All rights reserved.
    
!***************************************


subroutine test_density_balance_7pt(fin)

! Initialisation
! --------------

use shared_grid
use shared_plasma
use shared_beam

!use mod_build_ss
use shared_rf

use shared_timer
use shared_FPterms
 
use mod_ncint

use func_index
use mod_apply_operator

INTEGER, PARAMETER :: dp = KIND(1.0D0)


double PRECISION:: fin(nbig)!, intent(in) 

 REAL(dp), dimension(nbig) :: Lf_col_sp,Lf_col,Lf_col_sc,Lf_RF                  ! various operators applied to f

double PRECISION, dimension(nperp,npar) :: f1,f2,f3,f4,f5
double precision f1_sp(nperp,npar,nbulk)
double PRECISION n1,n2,n3,n4,n1_sp(nbulk),n5

double PRECISION, allocatable,dimension(:,:):: rf00

integer ix,iv,imu,ib

external dgemv
    
    ! We compute the df/dt total collisonal term in matrix form
    
    taum = 0.d0 ! temporary 
    
    CALL apply_operator(colin20,colin02,colin11,colin10,colin01,colin00, &
                        fin, Lf_col_sp)
            
     ! We compute the df/dt total collisonal term per species in matrix form
    
    do ib = 1,nbulk
    
     ! We compute the df/dt collisonal term in matrix form
            
    CALL apply_operator(colin20_sp(:,:,ib), colin02_sp(:,:,ib), &
                        colin11_sp(:,:,ib), colin10_sp(:,:,ib), &
                        colin01_sp(:,:,ib), colin00_sp(:,:,ib), &
                        fin, Lf_col_sp)
    do iv = 1,nperp
        do imu = 1,npar
            
            ix = index_mat(iv,imu)

    f1_sp(iv,imu,ib) = Lf_col_sp(ix)*jacob(iv,imu)
        enddo
    enddo
    
    
    enddo
    
    ! Self-collisions term
    
    if (isc/=0) then
    
         CALL apply_operator(sc20,sc02,sc11,sc10,sc01,sc00, &
                        fin, Lf_col_sc)
        
    endif
    
    
! -----------------------
    
    ! We compute the df/dt RF term in matrix form
    
    if(irf == -1) then
        
        allocate(rf00(nperp,npar))
        rf00 = 0.d0
        
        CALL apply_operator(rf20,rf02,rf11,rf10,rf01,rf00, &
                        fin, Lf_RF)
                    
        deallocate(rf00)
        
        
    else
        
        Lf_RF = 0.d0
        
    endif
    
    if (isource == -1) then    
        taum = 1.d0/taus
    else
        taum = 0.d0
    endif

!   
    !We integrate the whole df/dt (incl. source and particle losses) over d3v
        
do iv = 1,nperp
        do imu = 1,npar
            
            ix = index_mat(iv,imu)
            f1(iv,imu) = Lf_col(ix)*jacob(iv,imu)

    
            f2(iv,imu) = source(iv,imu)*jacob(iv,imu)
            f3(iv,imu) = -fin(ix)*taum*jacob(iv,imu)
            f4(iv,imu) = Lf_rf(ix)*jacob(iv,imu)
            
            if (isc/=0) f5(iv,imu) = Lf_col_sc(ix)*jacob(iv,imu)
            
        enddo
enddo

call ncint_2D(f1,n1)

do ib=1,nbulk
    call ncint_2D(f1_sp(:,:,ib),n1_sp(ib))
enddo

call ncint_2D(f2,n2)
call ncint_2D(f3,n3)
call ncint_2D(f4,n4)
if (isc/=0) then
    call ncint_2D(f5,n5)
else
    n5=0.d0
endif

write(*,*) ''
write(*,*) 'Particle Balance:'
write(*,*) '-----------------'
write(*,*) 'Collisions (total):      ',(n1+n5)/npart, '/s'!
write(*,*) '    '

do ib=1,nbulk
    
    if (ib == 1) then
        write(*,*) 'Collisions with electrons:      ',n1_sp(ib)/npart, '/s'
    else
        write(*,*) 'Collisions with ions ',ib-1, ': ',n1_sp(ib)/npart, '/s'
    endif
    
enddo
if (isc/=0) write(*,*) 'Self-collisions: ',n5/npart, '/s'

write(*,*) '    '   
write(*,*) 'Beam source:     ',n2/npart, '/s'
write(*,*) 'Particle losses: ',n3/npart, '/s'
write(*,*) 'RF term:         ',n4/npart, '/s'
write(*,*) '--------------------------------------------'
write(*,*) 'Total          : ',(n1+n2+n3+n4+n5)/npart, '/s'

!

end subroutine test_density_balance_7pt
