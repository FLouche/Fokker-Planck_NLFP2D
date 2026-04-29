!   test_power_balance.f90
!
!   FPColl_2D
!   
!   Created by Fabrice Louche on 04/09/25.
!   Copyright 2025 LPP-ERM/KMS. All rights reserved.
    
!***************************************

    subroutine test_power_balance(fin)
    
    ! Initialisation
! --------------

use shared_grid
use shared_plasma
use shared_beam
use shared_rf

use mod_build_ss

use shared_timer
use shared_FPterms

use mod_ncint

use func_index

implicit none

double PRECISION, intent(in) :: fin(nbig)

double precision :: bigm(nbig,nbig),f_wrk(nbig)
double precision, dimension(nperp,npar):: f1,f2
double precision, dimension(nbulk) :: pcoll

integer:: ib,ix,iv,imu
double precision pmass
double precision, dimension(nperp,npar):: ekin
double precision :: plosses,psource,pRF,pSC

double PRECISION, allocatable, dimension(:,:) :: rf00

external dgemv

data pmass/1.6726d-27/ !proton mass in kg

! ***************************************************
! Collisional power density
! ***************************************************

! We loop over each background species and compute the respective matrix operator

taum = 0.d0 ! temporary 


do ib = 1,nbulk
    
     ! We compute the df/dt collisonal term in matrix form
            
    call build_ss(colin20_sp(:,:,ib),colin02_sp(:,:,ib),colin11_sp(:,:,ib),colin10_sp(:,:,ib),colin01_sp(:,:,ib),colin00_sp(:,:,ib),bigm)
    
     ! We apply the operator to the solution
    
    call dgemv('N',nbig,nbig,1.d0,bigm,nbig,fin,1,0.d0,f_wrk,1)
    
    do iv = 1,nperp
                
        do imu = 1,npar
            
            ekin(iv,imu) = 0.5*pmass*aa*(vperp(iv)**2+vpar(imu)**2)
            ix = index_mat(iv,imu)
            f1(iv,imu) = ekin(iv,imu)*f_wrk(ix)*jacob(iv,imu)
            
        enddo
    enddo
    
    call ncint_2D(f1,pcoll(ib))
    
enddo

! Self-collisions

if(isc /= 0) then
    
    call build_ss(sc20,sc02,sc11,sc10,sc01,sc00,bigm)
            
    ! We apply the operator to the solution
    
        call dgemv('N',nbig,nbig,1.d0,bigm,nbig,fin,1,0.d0,f_wrk,1)
        
        do iv = 1,nperp
                
            do imu = 1,npar
            
              ix = index_mat(iv,imu)
              f1(iv,imu) = ekin(iv,imu)*f_wrk(ix)*jacob(iv,imu)
            
            enddo
         enddo
    
    call ncint_2D(f1,pSC)
    
else
    
    pSC = 0.d0
        
endif
    


! ***************************************************
! Beam (+ losses) power density
! ***************************************************

    if (isource == -1) then    
        taum = 1.d0/taus
        
        do iv = 1,nperp
            do imu = 1,npar
            
            ix = index_mat(iv,imu)
            f1(iv,imu) = -ekin(iv,imu)*fin(ix)*taum*jacob(iv,imu)
            f2(iv,imu) = ekin(iv,imu)*source(iv,imu)*jacob(iv,imu)
            
            enddo
        enddo
        
        call ncint_2D(f1,plosses)
        call ncint_2D(f2,psource)
        
    else 
        
        taum = 0.d0
        plosses = 0.d0
        psource = 0.d0
        
    endif
    
! ***************************************************
! RF power density
! ***************************************************

    ! We compute the df/dt RF term in matrix form
    
    if(irf == -1) then
        
        allocate(rf00(nperp,npar))
        rf00 = 0.d0
    
        call build_ss(rf20,rf02,rf11,rf10,rf01,rf00,bigm)
        
        deallocate(rf00)
        
          ! We apply the operator to the solution
    
        call dgemv('N',nbig,nbig,1.d0,bigm,nbig,fin,1,0.d0,f_wrk,1)
        
        do iv = 1,nperp
                
            do imu = 1,npar
            
              ix = index_mat(iv,imu)
              f1(iv,imu) = ekin(iv,imu)*f_wrk(ix)*jacob(iv,imu)
            
            enddo
         enddo
    
    call ncint_2D(f1,pRF)
        
    else ! no RF
        
        pRF = 0.d0
        
    endif

    
    write(*,*) ''
    write(*,*) 'Power density balance:'
    write(*,*) '-----------------'
    write(*,*) 'Collisions:      '
    write(*,*) '-----------'
    do ib=1,nbulk
        if (ib == 1) then
            write(*,*) 'Collisions with electrons:      ',pcoll(ib)/1d6, 'MW/m**3'
        else
            write(*,*) 'Collisions with ions ',ib-1, ': ',pcoll(ib)/1d6, 'MW/m**3'
        endif
    enddo
    if(isc /= 0) write(*,*) 'Self-collisions:      ',pSC/1d6, 'MW/m**3'
    write(*,*) '-----------------'
    write(*,*) 'Total collisions:            ',(SUM(pcoll,dim=1)+pSC)/1d6, 'MW/m**3'
    write(*,*) 'Beam source:     ',psource/1d6, 'MW/m**3'
    write(*,*) 'Particle losses: ',plosses/1d6, 'MW/m**3'
    write(*,*) ''
    write(*,*) 'RF power:        ',pRF/1d6, 'MW/m**3'
    write(*,*) '--------------------------------------------'

    write(*,*) 'Total          : ',(SUM(pcoll,dim=1)+pSC+plosses+psource+pRF)/1d6, 'MW/m**3'

    

! ===========================================================

end subroutine test_power_balance