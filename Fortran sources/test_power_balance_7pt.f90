!   test_power_balance.f90
!
!   FPColl_2D
!   
!   Created by Fabrice Louche on 04/09/25.
!   Copyright 2025 LPP-ERM/KMS. All rights reserved.
    
!***************************************

    subroutine test_power_balance_7pt(fin)
    
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

double precision, dimension(nbulk) :: pcoll

integer:: ib,ix,iv,imu
double precision :: plosses,psource,pRF,pSC,pSC_perp,pSC_par

double PRECISION, allocatable, dimension(:,:) :: rf00

external dgemv

! ***************************************************
! Collisional power density
! ***************************************************

! We loop over each background species and compute the respective matrix operator

taum = 0.d0 ! temporary 


call time_power_7pt(fin, npart, pcoll, pRF, psource, plosses, pSC, pSC_perp, pSC_par)

    
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

end subroutine test_power_balance_7pt