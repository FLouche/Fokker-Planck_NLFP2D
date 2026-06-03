!   test_power_balance_7pt.f90
!
!   FPColl_2D
!
!   Created by Fabrice Louche on 04/09/25.
!   Copyright 2025 LPP-ERM/KMS. All rights reserved.

!***************************************

subroutine test_power_balance_7pt(fin)

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

double precision, intent(in) :: fin(nbig)

double precision, dimension(nbulk) :: pcoll
integer  :: ib
double precision :: plosses, psource, pRF, pSC, pSC_perp, pSC_par

taum = 0.d0

call time_power_7pt(fin, npart, pcoll, pRF, psource, plosses, pSC, pSC_perp, pSC_par)

write(*,'(/,A)')    'POWER DENSITY BALANCE:'
write(*,'(A)')      '----------------------'
write(*,'(A)')      ''
write(*,'(A)')      '==========================='
write(*,'(A)')      ' Collisions:'
write(*,'(A)')      ' -----------'
do ib = 1, nbulk
    if (ib == 1) then
        write(*,'(A,ES12.5,A)') &
            '  Collisions with electrons:     ', pcoll(1)/1d6, ' MW/m**3'
    else
        write(*,'(A,I2,A,ES12.5,A)') &
            '  Collisions with ions           ', ib-1, ' :  ', pcoll(ib)/1d6, ' MW/m**3'
    end if
end do
write(*,'(A)')      ' ------------------'
write(*,'(A)')      ''
if (isc /= 0) then
    write(*,'(A)')      ' Self-collisions:'
    write(*,'(A,ES12.5,A)') &
        '     - P_SC_perp :  ', pSC_perp/1d6, ' MW/m**3'
    write(*,'(A,ES12.5,A)') &
        '     -  P_SC_par :  ', pSC_par/1d6,  ' MW/m**3'
    write(*,'(A)')      ' ------------------'
    write(*,'(A,ES12.5,A)') &
        '  Total self-collisions:        ', pSC/1d6, ' MW/m**3'
    write(*,'(A)')      ''
end if
write(*,'(A)')      ' -----------------'
write(*,'(A,ES12.5,A)') &
    '  Total collisions:             ', (SUM(pcoll)+pSC)/1d6, ' MW/m**3'
write(*,'(A)')      '==========================='
write(*,'(A)')      ''
write(*,'(A,ES12.5,A)') &
    '  Beam source:                  ', psource/1d6, ' MW/m**3'
write(*,'(A)')      ''
write(*,'(A,ES12.5,A)') &
    '  Particle losses:              ', plosses/1d6, ' MW/m**3'
write(*,'(A)')      ''
write(*,'(A)')      '==========================='
write(*,'(A)')      ''
write(*,'(A,ES12.5,A)') &
    '  RF power:                     ', pRF/1d6, ' MW/m**3'
write(*,'(A)')      ''
write(*,'(A)')      '  --------------------------------------------'
write(*,'(A,ES12.5,A)') &
    '  Total balance:                ', &
    (SUM(pcoll)+pSC+plosses+psource+pRF)/1d6, ' MW/m**3'
write(*,'(A)')      ''

end subroutine test_power_balance_7pt
