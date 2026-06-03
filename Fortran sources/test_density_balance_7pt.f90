! TEST: DOES THE GLOBAL EQUATION RESPECT PARTICLE DENSITY CONSERVATION
!
!   FPColl_2D
!
!   Created by Fabrice Louche on 04/09/25.
!   Modified on 19/11/2025: cylindrical coordinates
!
!   Copyright 2025 LPP-ERM/KMS. All rights reserved.

!***************************************

subroutine test_density_balance_7pt(fin)

use shared_grid
use shared_plasma
use shared_beam
use shared_rf
use shared_timer
use shared_FPterms
use mod_ncint
use func_index
use mod_apply_operator

INTEGER, PARAMETER :: dp = KIND(1.0D0)

double precision :: fin(nbig)

real(dp), dimension(nbig)        :: Lf_col_sp, Lf_col, Lf_col_sc, Lf_RF
double precision, dimension(nperp,npar) :: f1, f2, f3, f4, f5
double precision :: f1_sp(nperp,npar,nbulk)
double precision :: n1, n2, n3, n4, n5, n1_sp(nbulk)
double precision, allocatable :: rf00(:,:)
integer :: ix, iv, imu, ib

taum = 0.d0

call apply_operator(colin20, colin02, colin11, colin10, colin01, colin00, fin, Lf_col)

do ib = 1, nbulk
    call apply_operator(gammab(ib)*colin20_sp(:,:,ib), &
                        gammab(ib)*colin02_sp(:,:,ib), &
                        gammab(ib)*colin11_sp(:,:,ib), &
                        gammab(ib)*colin10_sp(:,:,ib), &
                        gammab(ib)*colin01_sp(:,:,ib), &
                        gammab(ib)*colin00_sp(:,:,ib), &
                        fin, Lf_col_sp)
    do iv = 1, nperp
        do imu = 1, npar
            ix = index_mat(iv,imu)
            f1_sp(iv,imu,ib) = Lf_col_sp(ix)*jacob(iv,imu)
        end do
    end do
end do

if (isc /= 0) then
    call apply_operator(sc20, sc02, sc11, sc10, sc01, sc00, fin, Lf_col_sc)
end if

if (irf == -1) then
    allocate(rf00(nperp,npar))
    rf00 = 0.d0
    call apply_operator(rf20, rf02, rf11, rf10, rf01, rf00, fin, Lf_RF)
    deallocate(rf00)
else
    Lf_RF = 0.d0
end if

if (isource == -1) then
    taum = 1.d0/taus
else
    taum = 0.d0
end if

do iv = 1, nperp
    do imu = 1, npar
        ix = index_mat(iv,imu)
        f1(iv,imu) = Lf_col(ix)*jacob(iv,imu)
        f2(iv,imu) = source(iv,imu)*jacob(iv,imu)
        f3(iv,imu) = -fin(ix)*taum*jacob(iv,imu)
        f4(iv,imu) = Lf_RF(ix)*jacob(iv,imu)
        if (isc /= 0) f5(iv,imu) = Lf_col_sc(ix)*jacob(iv,imu)
    end do
end do

call ncint_2D(f1, n1)
do ib = 1, nbulk
    call ncint_2D(f1_sp(:,:,ib), n1_sp(ib))
end do
call ncint_2D(f2, n2)
call ncint_2D(f3, n3)
call ncint_2D(f4, n4)
if (isc /= 0) then
    call ncint_2D(f5, n5)
else
    n5 = 0.d0
end if

write(*,'(/,A)')    'PARTICLE DENSITY BALANCE:'
write(*,'(A)')      '-------------------------'
write(*,'(A)')      ''
write(*,'(A)')      '==========================='
write(*,'(A)')      ' Collisions:'
write(*,'(A)')      ' -----------'
do ib = 1, nbulk
    if (ib == 1) then
        write(*,'(A,ES12.5,A)') &
            '  Collisions with electrons:     ', n1_sp(1)/npart, ' /s'
    else
        write(*,'(A,I2,A,ES12.5,A)') &
            '  Collisions with ions           ', ib-1, ' :  ', n1_sp(ib)/npart, ' /s'
    end if
end do
write(*,'(A)')      ' ------------------'
write(*,'(A)')      ''
if (isc /= 0) then
    write(*,'(A,ES12.5,A)') &
        '  Self-collisions:              ', n5/npart, ' /s'
    write(*,'(A)')      ''
end if
write(*,'(A)')      ' -----------------'
write(*,'(A,ES12.5,A)') &
    '  Total collisions:             ', (n1+n5)/npart, ' /s'
write(*,'(A)')      '==========================='
write(*,'(A)')      ''
write(*,'(A,ES12.5,A)') &
    '  Beam source:                  ', n2/npart, ' /s'
write(*,'(A)')      ''
write(*,'(A,ES12.5,A)') &
    '  Particle losses:              ', n3/npart, ' /s'
write(*,'(A)')      ''
write(*,'(A)')      '==========================='
write(*,'(A)')      ''
write(*,'(A,ES12.5,A)') &
    '  RF term:                      ', n4/npart, ' /s'
write(*,'(A)')      ''
write(*,'(A)')      '  --------------------------------------------'
write(*,'(A,ES12.5,A)') &
    '  Total balance:                ', (n1+n2+n3+n4+n5)/npart, ' /s'
write(*,'(A)')      ''

end subroutine test_density_balance_7pt
