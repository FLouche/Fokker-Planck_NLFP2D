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
use mod_coul_fv, only: coll_apply_species, coll_apply_total, sc_apply, rf_apply

INTEGER, PARAMETER :: dp = KIND(1.0D0)

double precision :: fin(nbig)

real(dp), dimension(nbig)        :: Lf_col_sp, Lf_col, Lf_col_sc, Lf_RF
double precision, dimension(nperp,npar) :: f1, f2, f3, f4, f5
double precision :: f1_sp(nperp,npar,nbulk)
double precision :: n1, n2, n3, n4, n5, n1_sp(nbulk)
double precision, allocatable :: rf00(:,:)
integer :: ix, iv, imu, ib

taum = 0.d0

call coll_apply_total(fin, Lf_col)

do ib = 1, nbulk
    call coll_apply_species(ib, fin, Lf_col_sp)
    do iv = 1, nperp
        do imu = 1, npar
            ix = index_mat(iv,imu)
            f1_sp(iv,imu,ib) = Lf_col_sp(ix)*jacob(iv,imu)
        end do
    end do
end do

if (isc /= 0) then
    call sc_apply(fin, Lf_col_sc)
end if

if (irf == -1) then
    call rf_apply(fin, Lf_RF)   ! 7-point or flux form (i_cc)
!DEBUG-OFF!     ! DEBUG: the whole RF term L_RF[f] on the grid (vperp, vpar, L_RF, L_RF*jacob)
!DEBUG-OFF!     open(125,file='LF_RF_2D.txt',status='unknown')
!DEBUG-OFF!     do iv = 1, nperp
!DEBUG-OFF!         do imu = 1, npar
!DEBUG-OFF!             ix = index_mat(iv,imu)
!DEBUG-OFF!             write(125,'(4ES24.15E3)') vperp(iv), vpar(imu), Lf_RF(ix), Lf_RF(ix)*jacob(iv,imu)
!DEBUG-OFF!         end do
!DEBUG-OFF!     end do
!DEBUG-OFF!     close(125)
!DEBUG-OFF!     ! DEBUG: the RF term applied to the initial Maxwellian fstix (Teff), i.e. the
!DEBUG-OFF!     ! RF part of the first time step for istart = 1 (rf** are fixed in time)
!DEBUG-OFF!     block
!DEBUG-OFF!         real(dp) :: fs(nbig), Lf_s(nbig)
!DEBUG-OFF!         do iv = 1, nperp
!DEBUG-OFF!             do imu = 1, npar
!DEBUG-OFF!                 fs(index_mat(iv,imu)) = fstix(iv,imu)
!DEBUG-OFF!             end do
!DEBUG-OFF!         end do
!DEBUG-OFF!         allocate(rf00(nperp,npar))
!DEBUG-OFF!         rf00 = 0.d0
!DEBUG-OFF!         call apply_operator(rf20, rf02, rf11, rf10, rf01, rf00, fs, Lf_s)
!DEBUG-OFF!         deallocate(rf00)
!DEBUG-OFF!         open(126,file='LF_RF_fstix_2D.txt',status='unknown')
!DEBUG-OFF!         do iv = 1, nperp
!DEBUG-OFF!             do imu = 1, npar
!DEBUG-OFF!                 ix = index_mat(iv,imu)
!DEBUG-OFF!                 write(126,'(5ES24.15E3)') vperp(iv), vpar(imu), Lf_s(ix), &
!DEBUG-OFF!                                           Lf_s(ix)*jacob(iv,imu), fs(ix)
!DEBUG-OFF!             end do
!DEBUG-OFF!         end do
!DEBUG-OFF!         close(126)
!DEBUG-OFF!         ! DEBUG: the collision operator on the same Maxwellian: total (all
!DEBUG-OFF!         ! backgrounds, colin**) and per background ib = 1 (e), 2 (first ion)
!DEBUG-OFF!         ! columns: vperp, vpar, L_col, L_col*jacob, L_e, L_ion1
!DEBUG-OFF!         block
!DEBUG-OFF!             real(dp) :: Lc(nbig), Lb(nbig, 2)
!DEBUG-OFF!             integer :: kb
!DEBUG-OFF!             call apply_operator(colin20, colin02, colin11, colin10, colin01, colin00, fs, Lc)
!DEBUG-OFF!             Lb = 0.d0
!DEBUG-OFF!             do kb = 1, min(nbulk, 2)
!DEBUG-OFF!                 call apply_operator(gammab(kb)*colin20_sp(:,:,kb), gammab(kb)*colin02_sp(:,:,kb), &
!DEBUG-OFF!                                     gammab(kb)*colin11_sp(:,:,kb), gammab(kb)*colin10_sp(:,:,kb), &
!DEBUG-OFF!                                     gammab(kb)*colin01_sp(:,:,kb), gammab(kb)*colin00_sp(:,:,kb), &
!DEBUG-OFF!                                     fs, Lb(:,kb))
!DEBUG-OFF!             end do
!DEBUG-OFF!             open(127,file='LF_col_fstix_2D.txt',status='unknown')
!DEBUG-OFF!             do iv = 1, nperp
!DEBUG-OFF!                 do imu = 1, npar
!DEBUG-OFF!                     ix = index_mat(iv,imu)
!DEBUG-OFF!                     write(127,'(6ES24.15E3)') vperp(iv), vpar(imu), Lc(ix), &
!DEBUG-OFF!                                               Lc(ix)*jacob(iv,imu), Lb(ix,1), Lb(ix,2)
!DEBUG-OFF!                 end do
!DEBUG-OFF!             end do
!DEBUG-OFF!             close(127)
!DEBUG-OFF!         end block
!DEBUG-OFF!     end block
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

!DEBUG-OFF! ! DEBUG: the collision operator on the solution fin, total and per background
!DEBUG-OFF! ! (ib = 1 e, 2 first ion; f1_sp holds L_ib*jacob) -- columns:
!DEBUG-OFF! ! vperp, vpar, L_col, L_col*jacob, L_e, L_ion1, fin
!DEBUG-OFF! open(128,file='LF_col_2D.txt',status='unknown')
!DEBUG-OFF! do iv = 1, nperp
!DEBUG-OFF!     do imu = 1, npar
!DEBUG-OFF!         ix = index_mat(iv,imu)
!DEBUG-OFF!         write(128,'(7ES24.15E3)') vperp(iv), vpar(imu), Lf_col(ix), f1(iv,imu), &
!DEBUG-OFF!             f1_sp(iv,imu,1)/jacob(iv,imu), f1_sp(iv,imu,min(2,nbulk))/jacob(iv,imu), fin(ix)
!DEBUG-OFF!     end do
!DEBUG-OFF! end do
!DEBUG-OFF! close(128)

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
