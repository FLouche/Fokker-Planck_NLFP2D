!   test_momentum_balance_7pt.f90
!
!   FPColl_2D
!
!   Created by Fabrice Louche on 29/05/2026.
!   Copyright 2026 LPP-ERM/KMS. All rights reserved.

!***************************************
!  Computes the perpendicular and parallel momentum transfer rate densities
!  (N/m**3) for each collision operator and the RF operator.
!
!  dPperp/dt = integral( ma * vperp * [df/dt]_op ) dV
!  dPpar/dt  = integral( ma * vpar  * [df/dt]_op ) dV
!***************************************

subroutine test_momentum_balance_7pt(fin)

use shared_grid
use shared_plasma
use shared_beam
use shared_rf
use shared_timer
use shared_FPterms
use mod_ncint
use func_index
use mod_apply_operator

implicit none

INTEGER, PARAMETER :: dp = KIND(1.0D0)

double precision, intent(in) :: fin(nbig)

real(dp), dimension(nbig)      :: Lf
real(dp), dimension(nperp,npar) :: fint_perp, fint_par
real(dp), allocatable          :: rf00(:,:)

real(dp), dimension(nbulk) :: mcoll_perp, mcoll_par
real(dp) :: mRF_perp,  mRF_par
real(dp) :: msrc_perp, msrc_par
real(dp) :: mloss_perp, mloss_par
real(dp) :: mSC_perp,  mSC_par

real(dp), parameter :: pmass = 1.6726d-27   ! proton mass [kg]
real(dp) :: taum_save, mom_fac

integer :: iv, ip, ix, ib

! Save taum; collisional operators do not include the loss term
taum_save = taum
taum      = 0.0_dp

! Pre-compute the mass factor (constant over the grid)
mom_fac = pmass * aa

!================================================================
! 1.  Bulk-species collisional momentum transfer
!================================================================
do ib = 1, nbulk

    call apply_operator(gammab(ib)*colin20_sp(:,:,ib), &
                        gammab(ib)*colin02_sp(:,:,ib), &
                        gammab(ib)*colin11_sp(:,:,ib), &
                        gammab(ib)*colin10_sp(:,:,ib), &
                        gammab(ib)*colin01_sp(:,:,ib), &
                        gammab(ib)*colin00_sp(:,:,ib), &
                        fin, Lf)

    do iv = 1, nperp
        do ip = 1, npar
            ix = index_mat(iv,ip)
            fint_perp(iv,ip) = mom_fac * vperp(iv) * Lf(ix) * jacob(iv,ip)
            fint_par (iv,ip) = mom_fac * vpar (ip) * Lf(ix) * jacob(iv,ip)
        end do
    end do

    call ncint_2d(fint_perp, mcoll_perp(ib))
    call ncint_2d(fint_par,  mcoll_par (ib))

end do

!================================================================
! 2.  Self-collision momentum transfer
!================================================================
if (isc /= 0) then

    call apply_operator(sc20, sc02, sc11, sc10, sc01, sc00, fin, Lf)

    do iv = 1, nperp
        do ip = 1, npar
            ix = index_mat(iv,ip)
            fint_perp(iv,ip) = mom_fac * vperp(iv) * Lf(ix) * jacob(iv,ip)
            fint_par (iv,ip) = mom_fac * vpar (ip) * Lf(ix) * jacob(iv,ip)
        end do
    end do

    call ncint_2d(fint_perp, mSC_perp)
    call ncint_2d(fint_par,  mSC_par)

else
    mSC_perp = 0.0_dp
    mSC_par  = 0.0_dp
end if

!================================================================
! 3.  RF momentum transfer
!================================================================
if (irf == -1) then

    allocate(rf00(nperp,npar))
    rf00 = 0.0_dp
    call apply_operator(rf20, rf02, rf11, rf10, rf01, rf00, fin, Lf)
    deallocate(rf00)

    do iv = 1, nperp
        do ip = 1, npar
            ix = index_mat(iv,ip)
            fint_perp(iv,ip) = mom_fac * vperp(iv) * Lf(ix) * jacob(iv,ip)
            fint_par (iv,ip) = mom_fac * vpar (ip) * Lf(ix) * jacob(iv,ip)
        end do
    end do

    call ncint_2d(fint_perp, mRF_perp)
    call ncint_2d(fint_par,  mRF_par)

else
    mRF_perp = 0.0_dp
    mRF_par  = 0.0_dp
end if

!================================================================
! 4.  Beam source and particle losses
!================================================================
if (isource == -1) then

    taum = 1.0_dp / taus

    do iv = 1, nperp
        do ip = 1, npar
            ix = index_mat(iv,ip)
            fint_perp(iv,ip) = mom_fac * vperp(iv) * (-fin(ix)*taum) * jacob(iv,ip)
            fint_par (iv,ip) = mom_fac * vpar (ip) * (-fin(ix)*taum) * jacob(iv,ip)
        end do
    end do
    call ncint_2d(fint_perp, mloss_perp)
    call ncint_2d(fint_par,  mloss_par)

    do iv = 1, nperp
        do ip = 1, npar
            fint_perp(iv,ip) = mom_fac * vperp(iv) * source(iv,ip) * jacob(iv,ip)
            fint_par (iv,ip) = mom_fac * vpar (ip) * source(iv,ip) * jacob(iv,ip)
        end do
    end do
    call ncint_2d(fint_perp, msrc_perp)
    call ncint_2d(fint_par,  msrc_par)

else
    mloss_perp = 0.0_dp;  mloss_par = 0.0_dp
    msrc_perp  = 0.0_dp;  msrc_par  = 0.0_dp
end if

taum = taum_save

!================================================================
! 5.  Print results
!================================================================

write(*,*) ''
write(*,*) 'Perpendicular momentum balance:'
write(*,*) '-----------------'
write(*,*) 'Collisions:'
write(*,*) '-----------'
do ib = 1, nbulk
    if (ib == 1) then
        write(*,*) 'Collisions with electrons:      ', mcoll_perp(ib), 'N/m**3'
    else
        write(*,*) 'Collisions with ions ', ib-1, ': ', mcoll_perp(ib), 'N/m**3'
    end if
end do
if (isc /= 0) write(*,*) 'Self-collisions:                ', mSC_perp, 'N/m**3'
write(*,*) '-----------------'
write(*,*) 'Total collisions:               ', (SUM(mcoll_perp)+mSC_perp), 'N/m**3'
write(*,*) 'Beam source:                    ', msrc_perp,  'N/m**3'
write(*,*) 'Particle losses:                ', mloss_perp, 'N/m**3'
write(*,*) ''
write(*,*) 'RF term:                        ', mRF_perp,   'N/m**3'
write(*,*) '--------------------------------------------'
write(*,*) 'Total          :                ', &
           (SUM(mcoll_perp)+mSC_perp+mloss_perp+msrc_perp+mRF_perp), 'N/m**3'

write(*,*) ''
write(*,*) 'Parallel momentum balance:'
write(*,*) '-----------------'
write(*,*) 'Collisions:'
write(*,*) '-----------'
do ib = 1, nbulk
    if (ib == 1) then
        write(*,*) 'Collisions with electrons:      ', mcoll_par(ib), 'N/m**3'
    else
        write(*,*) 'Collisions with ions ', ib-1, ': ', mcoll_par(ib), 'N/m**3'
    end if
end do
if (isc /= 0) write(*,*) 'Self-collisions:                ', mSC_par, 'N/m**3'
write(*,*) '-----------------'
write(*,*) 'Total collisions:               ', (SUM(mcoll_par)+mSC_par), 'N/m**3'
write(*,*) 'Beam source:                    ', msrc_par,  'N/m**3'
write(*,*) 'Particle losses:                ', mloss_par, 'N/m**3'
write(*,*) ''
write(*,*) 'RF term:                        ', mRF_par,   'N/m**3'
write(*,*) '--------------------------------------------'
write(*,*) 'Total          :                ', &
           (SUM(mcoll_par)+mSC_par+mloss_par+msrc_par+mRF_par), 'N/m**3'

end subroutine test_momentum_balance_7pt
