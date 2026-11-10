!**********************************************************
!
! This routine builds the matrix containing the terms
!  of the FP equation written on a 2D grid using a
!   finite-differences differentiation scheme
!
!   This version includes:
!
!    - arbitrary (non-uniform) grid in vperp
!    - uniform grid in vpar
!    - beam source
!    - steady-state equation
!
!   Finite-difference scheme (3-point, 2nd order):
!
!   v⊥ direction (non-uniform, h_m = vp(i)-vp(i-1),
!                               h_p = vp(i+1)-vp(i)):
!
!     df/dvp:   w_l = -h_p/(h_m*(h_m+h_p))
!               w_n = (h_p-h_m)/(h_m*h_p)
!               w_r =  h_m/(h_p*(h_m+h_p))
!
!     d2f/dvp2: w_l =  2/(h_m*(h_m+h_p))
!               w_n = -2/(h_m*h_p)
!               w_r =  2/(h_p*(h_m+h_p))
!
!   v∥ direction (uniform spacing dvpar = dmu):
!
!     df/dva:   [-f(j-1) + f(j+1)] / (2*dmu)
!     d2f/dva2: [f(j-1) - 2*f(j) + f(j+1)] / dmu^2
!
!   Mixed deriv (outer product of 3-pt 1st-deriv weights):
!     d2f/dvp_dva = w_vp(di) * w_vpa(dj)
!
!   Boundary conditions:
!     iv=nperp, ip=1, ip=npar : f = 0  (Dirichlet)
!     iv=1                    : df/dvp = 0  (Neumann, 2nd-order one-sided)
!
!   Solvability constraint (isource==0):
!     f(imid,jmid) = 1  (fixes normalisation of null-space solution)
!
!    by Fabrice Louche, extended for arbitrary vperp grid
!
!**********************************************************

module mod_build_ss

    contains

subroutine build_ss(all20,all02,all11,all10,all01,all00,bigm,bigv)

! Initialisation
! --------------

use shared_grid
use shared_plasma
use shared_beam

use shared_timer
use shared_rf

use func_index

implicit none

double precision, dimension(nperp,npar), intent(in) :: &
    all20, all02, all11, all10, all01, all00
double precision, dimension(nbig,nbig),  intent(out) :: bigm
double precision, dimension(nbig),       intent(out), optional :: bigv

! Local variables
! ---------------

! v⊥ finite-difference weights for current point iv
! offsets: _l = left (iv-1), _n = centre (iv), _r = right (iv+1)
double precision :: alpha_l, alpha_r, alpha_n   ! 1st derivative weights
double precision :: beta_l,  beta_r,  beta_n    ! 2nd derivative weights
double precision :: h_m, h_p                    ! left and right spacings

! v∥ weights (uniform grid, stored as scalars for clarity)
double precision :: wpa1_l, wpa1_r              ! df/dvpa: -1/(2dmu), +1/(2dmu)
double precision :: wpa2_l, wpa2_n, wpa2_r      ! d2f/dvpa2

integer :: iv, ip, ix1, ix2

!====================================================================
!
! Pre-compute uniform v∥ weights (same for all interior points)
!
wpa1_l =  -1.d0 / (2.d0*dvpar)
wpa1_r =  +1.d0 / (2.d0*dvpar)

wpa2_l =   1.d0 / dmu2
wpa2_n =  -2.d0 / dmu2
wpa2_r =   1.d0 / dmu2

!====================================================================
!
! Initialise full matrix to zero
! (sparse assembly: only non-zero entries are set below)
!
bigm = 0.d0
if (present(bigv)) bigv = 0.d0

!====================================================================

v_loop: do iv = 1, nperp

    !----------------------------------------------------------------
    ! Compute v⊥ finite-difference weights for point iv
    ! using the local (possibly non-uniform) grid spacings.
    !
    ! For interior points: use the 3-point non-uniform formula.
    ! For iv=1 and iv=nperp: these rows are handled as BCs below,
    ! but we still compute weights to avoid uninitialized values.
    !----------------------------------------------------------------

    if (iv > 1 .and. iv < nperp) then

        h_m = vperp(iv)   - vperp(iv-1)
        h_p = vperp(iv+1) - vperp(iv)

        ! 1st derivative weights
        alpha_l = -h_p          / (h_m*(h_m+h_p))   ! w(iv-1)
        alpha_n =  (h_p - h_m)  / (h_m*h_p)          ! w(iv)
        alpha_r =  h_m          / (h_p*(h_m+h_p))   ! w(iv+1)

        ! 2nd derivative weights
        beta_l  =  2.d0 / (h_m*(h_m+h_p))            ! w(iv-1)
        beta_n  = -2.d0 / (h_m*h_p)                  ! w(iv)
        beta_r  =  2.d0 / (h_p*(h_m+h_p))            ! w(iv+1)

    else if (iv == 1) then

        ! One-sided weights at iv=1 (used only for the Neumann row)
        h_p     = vperp(2) - vperp(1)
        alpha_l =  0.d0
        alpha_n = -3.d0 / (2.d0*h_p)
        alpha_r =  4.d0 / (2.d0*h_p)
        ! Note: the 3-point 2nd-order forward 1st-deriv also needs iv+2;
        ! since the Neumann row only uses alpha_n, alpha_r and
        ! the coefficient of iv+2 (= -1/(2*h_p)), we handle it
        ! explicitly in the BC block below.
        beta_l  =  0.d0
        beta_n  =  0.d0
        beta_r  =  0.d0

    else  ! iv == nperp: Dirichlet, weights not used

        alpha_l = 0.d0;  alpha_n = 0.d0;  alpha_r = 0.d0
        beta_l  = 0.d0;  beta_n  = 0.d0;  beta_r  = 0.d0

    end if

    mu_loop: do ip = 1, npar

        ix1 = index_mat(iv, ip)

        !============================================================
        ! BC: Dirichlet at outer v⊥ boundary and v∥ boundaries
        !     f = 0  -> row is simply  1 * f(iv,ip) = 0
        !============================================================
        BC_vmax: if (iv == nperp .or. ip == 1 .or. ip == npar) then

            bigm(ix1,ix1) = 1.d0
            if (present(bigv)) bigv(ix1) = 0.d0

        !============================================================
        ! BC: Neumann at vperp=0 (iv=1)
        !     df/dvp = 0  -> 2nd-order one-sided stencil:
        !     (-3*f(1) + 4*f(2) - f(3)) / (2*h_p) = 0
        !============================================================
        else BC_vmax

            BC_vmin: if (iv == 1) then

                h_p = vperp(2) - vperp(1)
                
                bigm(ix1, ix1) = -1.5d0   ! f(1,ip)
                
                ix2 = index_mat(iv+1, ip)
                bigm(ix1, ix2) =  2.d0   ! f(2,ip)
                
                ix2 = index_mat(iv+2, ip)
                bigm(ix1, ix2) = -0.5d0   ! f(3,ip)
                
                if (present(bigv)) bigv(ix1) = 0.d0

    

            !========================================================
            ! Interior points: assemble the full PDE stencil
            !
            !   A*f + B*df/dvp + C*df/dva
            !       + D*d2f/dvp2 + E*d2f/dvp_dva + F*d2f/dva2 = rhs
            !
            ! The 9 possible stencil entries are:
            !
            !   (iv-1, ip-1)  (iv-1, ip)  (iv-1, ip+1)
            !   (iv  , ip-1)  (iv  , ip)  (iv  , ip+1)
            !   (iv+1, ip-1)  (iv+1, ip)  (iv+1, ip+1)
            !
            ! Contributions per term:
            !
            !   A:    -> (iv,   ip  ) only
            !   B:    -> (iv-1, ip  ), (iv,   ip  ), (iv+1, ip  )
            !   C:    -> (iv,   ip-1), (iv,   ip  ), (iv,   ip+1)  [centre=0]
            !   D:    -> (iv-1, ip  ), (iv,   ip  ), (iv+1, ip  )
            !   E:    -> outer product of B-weights and C-weights
            !            (all 4 off-diagonal corners, centre=0)
            !   F:    -> (iv,   ip-1), (iv,   ip  ), (iv,   ip+1)
            !========================================================
            else BC_vmin

                !--- (iv, ip): centre ---
                bigm(ix1,ix1) = all00(iv,ip)          &   ! A
                              + alpha_n*all10(iv,ip)   &   ! B, centre
                              + beta_n *all20(iv,ip)   &   ! D, centre
                              + wpa2_n *all02(iv,ip)   &   ! F, centre
                              - taum
                ! Note: C and E have zero weight at centre offset (j=0)

                !--- (iv+1, ip): right in v⊥, same v∥ ---
                ix2 = index_mat(iv+1, ip)
                bigm(ix1,ix2) = alpha_r*all10(iv,ip)  &   ! B
                              + beta_r *all20(iv,ip)       ! D

                !--- (iv-1, ip): left in v⊥, same v∥ ---
                ix2 = index_mat(iv-1, ip)
                bigm(ix1,ix2) = alpha_l*all10(iv,ip)  &   ! B
                              + beta_l *all20(iv,ip)       ! D

                !--- (iv, ip+1): same v⊥, right in v∥ ---
                ix2 = index_mat(iv, ip+1)
                bigm(ix1,ix2) = wpa1_r*all01(iv,ip)   &   ! C
                              + wpa2_r*all02(iv,ip)        ! F

                !--- (iv, ip-1): same v⊥, left in v∥ ---
                ix2 = index_mat(iv, ip-1)
                bigm(ix1,ix2) = wpa1_l*all01(iv,ip)   &   ! C
                              + wpa2_l*all02(iv,ip)        ! F

                !--- Mixed derivative E: 4 corner points ---
                !    d2f/dvp_dva = alpha(di) * wpa1(dj)

                ix2 = index_mat(iv+1, ip+1)            ! (+1,+1)
                bigm(ix1,ix2) = alpha_r*wpa1_r*all11(iv,ip)

                ix2 = index_mat(iv+1, ip-1)            ! (+1,-1)
                bigm(ix1,ix2) = alpha_r*wpa1_l*all11(iv,ip)

                ix2 = index_mat(iv-1, ip+1)            ! (-1,+1)
                bigm(ix1,ix2) = alpha_l*wpa1_r*all11(iv,ip)

                ix2 = index_mat(iv-1, ip-1)            ! (-1,-1)
                bigm(ix1,ix2) = alpha_l*wpa1_l*all11(iv,ip)

                ! Note: E contributes zero at (iv,ip+/-1) and (iv+/-1,ip)
                ! because either alpha_n=0 or wpa1 centre=0 would be needed;
                ! alpha_n is NOT zero in general (non-uniform grid), but
                ! the v∥ 1st-deriv weight at dj=0 IS zero (centred scheme).
                ! Similarly, the v⊥ weight at di=0 multiplies wpa1 which
                ! is non-zero only at dj=+/-1, so (iv, ip+/-1) already
                ! covered above; but alpha_n*wpa1(dj=/=0) gives a
                ! contribution to (iv, ip+/-1) from the E term:
                ix2 = index_mat(iv, ip+1)
                bigm(ix1,ix2) = bigm(ix1,ix2) + alpha_n*wpa1_r*all11(iv,ip)

                ix2 = index_mat(iv, ip-1)
                bigm(ix1,ix2) = bigm(ix1,ix2) + alpha_n*wpa1_l*all11(iv,ip)

                !--- RHS ---
                if (present(bigv)) bigv(ix1) = -source(iv,ip)

            end if BC_vmin

        end if BC_vmax

    end do mu_loop

end do v_loop

!====================================================================
!
! Solvability constraint: sourceless steady-state case
! Impose f(imid,jmid) = 1 to fix the normalisation
!
!====================================================================
source_term: if (isource == 0 .AND. PRESENT(bigv)) then

    ix1 = index_mat(imid, jmid)

    bigm(ix1,:) = 0.d0          ! clear entire row
    bigm(ix1,ix1) = 1.d0        ! f(imid,jmid) = 1
    bigv(ix1) = 1.d0

    write(*,*) 'build_ss: imposing f=1 at (imid,jmid) = ', imid, jmid

end if source_term

! FOR DIAGNOSTIC: save non-zero elements

!open(400,file='bigm_nz.txt',status='unknown')
!do ix1=1,nbig
!    do ix2=1,nbig
!        if(bigm(ix1,ix2) /= 0.d0) write(400,*) ix1,ix2,bigm(ix1,ix2)
!    enddo
!enddo
!close(400)

!====================================================================

end subroutine build_ss

!**********************************************************

end module mod_build_ss