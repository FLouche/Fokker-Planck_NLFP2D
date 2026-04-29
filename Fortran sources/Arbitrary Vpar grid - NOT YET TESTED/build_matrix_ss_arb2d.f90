!**********************************************************
!
! This routine builds the matrix containing the terms
!  of the FP equation written on a 2D grid using a
!   finite-differences differentiation scheme
!
!   This version includes:
!
!    - arbitrary (non-uniform) grid in vperp AND vpar
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
!   v∥ direction (non-uniform, g_m = vpa(j)-vpa(j-1),
!                               g_p = vpa(j+1)-vpa(j)):
!
!     df/dva:   w_l = -g_p/(g_m*(g_m+g_p))
!               w_n = (g_p-g_m)/(g_m*g_p)
!               w_r =  g_m/(g_p*(g_m+g_p))
!
!     d2f/dva2: w_l =  2/(g_m*(g_m+g_p))
!               w_n = -2/(g_m*g_p)
!               w_r =  2/(g_p*(g_m+g_p))
!
!   Mixed deriv: outer product of v⊥ and v∥ 1st-deriv weights
!
!   Boundary conditions:
!     iv=nperp, ip=1, ip=npar : f = 0  (Dirichlet)
!     iv=1                    : df/dvp = 0  (Neumann, 2nd-order one-sided)
!
!   Solvability constraint (isource==0):
!     f(imid,jmid) = 1
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

! v⊥ finite-difference weights (non-uniform)
double precision :: alpha_l, alpha_r, alpha_n   ! 1st derivative
double precision :: beta_l,  beta_r,  beta_n    ! 2nd derivative
double precision :: h_m, h_p                    ! left/right v⊥ spacings

! v∥ finite-difference weights (non-uniform, computed per point)
double precision :: gam_l, gam_r, gam_n         ! 1st derivative
double precision :: del_l, del_r, del_n          ! 2nd derivative
double precision :: g_m, g_p                    ! left/right v∥ spacings

integer :: iv, ip, ix1, ix2

!====================================================================
!
! Initialise full matrix to zero
!
bigm = 0.d0
if (present(bigv)) bigv = 0.d0

!====================================================================

imid = int(2*nperp/5)+1

v_loop: do iv = 1, nperp

    !----------------------------------------------------------------
    ! v⊥ weights at point iv
    !----------------------------------------------------------------
    if (iv > 1 .and. iv < nperp) then

        h_m = vperp(iv)   - vperp(iv-1)
        h_p = vperp(iv+1) - vperp(iv)

        alpha_l = -h_p         / (h_m*(h_m+h_p))
        alpha_n =  (h_p - h_m) / (h_m*h_p)
        alpha_r =  h_m         / (h_p*(h_m+h_p))

        beta_l  =  2.d0 / (h_m*(h_m+h_p))
        beta_n  = -2.d0 / (h_m*h_p)
        beta_r  =  2.d0 / (h_p*(h_m+h_p))

    else if (iv == 1) then

        ! One-sided (used only in the Neumann BC row)
        h_p     = vperp(2) - vperp(1)
        alpha_l =  0.d0
        alpha_n = -3.d0 / (2.d0*h_p)
        alpha_r =  4.d0 / (2.d0*h_p)
        beta_l  =  0.d0;  beta_n = 0.d0;  beta_r = 0.d0

    else   ! iv == nperp: Dirichlet, not used

        alpha_l = 0.d0;  alpha_n = 0.d0;  alpha_r = 0.d0
        beta_l  = 0.d0;  beta_n  = 0.d0;  beta_r  = 0.d0

    end if

    mu_loop: do ip = 1, npar

        ix1 = index_mat(iv, ip)

        !============================================================
        ! BC: Dirichlet at outer v⊥ and v∥ boundaries
        !============================================================
        BC_vmax: if (iv == nperp .or. ip == 1 .or. ip == npar) then

            bigm(ix1,ix1) = 1.d0
            if (present(bigv)) bigv(ix1) = 0.d0

        !============================================================
        ! BC: Neumann at vperp=0 (iv=1): df/dvp = 0
        !     (-3*f(1) + 4*f(2) - f(3)) / (2*h_p) = 0
        !     (hardcoded scaling to match build_ss convention)
        !============================================================
        else BC_vmax

            BC_vmin: if (iv == 1) then

                bigm(ix1, ix1)              = -1.5d0
                ix2 = index_mat(iv+1, ip);  bigm(ix1, ix2) =  2.0d0
                ix2 = index_mat(iv+2, ip);  bigm(ix1, ix2) = -0.5d0
                if (present(bigv)) bigv(ix1) = 0.d0

            !========================================================
            ! Interior points: full PDE stencil
            !
            ! v∥ weights computed locally from vpar(ip-1), vpar(ip),
            ! vpar(ip+1) — valid for any non-uniform vpar grid.
            !========================================================
            else BC_vmin

                !--- v∥ weights at point ip (non-uniform) -----------
                g_m = vpar(ip)   - vpar(ip-1)
                g_p = vpar(ip+1) - vpar(ip)

                gam_l = -g_p         / (g_m*(g_m+g_p))   ! w(ip-1), df/dva
                gam_n =  (g_p - g_m) / (g_m*g_p)          ! w(ip),   df/dva
                gam_r =  g_m         / (g_p*(g_m+g_p))   ! w(ip+1), df/dva

                del_l =  2.d0 / (g_m*(g_m+g_p))           ! w(ip-1), d2f/dva2
                del_n = -2.d0 / (g_m*g_p)                 ! w(ip),   d2f/dva2
                del_r =  2.d0 / (g_p*(g_m+g_p))           ! w(ip+1), d2f/dva2

                !--- (iv, ip): centre ---
                bigm(ix1,ix1) = all00(iv,ip)            &   ! A
                              + alpha_n*all10(iv,ip)     &   ! B, centre
                              + beta_n *all20(iv,ip)     &   ! D, centre
                              + gam_n  *all01(iv,ip)     &   ! C, centre
                              + del_n  *all02(iv,ip)     &   ! F, centre
                              - taum                         ! damping

                !--- (iv+1, ip): right v⊥, same v∥ ---
                ix2 = index_mat(iv+1, ip)
                bigm(ix1,ix2) = alpha_r*all10(iv,ip)    &   ! B
                              + beta_r *all20(iv,ip)         ! D

                !--- (iv-1, ip): left v⊥, same v∥ ---
                ix2 = index_mat(iv-1, ip)
                bigm(ix1,ix2) = alpha_l*all10(iv,ip)    &   ! B
                              + beta_l *all20(iv,ip)         ! D

                !--- (iv, ip+1): same v⊥, right v∥ ---
                ix2 = index_mat(iv, ip+1)
                bigm(ix1,ix2) = gam_r*all01(iv,ip)      &   ! C
                              + del_r*all02(iv,ip)           ! F

                !--- (iv, ip-1): same v⊥, left v∥ ---
                ix2 = index_mat(iv, ip-1)
                bigm(ix1,ix2) = gam_l*all01(iv,ip)      &   ! C
                              + del_l*all02(iv,ip)           ! F

                !--- Mixed derivative E: corners ---
                ! d2f/dvp_dva = alpha(di) * gam(dj)
                ! Corners (di!=0, dj!=0):
                ix2 = index_mat(iv+1, ip+1)
                bigm(ix1,ix2) = alpha_r*gam_r*all11(iv,ip)

                ix2 = index_mat(iv+1, ip-1)
                bigm(ix1,ix2) = alpha_r*gam_l*all11(iv,ip)

                ix2 = index_mat(iv-1, ip+1)
                bigm(ix1,ix2) = alpha_l*gam_r*all11(iv,ip)

                ix2 = index_mat(iv-1, ip-1)
                bigm(ix1,ix2) = alpha_l*gam_l*all11(iv,ip)

                ! E at (iv, ip+/-1): alpha_n * gam_r/gam_l
                ! (non-zero when alpha_n != 0, i.e. non-uniform v⊥ grid)
                ix2 = index_mat(iv, ip+1)
                bigm(ix1,ix2) = bigm(ix1,ix2) + alpha_n*gam_r*all11(iv,ip)

                ix2 = index_mat(iv, ip-1)
                bigm(ix1,ix2) = bigm(ix1,ix2) + alpha_n*gam_l*all11(iv,ip)

                ! E at (iv+/-1, ip): alpha_r/alpha_l * gam_n
                ! (non-zero when gam_n != 0, i.e. non-uniform v∥ grid)
                ix2 = index_mat(iv+1, ip)
                bigm(ix1,ix2) = bigm(ix1,ix2) + alpha_r*gam_n*all11(iv,ip)

                ix2 = index_mat(iv-1, ip)
                bigm(ix1,ix2) = bigm(ix1,ix2) + alpha_l*gam_n*all11(iv,ip)

                !--- RHS ---
                if (present(bigv)) bigv(ix1) = -source(iv,ip)

            end if BC_vmin

        end if BC_vmax

    end do mu_loop

end do v_loop

!====================================================================
! Solvability constraint
!====================================================================
source_term: if (isource == 0 .AND. PRESENT(bigv)) then

    ix1 = index_mat(imid, jmid)
    bigm(ix1,:)   = 0.d0
    bigm(ix1,ix1) = 1.d0
    bigv(ix1)     = 1.d0
    write(*,*) 'build_ss: imposing f=1 at (imid,jmid) = ', imid, jmid

end if source_term

!====================================================================

end subroutine build_ss

end module mod_build_ss
