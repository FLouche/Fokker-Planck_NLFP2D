!***********************************************
!* Computation of the density at each timestep *
!***********************************************
!
! Version 2: adapted for general grids 
! 25/03/2026: F. Louche

subroutine time_density(f,dens)

use shared_grid
!use integrate_2d_module

implicit none

double precision, intent(in) :: f(nperp,npar)
double precision, intent(out) :: dens

common/mathcons/pi,twopi

double precision:: dvp_local, pi,twopi
integer iv,imu
 
!dens = integrate_2d(f*jacob, vperp, vpar, nperp, npar)

! Compute 2π ∫∫ f * vperp * dvperp * dvpar
dens = 0.d0
DO iv = 1, nperp-1
    dvp_local = vperp(iv+1) - vperp(iv)   ! works for any grid
    DO imu = 2, npar-1                                ! skip Dirichlet boundaries
        dens = dens + twopi * f(iv,imu) * vperp(iv) &
                            * dvp_local * dvpar
    END DO
END DO


end subroutine time_density

!==============================================