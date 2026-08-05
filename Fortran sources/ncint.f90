!
!   ncint.f
!   nlfp2011
!   
!   Created by Fabrice Louche on 20/01/11.
!   Copyright 2011 LPP-ERM/KMS. All rights reserved.
!   
!
!***********************************************
!*          NEWTON-COTES FORMULAS              *
!* on a grid of arbitrary distributed points   *
!*                             *
!***********************************************

module mod_ncint

contains

subroutine ncint_2D(f,s)

! Initialisation
! -------------- 

use shared_grid
       
implicit none

double precision, intent(in) :: f(nperp,npar)
double precision, intent(out):: s

double precision dvp_local
integer i,j

! ===========================================================================

! F is the values of the function to be integrated on the complete grid
! s  : output value of the integral

s=0.d0

! Compute 2π ∫∫ f * vperp * dvperp * dvpar
!
! v⊥ : trapezoidal rule on an arbitrary (possibly non-uniform) grid — each node
!      carries half the span of its two adjacent intervals.
!
!      This replaces a left-endpoint RECTANGLE rule,
!          do i = 1, nperp-1;  s = s + f(i,j)*(vperp(i+1)-vperp(i))*dvpar
!      which is only 1st-order accurate when the v⊥ grid is non-uniform.  On a
!      UNIFORM grid, with the integrand vanishing at both ends, the rectangle
!      sum coincides with the trapezoidal sum, so the defect was invisible for
!      ising=0 and appeared only for ising=±1.
!
!      NOTE: time_density (time_comps_mod.f90) does NOT call this routine — it
!      carries its own copy of the same quadrature, and that copy is the one in
!      the renormalisation path (fout = fout*npart/dens_tmp).  Both were fixed
!      together and must be kept in step.
!
!      The bias this removes was measured against an accurate integration of
!      the same fstix: 3.8% (ising=1) and 5.9% (ising=-1) at nperp=31, falling
!      only as O(h) (1.5% at nperp=121 on the two-domain grid); ising=0 is
!      unaffected at the 1e-6 level.
!
!      This routine feeds the analysis moments, the beam-source density, the
!      fstix renormalisation in consts.f90, and the power/momentum diagnostics.
!
! v∥ : uniform grid; the j=2..npar-1 sum already IS the trapezoidal rule,
!      because f vanishes on the two Dirichlet boundaries j=1 and j=npar.

DO i = 1, nperp

    IF (i == 1) THEN
        dvp_local = 0.5d0*(vperp(2) - vperp(1))
    ELSE IF (i == nperp) THEN
        dvp_local = 0.5d0*(vperp(nperp) - vperp(nperp-1))
    ELSE
        dvp_local = 0.5d0*(vperp(i+1) - vperp(i-1))
    END IF

    DO j = 2, npar-1                                ! skip Dirichlet boundaries
        s = s + f(i,j) * dvp_local * dvpar
    END DO

END DO

!
!============================================================================
!
end subroutine ncint_2D
	   


!
!============================================================================

end module mod_ncint