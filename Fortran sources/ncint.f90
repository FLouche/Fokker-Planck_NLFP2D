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

!===============================================================================
! SIMPSON_WEIGHTS
!
! Composite Simpson weights w(1:n) for an arbitrary strictly-increasing grid
! x(1:n), such that  integral f dx  ~=  sum_i w(i)*f(i).
!
! Each consecutive triple (x_i, x_i+1, x_i+2) is integrated exactly for the
! quadratic through those three points:
!
!   int_{x0}^{x2} p dx = (h0+h1)/6 * [ (2 - h1/h0) f0
!                                    + (h0+h1)^2/(h0*h1) f1
!                                    + (2 - h0/h1) f2 ]
!
! which reduces to the familiar h/3*(f0 + 4 f1 + f2) when h0 = h1.  The rule is
! 4th-order accurate, against 2nd order for the trapezoidal rule it replaces.
!
! If n is even one interval is left over; it is closed with a trapezoidal panel,
! which costs 2nd-order accuracy on that single interval only.  Use an odd
! number of grid points to keep the rule 4th order everywhere.
!===============================================================================
  subroutine simpson_weights(x, n, w)

    implicit none

    integer, intent(in)           :: n
    double precision, intent(in)  :: x(n)
    double precision, intent(out) :: w(n)

    double precision :: h0, h1, hs
    integer :: i

    w = 0.d0

    if (n < 2) return

    if (n == 2) then                       ! nothing better available
      w(1) = 0.5d0*(x(2)-x(1))
      w(2) = w(1)
      return
    end if

    i = 1
    do while (i + 2 <= n)
      h0 = x(i+1) - x(i)
      h1 = x(i+2) - x(i+1)
      hs = h0 + h1
      w(i)   = w(i)   + hs/6.d0 * (2.d0 - h1/h0)
      w(i+1) = w(i+1) + hs**3 / (6.d0*h0*h1)
      w(i+2) = w(i+2) + hs/6.d0 * (2.d0 - h0/h1)
      i = i + 2
    end do

    if (i < n) then                        ! even n: one interval left over
      h0 = x(i+1) - x(i)
      w(i)   = w(i)   + 0.5d0*h0
      w(i+1) = w(i+1) + 0.5d0*h0
    end if

  end subroutine simpson_weights

!===============================================================================
! NCINT_2D
!
! s = integral integral f dvperp dvpar, with f already carrying the Jacobian.
!
! Composite Simpson in BOTH directions (see simpson_weights).  This replaced a
! left-endpoint RECTANGLE rule in v_perp, which was only 1st-order accurate on a
! non-uniform grid -- on a uniform grid, with the integrand vanishing at both
! ends, it coincided with the trapezoidal rule, which is why the defect was
! invisible for ising=0 and appeared only for ising=+-1.
!
! Accuracy here matters beyond diagnostics: every solver renormalises its
! solution through this integral (fout = fout*npart/dens_tmp), so the quadrature
! order caps the ABSOLUTE accuracy of f.  With the trapezoidal rule that cap was
! 2nd order, which the 7-point stencil (6th order) could not reach.
!
! The v_par boundaries j=1 and j=npar are Dirichlet (f=0) and contribute
! nothing, but they are included so the rule is the plain composite one and no
! special-casing is needed.
!
! The weights depend only on the grid, so they are built once and cached; the
! cache is rebuilt automatically if nperp or npar changes.
!===============================================================================
subroutine ncint_2D(f,s)

! Initialisation
! --------------

use shared_grid

implicit none

double precision, intent(in) :: f(nperp,npar)
double precision, intent(out):: s

double precision, allocatable, save :: wperp(:), wpar(:)
integer, save :: n_perp_cached = -1, n_par_cached = -1

double precision :: col
integer i,j

! ===========================================================================

if (n_perp_cached /= nperp .or. n_par_cached /= npar) then
    if (allocated(wperp)) deallocate(wperp)
    if (allocated(wpar))  deallocate(wpar)
    allocate(wperp(nperp), wpar(npar))
    call simpson_weights(vperp, nperp, wperp)
    call simpson_weights(vpar,  npar,  wpar)
    n_perp_cached = nperp
    n_par_cached  = npar
end if

s = 0.d0

do i = 1, nperp
    col = 0.d0
    do j = 1, npar
        col = col + f(i,j) * wpar(j)
    end do
    s = s + col * wperp(i)
end do

!
!============================================================================
!
end subroutine ncint_2D



!
!============================================================================

end module mod_ncint
