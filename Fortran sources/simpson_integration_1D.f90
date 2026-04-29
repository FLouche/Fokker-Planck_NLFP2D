module integration_mod
  use iso_fortran_env, only: dp => real64
  implicit none
  private
  public :: integrate_grid, integrate_grid_accurate, dp
  
contains

  !> Integrate function F on grid from x(j) to x(k) using adaptive high-order methods
  function integrate_grid(x, F, j, k) result(integral)
    real(dp), intent(in) :: x(:)    ! Grid points (must be sorted)
    real(dp), intent(in) :: F(:)    ! Function values at grid points
    integer, intent(in) :: j, k     ! Integration limits (indices)
    real(dp) :: integral
    
    integer :: n_intervals
    
    ! Check inputs
    if (j < 1 .or. k > size(x) .or. j >= k) then
     ! write(*,*) 'Error: Invalid integration limits'
      integral = 0.0_dp
      return
    end if
    
    n_intervals = k - j
    
    ! Special case: single interval
    if (n_intervals == 1) then
      integral = trapezoidal_single(x(j), x(k), F(j), F(k))
      return
    end if
    
    ! Special case: two intervals - use Simpson's rule
    if (n_intervals == 2) then
      integral = simpson_three_point(x(j), x(j+1), x(k), F(j), F(j+1), F(k))
      return
    end if
    
    ! For longer intervals, use adaptive method based on grid uniformity
    if (is_uniform_grid(x, j, k)) then
      integral = integrate_uniform_grid(x, F, j, k)
    else
      integral = integrate_nonuniform_grid(x, F, j, k)
    end if
    
  end function integrate_grid
  
  !> Accurate integration using cubic spline interpolation for better accuracy with few points
  function integrate_grid_accurate(x, F, j, k, n_refine) result(integral)
    real(dp), intent(in) :: x(:)    ! Grid points (must be sorted)
    real(dp), intent(in) :: F(:)    ! Function values at grid points
    integer, intent(in) :: j, k     ! Integration limits (indices)
    integer, intent(in), optional :: n_refine  ! Refinement factor (default 4)
    real(dp) :: integral
    
    integer :: n_ref, i_start, i_end, n_base
    integer :: i
    real(dp), allocatable :: d2y(:)
    
    ! Check inputs
    if (j < 1 .or. k > size(x) .or. j >= k) then
     ! write(*,*) 'Error: Invalid integration limits'
      integral = 0.0_dp
      return
    end if
    
    ! Set refinement factor
    if (present(n_refine)) then
      n_ref = max(1, n_refine)
    else
      n_ref = 4
    end if
    
    n_base = k - j + 1
    
    ! For very few points, use simple Hermite-based integration
    if (n_base <= 2) then
      integral = integrate_with_hermite(x, F, j, k)
      return
    end if
    
    ! Extend the range to include neighboring points for better spline accuracy
    i_start = max(1, j - 1)
    i_end = min(size(x), k + 1)
    
    ! Compute spline coefficients on extended range
    allocate(d2y(i_start:i_end))
    call compute_natural_spline_coeffs(x(i_start:i_end), F(i_start:i_end), d2y)
    
    ! Integrate using the spline analytically on each original interval
    integral = 0.0_dp
    do i = j, k - 1
      integral = integral + integrate_spline_segment(x(i_start:i_end), F(i_start:i_end), &
                                                      d2y, x(i), x(i+1), i_start)
    end do
    
    deallocate(d2y)
    
  end function integrate_grid_accurate
  
  !> Integrate using Hermite interpolation for very few points
  function integrate_with_hermite(x, F, j, k) result(integral)
    real(dp), intent(in) :: x(:), F(:)
    integer, intent(in) :: j, k
    real(dp) :: integral
    real(dp) :: dfa, dfb, h
    
    ! Estimate derivatives using finite differences
    if (j > 1 .and. j < size(x)) then
      dfa = (F(j+1) - F(j-1)) / (x(j+1) - x(j-1))
    else if (j == 1 .and. size(x) >= 2) then
      dfa = (F(j+1) - F(j)) / (x(j+1) - x(j))
    else
      dfa = 0.0_dp
    end if
    
    if (k > 1 .and. k < size(x)) then
      dfb = (F(k+1) - F(k-1)) / (x(k+1) - x(k-1))
    else if (k == size(x) .and. k > 1) then
      dfb = (F(k) - F(k-1)) / (x(k) - x(k-1))
    else
      dfb = 0.0_dp
    end if
    
    h = x(k) - x(j)
    
    ! Analytical integral of cubic Hermite polynomial
    integral = h * (0.5_dp * (F(j) + F(k)) + h * (dfa - dfb) / 12.0_dp)
    
  end function integrate_with_hermite
  
  !> Integrate a cubic spline segment analytically
  function integrate_spline_segment(x, y, d2y, xa, xb, offset) result(integral)
    real(dp), intent(in) :: x(:), y(:), d2y(:)
    real(dp), intent(in) :: xa, xb
    integer, intent(in) :: offset
    real(dp) :: integral
    
    integer :: n, i_seg
    real(dp) :: h, dx
    real(dp) :: c0, c1, c2, c3
    
    n = size(x)
    
    ! Find which spline segment contains [xa, xb]
    i_seg = 1
    do i_seg = 1, n - 1
      if (xa >= x(i_seg) - 1.0e-10_dp .and. xb <= x(i_seg+1) + 1.0e-10_dp) then
        exit
      end if
    end do
    
    ! For a cubic spline on interval [x(i), x(i+1)], we have:
    ! S(t) = a + b*t + c*t^2 + d*t^3, where t = (x - x(i)) / h
    ! 
    ! The spline can be written as:
    ! S(x) = y(i) * (1-t) + y(i+1) * t 
    !        + [(1-t)^3 - (1-t)] * h^2/6 * d2y(i)
    !        + [t^3 - t] * h^2/6 * d2y(i+1)
    
    h = x(i_seg+1) - x(i_seg)
    
    ! Get coefficients for the polynomial representation
    ! S(x) = c0 + c1*(x-x(i)) + c2*(x-x(i))^2 + c3*(x-x(i))^3
    call get_spline_poly_coeffs(x(i_seg), x(i_seg+1), y(i_seg), y(i_seg+1), &
                                 d2y(i_seg), d2y(i_seg+1), c0, c1, c2, c3)
    
    ! Integrate polynomial from xa to xb
    dx = xb - xa
    integral = integrate_cubic_poly(c0, c1, c2, c3, xa - x(i_seg), xb - x(i_seg))
    
  end function integrate_spline_segment
  
  !> Get polynomial coefficients for cubic spline segment
  subroutine get_spline_poly_coeffs(x1, x2, y1, y2, d2y1, d2y2, c0, c1, c2, c3)
    real(dp), intent(in) :: x1, x2, y1, y2, d2y1, d2y2
    real(dp), intent(out) :: c0, c1, c2, c3
    real(dp) :: h
    
    h = x2 - x1
    
    ! Coefficients for S(x) = c0 + c1*t + c2*t^2 + c3*t^3, where t = x - x1
    c0 = y1
    c1 = (y2 - y1) / h - h * (2.0_dp * d2y1 + d2y2) / 6.0_dp
    c2 = d2y1 / 2.0_dp
    c3 = (d2y2 - d2y1) / (6.0_dp * h)
    
  end subroutine get_spline_poly_coeffs
  
  !> Integrate cubic polynomial c0 + c1*t + c2*t^2 + c3*t^3 from ta to tb
  function integrate_cubic_poly(c0, c1, c2, c3, ta, tb) result(integral)
    real(dp), intent(in) :: c0, c1, c2, c3, ta, tb
    real(dp) :: integral
    
    integral = c0 * (tb - ta) + &
               c1 * (tb**2 - ta**2) / 2.0_dp + &
               c2 * (tb**3 - ta**3) / 3.0_dp + &
               c3 * (tb**4 - ta**4) / 4.0_dp
    
  end function integrate_cubic_poly
  
  !> Compute natural cubic spline second derivatives
  subroutine compute_natural_spline_coeffs(x, y, d2y)
    real(dp), intent(in) :: x(:), y(:)
    real(dp), intent(out) :: d2y(:)
    
    integer :: n, i, i_base
    real(dp), allocatable :: a(:), b(:), c(:), r(:), soln(:)
    real(dp) :: h1, h2
    
    n = size(x)
    i_base = lbound(d2y, 1)
    
    ! Natural spline: d2y = 0 at endpoints
    d2y(i_base) = 0.0_dp
    d2y(i_base + n - 1) = 0.0_dp
    
    if (n == 2) then
      return
    end if
    
    allocate(a(n-2), b(n-2), c(n-2), r(n-2), soln(n-2))
    
    ! Build tridiagonal system
    do i = 2, n - 1
      h1 = x(i) - x(i-1)
      h2 = x(i+1) - x(i)
      
      a(i-1) = h1
      b(i-1) = 2.0_dp * (h1 + h2)
      c(i-1) = h2
      r(i-1) = 6.0_dp * ((y(i+1) - y(i)) / h2 - (y(i) - y(i-1)) / h1)
    end do
    
    ! Solve tridiagonal system
    call solve_tridiagonal(a, b, c, r, soln)
    
    ! Copy solution
    do i = 1, n - 2
      d2y(i_base + i) = soln(i)
    end do
    
    deallocate(a, b, c, r, soln)
    
  end subroutine compute_natural_spline_coeffs
  
  !> Solve tridiagonal system
  subroutine solve_tridiagonal(a, b, c, r, x)
    real(dp), intent(in) :: a(:), b(:), c(:), r(:)
    real(dp), intent(out) :: x(:)
    
    integer :: n, i
    real(dp), allocatable :: cp(:), rp(:)
    real(dp) :: m
    
    n = size(b)
    allocate(cp(n), rp(n))
    
    ! Forward elimination
    cp(1) = c(1) / b(1)
    rp(1) = r(1) / b(1)
    
    do i = 2, n
      m = b(i) - a(i) * cp(i-1)
      if (i < n) then
        cp(i) = c(i) / m
      end if
      rp(i) = (r(i) - a(i) * rp(i-1)) / m
    end do
    
    ! Back substitution
    x(n) = rp(n)
    do i = n - 1, 1, -1
      x(i) = rp(i) - cp(i) * x(i+1)
    end do
    
    deallocate(cp, rp)
    
  end subroutine solve_tridiagonal
  
  !> Single interval trapezoidal rule
  function trapezoidal_single(x1, x2, f1, f2) result(integral)
    real(dp), intent(in) :: x1, x2, f1, f2
    real(dp) :: integral
    integral = 0.5_dp * (x2 - x1) * (f1 + f2)
  end function trapezoidal_single
  
  !> Three-point Simpson's rule for non-uniform grid
  function simpson_three_point(x1, x2, x3, f1, f2, f3) result(integral)
    real(dp), intent(in) :: x1, x2, x3, f1, f2, f3
    real(dp) :: integral
    real(dp) :: h1, h2
    
    h1 = x2 - x1
    h2 = x3 - x2
    
    if (abs(h1 - h2) < 1.0e-12_dp * (h1 + h2)) then
      ! Uniform spacing - standard Simpson's rule
      integral = h1 / 3.0_dp * (f1 + 4.0_dp*f2 + f3)
    else
      ! Non-uniform spacing - exact integration of quadratic through 3 points
      integral = (h1 + h2) / 6.0_dp * ((2.0_dp - h2/h1) * f1 + &
                 (h1 + h2)**2 / (h1 * h2) * f2 + (2.0_dp - h1/h2) * f3)
    end if
  end function simpson_three_point
  
  !> Check if grid segment is uniform
  function is_uniform_grid(x, j, k) result(uniform)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: j, k
    logical :: uniform
    real(dp) :: h_avg, tol
    integer :: i
    
    h_avg = (x(k) - x(j)) / (k - j)
    tol = 1.0e-10_dp * h_avg
    uniform = .true.
    
    do i = j, k-1
      if (abs((x(i+1) - x(i)) - h_avg) > tol) then
        uniform = .false.
        return
      end if
    end do
  end function is_uniform_grid
  
  !> Integration for uniform grids using Boole's rule (5-point, 6th order)
  function integrate_uniform_grid(x, F, j, k) result(integral)
    real(dp), intent(in) :: x(:), F(:)
    integer, intent(in) :: j, k
    real(dp) :: integral
    real(dp) :: h
    integer :: i, n, n_boole, n_remaining
    
    n = k - j
    h = x(j+1) - x(j)
    integral = 0.0_dp
    
    ! Use Boole's rule (5-point) for groups of 4 intervals
    n_boole = n / 4
    
    do i = 0, n_boole - 1
      integral = integral + boole_rule(h, F(j + 4*i:j + 4*i + 4))
    end do
    
    ! Handle remaining intervals
    n_remaining = n - 4 * n_boole
    
    if (n_remaining == 1) then
      integral = integral + trapezoidal_single(x(j+4*n_boole), x(j+4*n_boole+1), &
                                                F(j+4*n_boole), F(j+4*n_boole+1))
    else if (n_remaining == 2) then
      integral = integral + simpson_uniform(h, F(j+4*n_boole), F(j+4*n_boole+1), &
                                             F(j+4*n_boole+2))
    else if (n_remaining == 3) then
      integral = integral + simpson_3_8_rule(h, F(j+4*n_boole:j+4*n_boole+3))
    end if
    
  end function integrate_uniform_grid
  
  !> Boole's rule (5-point Newton-Cotes, 6th order accurate)
  function boole_rule(h, f) result(integral)
    real(dp), intent(in) :: h
    real(dp), intent(in) :: f(0:4)
    real(dp) :: integral
    
    integral = 2.0_dp * h / 45.0_dp * &
               (7.0_dp*f(0) + 32.0_dp*f(1) + 12.0_dp*f(2) + 32.0_dp*f(3) + 7.0_dp*f(4))
  end function boole_rule
  
  !> Simpson's rule for uniform spacing
  function simpson_uniform(h, f0, f1, f2) result(integral)
    real(dp), intent(in) :: h, f0, f1, f2
    real(dp) :: integral
    
    integral = h / 3.0_dp * (f0 + 4.0_dp*f1 + f2)
  end function simpson_uniform
  
  !> Simpson's 3/8 rule for 3 uniform intervals
  function simpson_3_8_rule(h, f) result(integral)
    real(dp), intent(in) :: h
    real(dp), intent(in) :: f(0:3)
    real(dp) :: integral
    
    integral = 3.0_dp * h / 8.0_dp * (f(0) + 3.0_dp*f(1) + 3.0_dp*f(2) + f(3))
  end function simpson_3_8_rule
  
  !> Integration for non-uniform grids
  function integrate_nonuniform_grid(x, F, j, k) result(integral)
    real(dp), intent(in) :: x(:), F(:)
    integer, intent(in) :: j, k
    real(dp) :: integral
    integer :: i
    
    integral = 0.0_dp
    i = j
    
    ! Process intervals in pairs when possible (for Simpson's rule)
    do while (i + 2 <= k)
      integral = integral + simpson_three_point(x(i), x(i+1), x(i+2), &
                                                 F(i), F(i+1), F(i+2))
      i = i + 2
    end do
    
    ! Handle last interval if odd number
    if (i < k) then
      integral = integral + trapezoidal_single(x(i), x(i+1), F(i), F(i+1))
    end if
    
  end function integrate_nonuniform_grid

end module integration_mod