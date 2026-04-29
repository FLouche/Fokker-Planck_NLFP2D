module integrate_2d_module
    implicit none
    private
    public :: integrate_2d, test_integration
    
contains

    !===========================================================================
    ! Function: integrate_2d
    !
    ! Purpose: Compute the 2D integral of F(x,y) over the domain [0,xmax] x [-ymax,ymax]
    !          using the trapezoidal rule with arbitrary point distributions
    !
    ! Inputs:
    !   F      - 2D array containing function values: F(i,j) = F(x(i), y(j))
    !   x      - Array of x coordinates (must be monotonically increasing, x(1)=0, x(nx)=xmax)
    !   y      - Array of y coordinates (must be monotonically increasing, y(1)=-ymax, y(ny)=ymax)
    !   nx     - Number of points in x direction
    !   ny     - Number of points in y direction
    !
    ! Returns:
    !   result - The computed integral value
    !
    ! Notes:
    !   - Uses composite trapezoidal rule in both dimensions
    !   - Handles arbitrary (non-uniform) spacing in x and y
    !   - For uniform grids, this reduces to standard trapezoidal rule
    !   - Accuracy is O(h^2) for smooth functions
    !   - F array is indexed as F(i,j) where i corresponds to x and j to y
    !===========================================================================
    function integrate_2d(F, x, y, nx, ny) result(integral)
        implicit none
        
        ! Arguments
        integer, intent(in) :: nx, ny
        real(8), dimension(nx, ny), intent(in) :: F
        real(8), dimension(nx), intent(in) :: x
        real(8), dimension(ny), intent(in) :: y
        real(8) :: integral
        
        ! Local variables
        integer :: i, j
        real(8) :: dx_left, dy_bottom
        real(8) :: area, f_val
        real(8) :: sum_integral
        
        ! Initialize
        sum_integral = 0.0d0
        
        ! Check inputs
        if (nx < 2 .or. ny < 2) then
            write(*,*) 'Error: Need at least 2 points in each direction'
            integral = 0.0d0
            return
        end if
        
        ! Integrate using composite trapezoidal rule
        ! For arbitrary spacing, we treat each cell as a trapezoid in both directions
        
        do i = 1, nx-1
            do j = 1, ny-1
                ! Get spacing for this cell
                dx_left = x(i+1) - x(i)
                dy_bottom = y(j+1) - y(j)
                
                ! Area of this cell
                area = dx_left * dy_bottom
                
                ! Trapezoidal rule: average of function values at 4 corners
                ! multiplied by the cell area
                f_val = 0.25d0 * (F(i,   j) + &
                                  F(i+1, j) + &
                                  F(i,   j+1) + &
                                  F(i+1, j+1))
                
                sum_integral = sum_integral + area * f_val
            end do
        end do
        
        integral = sum_integral
        
    end function integrate_2d


    !===========================================================================
    ! Subroutine: test_integration
    !
    ! Purpose: Test the integration routine with known analytical results
    !===========================================================================
    subroutine test_integration()
        implicit none
        
        integer :: nx, ny, i, j
        real(8), dimension(:), allocatable :: x, y
        real(8), dimension(:,:), allocatable :: F_array
        real(8) :: xmax, ymax, result, exact, error
        
        write(*,*) '========================================='
        write(*,*) 'Testing 2D Integration Routine'
        write(*,*) '========================================='
        write(*,*)
        
        ! Test 1: Constant function F(x,y) = 1
        ! Exact integral = xmax * 2*ymax
        write(*,*) 'Test 1: F(x,y) = 1'
        xmax = 2.0d0
        ymax = 3.0d0
        nx = 50
        ny = 50
        
        allocate(x(nx), y(ny), F_array(nx, ny))
        
        ! Create uniform grid
        do i = 1, nx
            x(i) = xmax * dble(i-1) / dble(nx-1)
        end do
        do j = 1, ny
            y(j) = -ymax + 2.0d0*ymax * dble(j-1) / dble(ny-1)
        end do
        
        ! Compute function values on grid
        do i = 1, nx
            do j = 1, ny
                F_array(i,j) = test_func_1(x(i), y(j))
            end do
        end do
        
        result = integrate_2d(F_array, x, y, nx, ny)
        exact = xmax * 2.0d0 * ymax
        error = abs(result - exact) / abs(exact)
        
        write(*,'(A,F15.8)') '  Computed: ', result
        write(*,'(A,F15.8)') '  Exact:    ', exact
        write(*,'(A,E12.4)') '  Rel Error:', error
        write(*,*)
        
        ! Test 2: Linear function F(x,y) = x + y
        ! Exact integral = xmax^2/2 * 2*ymax + xmax * (ymax^2 - ymax^2) = xmax^2 * ymax
        write(*,*) 'Test 2: F(x,y) = x + y'
        do i = 1, nx
            do j = 1, ny
                F_array(i,j) = test_func_2(x(i), y(j))
            end do
        end do
        
        result = integrate_2d(F_array, x, y, nx, ny)
        exact = xmax**2 * ymax
        error = abs(result - exact) / abs(exact)
        
        write(*,'(A,F15.8)') '  Computed: ', result
        write(*,'(A,F15.8)') '  Exact:    ', exact
        write(*,'(A,E12.4)') '  Rel Error:', error
        write(*,*)
        
        ! Test 3: Quadratic function F(x,y) = x^2 + y^2
        ! Exact integral = xmax^3/3 * 2*ymax + xmax * 2*ymax^3/3
        write(*,*) 'Test 3: F(x,y) = x^2 + y^2'
        do i = 1, nx
            do j = 1, ny
                F_array(i,j) = test_func_3(x(i), y(j))
            end do
        end do
        
        result = integrate_2d(F_array, x, y, nx, ny)
        exact = (xmax**3 / 3.0d0) * 2.0d0*ymax + xmax * 2.0d0*ymax**3/3.0d0
        error = abs(result - exact) / abs(exact)
        
        write(*,'(A,F15.8)') '  Computed: ', result
        write(*,'(A,F15.8)') '  Exact:    ', exact
        write(*,'(A,E12.4)') '  Rel Error:', error
        write(*,*)
        
        ! Test 4: Non-uniform grid
        write(*,*) 'Test 4: F(x,y) = 1 with non-uniform grid'
        ! Create non-uniform grid (clustered near boundaries)
        do i = 1, nx
            x(i) = xmax * (dble(i-1) / dble(nx-1))**2
        end do
        do j = 1, ny
            y(j) = -ymax + 2.0d0*ymax * (dble(j-1) / dble(ny-1))**2
        end do
        
        ! Recompute function values on new grid
        do i = 1, nx
            do j = 1, ny
                F_array(i,j) = test_func_1(x(i), y(j))
            end do
        end do
        
        result = integrate_2d(F_array, x, y, nx, ny)
        exact = xmax * 2.0d0 * ymax
        error = abs(result - exact) / abs(exact)
        
        write(*,'(A,F15.8)') '  Computed: ', result
        write(*,'(A,F15.8)') '  Exact:    ', exact
        write(*,'(A,E12.4)') '  Rel Error:', error
        write(*,*)
        
        deallocate(x, y, F_array)
        
        write(*,*) '========================================='
        write(*,*) 'Testing Complete'
        write(*,*) '========================================='
        
    end subroutine test_integration
    
    
    ! Test functions
    function test_func_1(x, y) result(f)
        real(8), intent(in) :: x, y
        real(8) :: f
        f = 1.0d0
    end function test_func_1
    
    function test_func_2(x, y) result(f)
        real(8), intent(in) :: x, y
        real(8) :: f
        f = x + y
    end function test_func_2
    
    function test_func_3(x, y) result(f)
        real(8), intent(in) :: x, y
        real(8) :: f
        f = x**2 + y**2
    end function test_func_3

end module integrate_2d_module



