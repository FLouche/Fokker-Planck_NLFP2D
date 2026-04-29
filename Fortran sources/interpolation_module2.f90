module interpolation_module
    
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)

    
contains

    
    subroutine cylindrical_to_spherical_interp(nr_cyl, nz_cyl, rcyl, z, F_cyl, &
                                             nr_sph, ntheta_sph, rsph, mu, F_sph)
    
! Interpolate function F from cylindrical (r, z) to spherical (r, mu) coordinates 
! where mu is cosine of the pitch-angle
    
  
  ! Cylindrical grid parameters
  integer :: nr_cyl, nz_cyl
  real(dp) :: rcyl(nr_cyl), z(nz_cyl)
  real(dp), dimension(nr_cyl,nz_cyl) :: F_cyl
  
  ! Spherical grid parameters
  integer :: nr_sph, ntheta_sph
  real(dp), dimension(nr_sph) :: rsph
  real(dp), dimension(ntheta_sph) :: mu
  real(dp), dimension(nr_sph, ntheta_sph) :: F_sph
  
  ! Local variables
  integer :: i, j
  real(dp) :: rcyl_interp, z_interp
    
  ! Perform interpolation from cylindrical to spherical
  print *, 'Performing interpolation...'
  do j = 1, ntheta_sph
    do i = 1, nr_sph
      ! Convert spherical coordinates to cylindrical
      rcyl_interp = rsph(i) * dsqrt(1-mu(j)**2)
      z_interp = rsph(i) * mu(j)
      
      ! Interpolate F at (rcyl_interp, z_interp)
      F_sph(i, j) = bilinear_interp(rcyl, z, F_cyl, nr_cyl, nz_cyl, &
                                     rcyl_interp, z_interp)
    end do
  end do
    
  print *, 'Interpolation complete.'
    
end subroutine cylindrical_to_spherical_interp


  function bilinear_interp(x, y, f, nx, ny, x0, y0) result(f0)
    ! Bilinear interpolation on a 2D grid
    implicit none
    integer, intent(in) :: nx, ny
    real(dp), dimension(nx), intent(in) :: x
    real(dp), dimension(ny), intent(in) :: y
    real(dp), dimension(nx, ny), intent(in) :: f
    real(dp), intent(in) :: x0, y0
    real(dp) :: f0
    
    integer :: i1, i2, j1, j2
    real(dp) :: t, u
    
    ! Find indices for x0
    if (x0 <= x(1)) then
      i1 = 1; i2 = 1; t = 0.0_dp
    else if (x0 >= x(nx)) then
      i1 = nx; i2 = nx; t = 0.0_dp
    else
      do i1 = 1, nx-1
        if (x0 >= x(i1) .and. x0 <= x(i1+1)) exit
      end do
      i2 = i1 + 1
      t = (x0 - x(i1)) / (x(i2) - x(i1))
    end if
    
    ! Find indices for y0
    if (y0 <= y(1)) then
      j1 = 1; j2 = 1; u = 0.0_dp
    else if (y0 >= y(ny)) then
      j1 = ny; j2 = ny; u = 0.0_dp
    else
      do j1 = 1, ny-1
        if (y0 >= y(j1) .and. y0 <= y(j1+1)) exit
      end do
      j2 = j1 + 1
      u = (y0 - y(j1)) / (y(j2) - y(j1))
    end if
    
    ! Bilinear interpolation
    f0 = (1.0_dp - t) * (1.0_dp - u) * f(i1, j1) + &
         t * (1.0_dp - u) * f(i2, j1) + &
         (1.0_dp - t) * u * f(i1, j2) + &
         t * u * f(i2, j2)
    
  end function bilinear_interp


end module interpolation_module