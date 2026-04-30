!
! *****************************************
! * Grid construction and related objects *
! *****************************************
!

module mod_grid

contains


subroutine make_grid
!
! Initialisation
! -----------------

use shared_grid

implicit none

integer, parameter :: dp = kind(1.0d0)

common/mathcons/pi,twopi

double precision pi,twopi

integer i,j

! Parameters
! ----------

nbig = nperp*npar

imid = int(2*nperp/5)+1   ! default; override in calling code if needed

jmid = int(npar+1)/2

! Parallel Velocity grid construction
! ----------------------------------------
!

dvpar = (vpar_max-vpar_min)/(npar-1)

do i = 1,npar

    vpar(i) = vpar_min+(i-1)*dvpar
    
enddo

! Perpendicular Velocity grid construction
! ----------------------------------------
!

! --> Case of homogeneous grid

dvperp=(vperp_max-vperp_min)/(nperp-1)

dv2 = dvperp**2
dmu2 = dvpar**2

! --> Case of inhomogeneous grid

!  First derivative constant

dvleft=(vbound-vperp_min)/(nsing-1)
dvright=(vperp_max-vbound)/(nperp-nsing)


do i=1,nperp

non_uniform_grid:    if(ising == 0) then
    
    vperp(i)=vperp_min+(i-1)*dvperp
    
else non_uniform_grid
    
    quadratic: if (ising == -1) then

   left:		if(i <= nsing) then
				  vperp(i)=vperp_min+(i-1)*dvleft
					           else left
 	              vperp(i)=vbound+(i-nsing)*dvright
                               endif left                          
else quadratic
    
    vperp(i) = vperp_min + vperp_max * (real(i-0.5d0, dp)/real(nperp-1, dp))**2!-0.5
    
endif quadratic
		         

	                 endif non_uniform_grid

enddo
 


!  -> Total jacobian calculation
!
! Jacob(v, mu) = 2 pi vperp 

do i =1,npar
    do j=1,nperp
    
        jacob(j,i) = twopi*vperp(j)
        
    enddo
enddo

!close(40)

end subroutine make_grid


!!=======================================================================
!! Computation of the distance |x-x'| between source point and observation
!!  point in cylindrical coordinates with axisymmetry,
!!  necessary for the evaluation of the Rosenbluth potential Psi,
!!  and φ integration using uniform quadrature
!!=======================================================================
!
!subroutine distance_v_uniform
!
!
!use shared_grid
!use func_index
!
!implicit none
!
!integer, parameter :: dp = kind(1.0d0)
!  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
!
!integer, parameter :: nphi = 64  ! Number of phi integration points
!
!
!integer :: i, j, ip, jp, kphi
!integer ix1,ix2
!    real(dp) :: r, z, rp, zp, phi, dphi
!    real(dp) :: distance
!    
!        ! Phi integration step
!    dphi = 2.0_dp * pi / real(nphi, dp)
!
!
!    ! Loop over observation points
!    do i = 1, nperp
!      r = vperp(i)
!      do j = 1, npar
!        z = vpar(j)
!        
!        ix1=index_mat(i,j)
!       
!        ! Sum over source points
!        do ip = 1, nperp
!          rp = vperp(ip)
!          do jp = 1, npar
!            zp = vpar(jp)
!            
!            ix2=index_mat(ip,jp)
!                      
!            ! Integration over phi
!            sum_phi(ix1,ix2) = 0.0_dp
!           
!            do kphi = 1, nphi
!              phi = (real(kphi, dp) - 0.5_dp) * dphi
!             
!              ! Distance between (r,z,0) and (rp,zp,phi)
!              distance = sqrt(r**2 + rp**2 - 2.0_dp*r*rp*cos(phi) + (z - zp)**2)
!             
!              if (distance > 0.0_dp) then
!                sum_phi(ix1,ix2) = sum_phi(ix1,ix2) + distance * dphi
!              end if
!            end do
!            
!          enddo
!        enddo
!        
!      enddo
!    enddo
!    
!
!end subroutine distance_v_uniform

!=======================================================================
! Computation of the distance |x-x'| between source point and observation
!  point in cylindrical coordinates with axisymmetry,
!  necessary for the evaluation of the Rosenbluth potential Psi,
!  and φ integration using  Gauss-Legendre quadrature
!=======================================================================

subroutine distance_v_gauss_legendre


use shared_grid
use func_index

implicit none

integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)

integer, parameter :: nphi = 16  ! Number of phi integration points


integer :: i, j, ip, jp, kphi
integer ix1,ix2
    real(dp) :: r, z, rp, zp, phi
    real(dp) :: distance
    
double precision, dimension(nphi) :: gauss_phi, gauss_weight
    
call initialize_gauss_legendre(nphi,gauss_phi, gauss_weight)


    ! Loop over observation points
    do i = 1, nperp
      r = vperp(i)
      do j = 1, npar
        z = vpar(j)
        
        ix1=index_mat(i,j)
       
        ! Sum over source points
        do ip = 1, nperp
          rp = vperp(ip)
          do jp = 1, npar
            zp = vpar(jp)
            
            ix2=index_mat(ip,jp)
                      
            ! Integration over phi
            sum_phi(ix1,ix2) = 0.0_dp
           
            do kphi = 1, nphi
              phi = gauss_phi(kphi)
             
              ! Distance between (r,z,0) and (rp,zp,phi)
              distance = sqrt(r**2 + rp**2 - 2.0_dp*r*rp*cos(phi) + (z - zp)**2)
             
              if (distance > 0.0_dp) then
                sum_phi(ix1,ix2) = sum_phi(ix1,ix2) + distance * gauss_weight(kphi)
              end if
            end do
            
          enddo
        enddo
        
      enddo
    enddo
    

end subroutine distance_v_gauss_legendre

  !=====================================================================
  ! Initialize Gauss-Legendre nodes and weights
  !=====================================================================
  subroutine initialize_gauss_legendre(nphi_gauss_fixed,gauss_phi, gauss_weight)
    implicit none
    
    integer, intent(in):: nphi_gauss_fixed
    double precision, dimension(nphi_gauss_fixed), intent(out) :: gauss_phi, gauss_weight
    
    call setup_gauss_legendre(nphi_gauss_fixed, gauss_phi, gauss_weight)
    
    print *, "Gauss-Legendre quadrature initialized with", nphi_gauss_fixed, "points"
    print *
    
  end subroutine initialize_gauss_legendre
  
  !=====================================================================
  ! Setup Gauss-Legendre quadrature for integration over [0, 2π]
  ! Returns nodes in [0, 2π] and weights
  !=====================================================================
  subroutine setup_gauss_legendre(n, x, w)
    implicit none
    integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)

    integer, intent(in) :: n
    double precision, intent(out) :: x(n), w(n)
    
    double precision :: x_std(n), w_std(n)
    integer :: i
    
    ! Get standard Gauss-Legendre nodes/weights for [-1, 1]
    call gauss_legendre_std(n, x_std, w_std)
    
    ! Transform from [-1, 1] to [0, 2π]
    ! x_new = (x_std + 1) * π
    ! w_new = w_std * π
    do i = 1, n
      x(i) = (x_std(i) + 1.0_dp) * pi
      w(i) = w_std(i) * pi
    end do
    
  end subroutine setup_gauss_legendre
  
  !=====================================================================
  ! Standard Gauss-Legendre quadrature on [-1, 1]
  !=====================================================================
  subroutine gauss_legendre_std(n, x, w)
    implicit none
    integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)

    integer, intent(in) :: n
    double precision, intent(out) :: x(n), w(n)
    
    integer :: i, j, k, m
    double precision :: z, z1, p1, p2, p3, pp
    double precision, parameter :: eps = 1.0d-15
    
    m = (n + 1) / 2
    
    do i = 1, m
      ! Initial approximation to the i-th root
      z = cos(pi * (real(i, dp) - 0.25_dp) / (real(n, dp) + 0.5_dp))
      
      ! Newton's method iteration
      do j = 1, 50
        p1 = 1.0_dp
        p2 = 0.0_dp
        
        ! Compute Legendre polynomial P_n(z) using recurrence
        do k = 1, n
          p3 = p2
          p2 = p1
          p1 = ((2.0_dp*real(k, dp) - 1.0_dp) * z * p2 - (real(k, dp) - 1.0_dp) * p3) / real(k, dp)
        end do
        
        ! pp = P_n'(z)
        pp = real(n, dp) * (z * p1 - p2) / (z*z - 1.0_dp)
        
        ! Newton update
        z1 = z
        z = z1 - p1/pp
        
        if (abs(z - z1) < eps) exit
      end do
      
      ! Store the root and its corresponding weight
      x(i) = -z
      x(n + 1 - i) = z
      w(i) = 2.0_dp / ((1.0_dp - z*z) * pp*pp)
      w(n + 1 - i) = w(i)
    end do
    
  end subroutine gauss_legendre_std
  
!!=======================================================================
!! Computation of the distance |x-x'| between source point and observation
!!  point in cylindrical coordinates with axisymmetry,
!!  necessary for the evaluation of the Rosenbluth potential Psi,
!!  and φ integration using analytical solution (elliptic integral of 2nd kind)
!!=======================================================================
!
!subroutine distance_v_elliptic_E
!
!
!use shared_grid
!use func_index
!!use elliptic_integrals
!use xcei
!
!implicit none
!
!integer, parameter :: dp = kind(1.0d0)
!
!integer :: i, j, ip, jp
!integer ix1,ix2
!real(dp) :: r, z, rp, zp
!real(dp) :: distance
!
!real(dp) :: m1,m2,d1,d2,e1,e2
!integer :: ierr
!
!
!
!    ! Loop over observation points
!    do i = 1, nperp
!      r = vperp(i)
!      do j = 1, npar
!        z = vpar(j)
!        
!        ix1=index_mat(i,j)
!       
!        ! Sum over source points
!        do ip = 1, nperp
!          rp = vperp(ip)
!          do jp = 1, npar
!            zp = vpar(jp)
!            
!            ix2=index_mat(ip,jp)
!                      
!            ! Argument of the elliptic integrals
!            
!            d1 = (r-rp)**2+(z-zp)**2                
!                
!            m1 = -4.d0*r*rp/d1
!            
!            d2 =  (r+rp)**2+(z-zp)**2  
!            
!             m2 = 4.d0*r*rp/d2
!            
!            !e1 = elliptic_e_parameter(m1, ierr)
!
!             e1 = ceie(m1)
!             
!           ! e2 = elliptic_e_parameter(m2, ierr)
!             
!             e2 = ceie(m2)
!                             
!            sum_phi(ix1,ix2) = 2.d0*(dsqrt(d1)*e1+dsqrt(d2)*e2)
!            if (isnan(sum_phi(ix1,ix2))) then
!                write(*,*) i,j,ip,jp, 'sum is NaN'
!                write(*,*) dsqrt(d1),e1,dsqrt(d2),e2
!                stop 
!            endif
!            
!          enddo
!        enddo
!        
!      enddo
!    enddo
!    
!end subroutine distance_v_elliptic_E
  

end module mod_grid