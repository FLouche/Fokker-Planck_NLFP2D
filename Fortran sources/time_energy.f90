!***************************************************
!* Computation of the temperature at each timestep *
!*                                                 *
!*  Version 2.0: adapted to general spacing grids  *
!*   30/03/2026; F.Louche                          *
!***************************************************

subroutine time_energy(f,dens,temp,tperp,tpar)

use shared_grid
use mod_ncint
use shared_plasma

implicit none

double precision, intent(in) :: f(nperp,npar),dens
double precision, intent(out) :: temp,tpar,tperp

double precision pmass, kev_in_J, mod

data pmass/1.6726d-27/ !proton mass in kg
data kev_in_J/1.60218d-16/ !convert keV to Joule

double precision,allocatable, dimension(:,:) :: fint
integer iv,ip

allocate(fint(nperp,npar))

    ! Particle kinetic energy
    
    ! -> total
    
    do iv = 1,nperp
        do ip = 1,npar
            fint(iv,ip) = f(iv,ip)*(vperp(iv)**2+vpar(ip)**2)*jacob(iv,ip)/dens
        enddo
    enddo
    
    call ncint_2d(fint,mod)
    temp = 0.5*pmass*aa*mod/kev_in_J
    
    ! -> Perpendicular energy
    
     do iv = 1,nperp
        do ip = 1,npar
            fint(iv,ip) = f(iv,ip)*vperp(iv)**2*jacob(iv,ip)/dens
        enddo
     enddo
    
    call ncint_2d(fint,mod)
    tperp = 0.5*pmass*aa*mod/kev_in_J
!    
     ! -> Parallel energy
    
     do iv = 1,nperp
        do ip = 1,npar
            fint(iv,ip) = f(iv,ip)*vpar(ip)**2*jacob(iv,ip)/dens
        enddo
     enddo
     
     call ncint_2d(fint,mod)
     tpar = 0.5*pmass*aa*mod/kev_in_J

deallocate(fint)
     
end subroutine time_energy

!==============================================
