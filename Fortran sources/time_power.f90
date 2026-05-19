!*******************************************************
!* Computation of the power densities at each timestep *
!*******************************************************

subroutine time_power(f,dens,pcoll,pRF,psource,plosses,pcoll_self)

use shared_grid
use mod_ncint
use shared_plasma
use shared_beam
use shared_rf
use shared_FPterms
use mod_build_ss
use shared_timer



use func_index


implicit none

double precision, intent(in) :: f(nbig),dens
double precision, intent(out) :: pcoll(nbulk),pRF,psource,plosses,pcoll_self
double precision ekin(nperp,npar),normfac
double precision, dimension(nperp,npar):: f1,f2

double precision :: bigm(nbig,nbig),f_wrk(nbig)
double PRECISION, allocatable, dimension(:,:) :: rf00

double precision pmass

data pmass/1.6726d-27/ !proton mass in kg

double precision,allocatable, dimension(:,:) :: fint
double precision :: taum_save     ! saved taum; restored on exit
integer iv,ip,ix,imu,ib

external dgemv

allocate(fint(nperp,npar))


! ***************************************************
! Collisional power density
! ***************************************************

    taum_save = taum
    taum = 0.d0 ! temporary: build_ss reads taum; zero it for diagnostic calls
    
    if(isource == 0)then
        normfac = npart/dens
    else
        normfac=1.d0
    endif
    

do ib = 1,nbulk
    
     ! We compute the df/dt collisonal term in matrix form
            
    call build_ss(colin20_sp(:,:,ib),colin02_sp(:,:,ib),colin11_sp(:,:,ib),colin10_sp(:,:,ib),colin01_sp(:,:,ib),colin00_sp(:,:,ib),bigm)
    
     ! We apply the operator to the solution
    
    call dgemv('N',nbig,nbig,1.d0,bigm,nbig,f,1,0.d0,f_wrk,1)
    
    ! -> total
    
    do iv = 1,nperp
        do ip = 1,npar
            ekin(iv,ip) = 0.5*pmass*aa*(vperp(iv)**2+vpar(ip)**2)*normfac
            ix = index_mat(iv,ip)
            fint(iv,ip) = ekin(iv,ip)*f_wrk(ix)*jacob(iv,ip)
        enddo
    enddo
    
    call ncint_2d(fint,pcoll(ib))

enddo

!  --> self-collisions

if (isc == -1) then
    
         ! We compute the df/dt collisonal term in matrix form
            
    call build_ss(sc20,sc02,sc11,sc10,sc01,sc00,bigm)
    
     ! We apply the operator to the solution
    
    call dgemv('N',nbig,nbig,1.d0,bigm,nbig,f,1,0.d0,f_wrk,1)
    
    ! -> total
    
    do iv = 1,nperp
        do ip = 1,npar
            ix = index_mat(iv,ip)
            fint(iv,ip) = ekin(iv,ip)*f_wrk(ix)*jacob(iv,ip)
        enddo
    enddo
    
    call ncint_2d(fint,pcoll_self)
    
else
    
    pcoll_self=0.d0

endif

! ***************************************************
! RF power density
! ***************************************************

    ! We compute the df/dt RF term in matrix form
    
    if(irf == -1) then
        
        allocate(rf00(nperp,npar))
        rf00 = 0.d0
    
        call build_ss(rf20,rf02,rf11,rf10,rf01,rf00,bigm)
        
        deallocate(rf00)
        
          ! We apply the operator to the solution
    
        call dgemv('N',nbig,nbig,1.d0,bigm,nbig,f,1,0.d0,f_wrk,1)
        
        do iv = 1,nperp
                
            do imu = 1,npar
            
              ix = index_mat(iv,imu)
              fint(iv,imu) = ekin(iv,imu)*f_wrk(ix)*jacob(iv,imu)
            
            enddo
         enddo
    
    call ncint_2D(fint,pRF)
        
    else ! no RF
        
        pRF = 0.d0
        
    endif    
    
! ***************************************************
! Beam (+ losses) power density
! ***************************************************

    if (isource == -1) then    
        taum = 1.d0/taus
        
        do iv = 1,nperp
            do imu = 1,npar
            
            ix = index_mat(iv,imu)
            f1(iv,imu) = -ekin(iv,imu)*f(ix)*taum*jacob(iv,imu)
            f2(iv,imu) = ekin(iv,imu)*source(iv,imu)*jacob(iv,imu)
            
            enddo
        enddo
        
        call ncint_2D(f1,plosses)
        call ncint_2D(f2,psource)
        
    else

        plosses = 0.d0
        psource = 0.d0

    endif

    taum = taum_save   ! restore original value for the caller

deallocate(fint)
     
end subroutine time_power

!==============================================
