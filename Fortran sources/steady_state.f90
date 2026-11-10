!*******************************************************************
!*          Resolution of the Fokker-Planck Equation               *
!*               with the linear collision term.   
!*               Steady-state solution                             * 
!*******************************************************************

module mod_linear

contains

subroutine FP_steady_state(all20,all02,all11,all10,all01,all00,fout)

use shared_grid
use shared_plasma
use shared_beam
use shared_timer
use time_comps_mod

!use func_index

implicit none

double precision, dimension(nperp,npar),intent(in) :: all20,all02,all11,all10,all01,all00
double precision, dimension(nperp,npar),intent(out) :: fout

logical no_beam

double precision dens_tmp

integer:: iv,imu

!====================================================================


if(isource == 0) then
    no_beam = .true.
else
    no_beam = .false.
endif


call solve_fp_pardiso(all00,all10,all01,all20,all11,all02,no_beam, fout)

!====================================================================

! We compute the density associated with our solution and renormalize
! the distribution function to the correct value

call time_density(fout,dens_tmp)
write(*,*) 'Unnormalized density is ',dens_tmp

write(*,*)''
write(*,*) 'Renormalizing...'

write(*,*)''


fout = fout*npart/dens_tmp

! The properly renormalized solution is stored on disk

open(40,file=TRIM(outfile('fout.txt')), status='unknown')

do iv=1,nperp
        do imu=1,npar
        write(40,*) vperp(iv), vpar(imu), fout(iv,imu)
        enddo
enddo

close(40)	

! We also store the solution (in column vector form) 

 open(42,file=TRIM(outfile('xout.dat')),status='unknown')
 write(42,*) 0.d0 !First timestep
 do iv=1,nperp
        do imu=1,npar
  !       ix = index_mat(iv,imu)
         write(42,*) fout(iv,imu)*npart/dens_tmp
        enddo
 enddo
 close(42)
                            
               
end subroutine  FP_steady_state

!***********************************************************************

end module mod_linear
