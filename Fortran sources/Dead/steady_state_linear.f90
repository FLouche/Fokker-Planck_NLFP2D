! ---------------------------------------------------------------------------
! DEAD CODE, kept for reference only. Not in the Visual Studio project.
! Moved here on 2026-10-27 from 'steady_state.f90': subroutine linear, the original steady-state path
! (build_ss -> sparse_solve_ss -> pardiso_solve_steady). Never called: the
! steady state goes through FP_steady_state -> solve_fp_pardiso. It was a
! contained routine of module mod_linear (steady_state.f90).
! ---------------------------------------------------------------------------

subroutine linear(all20,all02,all11,all10,all01,all00,fout)

! Initialisation
! --------------

use shared_grid
use shared_plasma

use mod_build_ss

use mod_sparse_solve

use shared_timer

use func_index

use time_comps_mod

implicit none

intrinsic SIGN,DLOG10,DABS

double precision, dimension(nperp,npar),intent(in) :: all20,all02,all11,all10,all01,all00
double precision, dimension(nperp,npar),intent(out) :: fout
double precision, dimension(nbig) :: xout

double precision :: bigm(nbig,nbig),bigv(nbig)

double precision dens_tmp

integer iv,imu,ix


!======================================================================
!
! Construction of the Coefficients Matrix & the Independant Term Vector
! ---------------------------------------------------------------------

call build_ss(all20,all02,all11,all10,all01,all00,bigm,bigv)


!======================================================================
!
! Resolution of the Steady-State Equation
! ---------------------------------------

call sparse_solve_ss(-1,bigm,bigv,xout)

do iv=1,nperp
        do imu=1,npar
        ix = index_mat(iv,imu)
        fout(iv,imu) = xout(ix)
    enddo
enddo

!====================================================================

! We compute the density associated with our solution and renormalize
! the distribution function to the correct value


call time_density(fout,dens_tmp)


! Rescale to physical density n_species

write(*,*) 'Unnormalized density is ',dens_tmp

write(*,*)''
write(*,*) 'Renormalizing...'

write(*,*)''

! Rescale to physical density n_species


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
         ix = index_mat(iv,imu)
         write(42,*) xout(ix)*npart/dens_tmp
        enddo
 enddo
 close(42)
 
end subroutine linear
