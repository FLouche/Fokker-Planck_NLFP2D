!*******************************************************************
!*     Resolution of the time-dependant Fokker-Planck Equation     *
!*               with Cranck-Nicholson Method                      * 
!*                                                                 * 
!*******************************************************************
!    
! ***********************************************************************************
! 
!  The system is solved taking into account the sparsity of the matrix.
!
!  The PARDISO (Parallel Direct Sparse Solver) interface is used. It's provided 
!   by Intel(c) Math Kernel Libray (MKL). 
!   See "Intel(R) Math Kernel Library Reference Manual", ch. 8 p. 2495
!    or http://www.pardiso.org 
    
  !  Version 3.0 - 30/03/2026 - F. Louche
!
! ***********************************************************************************

module mod_timefp3

    contains
    
    subroutine timefp(all00_lin,all10_lin,all01_lin,all11_lin,all20_lin,all02_lin,fstart,fout,otime)
    
    ! NO NON-LINEAR CONTRIBUTIONS
    
    use mod_build_td
    use mod_build_ss
    
    use nlterm
    
    use shared_plasma
    use shared_grid
    use shared_timer
    use shared_beam
    use shared_RF
    
    use shared_FPterms
    
    use pardiso_solver
        
    use func_index
    

    !USE mod_eigen_bigm
    
    implicit none
    
    TYPE(pardiso_handle_t) :: handle
    
    external time_density, time_energy, time_power
    external sparse_matrix_vect_mult
        
    double PRECISION, dimension(nperp,npar),intent(in) :: all00_lin,all10_lin,all01_lin,all20_lin,all02_lin,all11_lin
    double PRECISION, dimension(nperp,npar):: all00,all10,all01,all20,all02,all11
    double PRECISION, dimension(nbig), intent(inout):: fstart
    double precision, dimension(nperp,npar),intent(out) :: fout

    
    double PRECISION, dimension(:), allocatable :: f_one,f_wrk
    
    double precision,dimension(nbig) :: bigv,xout
    double precision :: bigm_ss(nbig,nbig), bigm(nbig,nbig)
    double PRECISION dens_tmp, time,otime, tk,tkperp,tkpar,vtk
    
    double PRECISION pcoll(nbulk),pRF,psource,plosses,pcoll_self
    
    character (len=2) :: ibString

    integer :: flag1
    integer :: nz

    
    integer itime,ix,iv,imu,ib
    integer i,j,index
    
    
    !  PARDISO input parameters
!  ------------------------  
! Internal solver memory pointer for 64-bit architectures
 !   integer*8 pt(64)
! Other variables
    integer maxfct, mnum, mtype, n, nrhs, error, msglvl!, phase
 !   integer, dimension(64) :: iparm
    integer, dimension(:), allocatable :: ia,ja
    double precision, dimension(:), allocatable :: a ! b is actually bigv
 !   integer :: idum
!    double precision :: ddum		

        
        error = 0 ! initialize error flag
        msglvl = 0 ! print statistical information
        mtype = 11 ! real unsymmetric
        maxfct = 1
        mnum = 1
        nrhs = 1

    ! =======================================================================    
    !  Output files
    !
    !   In case of new run, the files are created, if not results will be 
    !   appended to old files
    !
        
    if (otime == 0.d0) then
        open(45,file='density_vs_time.txt',status='unknown')
        open(46,file='energy_vs_time.txt',status='unknown')
        open(470,file='power_coll_tot_vs_time.txt',status='unknown')
        do ib=1,nbulk
            if (ib==1) then
                open(471,file='power_coll_e_vs_time.txt',status='unknown')
            else
                write(ibString,'(i2)') ib-1
                open(470+ib,file='power_coll_ion'//ibString//'_vs_time.txt',status='unknown')
            endif
            
            
        enddo
        
        if (irf==-1) open(480,file='power_RF_vs_time.txt',status='unknown')
         
        if(isource==-1) open(490,file='power_NBI_vs_time.txt',status='unknown')
         
        if(isc /= 0) open(500,file='power_coll_self_vs_time.txt',status='unknown')
        
    else
        
        open(45,file='density_vs_time.txt',status='old', access = 'append')
        open(46,file='energy_vs_time.txt',status='old', access = 'append')
        open(470,file='power_coll_tot_vs_time.txt',status='old', access = 'append')
        do ib=1,nbulk
            if (ib==1) then
                open(471,file='power_coll_e_vs_time.txt',status='old', access = 'append')
            else
                write(ibString,'(i2)') ib-1
                open(470+ib,file='power_coll_ion'//ibString//'_vs_time.txt',status='old', access = 'append')
            endif
        enddo
        
        if (irf==-1) open(480,file='power_RF_vs_time.txt',status='old', access = 'append')
        if(isource==-1) open(490,file='power_NBI_vs_time.txt',status='old', access = 'append')
        if(isc /= 0) open(500,file='power_coll_self_vs_time.txt',status='old', access = 'append')
        
    endif
    
! =======================================================================    
 ! We start by computing the "steady-state" matrix
    
            
   call build_ss(all20_lin,all02_lin,all11_lin,all10_lin,all01_lin,all00_lin,bigm_ss)! we do not need bigv  
   
    
    ! When the self-collisions are not taken into account, this matrix is constant through the whole computation
    ! Only the rhs will change
    
    allocate(f_one(nbig),f_wrk(nbig))
    
    !	 

    nrhs = 1
    n = nbig
    
    
    time_loop: do itime= 1, ntimes
                         
        time = otime+itime*timestep
                
        write(*,*) 'Time is ',time, ' s'
        
    ! If self-collisions with non-Maxwellian are allowed (icn =-1) 
    !  the self-collision diffusion/fricton ion must be evaluated with the previous timestep vdf
    ! and a self-collision "BigM" must be computed and added to "BigM_SS"

        !if (isc == -1) then
        !    call main_nlterm(fstart,sc00,sc10,sc01,sc20,sc11,sc02)
        !                
        !    all00 = all00_lin+sc00
        !    all10 = all10_lin+sc10
        !    all01 = all01_lin+sc01
        !    all20 = all20_lin+sc20
        !    all02 = all02_lin+sc02
        !    all11 = all11_lin+sc11
        !     
        !    call build_ss(all20,all02,all11,all10,all01,all00,bigm_ss)
        !    
        !
        !endif
        
        
       bigm = 0.d0
       bigv = 0.d0       
            
    if(icn == -1) then
        call build_cn(fstart,bigm_ss,bigm,bigv)! Crank-Nicholson implicit scheme
    else
        call build_imp(fstart,bigm_ss,bigm,bigv)! fully implicit scheme
    endif
    
    
    first_timestep: if(itime == 1) then
        
! Settings
! --------       
      
! --> computes the number of non-zeros elements of Bigm
!      
nz=0         		
do i=1,nbig
	do j=1,nbig
	 if(bigm(i,j)/=0.d0)  then
		       nz=nz+1
	 endif
	enddo
enddo

write(*,*) 'Number of non-zeros in Bigm',nz, ' = ',100.d0*nz/nbig**2, ' %'
	   
! --> allocates memory space to working arrays
!
allocate(a(1:nz))
!allocate(x(1:nbig))
allocate(ia(1:(nbig+1)))
allocate(ja(1:nz))

! --> non-zeros elements of sparse Bigm matrix are stored
!
index=0
do i=1,nbig
	flag1 = 0
	do j=1,nbig
		if(bigm(i,j) /= 0.d0) then
			index=index+1
			if(flag1 == 0) then
				flag1 = 1
				ia(i)= index
			endif
		    a(index)=bigm(i,j)
			ja(index)=j
		endif
	enddo
enddo 
	   
ia(nbig+1) = nz+1

  !============================================================================
  !  Time-dependent, A CONSTANT  (b changes every step)
  !           Phase 11 + 22  done once in _init
  !           Phase 33        done each step in _step (a_changed=.FALSE.)
  !============================================================================


CALL pardiso_solve_init(handle, N, a, ia, ja, &
                          a_constant=.TRUE., mtype=11, msglvl=0, error=error)
  IF (error /= 0) STOP 'Init failed'
  
 endif first_timestep

   !-- a_changed = .FALSE. => skip phase 22, reuse factorisation
  
    CALL pardiso_solve_step(handle, a, ia, ja, bigv, xout, &
                            a_changed=.FALSE., error=error)
            
    ! =======================================================================    


do iv=1,nperp
        do imu=1,npar
        ix = index_mat(iv,imu)
        fout(iv,imu) = xout(ix)
        enddo
enddo

call time_density(fout,dens_tmp)
write(*,*) 'Unnormalised density is ',dens_tmp
write(45,*) time,dens_tmp

call time_energy(fout,dens_tmp,tk,tkperp,tkpar)
write(46,*) time,tk,tkperp

! Effective thermal velocity is computed (only uselful for self-collisions)
    
    vtk=9.79d3*dsqrt(tk*1.d3/aa)! tk is evaluated in keV --> convert to eV
    !write(*,*) 'Temperature is ',tk
    !write(*,*) 'Effective thermal velocity is ',vtk
    


call time_power(xout,dens_tmp,pcoll,pRF,psource,plosses,pcoll_self)

write(470,*) time,dabs(sum(pcoll)+pcoll_self)/1.d6
do ib=1,nbulk
    write(470+ib,*) time,(pcoll(ib)/1.d6)
enddo

if(irf==-1) write(480,*) time,pRF/1.d6
    
if(isource==-1) write(490,*) time,psource/1.d6,plosses/1.d6

if(isc==-1) write(500,*) time,pcoll_self/1.d6


fstart = xout*npart/dens_tmp
   
    enddo time_loop
    
    
    ! ======================================================
    
    if(isource==-1) close(490)
    
    if(irf==-1) close(480)
    
    if(isc==-1) close(500)
    
    close(470)
    close(46)
    close(45)
    
    ! Termination and release of memory
! =================================
  
    CALL pardiso_solve_finalize(handle, ia, ja, error)
    WRITE(*,*) 'Solve completed ... '


deallocate(a,ia)
    
! =======================================================================    
! For sourceless case the solution is re-normalised
    
!    if (isource == 0) then
!        
!write(*,*)''
!write(*,*) 'Renormalizing...'
!
!write(*,*)''
!
!
!fout = fout*npart/dens_tmp
!
!    endif
    
! =======================================================================    

    
    open(40,file='fout.txt', status='unknown')
    open(42,file='xout.dat',status='unknown')

    write(42,*) time

do iv=1,nperp
        do imu=1,npar
        write(40,*) vperp(iv), vpar(imu), fout(iv,imu)
         ix = index_mat(iv,imu)
         write(42,*) xout(ix)
        enddo
enddo

close(42)
close(41)
close(40)	

    
    end subroutine timefp
        
 ! =======================================================================    
    
end module mod_timefp3
    
    
