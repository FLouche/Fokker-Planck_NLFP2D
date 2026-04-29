!**********************************************************
!
! This routine builds the matrix containing the terms
!  of the FP equation written on a 2D grid using a
!   finite-differences differentiation scheme
!
!   This version includes:
!
!    - homogeneous grid only (ising = 0)
!    - only Coulomb collisions (no df/dvdmu term)
!   - beam source
!    - steady-state equation 
!
!    by Fabrice Louche on 20/08/2025
!
!   Update 14/01/26:
!
!  - inhomogeneous grids in vperp
!  - RF term included
!
!**********************************************************

module mod_build_ss_1

    contains
    
subroutine build_ss_1(all20,all02,all11,all10,all01,all00,bigm,bigv)

! Initialisation
! --------------

use shared_grid
use shared_plasma
use shared_beam 

use shared_timer
use shared_rf

use func_index

implicit none

double precision, dimension(nperp,npar),intent(in) :: all20,all02,all11,all10,all01,all00
double precision,  dimension (nbig,nbig), intent(out) :: bigm
double precision, dimension (nbig),intent(out), optional :: bigv

double precision alpha_l,alpha_r,alpha_n,beta_l,beta_r,beta_n

integer iv, ip, ix1, ix2

!====================================================================
!
! We write the equation for dfij/dt (i -> velocity, j -> pitch-angle)
!   
!  The equations are ordered according to the following rule:
!
!     dfij/dt -> line (i-1)*npar + j
!
!  The column are filled following the rule:
!
!    factor in front fij -> column (i-1)*npar + j for a given equation 
    

imid = 1!int(4*nperp/5)+1

v_loop: do iv = 1,nperp
    
    if (ising == -1) then
    
        if(iv == nsing) then
            
            alpha_l = dvleft/dvright/(dvleft+dvright)
            alpha_r = dvright/dvleft/(dvleft+dvright)
            alpha_n = (dvleft-dvright)/dvleft/dvright
    
            beta_r = 2.d0/dvright/(dvleft+dvright)
            beta_l = 2.d0/dvleft/(dvleft+dvright)
            beta_n = 2.d0/dvleft/dvright
            
        else
            
            if(iv < nsing) then
                
                alpha_l = 1.d0/2.d0/dvleft
                alpha_r = 1.d0/2.d0/dvleft
                alpha_n = 0.d0
                
                beta_r = 1.d0/dvleft**2
                beta_l = 1.d0/dvleft**2
                beta_n = 2.d0/dvleft**2
                
            else
                
                alpha_l = 1.d0/2.d0/dvright
                alpha_r = 1.d0/2.d0/dvright
                alpha_n = 0.d0
                
                beta_r = 1.d0/dvright**2
                beta_l = 1.d0/dvright**2
                beta_n = 2.d0/dvright**2
                
            endif
            
        endif
        
    else
        
        alpha_l = 1.d0/2.d0/dvperp
        alpha_r = 1.d0/2.d0/dvperp
        alpha_n = 0.d0
        
        beta_r = 1.d0/dv2
        beta_l = 1.d0/dv2
        beta_n = 2.d0/dv2
        
    endif    
    
    mu_loop: do ip = 1,npar
            
        ix1 = index_mat(iv,ip)
        
        BC_vmax: if (iv==nperp.or.ip==1.or.ip==npar) then ! at v = vmax we impose f = 0
            
                    bigm(ix1,ix1) = 1.d0!all00(iv,ip)
                    if(present(bigv)) bigv(ix1) = 0
                    
        else BC_vmax
            
             BC_vmin: if(iv == 1) then
                 
                 ! at v=0 we impose df/dvperp=0
                                
                    bigm(ix1,ix1) = -1.5d0
                                    
                    ix2 = index_mat(iv+1,ip)
                    bigm(ix1,ix2) = 2.d0
                                
                    ix2 = index_mat(iv+2,ip)
                    bigm(ix1,ix2) = -0.5d0
                                    
                    if(present(bigv)) bigv(ix1) = 0
                                    
                                
                else BC_vmin
            
            !                    
            !                    ! General formulation in the inside domain
            !                
                                bigm(ix1,ix1) = all00(iv,ip)-alpha_n*all10(iv,ip)-beta_n*all20(iv,ip)-2.d0*all02(iv,ip)/dmu2-taum
                                                                
                                ix2 = index_mat(iv+1,ip)
                                bigm(ix1,ix2) = alpha_l*all10(iv,ip)+beta_r*all20(iv,ip)
                                                                
                                ix2 = index_mat(iv-1,ip)
                                bigm(ix1,ix2) = beta_l*all20(iv,ip)-alpha_r*all10(iv,ip)
                                
                                ix2 = index_mat(iv,ip+1)
                                bigm(ix1,ix2) = all02(iv,ip)/dmu2+all01(iv,ip)/2.d0/dvpar-alpha_n*all11(iv,ip)/2.d0/dvpar                                
                                
                                ix2 = index_mat(iv,ip-1)
                                bigm(ix1,ix2) = all02(iv,ip)/dmu2-all01(iv,ip)/2.d0/dvpar+alpha_n*all11(iv,ip)/2.d0/dvpar
                                
                                ix2 = index_mat(iv-1,ip-1)
                                bigm(ix1,ix2) = alpha_r*all11(iv,ip)/2.d0/dvpar
                                
                                ix2 = index_mat(iv-1,ip+1)
                                bigm(ix1,ix2) = -alpha_r*all11(iv,ip)/2.d0/dvpar
                                
                                ix2 = index_mat(iv+1,ip-1)
                                bigm(ix1,ix2) = -alpha_l*all11(iv,ip)/2.d0/dvpar
                                
                                ix2 = index_mat(iv+1,ip+1)
                                bigm(ix1,ix2) = alpha_l*all11(iv,ip)/2.d0/dvpar
                                
                                if(present(bigv))bigv(ix1) = -source(iv,ip)
                                
            !                    
                
                endif BC_vmin
        
        endif BC_vmax      
        
        

    enddo mu_loop
    
                
enddo v_loop
          
! Sourceless steady-state case -> we impose f = 1 somewhere (imid, jmid)

source_term: if (isource == 0.AND.PRESENT(bigv)) then  !.ntimes==0.and

    ix1 = index_mat(imid,jmid)
    
    do ix2= 1,nbig
        
       bigm(ix1,ix2) = 0.d0
       
    enddo
                        
    bigm(ix1,ix1) = 1.d0
                            
    bigv(ix1) = 1.d0
    
    write(*,*) 'Imposing f=1 somewhere...'
    
        
endif source_term
        
                    
!====================================================================

end subroutine build_ss_1

!**********************************************************


end module mod_build_ss_1
