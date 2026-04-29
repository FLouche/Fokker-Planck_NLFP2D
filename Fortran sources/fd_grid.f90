module fd_grid
    
    contains
    
! ====================================================================
 !
 ! Evaluation of the 1st derivative along vperp
 ! of a function given on 2D grid
 
 subroutine finite_diff_1D_vperp(f,dfpe)
 
 use shared_grid

 implicit none
 
 double precision, dimension(nperp,npar), intent(in):: f
 double precision, dimension(nperp,npar), intent(out):: dfpe

 integer i
 
 double precision alpha_l,alpha_r,alpha_n
 
 do i=1,nperp
     
         if (ising == -1) then
    
          if(i == nsing) then
            
            alpha_l = dvleft/dvright/(dvleft+dvright)
            alpha_r = dvright/dvleft/(dvleft+dvright)
            alpha_n = (dvleft-dvright)/dvleft/dvright
                
         else
            
            if(i < nsing) then
                
                alpha_l = 1.d0/2.d0/dvleft
                alpha_r = 1.d0/2.d0/dvleft
                alpha_n = 0.d0
                
            else
                
                alpha_l = 1.d0/2.d0/dvright
                alpha_r = 1.d0/2.d0/dvright
                alpha_n = 0.d0
                                
            endif
            
        endif
        
    else
        
        alpha_l = 1.d0/2.d0/dvperp
        alpha_r = 1.d0/2.d0/dvperp
        alpha_n = 0.d0

    endif
    
    if (i == 1) then
        
        dfpe(i,:) = alpha_l*(4.d0*f(i+1,:)-3.d0*f(i,:)-f(i+2,:))
        
    else 
        if (i == nperp) then
            
            dfpe(i,:) = alpha_l*(-4.d0*f(i-1,:)+3.d0*f(i,:)+f(i-2,:))
            
        else
            
            dfpe(i,:) = alpha_r*f(i+1,:)-alpha_l*f(i-1,:)-alpha_n*f(i,:)
            
        endif
    endif
    
 enddo
 
 
 end subroutine finite_diff_1D_vperp
 
  ! ====================================================================
 !
 ! Evaluation of the 1st derivative along vpar
 ! of a function given on 2D grid
 
 subroutine finite_diff_1D_vpar(f,dfpa)
 
 use shared_grid

 implicit none
 
 double precision, dimension(nperp,npar), intent(in):: f
 double precision, dimension(nperp,npar), intent(out):: dfpa

 integer i
 
 double precision alpha_l
 
 
  do i =1,npar
     
     alpha_l=1.d0/2.d0/dvpar
     
         if (i == 1) then
        
        dfpa(i,:) = alpha_l*(4.d0*f(:,i+1)-3.d0*f(:,i)-f(:,i+2))
        
    else 
        if (i == npar) then
            
            dfpa(i,:) = alpha_l*(-4.d0*f(:,i-1)+3.d0*f(:,i)+f(:,i-2))
            
        else
            
            dfpa(i,:) = alpha_l*(f(:,i+1)-f(:,i-1))
            
        endif
    endif

 enddo

 end subroutine finite_diff_1D_vpar
 
   ! ====================================================================
 !
 !! Evaluation of the 1st derivative along v
 !! of a function given on 1D grid
 !
 !subroutine finite_diff_1D_v(f,dfv)
 !
 !use spherical_grid
 !
 !implicit none
 !
 !double precision, dimension(nv), intent(in):: f
 !double precision, dimension(nv), intent(out):: dfv
 !
 !integer i
 !
 !double precision alpha_l
 !
 !
 ! do i =1,nv
 !    
 !    alpha_l=1.d0/2.d0/delta_v
 !    
 !        if (i == 1) then
 !       
 !       dfv(i) = alpha_l*(4.d0*f(i+1)-3.d0*f(i)-f(i+2))
 !       
 !   else 
 !       if (i == nv) then
 !           
 !           dfv(i) = alpha_l*(-4.d0*f(i-1)+3.d0*f(i)+f(i-2))
 !           
 !       else
 !           
 !           dfv(i) = alpha_l*(f(i+1)-f(i-1))
 !           
 !       endif
 !   endif
 !
 !enddo
 !
 !end subroutine finite_diff_1D_v
 
 end module fd_grid
