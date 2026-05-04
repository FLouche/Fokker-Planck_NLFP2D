
!   analysis.f90
!
!   FPColl_2D
!   
!   Created by Fabrice Louche on 18/08/25.
!
!    Modified on 19/11/2025: cylindrical coordinates
!
!   Copyright 2025 LPP-ERM/KMS. All rights reserved.
    
!***************************************

module mod_anal

    contains
    
    subroutine analysis(f)
    
    use shared_grid
    use shared_plasma
    use shared_beam
    use shared_timer
        
    use mod_ncint
    use func_index
    
    implicit none
    
    double PRECISION, dimension(nperp,npar), intent(in) :: f
    double PRECISION, dimension(nperp,npar) :: fint
    
    double precision, allocatable, dimension(:) :: xout

    
    double precision pmass, kev_in_J
    double PRECISION mod0,mod2,mod2_perp,mod2_par
    
    integer iv,ip
    
    external test_density_balance, test_power_balance
    
    data pmass/1.6726d-27/ !proton mass in kg
    data kev_in_J/1.60218d-16/ !conperpert keV to Joule
    
    ! =====================================================
    
    ! Particle Density
    
    do iv = 1,nperp
        do ip = 1,npar
            fint(iv,ip) = f(iv,ip)*jacob(iv,ip)
        enddo
    enddo
    
    call ncint_2d(fint,mod0)
    
    write(*,*) 'Density is ', mod0, ' /m3'
    
     ! =====================================================
    
    ! Particle kinetic energy
    
    ! -> total
    
    do iv = 1,nperp
        do ip = 1,npar
            fint(iv,ip) = f(iv,ip)*(vperp(iv)**2+vpar(ip)**2)*jacob(iv,ip)/mod0
        enddo
    enddo
    
    call ncint_2d(fint,mod2)
    
    write(*,*) 'Total Kinetic Energy is ', 0.5*pmass*aa*mod2/kev_in_J, 'keV'
    
    open(40,file=TRIM(outfile('Ekin.txt')), status='unknown')


do iv=1,nperp
        do ip=1,npar
        write(40,*) vperp(iv), vpar(ip), fint(iv,ip)*0.5*pmass*aa/kev_in_J
        enddo
enddo

close(40)


    
    ! -> Perpendicular energy
    
     do iv = 1,nperp
        do ip = 1,npar
            fint(iv,ip) = f(iv,ip)*vperp(iv)**2*jacob(iv,ip)/mod0
        enddo
     enddo
     
    
    call ncint_2d(fint,mod2_perp)
    
    write(*,*) 'Perpendicular Kinetic Energy is ', 0.5*pmass*aa*mod2_perp/kev_in_J, 'keV'
    
    open(40,file=TRIM(outfile('Ekin_perp.txt')), status='unknown')
    open(41,file=TRIM(outfile('Ekin_perp_at_vpar0.txt')),status='unknown')


do iv=1,nperp
    write(41,*) vperp(iv), fint(iv,jmid)*0.5*pmass*aa/kev_in_J
        do ip=1,npar
        write(40,*) vperp(iv), vpar(ip), fint(iv,ip)*0.5*pmass*aa/kev_in_J
        enddo
enddo

close(41)
close(40)
!    
     ! -> Parallel energy
    
     do iv = 1,nperp
        do ip = 1,npar
            fint(iv,ip) = f(iv,ip)*vpar(ip)**2*jacob(iv,ip)/mod0
        enddo
     enddo
     
     

    call ncint_2d(fint,mod2_par)
    
    write(*,*) 'Parallel Kinetic Energy is ', 0.5*pmass*aa*mod2_par/kev_in_J, 'keV'
    
!    open(40,file='Ekin_par.txt', status='unknown')
!
!
!do iv=1,nperp
!        do ip=1,npar
!        write(40,*) v(iv), mu(ip), fint(iv,ip)*0.5*pmass*aa/kev_in_J
!        enddo
!enddo
!
!close(40)
    
! We test the balance of particle density and collisional distributed power on background species

allocate(xout(nbig))
xout=0.d0

do iv = 1,nperp
        do ip = 1,npar
            
            xout(index_mat(iv,ip)) = f(iv,ip)
            
        enddo
enddo

if(ifd7 == -1) then
    call test_density_balance_7pt(xout)
    call test_power_balance_7pt(xout)
else
    call test_density_balance(xout)
    call test_power_balance(xout)
endif



deallocate(xout)



! -------------------------------------------------------------------------------------
! TEST: shape of the vdf along boundaries
!
    

! Curves of F along boundaries

!open(40,file='F_at_vperp0.txt',status='unknown')
!
!do ip=1,npar
!    write(40,*) vpar(ip),f(1,ip)
!enddo
!
!close(40)
!
!open(40,file='F_at_vparmin.txt',status='unknown')
!
!do iv=1,nperp
!    write(40,*) vperp(iv),f(iv,1)
!enddo
!
!close(40)
!
!open(40,file='F_at_vparmax.txt',status='unknown')
!
!do iv=1,nperp
!    write(40,*) vperp(iv),f(iv,npar)
!enddo
!
!close(40)

    
    !***************************************

    
    end subroutine analysis
    
end module mod_anal    