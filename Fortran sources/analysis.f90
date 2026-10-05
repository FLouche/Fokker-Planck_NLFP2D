
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
    double precision teff,Tperp,Tpar,anisotropy_perp
    
    integer iv,ip
    
    
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
    
    write(*,*) 'Average Kinetic Energy is ', 0.5*pmass*aa*mod2/kev_in_J, 'keV'
    
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
    
    Tperp = 0.5*pmass*aa*mod2_perp/kev_in_J
    write(*,*) 'Perpendicular Kinetic Energy is ',Tperp , 'keV'
    
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
    
    ! Tpar = m<v_par^2> (plasma convention: no factor ½).
    ! True parallel kinetic energy E_par = ½m<v_par^2> = 0.5*Tpar.
    Tpar = pmass*aa*mod2_par/kev_in_J

    write(*,*) 'Parallel Kinetic Energy is ', 0.5d0*Tpar, 'keV'

    ! Teff = (2/3)*E_total = (2*Tperp + Tpar)/3 is correct given the above.
    Teff =(2.d0*Tperp+Tpar)/3.d0

    write(*,*) 'Effective temperature is ',Teff, 'keV'

    ! Perpendicular anisotropy: 0% (all parallel) -> 50% (Maxwellian) -> 100% (all perp).
    ! = 100 * E_perp / (E_perp + 2*E_par) = 100 * Tperp / (Tperp + Tpar).
    anisotropy_perp = merge(100.d0*Tperp/(Tperp + Tpar), 0.d0, (Tperp + Tpar) > 0.d0)
    write(*,*) 'Perpendicular anisotropy is    ', anisotropy_perp, '%  (50% = Maxwellian)'


    open(40,file=TRIM(outfile('Ekin_par.txt')), status='unknown')


do iv=1,nperp
        do ip=1,npar
        write(40,*) vperp(iv), vpar(ip), fint(iv,ip)*pmass*aa/kev_in_J
        enddo
enddo

close(40)
    
! We test the balance of particle density and collisional distributed power on background species

allocate(xout(nbig))
xout=0.d0

do iv = 1,nperp
        do ip = 1,npar
            
            xout(index_mat(iv,ip)) = f(iv,ip)
            
        enddo
enddo

call test_power_balance_7pt(xout)
call test_density_balance_7pt(xout)
call test_momentum_balance_7pt(xout)



deallocate(xout)

    
    !***************************************

    
    end subroutine analysis
    
end module mod_anal    