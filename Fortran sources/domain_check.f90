!====================================================================
!
!   domain_check.f90
!
!   FPColl_2D
!
!   Domain adequacy check, called once at the end of the run (moved out
!   of main-FP_Coll_2D.f90, where it was the block domain_check).
!
!====================================================================

module mod_domain_check

implicit none

contains

!--------------------------------------------------------------------
! Domain adequacy check
! ---------------------
!
! The Dirichlet walls impose f=0, while the true solution there is
!   f(wall)/f(peak) = exp(-(V/vth)^2/2).
! Whatever that ratio is, it is the floor on the solution error, and it does
! NOT improve with grid refinement.  The check must use the FINAL (heated)
! temperature, not the cold background one: an RF or beam case can be perfectly
! well sized for its initial Maxwellian and badly under-sized for the
! distribution it actually evolves into.
!
! Target ratio 1e-12 (i.e. |v|max ~ 7.4 vth) was measured to give a shape error
! ~2e-8, against ~3e-5 at 5 vth.
!--------------------------------------------------------------------
subroutine domain_adequacy_check(fout)

    use shared_grid,    only: nperp, npar, vperp_max, vpar_min, vpar_max
    use shared_plasma,  only: aa
    use time_comps_mod, only: time_density, time_energy

    double precision, intent(in) :: fout(nperp, npar)

    double precision :: dens_chk, teff_chk, vth_chk, v_req
    double precision :: r_perp, r_par
    double precision, parameter :: f_target = 1.d-12

    call time_density(fout, dens_chk)
    call time_energy(fout, dens_chk, teff=teff_chk)

    if (teff_chk > 0.d0) then

        vth_chk = 9.79d3*dsqrt(teff_chk*1.d3/aa)
        v_req   = vth_chk*dsqrt(2.d0*dlog(1.d0/f_target))

        r_perp = vperp_max/vth_chk
        r_par  = min(dabs(vpar_min), dabs(vpar_max))/vth_chk

        write(*,*) ' '
        write(*,*) 'DOMAIN ADEQUACY CHECK (Dirichlet wall)'
        write(*,*) '--------------------------------------'
        write(*,'(A,F9.3,A,ES10.3,A)') '  final Teff  = ', teff_chk, &
             ' keV   ->  vth = ', vth_chk, ' m/s'
        write(*,'(A,ES9.2,A,ES10.3,A)') '  for f(wall)/f(peak) < ', f_target, &
             '  recommend |v|max > ', v_req, ' m/s'
        write(*,'(A,ES10.3,A,F6.2,A,ES9.2,A)') '  vperp_max   = ', vperp_max, &
             '  (', r_perp, ' vth)  f(wall)/f(peak) = ', &
             dexp(-0.5d0*r_perp**2), '  '//trim(merge('OK       ', &
             'TOO SMALL', vperp_max >= v_req))
        write(*,'(A,ES10.3,A,F6.2,A,ES9.2,A)') '  |vpar|max   = ', &
             min(dabs(vpar_min),dabs(vpar_max)), &
             '  (', r_par, ' vth)  f(wall)/f(peak) = ', &
             dexp(-0.5d0*r_par**2), '  '//trim(merge('OK       ', &
             'TOO SMALL', min(dabs(vpar_min),dabs(vpar_max)) >= v_req))

        if (vperp_max < v_req .or. &
            min(dabs(vpar_min),dabs(vpar_max)) < v_req) then
            write(*,*) ' '
            write(*,'(A,ES10.3,A)') '  ==> the domain truncates the tail. '// &
                 'Re-run with |v|max >= ', v_req, ' m/s;'
            write(*,*) '      refining the grid will NOT reduce this error.'
        end if

    end if

end subroutine domain_adequacy_check

end module mod_domain_check
