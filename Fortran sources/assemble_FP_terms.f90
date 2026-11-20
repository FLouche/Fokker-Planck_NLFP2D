! ***********************************************************
! * Assemble the terms of the linear Fokker-Planck equation *
! *                                                         *
! *   Version 1: 28/05/2026 - F.Louche                      *
! *                                                         *
! ***********************************************************

module assemble_FP_lin

    contains
    
    subroutine assemble_FP_terms(all00,all10,all01,all20,all11,all02)
    
    use shared_grid
    use shared_plasma
    use shared_rf
    use shared_FPterms
    
    implicit none
    
    double precision, dimension(nperp,npar), intent(out) :: all00,all10,all01,all20,all11,all02
    
    integer ib,ipe
    
    ! =========================================================
    ! Linear collision terms: accumulate gammab*colin**_sp without
    ! modifying colin**_sp so this routine can be called repeatedly.

    colin20 = 0.d0; colin02 = 0.d0; colin11 = 0.d0
    colin10 = 0.d0; colin01 = 0.d0; colin00 = 0.d0
    
  !      open(444,file='colin_D.txt',status='unknown')


    do ib = 1, nbulk
        colin20 = colin20 + gammab(ib)*colin20_sp(:,:,ib)
        colin02 = colin02 + gammab(ib)*colin02_sp(:,:,ib)
        colin11 = colin11 + gammab(ib)*colin11_sp(:,:,ib)
        colin10 = colin10 + gammab(ib)*colin10_sp(:,:,ib)
        colin01 = colin01 + gammab(ib)*colin01_sp(:,:,ib)
        colin00 = colin00 + gammab(ib)*colin00_sp(:,:,ib)
                
    enddo
    
   !         do ipe=1,nperp
  !          write(444,*) vperp(ipe),gammab(2)*colin20_sp(ipe,imid,2)
  !      enddo

    
  !  close(444)
    
 !   write(*,*) 'Gamma for D ',gammab(2)
    
    ! =========================================================
    ! RF heating term is added 
    
    if (irf == -1) then
    
    all00 = colin00
    all10 = colin10+rf10
    all01 = colin01+rf01
    all20 = colin20+rf20
    all02 = colin02+rf02
    all11 = colin11+rf11

    else
        
    all00 = colin00
    all10 = colin10
    all01 = colin01
    all20 = colin20
    all02 = colin02
    all11 = colin11

    endif

    ! =========================================================
    ! Energy-confinement loss L_E (multispecies extension, namelist tau_E;
    ! main sets tau_E = 0 when there is no source, so this is skipped then).
    call add_lossE_terms(all00, all10, all01)

    end subroutine assemble_FP_terms

    !-----------------------------------------------------------------
    ! L_E f = (1/2 tau_E) [ (1/vperp) d(vperp^2 f)/dvperp + d(vpar f)/dvpar ]
    !       = (1/2 tau_E) [ 3 f + vperp df/dvperp + vpar df/dvpar ]
    ! i.e.  A += 3/(2 tau_E),  B += vperp/(2 tau_E),  C += vpar/(2 tau_E).
    ! Conserves particles; dW/dt = -W/tau_E.  Nothing is added when
    ! tau_E = 0, so the default reproduces v3.2 to the bit.
    !-----------------------------------------------------------------
    subroutine add_lossE_terms(a00, a10, a01)

    use shared_grid
    use shared_multi, only: tau_E

    implicit none

    double precision, dimension(nperp,npar), intent(inout) :: a00, a10, a01
    double precision :: r
    integer :: i, j

    if (tau_E <= 0.d0) return
    r = 0.5d0 / tau_E
    do j = 1, npar
        do i = 1, nperp
            a00(i,j) = a00(i,j) + 3.d0 * r
            a10(i,j) = a10(i,j) + vperp(i) * r
            a01(i,j) = a01(i,j) + vpar(j) * r
        enddo
    enddo

    end subroutine add_lossE_terms
    
end module assemble_FP_lin

    
 
        