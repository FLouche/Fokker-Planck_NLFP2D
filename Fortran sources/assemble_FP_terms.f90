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
    
    integer ib
    
    ! =========================================================
    ! Linear collision terms: accumulate gammab*colin**_sp without
    ! modifying colin**_sp so this routine can be called repeatedly.

    colin20 = 0.d0; colin02 = 0.d0; colin11 = 0.d0
    colin10 = 0.d0; colin01 = 0.d0; colin00 = 0.d0

    do ib = 1, nbulk
        colin20 = colin20 + gammab(ib)*colin20_sp(:,:,ib)
        colin02 = colin02 + gammab(ib)*colin02_sp(:,:,ib)
        colin11 = colin11 + gammab(ib)*colin11_sp(:,:,ib)
        colin10 = colin10 + gammab(ib)*colin10_sp(:,:,ib)
        colin01 = colin01 + gammab(ib)*colin01_sp(:,:,ib)
        colin00 = colin00 + gammab(ib)*colin00_sp(:,:,ib)
    enddo
    
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


    end subroutine assemble_FP_terms
    
end module assemble_FP_lin

    
 
        