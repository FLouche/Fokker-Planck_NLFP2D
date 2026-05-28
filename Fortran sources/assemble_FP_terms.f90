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
    ! Linear collision terms
    
    do ib=1,nbulk
        
        colin20_sp(:,:,ib) = gammab(ib)*colin20_sp(:,:,ib)
        colin02_sp(:,:,ib) = gammab(ib)*colin02_sp(:,:,ib)
        colin11_sp(:,:,ib) = gammab(ib)*colin11_sp(:,:,ib)
        colin10_sp(:,:,ib) = gammab(ib)*colin10_sp(:,:,ib)
        colin01_sp(:,:,ib) = gammab(ib)*colin01_sp(:,:,ib)
        colin00_sp(:,:,ib) = gammab(ib)*colin00_sp(:,:,ib)
        
    enddo
    
    colin20 = sum(colin20_sp, DIM = 3)

    colin02 = sum(colin02_sp, DIM = 3)

    colin11 = sum(colin11_sp, DIM = 3)

    colin10 = sum(colin10_sp, DIM = 3)

    colin01 = sum(colin01_sp, DIM = 3)

    colin00 = sum(colin00_sp, DIM = 3)
    
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

    
 
        