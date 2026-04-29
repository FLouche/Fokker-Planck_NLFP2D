include 'mkl_spblas.f90'

    subroutine sparse_matrix_vect_mult(mat,vec,y)
    
    use shared_grid
    use mkl_spblas

    implicit none
        
    double PRECISION, dimension(nbig,nbig), intent(in) :: mat
    double PRECISION, dimension(nbig), intent(in) :: vec
    double PRECISION, dimension(nbig), intent(out) :: y
    
    integer, dimension(nbig*nbig) :: rows_nz_temp,cols_nz_temp
    double PRECISION, dimension(nbig*nbig) :: values_temp
    integer, allocatable,dimension(:) :: rows_nz,cols_nz
    double PRECISION, allocatable, dimension(:) :: values
    
    type(sparse_matrix_t) :: b
    type(matrix_descr) :: descr

    
    integer:: i1,i2, nnz, stat
        
    ! -------------------------------------------------------
    ! Non-zero elements of array are counted and localized
    
    nnz = 0
    
    do i1=1,nbig
        do i2=1,nbig
            if(mat(i1,i2) /= 0.d0) then
                nnz = nnz+1
                rows_nz_temp(nnz) = i1
                cols_nz_temp(nnz) = i2
                values_temp(nnz) = mat(i1,i2)
            endif
        enddo
    enddo
    
  !  write(*,*) 'Number of non-zero elements = ',nnz
    
    allocate(rows_nz(nnz),cols_nz(nnz),values(nnz))
    
    do i1=1,nnz
        rows_nz(i1) = rows_nz_temp(i1)
        cols_nz(i1) = cols_nz_temp(i1)
        values(i1)  = values_temp(i1)
    enddo
    
    ! We create the sparse matrix handle in coordinates storage
    
    stat = mkl_sparse_d_create_coo(b,SPARSE_INDEX_BASE_ONE,nbig,nbig,nnz,rows_nz,cols_nz,values)
   ! print *, "stat create = ", stat

    descr%type = SPARSE_MATRIX_TYPE_GENERAL
    
    stat = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE,1.d0,b,descr,vec,0.d0,y)
   ! print *, "stat mv = ", stat
        
    end subroutine sparse_matrix_vect_mult
    