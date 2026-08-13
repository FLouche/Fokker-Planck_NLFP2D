!***********************************************************************
!
! Test module to derive the self-collisions FP coefficients
! in terms of the Rosenbluth-Trubnikov potentials 
! 
! Method
!
! - elements of diffusion and frictions evaluated from Karney expressions
!       (23 a, b + 24 b, c) for a general background by numerical 2D integration
! - derivatives evaluated numerically with finite-differences
!
!***********************************************************************

    module nlterm
    
    implicit none
            
    contains
    
    subroutine main_nlterm(xout,teff,sc00,sc10,sc01,sc20,sc11,sc02)
    
    use shared_grid
    use mod_grid
    use shared_timer
    use mod_ncint
    
    use func_index
    
    use shared_plasma
    
    use derivatives_2d
    
    use coulomb_log_mod
    
    double precision, intent(in), dimension (nbig):: xout
    double precision, intent(in) :: teff
    double precision, dimension (nperp,npar):: fout,psi,dpsidpe,d2psidpe2,d2psidpepa
    double precision, dimension (nperp,npar):: dpsidpa,d2psidpa2,phi,dphidpe,dphidpa,d2phidpe2,d2phidpa2
    double precision, dimension (nperp,npar):: d3psidpe3, d3psidpepa2,d3psidpe2pa,d3psidpa3
    double precision, dimension (nperp,npar):: Lphi_radial
    double precision, intent(out), dimension (nperp,npar):: sc00,sc10,sc01,sc20,sc11,sc02
        
!integer, parameter :: dp = selected_real_kind(15, 307)
    integer ipe,ipa,ix
        
    !double precision, allocatable, dimension(:,:) :: dpepe,dpepa,dpapa,fpe,fpa
    !double precision, allocatable, dimension(:,:) :: dfpe,dfpa,d_dpepe_1,d_dpapa_2,d_dpepa_1,d_dpepa_2
        
    common/mathcons/pi,twopi
    
    double precision pi,twopi,gamma0,cte0
    double precision coef,ta_ev,lnaa

    data gamma0/2.390775d-1/ ! This is e^4/(4pi Eps0^2 mp^2) 

    external cblin

    do ipe=1,nperp
        do ipa=1,npar
         ix = index_mat(ipe,ipa)
         fout(ipe,ipa) = xout(ix)
        enddo
    enddo
    
    ! The first Rosenbluth potential is computed on the grid
    
  call compute_psi(fout,psi)      

  
    ! Various derivatives are evaluated
    
!! dPsi/Dvperp
!    

  call deriv_x1(psi,vperp,nperp,npar,dpsidpe)
       
!
!! d2Psi/Dvperp2
!

  call deriv_x2(psi,vperp,nperp,npar,d2psidpe2)
  
  !
!! d3Psi/Dvperp3
!

  call deriv_x3(psi,vperp,nperp,npar,d3psidpe3)

   
!
!! d2Psi/DvperpDvpar
!    
 call deriv_xy(psi, vperp, nperp, npar, dvpar, d2psidpepa)!  
 
 
 !
!! d3Psi/DvperpDvpar2
!    
 call deriv_y1(d2psidpepa,  nperp, npar, dvpar, d3psidpepa2)!  
 
  !
!! d3Psi/Dvperp2Dvpar
!     
 call deriv_x1(d2psidpepa, vperp, nperp, npar, d3psidpe2pa)!  


 
!! dPsi/Dvpar
!    

 call  deriv_y1(psi, nperp, npar, dvpar, dpsidpa)
  
  !    
!! d2Psi/Dvpar2
!    
 call  deriv_y2(psi, nperp, npar, dvpar, d2psidpa2)
 
 !! d3Psi/Dvpar3
!    
 call  deriv_y1(d2psidpa2, nperp, npar, dvpar, d3psidpa3)

  
!     
!! The second potential phi is obtained from Poisson's equation
!    
! OLD (singular at vperp=0):
!do ipe = 1, nperp
!  phi(ipe,:) = d2psidpe2(ipe,:) + dpsidpe(ipe,:)/vperp(ipe) + d2psidpa2(ipe,:)
!enddo

! NEW (regularised):
call laplacian_radial(psi, vperp, nperp, npar, Lphi_radial)
phi(:,:) = Lphi_radial(:,:) + d2psidpa2(:,:)

!n_bad = 10
!n_fit = 14
!call regularise_axis(phi, vperp, nperp, npar, n_bad, n_fit)

! First pass: raw d2phi used only for detection
CALL deriv_x2(phi, vperp, nperp, npar, d2phidpe2)

! Detect n_bad and correct phi in-place
CALL regularise_axis_3(phi, d2phidpe2, vperp, nperp, npar)

! Second pass: recompute derivatives on the corrected phi
CALL deriv_x1(phi, vperp, nperp, npar, dphidpe)
CALL deriv_x2(phi, vperp, nperp, npar, d2phidpe2)


!    
!    ! Dphi/Dvpar
!    
call deriv_y1(phi, nperp, npar, dvpar, dphidpa)    

!    
!    ! D2phi/Dvpar2
!    
call deriv_y2(phi, nperp, npar, dvpar, d2phidpa2)    


! -------------------------------------------------------------------------
! Components of the diffusion/friction for self-collisions
! --------------------------------------------------------
!

    !
    !! Factors in front of each term coming from self-collisions are
    !! evaluated
    !
    !

! The gammaa factor depends on the time varying Coulomb logarithm and should evaluated here

ta_eV = teff * 1.0d3

call coulomb_log_ab(za, aa, ta_eV, npart, za, aa, ta_eV, npart, lnaa)

cte0=gamma0*lnaa*(za/aa)**2
 
 gammaa = cte0*npart*za**2

 coef = -4.d0*pi*gammaa/npart

    do ipa=1,npar

        do ipe=1,nperp

            sc00(ipe,ipa) = -coef*(dphidpe(ipe,ipa)/vperp(ipe)+d2phidpe2(ipe,ipa)+d2phidpa2(ipe,ipa))
            sc10(ipe,ipa) = coef*(d2psidpe2(ipe,ipa)/vperp(ipe)+d3psidpe3(ipe,ipa)+d3psidpepa2(ipe,ipa)-dphidpe(ipe,ipa))
            sc01(ipe,ipa) = coef*(d2psidpepa(ipe,ipa)/vperp(ipe)+d3psidpe2pa(ipe,ipa)+d3psidpa3(ipe,ipa)-dphidpa(ipe,ipa))
           !!
            sc20(ipe,ipa) = coef*d2psidpe2(ipe,ipa)
            sc11(ipe,ipa) = coef*2.d0*d2psidpepa(ipe,ipa)
            sc02(ipe,ipa) = coef*d2psidpa2(ipe,ipa)


        enddo
    enddo


    end subroutine main_nlterm
    
   !***********************************************************************
 
!=======================================================================
! Compute the first Rosenbluth potential:
!     Psi(v) = -1/(8π) ∫ |v-v'| f(v') d^3v'
! in cylindrical velocity-space (v⊥, v∥) with axisymmetry.
!
! Replaces the original O(N^4) quadruple loop + per-point integrate_2d
! call with a single kernel apply: psi_vec = prefactor * (K g),
! where g(ix2) = f(ip,jp)*vperp(ip)*w_vperp(ip)*w_vpar(jp).
! Weights are the standard 2D trapezoidal rule, identical to integrate_2d.
!=======================================================================

subroutine compute_psi(f_values, psi_values)

  use shared_grid
  use func_index
  use mod_phi_kernel

  implicit none

  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: prefactor = -1.0_dp / (32.0_dp * atan(1.0_dp))

  double precision, dimension(nperp,npar), intent(in)  :: f_values
  double precision, dimension(nperp,npar), intent(out) :: psi_values

  double precision :: g(nbig), psi_vec(nbig)
  double precision :: w_vperp(nperp), w_vpar(npar)
  integer :: i, j

  ! Trapezoidal weights for v⊥ (non-uniform grid)
  w_vperp(1) = 0.5_dp * (vperp(2) - vperp(1))
  do i = 2, nperp-1
    w_vperp(i) = 0.5_dp * (vperp(i+1) - vperp(i-1))
  end do
  w_vperp(nperp) = 0.5_dp * (vperp(nperp) - vperp(nperp-1))

  ! Trapezoidal weights for v∥ (uniform grid)
  w_vpar(1)    = 0.5_dp * dvpar
  w_vpar(npar) = 0.5_dp * dvpar
  do j = 2, npar-1
    w_vpar(j) = dvpar
  end do

  ! Source vector: g(ix) = f * v⊥ * integration weights
  do i = 1, nperp
    do j = 1, npar
      g(index_mat(i,j)) = f_values(i,j) * vperp(i) * w_vperp(i) * w_vpar(j)
    end do
  end do

  ! psi_vec = prefactor * sum_phi * g
  !
  ! This was a single DGEMV against the dense sum_phi, which had to stream the
  ! whole matrix from RAM on every time step (11.9 GiB at 200x200).  The matrix
  ! is block-Toeplitz with Toeplitz blocks, so the same product is obtained from
  ! the compressed kernel by FFT, touching ~10^2 MiB instead.  Result identical
  ! to round-off; see mod_phi_kernel.
  call phi_kernel_matvec(g, psi_vec)
  psi_vec = prefactor * psi_vec

  ! Unpack to 2D
  do i = 1, nperp
    do j = 1, npar
      psi_values(i,j) = psi_vec(index_mat(i,j))
    end do
  end do

end subroutine compute_psi


SUBROUTINE regularise_axis_3(phi, d2phi_raw, vperp, nperp, npar)
  !-------------------------------------------------------------------
  ! Detects and corrects near-axis contamination in phi automatically.
  !
  ! Detection : uses d2phi_raw (DERIV_X2(phi), already computed) with
  !             a local smoothness criterion — no physical parameters.
  ! Correction: fits phi = a0 + a1*vp^2 + a2*vp^4 (npoly=3 coefficients)
  !             from phi(n_bad+1:n_bad+n_fit,:) and replaces phi(1:n_bad,:).
  !
  ! Uses LAPACK dgesv to solve the 3x3 normal equations.
  ! Requires IEEE_ARITHMETIC for IEEE_IS_FINITE.
  !-------------------------------------------------------------------
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
  INTEGER,  INTENT(IN)    :: nperp, npar
   double precision, INTENT(IN)    :: vperp(nperp)
   double precision, INTENT(IN)    :: d2phi_raw(nperp, npar)  ! for detection only
   double precision, INTENT(INOUT) :: phi(nperp, npar)

  INTEGER,  PARAMETER :: n_fit      = 8       ! points used for phi fit
  INTEGER,  PARAMETER :: npoly      = 3       ! a0 + a1*vp^2 + a2*vp^4
   double precision, PARAMETER :: thr_smooth = 0.01d0 ! 1%:  seed smoothness threshold
   double precision, PARAMETER :: thr_bad    = 0.00005d0! 0.005%: bad-point threshold

  INTEGER  :: i, j, i_seed, n_bad, info
   double precision :: r1, r2, dev, ref_i, vp2i
   double precision :: col(nperp)

  ! Detection: linear fit of d2phi in vp^2 (degree 1 sufficient for detection)
   double precision :: vp2_det(n_fit)
   double precision :: sum1, sumx, sumx2, sumy, sumxy, denom, c0_det, c1_det

  ! Correction: degree-2 fit of phi in vp^2 via normal equations + dgesv
   double precision :: A(n_fit, npoly), ATA(npoly,npoly), ATb(npoly,1)
   double precision :: ATA_copy(npoly,npoly)
   double precision :: vp2_fit(n_fit)
  INTEGER  :: ipiv(npoly)

  !--------------------------------------------------------------------
  ! Step 1: find i_seed — first index where d2phi is locally smooth
  !--------------------------------------------------------------------
  col    = d2phi_raw(:, 1)   ! use column 1 as representative
  i_seed = nperp - n_fit     ! safe fallback

  DO i = 1, nperp - 2
    IF (.NOT. (IEEE_IS_FINITE(col(i))   .AND. &
               IEEE_IS_FINITE(col(i+1)) .AND. &
               IEEE_IS_FINITE(col(i+2)))) CYCLE
    IF (col(i) == 0.d0 .OR. col(i+1) == 0.d0) CYCLE
    r1 = ABS(col(i+1)/col(i)   - 1.d0)
    r2 = ABS(col(i+2)/col(i+1) - 1.d0)
    IF (r1 < thr_smooth .AND. r2 < thr_smooth) THEN
      i_seed = i
      EXIT
    END IF
  END DO

  !--------------------------------------------------------------------
  ! Step 2: build local linear fit of d2phi vs vp^2 at i_seed
  !         (degree 1 is enough for detection purposes)
  !--------------------------------------------------------------------
  DO i = 1, n_fit
    vp2_det(i) = vperp(i_seed + i - 1)**2
  END DO
  sum1  = REAL(n_fit, 8)
  sumx  = SUM(vp2_det)
  sumx2 = SUM(vp2_det**2)
  sumy  = SUM(col(i_seed : i_seed+n_fit-1))
  sumxy = SUM(vp2_det * col(i_seed : i_seed+n_fit-1))
  denom = sum1*sumx2 - sumx**2
  c0_det = (sumy*sumx2  - sumxy*sumx) / denom
  c1_det = (sum1*sumxy  - sumx*sumy)  / denom

  !--------------------------------------------------------------------
  ! Step 3: scan backward from i_seed-1, track outermost bad index
  !--------------------------------------------------------------------
  n_bad = 0
  DO i = i_seed - 1, 1, -1
    IF (.NOT. IEEE_IS_FINITE(col(i))) THEN
      n_bad = i_seed - 1
      EXIT
    END IF
    ref_i = c0_det + c1_det * vperp(i)**2
    dev   = ABS(col(i) - ref_i) / (ABS(ref_i) + 1d-30)
    IF (dev > thr_bad) n_bad = MAX(n_bad, i)
  END DO

  !WRITE(*,'(A,I4,A,ES10.3)') &
  !  '  regularise_axis: n_bad=', n_bad, &
  !  '  vp_bad_max=', vperp(MAX(n_bad, 1))

  IF (n_bad < 1) RETURN

  !--------------------------------------------------------------------
  ! Step 4: fit phi = a0 + a1*vp^2 + a2*vp^4 from n_bad+1:n_bad+n_fit
  !         Solve normal equations ATA*x = ATb via dgesv, for all j
  !--------------------------------------------------------------------
  DO i = 1, n_fit
    vp2_fit(i)  = vperp(n_bad + i)**2
    A(i, 1)     = 1.d0
    A(i, 2)     = vp2_fit(i)
    A(i, 3)     = vp2_fit(i)**2
  END DO
  ATA = MATMUL(TRANSPOSE(A), A)   ! 3x3, same for all j

  DO j = 1, npar
    ATb(:, 1) = MATMUL(TRANSPOSE(A), phi(n_bad+1 : n_bad+n_fit, j))
    ATA_copy  = ATA                ! dgesv overwrites its A argument
    CALL dgesv(npoly, 1, ATA_copy, npoly, ipiv, ATb, npoly, info)
    IF (info /= 0) THEN
      WRITE(*,*) 'regularise_axis: dgesv failed, info=', info, ' j=', j
      CYCLE
    END IF
    DO i = 1, n_bad
      vp2i      = vperp(i)**2
      phi(i, j) = ATb(1,1) + ATb(2,1)*vp2i + ATb(3,1)*vp2i**2
    END DO
  END DO

END SUBROUTINE regularise_axis_3

    
!subroutine test_maxwell
!
!use shared_grid
!use shared_plasma
!use mod_ncint
!use mod_grid
!    use derivatives_2d
!
!
!implicit none
!
!integer ipe,ipa
!
!double precision sq2,v1,v2,arg,arg2,coef1,func1,derfarg,chandra,chandrap,dgonv
!double precision, dimension(nperp,npar) :: phi,dphidpe, d2phidpe2
!
!
!common/mathcons/pi,twopi
!
!double precision pi,twopi
!
!    ! ==============================
!    ! TEST : maxwellian background
!    ! ==============================
!    
!! We compare Dperperp with the analytical expression for a
!    !  Maxwellian background
!! Dperpperpis the factor in front of d2fdvperp2
!
!     open(41,file='D2PhiDpe2_max_0.txt',status='unknown')
!     open(42,file='Phi_max_0.txt',status='unknown')
!     open(43,file='DPhiDpe_max_0.txt',status='unknown')
!
!         ! Write header
!    write(41, '(A)')  "# Vperp D2Phi/Dvperp2(vperp,0)"
!    write(42, '(A)')  "# Vperp Phi(vperp,0)"
!    write(43, '(A)')  "# Vperp DPhi/Dvperp(vperp,0)"
!
!
!do ipa=1,npar
!        do ipe=1,nperp
!    
!!call cblin(ipe,ipa,vteff,gammaa,1.d0,c20,c02,c11,c10,c01,c00)
!
!
!
!
!    ! We compare Psi with the analytical expression for a
!    !  Maxwellian background
!!         
!            sq2 = dsqrt(2.d0)
!!
!!
!v1=dSQRT(vperp(ipe)**2+vpar(ipa)**2)
!v2=vperp(ipe)**2+vpar(ipa)**2
!    
!! Main mathematical functions
!
!arg=v1/sq2/vteff
!
!!write(*,*) 'Test 0',vperp(ipe),vpar(ipa),arg
!
!arg2=v2/2.d0/vteff**2
!
!coef1=2.d0/dsqrt(pi)
!
!derfarg = coef1*dexp(-arg2) ! Erf'[u]
!
!!if (arg < 1d-4) then
!!    func1 = coef1*arg
!!   chandra = (2*arg/3.d0-2*arg2/5.d0)/DSQRT(pi)
!!else
!
!    func1 = derf(arg) ! Erf[u]
!
!chandra = (func1-arg*derfarg)/(2*arg**2)
!
!!endif
!
!! First derivative of G wrt. its argument
!
!chandrap = derfarg-2*chandra/arg
!
!! First derivative of G/v wrt. v
!
!dGonv = (arg*chandrap-chandra)/v2
!
!
!
!!dpsidve_max(ipe,ipa)=-npart/16.d0/pi*vperp(ipe)/v1*(derfarg/arg+func1*(2.d0-1d0/arg2))
!!dpsidva_max(ipe,ipa)=-npart/16.d0/pi*vpar(ipa)/v1*(derfarg/arg+func1*(2.d0-1d0/arg2))
!
!!psi_max(ipe,ipa)=-npart/8.d0/pi*v1*(derfarg/2.d0/arg+func1*(1.d0+1d0/2.d0/arg2))
!!d3p(ipe,ipa) = -npart/4.d0/pi*vperp(ipe)/v1*dGonv
!!d3p(ipe,ipa) = -npart/4.d0/pi*vperp(ipe)/v1*(arg*derfarg-3.d0*chandra)/v2
!if (arg < 1d-4) then
!    phi(ipe,ipa) =  -npart/4.d0/pi*coef1/sq2/vteff*(1.d0-arg2/3.d0)
!else
! phi(ipe,ipa) =  -npart/4.d0/pi*func1/v1   
!endif
!
!!dphidpe(ipe,ipa) =  npart/4.d0/pi/vteff**2*chandra*vperp(ipe)/v1
!
!!
!        enddo
!enddo
!
!call deriv_x2(phi,vperp,nperp,npar,d2phidpe2)
!call deriv_x1(phi,vperp,nperp,npar,dphidpe)
!
!do ipe=1,nperp
!    
!!if(ipa == (npar+1)/2) 
!    write(41,*) vperp(ipe),d2phidpe2(ipe,(npar+1)/2)!then!(2*npar+1)/3)
!      write(42,*) vperp(ipe),phi(ipe,(npar+1)/2)  
!      write(43,*) vperp(ipe),dphidpe(ipe,(npar+1)/2)
!!endif
!!
!enddo
!!enddo
!    
!close(43)
!close(42)
!    close(41)
!
!!    enddo
!!enddo
!
!!open(61,file='d2PsiDpe2_max_z=0_3.txt',status='unknown')
!
!  !   call deriv_x1(dpsidve_max,vperp,nperp,npar,d2p)
!     
!!do ipe=1,nperp
!!      !  psi_max_0(ipe) = c20*(-npart/4.d0/pi/gammaa)
!!  !  c20_lim = 2/dsqrt(2.d0*pi)/vteff*(1.d0/3.d0-arg2/5.d0)
!!    psi_max_0(ipe) = d2p(ipe,(npar+1)/2)
!!    write(61,*) vperp(ipe),psi_max_0(ipe)!,c20_lim*(-npart/4.d0/pi)!,ipa
!!        enddo
!!
!!close(61)
!
!    
!    ! ==============================
!    ! END TEST : maxwellian background
!    ! ==============================
!    
!end subroutine test_maxwell

!=====================================================================
  ! Output results to file
  !=====================================================================
!  subroutine output_results(psi_values)!,psi_0
!  
!  use shared_grid
!
!    implicit none
!    integer :: i, j, ios
!    character(len=*), parameter :: filename = "D2PhiDvperp2_integral_6_6.txt"
! !   character(len=*), parameter :: filename0 = "d2psidpe2_results_z=0_1.txt"
!!    character(len=*), parameter :: filename_c = "d2psi_results_z=0-comp.txt"
!    
!double precision, dimension(nperp,npar), intent(in) :: psi_values
!!double precision, dimension(nperp), intent(in) :: psi_0
!
!    
!   
!    open(unit=10, file=filename, status='replace', action='write', iostat=ios)
!    if (ios /= 0) then
!      print *, "Error opening file for writing"
!      return
!    end if
!    !open(unit=11, file=filename0, status='replace', action='write', iostat=ios)
!    !if (ios /= 0) then
!    !  print *, "Error opening file for writing"
!    !  return
!    !end if
!    !    open(unit=12, file=filename_c, status='replace', action='write', iostat=ios)
!    !if (ios /= 0) then
!    !  print *, "Error opening file for writing"
!    !  return
!    !end if
!
!
!   write(*,*) 'Slice at vpar = ',vpar((npar+1)/2)
!   
!    ! Write header
!    write(10, '(A)') "# Vperp Phi(vperp,0)"
!    write(10, *)
!   
!    ! Write data
!    do i = 1, nperp
!      do j = 1, npar
!     !   write(10, '(4ES15.6)') vperp(i), vpar(j), psi_values(i,j)
!        if (j==(npar+1)/2) then!(npar+1)/2
!            write(10, '(4ES15.6)') vperp(i), psi_values(i,j)
!       !     write(12, '(4ES15.6)') vperp(i), dabs(psi_values(i,j)-psi_0(i))/psi_0(i)*100.d0
!        endif
!        
!      end do
!!      write(10, *)  ! Blank line between radial slices for gnuplot
!    end do
!    
!!    close(12)
!    close(11)
!!    close(10)
!!    print *, "Results written to ", filename
!   
!  end subroutine output_results

    end module nlterm
