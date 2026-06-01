!***************************************************************
!* Computation of the perpendicular and parallel momentum      *
!*   transfer rate densities (N/m**3) at each timestep        *
!*   Version using 7-point Fornberg stencil                    *
!*                                                             *
!* For each operator the pattern mirrors time_power_7pt:       *
!*   1. Apply operator to f → Lf                               *
!*   2. Integrate  ma*v_component * Lf * jacob  over v-space   *
!*                                                             *
!*   F. Louche – June 2026                                     *
!***************************************************************

SUBROUTINE time_momentum_7pt(f, dens, &
    mcoll_perp, mcoll_par, &
    mRF_perp,   mRF_par,   &
    msrc_perp,  msrc_par,  &
    mloss_perp, mloss_par, &
    mSC_perp,   mSC_par)

  USE shared_grid
  USE mod_ncint
  USE shared_plasma
  USE shared_beam
  USE shared_rf
  USE shared_FPterms
  USE shared_timer
  USE func_index
  USE mod_apply_operator

  IMPLICIT NONE

  INTEGER, PARAMETER :: dp = KIND(1.0D0)

  !--- Arguments ---------------------------------------------------
  REAL(dp), INTENT(IN)  :: f(nbig), dens
  REAL(dp), INTENT(OUT) :: mcoll_perp(nbulk), mcoll_par(nbulk)
  REAL(dp), INTENT(OUT) :: mRF_perp,   mRF_par
  REAL(dp), INTENT(OUT) :: msrc_perp,  msrc_par
  REAL(dp), INTENT(OUT) :: mloss_perp, mloss_par
  REAL(dp), INTENT(OUT) :: mSC_perp,   mSC_par

  !--- Local arrays ------------------------------------------------
  REAL(dp) :: fint_perp(nperp,npar), fint_par(nperp,npar)
  REAL(dp) :: Lf(nbig)
  REAL(dp), ALLOCATABLE :: rf00(:,:)

  !--- Scalars -----------------------------------------------------
  REAL(dp) :: normfac, mom_fac
  REAL(dp) :: taum_save
  REAL(dp), PARAMETER :: pmass = 1.6726d-27   ! proton mass [kg]

  INTEGER :: iv, ip, ix, ib

  !================================================================
  ! 0.  Preliminary
  !================================================================
  taum_save = taum
  taum      = 0.0_dp   ! collisional operators do not include the loss term

  IF (isource == 0) THEN
    normfac = npart / dens
  ELSE
    normfac = 1.0_dp
  END IF

  mom_fac = pmass * aa * normfac   ! [kg] — mass of the beam ion (with normfac)

  !================================================================
  ! 1.  Bulk-species collisional momentum transfer
  !================================================================
  DO ib = 1, nbulk

    CALL apply_operator(gammab(ib)*colin20_sp(:,:,ib), &
                        gammab(ib)*colin02_sp(:,:,ib), &
                        gammab(ib)*colin11_sp(:,:,ib), &
                        gammab(ib)*colin10_sp(:,:,ib), &
                        gammab(ib)*colin01_sp(:,:,ib), &
                        gammab(ib)*colin00_sp(:,:,ib), &
                        f, Lf)

    DO iv = 1, nperp
      DO ip = 1, npar
        ix = index_mat(iv, ip)
        fint_perp(iv,ip) = mom_fac * vperp(iv) * Lf(ix) * jacob(iv,ip)
        fint_par (iv,ip) = mom_fac * vpar (ip) * Lf(ix) * jacob(iv,ip)
      END DO
    END DO

    CALL ncint_2d(fint_perp, mcoll_perp(ib))
    CALL ncint_2d(fint_par,  mcoll_par (ib))

  END DO

  !================================================================
  ! 2.  Self-collision momentum transfer
  !================================================================
  IF (isc /= 0) THEN

    CALL apply_operator(sc20, sc02, sc11, sc10, sc01, sc00, f, Lf)

    DO iv = 1, nperp
      DO ip = 1, npar
        ix = index_mat(iv, ip)
        fint_perp(iv,ip) = mom_fac * vperp(iv) * Lf(ix) * jacob(iv,ip)
        fint_par (iv,ip) = mom_fac * vpar (ip) * Lf(ix) * jacob(iv,ip)
      END DO
    END DO

    CALL ncint_2d(fint_perp, mSC_perp)
    CALL ncint_2d(fint_par,  mSC_par)

  ELSE
    mSC_perp = 0.0_dp
    mSC_par  = 0.0_dp
  END IF

  !================================================================
  ! 3.  RF momentum transfer
  !================================================================
  IF (irf == -1) THEN

    ALLOCATE(rf00(nperp,npar))
    rf00 = 0.0_dp
    CALL apply_operator(rf20, rf02, rf11, rf10, rf01, rf00, f, Lf)
    DEALLOCATE(rf00)

    DO iv = 1, nperp
      DO ip = 1, npar
        ix = index_mat(iv, ip)
        fint_perp(iv,ip) = mom_fac * vperp(iv) * Lf(ix) * jacob(iv,ip)
        fint_par (iv,ip) = mom_fac * vpar (ip) * Lf(ix) * jacob(iv,ip)
      END DO
    END DO

    CALL ncint_2d(fint_perp, mRF_perp)
    CALL ncint_2d(fint_par,  mRF_par)

  ELSE
    mRF_perp = 0.0_dp
    mRF_par  = 0.0_dp
  END IF

  !================================================================
  ! 4.  Beam source and particle losses
  !================================================================
  IF (isource == -1) THEN

    taum = 1.0_dp / taus

    DO iv = 1, nperp
      DO ip = 1, npar
        ix = index_mat(iv, ip)
        fint_perp(iv,ip) = mom_fac * vperp(iv) * (-f(ix)*taum) * jacob(iv,ip)
        fint_par (iv,ip) = mom_fac * vpar (ip) * (-f(ix)*taum) * jacob(iv,ip)
      END DO
    END DO
    CALL ncint_2d(fint_perp, mloss_perp)
    CALL ncint_2d(fint_par,  mloss_par)

    DO iv = 1, nperp
      DO ip = 1, npar
        fint_perp(iv,ip) = mom_fac * vperp(iv) * source(iv,ip) * jacob(iv,ip)
        fint_par (iv,ip) = mom_fac * vpar (ip) * source(iv,ip) * jacob(iv,ip)
      END DO
    END DO
    CALL ncint_2d(fint_perp, msrc_perp)
    CALL ncint_2d(fint_par,  msrc_par)

  ELSE
    mloss_perp = 0.0_dp;  mloss_par = 0.0_dp
    msrc_perp  = 0.0_dp;  msrc_par  = 0.0_dp
  END IF

  taum = taum_save   ! restore for caller

END SUBROUTINE time_momentum_7pt
