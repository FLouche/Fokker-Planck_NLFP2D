!*******************************************************
!* Computation of the power densities at each timestep *
!*   Version using 7-point Fornberg stencil            *
!*   (fd_stencil_2d replaces build_ss + dgemv)         *
!*******************************************************
!
! For each power term, the pattern is:
!   1. Apply the relevant FP operator L to f using the
!      7-point stencil:   (Lf)(row) = sum_k coeff(k)*f(col_idx(k))
!      This replaces:  build_ss -> bigm,  dgemv -> bigm*f
!   2. Integrate  ekin * (Lf) * jacob  over velocity space
!
! The dense bigm(nbig,nbig) matrix and dgemv are eliminated entirely.
! Memory saving: O(nbig^2) -> O(49*nbig).
!
!*******************************************************
    

SUBROUTINE time_power_7pt(f, dens, pcoll, pRF, psource, plosses, &
                           pcoll_self, pcoll_self_perp, pcoll_self_par, &
                           tau_rf, tau_ii, tau_ie)

  USE shared_grid
  USE mod_ncint
  USE shared_plasma
  USE shared_beam
  USE shared_rf
  USE shared_FPterms
  USE mod_fd_stencil_2d     ! provides fd_stencil_2d
  USE shared_timer
  USE func_index

  use mod_apply_operator
  
  IMPLICIT NONE

  INTEGER, PARAMETER :: dp = KIND(1.0D0)

  !--- Arguments ---------------------------------------------------
  REAL(dp), INTENT(IN)  :: f(nbig), dens
  REAL(dp), INTENT(OUT) :: pcoll(nbulk), pRF, psource, plosses, pcoll_self
  REAL(dp), INTENT(OUT) :: pcoll_self_perp, pcoll_self_par
  REAL(dp), INTENT(OUT) :: tau_rf     ! RF tail formation time [s]
  REAL(dp), INTENT(OUT) :: tau_ii     ! effective ion-ion time [s]
  REAL(dp), INTENT(OUT) :: tau_ie     ! effective ion-electron time [s]

  !--- Local arrays ------------------------------------------------
  REAL(dp) :: ekin(nperp,npar)          ! kinetic energy at each node
  REAL(dp) :: fint(nperp,npar)          ! integrand array
  REAL(dp) :: Lf(nbig)                  ! operator applied to f
  REAL(dp), ALLOCATABLE :: rf00(:,:)

  !--- Scalars -----------------------------------------------------
  REAL(dp) :: normfac
  REAL(dp) :: ekin_dens          ! npart*Teff [J/m^3]
  REAL(dp) :: pcoll_i            ! summed background-ion collisional power
  REAL(dp) :: taum_save          ! saved taum; restored on exit
  REAL(dp), PARAMETER :: pmass = 1.6726d-27   ! proton mass [kg]

  INTEGER :: iv, ip, imu, ix, ib

  !================================================================
  ! 0.  Preliminary
  !================================================================
  ! fd_stencil_2d (via apply_operator) reads taum from shared_beam.
  ! Save it now so we restore the correct value on exit; otherwise the
  ! NL time-loop, which rebuilds aa_L every step via fd_stencil_2d
  ! *after* calling this routine, would see taum=0 and lose the -f/taus
  ! loss term from the operator matrix.
  taum_save = taum
  taum = 0.0_dp    ! collisional operators do not include the loss term

  IF (isource == 0) THEN
    normfac = npart / dens
  ELSE
    normfac = 1.0_dp
  END IF

  ! Kinetic energy at each grid node [J]
  DO iv = 1, nperp
    DO ip = 1, npar
      ekin(iv,ip) = 0.5_dp * pmass * aa &
                  * (vperp(iv)**2 + vpar(ip)**2) * normfac
    END DO
  END DO

  !================================================================
  ! 1.  Collisional power density  (one term per bulk species)
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
        fint(iv,ip) = ekin(iv,ip) * Lf(ix) * jacob(iv,ip)
      END DO
    END DO

    CALL ncint_2d(fint, pcoll(ib))

  END DO

  !================================================================
  ! 2.  Self-collision power density
  !================================================================
  IF (isc /= 0) THEN

    CALL apply_operator(sc20, sc02, sc11, sc10, sc01, sc00, f, Lf)

    DO iv = 1, nperp
      DO ip = 1, npar
        ix = index_mat(iv, ip)
        fint(iv,ip) = ekin(iv,ip) * Lf(ix) * jacob(iv,ip)
      END DO
    END DO

    CALL ncint_2d(fint, pcoll_self)

    ! Perpendicular component: weight by ½ m v⊥²
    DO iv = 1, nperp
      DO ip = 1, npar
        ix = index_mat(iv, ip)
        fint(iv,ip) = 0.5_dp*pmass*aa * vperp(iv)**2 * normfac * Lf(ix) * jacob(iv,ip)
      END DO
    END DO
    CALL ncint_2d(fint, pcoll_self_perp)

    ! Parallel component: weight by ½ m v∥²
    DO iv = 1, nperp
      DO ip = 1, npar
        ix = index_mat(iv, ip)
        fint(iv,ip) = 0.5_dp*pmass*aa * vpar(ip)**2  * normfac * Lf(ix) * jacob(iv,ip)
      END DO
    END DO
    CALL ncint_2d(fint, pcoll_self_par)

  ELSE
    pcoll_self      = 0.0_dp
    pcoll_self_perp = 0.0_dp
    pcoll_self_par  = 0.0_dp
  END IF

  !================================================================
  ! 3.  RF power density
  !================================================================
  IF (irf == -1) THEN

    ALLOCATE(rf00(nperp,npar))
    rf00 = 0.0_dp

    CALL apply_operator(rf20, rf02, rf11, rf10, rf01, rf00, f, Lf)

    DEALLOCATE(rf00)

    DO iv = 1, nperp
      DO imu = 1, npar
        ix = index_mat(iv, imu)
        fint(iv,imu) = ekin(iv,imu) * Lf(ix) * jacob(iv,imu)
      END DO
    END DO

    CALL ncint_2d(fint, pRF)

  ELSE
    pRF = 0.0_dp
  END IF

  !----------------------------------------------------------------
  ! 3b. Effective timescales:  tau = npart * Teff / P
  !----------------------------------------------------------------
  ! The energy stored in the distribution divided by the rate at which a given
  ! channel supplies or removes it.  Integrating ekin*f*jacob gives
  ! npart*<Ekin>, because ekin already carries the normfac = npart/dens factor,
  ! and Teff = (2/3)*<Ekin> is the same quantity written to Teff_vs_time (see
  ! time_energy in time_comps_mod.f90), so
  !
  !     npart * Teff[J]  =  (2/3) * INT( ekin * f * jacob )
  !
  ! and no separate temperature evaluation is needed here.  Numerator and
  ! denominator are both normalised to npart, so the ratios do not depend on
  ! the running density.
  !
  ! These are *effective* times built from the power actually flowing in each
  ! channel at the current Teff, not textbook collision times at a nominal
  ! temperature (those are printed once at startup by consts).  They therefore
  ! follow both the heating of the tail and its change of shape.  At steady
  ! state the RF input balances the collisional losses, so
  !
  !     1/tau_rf  =  1/tau_ii + 1/tau_ie  (+ self-collisions, which carry no
  !                                        net energy, + any source/losses)
  !
  ! which is a useful check on the three traces.
  DO iv = 1, nperp
    DO imu = 1, npar
      ix = index_mat(iv, imu)
      fint(iv,imu) = ekin(iv,imu) * f(ix) * jacob(iv,imu)
    END DO
  END DO
  CALL ncint_2d(fint, ekin_dens)              ! npart*<Ekin>  [J/m^3]
  ekin_dens = (2.0_dp / 3.0_dp) * ekin_dens   ! = npart*Teff  [J/m^3]

  !  -> RF tail formation time
  IF (pRF > 0.0_dp) THEN
    tau_rf = ekin_dens / pRF
  ELSE
    tau_rf = 0.0_dp
  END IF

  !  -> Effective ion-ion time: loss to the background ions.  pcoll is negative
  !     (energy leaving the resonant species), hence the ABS.
  pcoll_i = 0.0_dp
  DO ib = 2, nbulk
    pcoll_i = pcoll_i + pcoll(ib)
  END DO
  IF (ABS(pcoll_i) > 0.0_dp) THEN
    tau_ii = ekin_dens / ABS(pcoll_i)
  ELSE
    tau_ii = 0.0_dp
  END IF

  !  -> Effective ion-electron time: loss to the electrons
  IF (ABS(pcoll(1)) > 0.0_dp) THEN
    tau_ie = ekin_dens / ABS(pcoll(1))
  ELSE
    tau_ie = 0.0_dp
  END IF

  !================================================================
  ! 4.  Beam source power + losses
  !================================================================
  IF (isource == -1) THEN

    taum = 1.0_dp / taus

    DO iv = 1, nperp
      DO imu = 1, npar
        ix = index_mat(iv, imu)
        ! losses: -ekin * f * (1/taus)
        fint(iv,imu) = -ekin(iv,imu) * f(ix) * taum * jacob(iv,imu)
      END DO
    END DO
    CALL ncint_2d(fint, plosses)

    DO iv = 1, nperp
      DO imu = 1, npar
        ! source power: ekin * S
        fint(iv,imu) = ekin(iv,imu) * source(iv,imu) * jacob(iv,imu)
      END DO
    END DO
    CALL ncint_2d(fint, psource)

  ELSE
    plosses = 0.0_dp
    psource = 0.0_dp
  END IF

  taum = taum_save   ! restore original value for the caller

END SUBROUTINE time_power_7pt
