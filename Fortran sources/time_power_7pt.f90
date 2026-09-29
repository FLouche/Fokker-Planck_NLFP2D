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

!-----------------------------------------------------------------------
! Restored 2026-09-07 from commit c0180d0 (removed 17d31fe).  Analytic vs
! numerical Rosenbluth potentials at v_par=0; see
! CN_N2_CrankNicolson_instability/rosenbluth_retrieved/README.md
!-----------------------------------------------------------------------

SUBROUTINE sc_components_maxw_diag(vth)

  USE shared_grid           ! vperp, vpar, jmid, nperp, npar
  USE shared_plasma         ! aa, npart, gammaa
  USE shared_timer          ! outfile()

  IMPLICIT NONE
  INTEGER, PARAMETER :: dp = KIND(1.0D0)

  REAL(dp), INTENT(IN) :: vth

  REAL(dp), DIMENSION(nperp) :: Dpepe, Dpapa, Dpepa, Fpe, Fpa, psi_p, phi_p
  REAL(dp) :: pi, sq2, sqpi, vpe, vpa, v1, v2, arg, arg2
  REAL(dp) :: func1, derfarg, chandra, Theta, Phi, Psi
  INTEGER  :: iv

  pi   = 4.0_dp * ATAN(1.0_dp)
  sq2  = SQRT(2.0_dp)
  sqpi = SQRT(pi)
  vpa  = vpar(jmid)

  DO iv = 1, nperp
    vpe = vperp(iv)
    v2  = vpe*vpe + vpa*vpa
    v1  = SQRT(v2)
    arg = v1 / (sq2*vth)
    arg2= v2 / (2.0_dp*vth*vth)
    func1   = ERF(arg)
    derfarg = (2.0_dp/sqpi) * EXP(-arg2)
    chandra = (func1 - arg*derfarg) / (2.0_dp*arg*arg)
    Theta   = chandra / v1
    Phi     = (func1 - 3.0_dp*chandra) / (2.0_dp*v1*v2)
    Psi     = chandra / (v1*vth*vth)          ! maonmb = 1 (same species)
    Dpepe(iv) = gammaa * (Theta + vpa*vpa*Phi)
    Dpapa(iv) = gammaa * (Theta + vpe*vpe*Phi)
    Dpepa(iv) = gammaa * (-vpe*vpa*Phi)
    Fpe(iv)   = gammaa * (-vpe*Psi)
    Fpa(iv)   = gammaa * (-vpa*Psi)
    psi_p(iv) = -npart/(8.0_dp*pi) * v1 * (derfarg/(2.0_dp*arg) + func1*(1.0_dp + 1.0_dp/(2.0_dp*arg2)))
    phi_p(iv) = -npart/(4.0_dp*pi*v1) * func1
  END DO

  CALL wr1d('sc_Dpepe_at_vpar0.txt', Dpepe)
  CALL wr1d('sc_Dpapa_at_vpar0.txt', Dpapa)
  CALL wr1d('sc_Dpepa_at_vpar0.txt', Dpepa)
  CALL wr1d('sc_Fpe_at_vpar0.txt',   Fpe)
  CALL wr1d('sc_Fpa_at_vpar0.txt',   Fpa)
  CALL wr1d('sc_psi_at_vpar0.txt',   psi_p)
  CALL wr1d('sc_phi_at_vpar0.txt',   phi_p)
  WRITE(*,*) '  SC components at vpar=0 (D, F, psi, phi) written.'

CONTAINS
  SUBROUTINE wr1d(name, a)
    CHARACTER(*), INTENT(IN) :: name
    REAL(dp),     INTENT(IN) :: a(nperp)
    INTEGER :: i
    OPEN(521, file=TRIM(outfile(name)), status='unknown')
    DO i = 1, nperp
      WRITE(521,*) vperp(i), a(i)
    END DO
    CLOSE(521)
  END SUBROUTINE wr1d

END SUBROUTINE sc_components_maxw_diag

!*******************************************************
!* Density rate carried by each term of the operator   *
!*******************************************************
!
! Same pattern as time_power_7pt, with weight 1 in place of ekin:
!
!     dn_k/dt = INT( L_k f  d3v ) = INT( (L_k f) * jacob )
!
! The terms are linear in the coefficients, so on interior rows they add up
! to INT(L f), i.e. the dn/dt seen in density_vs_time (up to the O(dt) gap
! between f^{n+1} evaluated here and the theta-average the scheme advances).
! Every term is a divergence in the continuum, so each one should only move
! particles through the Dirichlet walls; a term that stays nonzero with the
! walls far from the tail is not conserving particles in its discrete form.
!
! Raw rates in m^-3 s^-1: no npart/dens normalisation, so the columns compare
! directly with the time derivative of density_vs_time.
!
! 29/09/2026: added to split the density drift of the ITER NLSC grid scan.
!*******************************************************

SUBROUTINE time_density_terms_7pt(f, ncoll, nsc, nRF, nsource, nlosses)

  USE shared_grid
  USE mod_ncint
  USE shared_plasma
  USE shared_beam
  USE shared_rf
  USE shared_FPterms
  USE shared_timer          ! isc
  USE func_index
  USE mod_apply_operator

  IMPLICIT NONE

  INTEGER, PARAMETER :: dp = KIND(1.0D0)

  REAL(dp), INTENT(IN)  :: f(nbig)
  REAL(dp), INTENT(OUT) :: ncoll(nbulk), nsc, nRF, nsource, nlosses

  REAL(dp) :: fint(nperp,npar), Lf(nbig)
  REAL(dp), ALLOCATABLE :: rf00(:,:)
  REAL(dp) :: taum_save
  INTEGER  :: iv, ip, ib

  ! apply_operator reads taum through fd_stencil_2d: keep the loss term out
  ! of the operator terms and restore it for the caller (see time_power_7pt).
  taum_save = taum
  taum = 0.0_dp

  DO ib = 1, nbulk
    CALL apply_operator(gammab(ib)*colin20_sp(:,:,ib), &
                        gammab(ib)*colin02_sp(:,:,ib), &
                        gammab(ib)*colin11_sp(:,:,ib), &
                        gammab(ib)*colin10_sp(:,:,ib), &
                        gammab(ib)*colin01_sp(:,:,ib), &
                        gammab(ib)*colin00_sp(:,:,ib), &
                        f, Lf)
    CALL integrate_Lf(Lf, ncoll(ib))
  END DO

  IF (isc /= 0) THEN
    CALL apply_operator(sc20, sc02, sc11, sc10, sc01, sc00, f, Lf)
    CALL integrate_Lf(Lf, nsc)
  ELSE
    nsc = 0.0_dp
  END IF

  IF (irf == -1) THEN
    ALLOCATE(rf00(nperp,npar))
    rf00 = 0.0_dp
    CALL apply_operator(rf20, rf02, rf11, rf10, rf01, rf00, f, Lf)
    DEALLOCATE(rf00)
    CALL integrate_Lf(Lf, nRF)
  ELSE
    nRF = 0.0_dp
  END IF

  IF (isource == -1) THEN
    DO iv = 1, nperp
      DO ip = 1, npar
        fint(iv,ip) = source(iv,ip) * jacob(iv,ip)
      END DO
    END DO
    CALL ncint_2d(fint, nsource)
    DO iv = 1, nperp
      DO ip = 1, npar
        fint(iv,ip) = -f(index_mat(iv,ip)) / taus * jacob(iv,ip)
      END DO
    END DO
    CALL ncint_2d(fint, nlosses)
  ELSE
    nsource = 0.0_dp
    nlosses = 0.0_dp
  END IF

  taum = taum_save

CONTAINS

  SUBROUTINE integrate_Lf(Lf_in, res)
    REAL(dp), INTENT(IN)  :: Lf_in(nbig)
    REAL(dp), INTENT(OUT) :: res
    DO iv = 1, nperp
      DO ip = 1, npar
        fint(iv,ip) = Lf_in(index_mat(iv,ip)) * jacob(iv,ip)
      END DO
    END DO
    CALL ncint_2d(fint, res)
  END SUBROUTINE integrate_Lf

END SUBROUTINE time_density_terms_7pt
