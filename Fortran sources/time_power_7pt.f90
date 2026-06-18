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
                           pcoll_self, pcoll_self_perp, pcoll_self_par)

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

  !--- Local arrays ------------------------------------------------
  REAL(dp) :: ekin(nperp,npar)          ! kinetic energy at each node
  REAL(dp) :: fint(nperp,npar)          ! integrand array
  REAL(dp) :: Lf(nbig)                  ! operator applied to f
  REAL(dp), ALLOCATABLE :: rf00(:,:)

  !--- Stencil workspace (reused for each operator call) -----------
  INTEGER  :: col_idx(49)
  REAL(dp) :: stencil_coeff(49)
  INTEGER  :: n_entries
  REAL(dp) :: rhs_ij

  !--- Scalars -----------------------------------------------------
  REAL(dp) :: normfac
  REAL(dp) :: taum_save          ! saved taum; restored on exit
  REAL(dp), PARAMETER :: pmass = 1.6726d-27   ! proton mass [kg]

  INTEGER :: iv, ip, imu, ix, ib, row, k

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


!*******************************************************************
! DIAGNOSTIC (SC_diagnostics branch only): velocity-space density of
! the self-collision power,
!
!     dP_SC/d3v = (1/2) m v^2 * C_SC[f]
!
! i.e. the integrand of P_SC *without* the jacobian (the d3v measure).
! Integrating dP_SC/d3v * jacob over the grid reproduces pcoll_self.
! Shows directly WHERE in velocity space the SC operator deposits
! (>0) or removes (<0) energy.  Writes a 2D map and a vpar=0 slice.
!*******************************************************************

SUBROUTINE sc_power_density_diag(f, dens)

  USE shared_grid
  USE shared_plasma
  USE shared_beam        ! isource, taum
  USE shared_FPterms     ! sc00..sc02
  USE shared_timer       ! isc, outfile()
  USE func_index
  USE mod_apply_operator

  IMPLICIT NONE

  INTEGER, PARAMETER :: dp = KIND(1.0D0)

  REAL(dp), INTENT(IN) :: f(nbig), dens

  REAL(dp) :: Lf(nbig)
  REAL(dp) :: pdens(nperp,npar)
  REAL(dp) :: normfac, taum_save
  REAL(dp), PARAMETER :: pmass = 1.6726d-27   ! proton mass [kg]
  INTEGER  :: iv, ip, ix

  IF (isc == 0) RETURN

  ! apply_operator (via fd_stencil_2d) reads taum from shared_beam;
  ! the collision operator must not include the -f/taus loss term.
  taum_save = taum
  taum = 0.0_dp

  IF (isource == 0) THEN
    normfac = npart / dens
  ELSE
    normfac = 1.0_dp
  END IF

  CALL apply_operator(sc20, sc02, sc11, sc10, sc01, sc00, f, Lf)

  DO iv = 1, nperp
    DO ip = 1, npar
      ix = index_mat(iv, ip)
      pdens(iv,ip) = 0.5_dp * pmass * aa &
                   * (vperp(iv)**2 + vpar(ip)**2) * normfac * Lf(ix)
    END DO
  END DO

  taum = taum_save

  ! 2D map
  OPEN(512, file=TRIM(outfile('sc_power_density.txt')), status='unknown')
  DO iv = 1, nperp
    DO ip = 1, npar
      WRITE(512,*) vperp(iv), vpar(ip), pdens(iv,ip)
    END DO
  END DO
  CLOSE(512)

  ! vpar = 0 slice (column jmid)
  OPEN(513, file=TRIM(outfile('sc_power_density_at_vpar0.txt')), status='unknown')
  DO iv = 1, nperp
    WRITE(513,*) vperp(iv), pdens(iv,jmid)
  END DO
  CLOSE(513)

  WRITE(*,*) '  SC power density (2D map + vpar=0 slice) written.'

END SUBROUTINE sc_power_density_diag


!*******************************************************************
! DIAGNOSTIC (SC_diagnostics branch): full SC friction/diffusion tensor
! and Rosenbluth potentials at v_par=0, for a MAXWELLIAN background of
! thermal speed vth (= sqrt(T/m)), using the cblin functions Theta, Phi, Psi.
! All D/F components carry the gammaa prefactor (the actual FP operator
! coefficients), so they are comparable with the nonlinear (isc=-1) case:
!   Dperperp = gammaa(Theta + vpar^2 Phi)   (= sc20)
!   Dparpar  = gammaa(Theta + vperp^2 Phi)  (= sc02)
!   Dperpar  = gammaa(-vperp vpar Phi)       (= sc11/2)
!   Fperp    = gammaa(-vperp Psi)
!   Fpar     = gammaa(-vpar  Psi)
! Rosenbluth potentials (Maxwellian analytic form):
!   psi = -npart/(8 pi) v [ derfarg/(2 arg) + func1 (1 + 1/(2 arg^2)) ]
!   phi = -npart/(4 pi v) func1
! Writes one v_par=0 slice file per quantity.
!*******************************************************************

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
