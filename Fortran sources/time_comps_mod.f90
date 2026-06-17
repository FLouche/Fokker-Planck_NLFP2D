! *****************************************
! *   Computes density and temperatures   *
! *         at each time step             *
! *       28/05/2026: F Louche            *
! *****************************************
    
    MODULE time_comps_mod

  IMPLICIT NONE

CONTAINS

  !***********************************************
  !* Computation of the density at each timestep *
  !***********************************************
  !
  ! Version 2: adapted for general grids
  ! 25/03/2026: F. Louche

  SUBROUTINE time_density(f, dens)

    USE shared_grid

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN)  :: f(nperp, npar)
    DOUBLE PRECISION, INTENT(OUT) :: dens

    COMMON/mathcons/ pi, twopi
    DOUBLE PRECISION :: dvp_local, pi, twopi
    INTEGER :: iv, imu

    dens = 0.d0
    DO iv = 1, nperp-1
      dvp_local = vperp(iv+1) - vperp(iv)
      DO imu = 2, npar-1
        dens = dens + twopi * f(iv,imu) * vperp(iv) &
                            * dvp_local * dvpar
      END DO
    END DO

  END SUBROUTINE time_density

  !***************************************************
  !* Computation of the temperature at each timestep *
  !*                                                 *
  !*  Version 2.0: adapted to general spacing grids  *
  !*   30/03/2026; F.Louche                          *
  !*   26/05/2026 (F. Louche):                       *
  !*      - new definition of parallel temperature   *
  !*      - effective temperature computed           *
  !*                                                 *
  !*   28/05/2026 (F.Louche)                         *
  !*     Version 2.2: adding optional arguments      *
  !*                                                 *
  !***************************************************

  SUBROUTINE time_energy(f, dens, temp, tperp, tpar, teff)

    USE shared_grid
    USE mod_ncint
    USE shared_plasma

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN)  :: f(nperp, npar), dens
    DOUBLE PRECISION, INTENT(OUT), optional :: temp, tpar, tperp, teff

    DOUBLE PRECISION :: pmass, kev_in_J, mod
    DATA pmass    / 1.6726d-27  /   ! proton mass in kg
    DATA kev_in_J / 1.60218d-16 /   ! convert keV to Joule

    DOUBLE PRECISION :: tperp_loc, tpar_loc
    DOUBLE PRECISION, ALLOCATABLE :: fint(:,:)
    INTEGER :: iv, ip

    ALLOCATE(fint(nperp, npar))

    ! Total kinetic energy (optional)
    IF (PRESENT(temp)) THEN
      DO iv = 1, nperp
        DO ip = 1, npar
          fint(iv,ip) = f(iv,ip) * (vperp(iv)**2 + vpar(ip)**2) &
                                  * jacob(iv,ip) / dens
        END DO
      END DO
      CALL ncint_2d(fint, mod)
      temp = 0.5d0 * pmass * aa * mod / kev_in_J
    END IF

    ! Perpendicular energy (always computed: needed for teff)
    DO iv = 1, nperp
      DO ip = 1, npar
        fint(iv,ip) = f(iv,ip) * vperp(iv)**2 * jacob(iv,ip) / dens
      END DO
    END DO
    CALL ncint_2d(fint, mod)
    tperp_loc = 0.5d0 * pmass * aa * mod / kev_in_J
    IF (PRESENT(tperp)) tperp = tperp_loc

    ! Parallel energy (always computed: needed for teff)
    DO iv = 1, nperp
      DO ip = 1, npar
        fint(iv,ip) = f(iv,ip) * vpar(ip)**2 * jacob(iv,ip) / dens
      END DO
    END DO
    CALL ncint_2d(fint, mod)
    tpar_loc = pmass * aa * mod / kev_in_J
    IF (PRESENT(tpar)) tpar = tpar_loc

    IF (PRESENT(teff)) teff = (2.d0 * tperp_loc + tpar_loc) / 3.d0

    DEALLOCATE(fint)

  END SUBROUTINE time_energy

  !***************************************************
  !* Density-characteristic (log-slope) temperature  *
  !*   Tn [keV] from the harmonic mean of v^2:        *
  !*     1/vth_n^2 = <1/v^2>,   Tn = m / <1/v^2>      *
  !*   Equivalent to the f-weighted local log-slope   *
  !*   temperature; dominated by the dense cold bulk  *
  !*   (cf. Tn_density_temperature note).             *
  !*   Code convention vth = 9.79e3*sqrt(T[eV]/A):    *
  !*     Tn[eV] = A / ( (9.79e3)^2 * <1/v^2> ).        *
  !*   05/2026: F. Louche                             *
  !***************************************************

  SUBROUTINE time_Tn(f, dens, Tn)

    USE shared_grid
    USE mod_ncint
    USE shared_plasma

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN)  :: f(nperp, npar), dens
    DOUBLE PRECISION, INTENT(OUT) :: Tn          ! [keV]

    DOUBLE PRECISION, PARAMETER :: cvth = 9.79d3 ! sqrt(e/m_p) [m/s per sqrt(eV/amu)]
    DOUBLE PRECISION, ALLOCATABLE :: fint(:,:)
    DOUBLE PRECISION :: v2, v2_floor, inv_v2_avg, mom
    INTEGER :: iv, ip

    ALLOCATE(fint(nperp, npar))

    ! f/v^2 integrand (jacob included).  The 1/v^2 weight is integrable;
    ! guard only the exact origin.
    v2_floor = MAX(vperp(1)**2, 1.0d0)
    DO iv = 1, nperp
      DO ip = 1, npar
        v2 = vperp(iv)**2 + vpar(ip)**2
        fint(iv,ip) = f(iv,ip) / MAX(v2, v2_floor) * jacob(iv,ip)
      END DO
    END DO
    CALL ncint_2d(fint, mom)          ! = n * <1/v^2>
    inv_v2_avg = mom / dens           ! <1/v^2>

    Tn = aa / (cvth**2 * inv_v2_avg) / 1.0d3   ! keV

    DEALLOCATE(fint)

  END SUBROUTINE time_Tn

END MODULE time_comps_mod
