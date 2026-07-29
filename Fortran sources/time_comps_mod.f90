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
  !*   Tn [keV] from the phase-space-weighted least-  *
  !*   squares slope of ln f vs v^2:                  *
  !*     slope = d(ln f)/d(v^2) ~ -1/(2 vth_n^2)      *
  !*     vth_n^2 = -1/(2 slope),  Tn = A vth_n^2/cvth^2*
  !*   For a Maxwellian ln f is exactly linear in v^2,*
  !*   so the slope -> -1/(2 vth^2) on ANY grid: this *
  !*   estimator has no 1/v^2 weight, no near-axis     *
  !*   singularity, is grid-convergent and stays > 0. *
  !*   Weighting by f*jacob lets the dense cold bulk  *
  !*   (the SC-drag population, isc=3 background)      *
  !*   dominate.  Code convention vth=9.79e3*sqrt(T/A)*
  !*   05/2026: F. Louche  (06/2026: log-slope form)  *
  !***************************************************

  SUBROUTINE time_Tn(f, dens, Tn)

    USE shared_grid
    USE shared_plasma

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN)  :: f(nperp, npar), dens   ! dens kept for interface
    DOUBLE PRECISION, INTENT(OUT) :: Tn          ! [keV]

    DOUBLE PRECISION, PARAMETER :: cvth = 9.79d3 ! sqrt(e/m_p) [m/s per sqrt(eV/amu)]
    ! core_frac (namelist, default 1e-2): fit only the bulk core, f > core_frac*max(f)
    DOUBLE PRECISION :: fmax, fcore, w, x, y, vth2, slope
    DOUBLE PRECISION :: sw, swx, swy, swxx, swxy, xbar, ybar, denom
    INTEGER :: iv, ip

    ! Restrict the ln f vs v^2 fit to the dense thermal CORE (f within
    ! core_frac of the peak).  A fit spanning the full bulk + RF tail returns
    ! a slope shallower than the bulk's, overestimating Tn (Tn > Teff) and
    ! making isc=3 overheat like isc=2.  Limiting to the core isolates the
    ! cold-bulk slope -> Tn < Teff, while staying grid-robust (the core is
    ! well resolved and free of the 1/v^2 singularity).  core_frac is the knob
    ! trading bulk-purity (smaller) against fit stability (larger).
    fmax = 0.0d0
    DO ip = 1, npar
      DO iv = 1, nperp
        IF (f(iv,ip) > fmax) fmax = f(iv,ip)
      END DO
    END DO
    fcore = core_frac * fmax

    sw = 0.0d0; swx = 0.0d0; swy = 0.0d0; swxx = 0.0d0; swxy = 0.0d0
    DO ip = 1, npar
      DO iv = 1, nperp
        IF (f(iv,ip) <= fcore) CYCLE
        w = f(iv,ip) * jacob(iv,ip)          ! phase-space (density) weight
        x = vperp(iv)**2 + vpar(ip)**2       ! v^2
        y = LOG(f(iv,ip))                    ! ln f
        sw   = sw   + w
        swx  = swx  + w*x
        swy  = swy  + w*y
        swxx = swxx + w*x*x
        swxy = swxy + w*x*y
      END DO
    END DO

    IF (sw <= 0.0d0) THEN                     ! no valid cells
      Tn = 0.0d0
      RETURN
    END IF

    xbar  = swx / sw
    ybar  = swy / sw
    denom = swxx - sw*xbar*xbar               ! Sum w (x-xbar)^2
    IF (denom <= 0.0d0) THEN                   ! degenerate (all f at one v^2)
      Tn = 0.0d0
      RETURN
    END IF

    slope = (swxy - sw*xbar*ybar) / denom      ! d(ln f)/d(v^2) ~ -1/(2 vth^2)
    IF (slope >= 0.0d0) THEN                    ! f not decreasing with v^2: no T
      Tn = 0.0d0
      RETURN
    END IF

    vth2 = -1.0d0 / (2.0d0 * slope)            ! [m^2/s^2]
    Tn   = aa * vth2 / cvth**2 / 1.0d3         ! keV  (T[eV] = A vth^2 / cvth^2)

  END SUBROUTINE time_Tn

END MODULE time_comps_mod
