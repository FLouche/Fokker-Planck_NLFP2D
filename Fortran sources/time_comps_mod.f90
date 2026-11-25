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
    USE mod_ncint

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN)  :: f(nperp, npar)
    DOUBLE PRECISION, INTENT(OUT) :: dens

    DOUBLE PRECISION, ALLOCATABLE :: fint(:,:)

    ! This routine used to carry its own inlined copy of ncint_2D's quadrature.
    ! The two then had to be kept in step by hand, and were not: both started as
    ! a left-endpoint rectangle rule, and a fix applied to ncint_2D alone
    ! changed nothing, because THIS copy is the one in the renormalisation path
    ! (fout = fout*npart/dens_tmp).  It now calls the shared routine, so the
    ! quadrature has exactly one definition and one order.
    ALLOCATE(fint(nperp, npar))
    fint = f * jacob                      ! jacob = 2*pi*vperp
    CALL ncint_2d(fint, dens)
    DEALLOCATE(fint)

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


  !***********************************************************
  !* VDF snapshot (namelist n_snap > 0; off by default)      *
  !*                                                         *
  !* Called by timefp_7pt and timefp_7pt_nl after each step's *
  !* solve.  Writes f every n_snap steps to vdf_snap_<step>  *
  !* .txt (through outfile, so -<casename> is appended), in   *
  !* the 3-column layout of fout.txt with the time on a       *
  !* header line.  <step> is the total step count            *
  !* nstep_restart + itime: a restart (iold=-1) continues the *
  !* numbering and the cadence of the run it restarts from,   *
  !* instead of overwriting its snapshots.  f is the raw      *
  !* solution, without the end-of-run renormalisation applied *
  !* to fout.txt in sourceless runs.  Unit 530.               *
  !***********************************************************
  SUBROUTINE write_vdf_snapshot(itime, time, f)

    USE shared_grid      ! n_snap, nstep_restart, nperp, npar, vperp, vpar
    USE shared_timer     ! outfile

    IMPLICIT NONE

    INTEGER,          INTENT(IN) :: itime            ! steps done in this run
    DOUBLE PRECISION, INTENT(IN) :: time
    DOUBLE PRECISION, INTENT(IN) :: f(nperp, npar)

    CHARACTER(len=64) :: snapname
    INTEGER :: istep, iv, ip

    IF (n_snap <= 0) RETURN
    istep = nstep_restart + itime
    IF (MOD(istep, n_snap) /= 0) RETURN

    WRITE(snapname,'(A,I6.6,A)') 'vdf_snap_', istep, '.txt'
    OPEN(530, file=TRIM(outfile(snapname)), status='unknown')
    WRITE(530,'(A,ES16.8)') '# time = ', time
    DO iv = 1, nperp
      DO ip = 1, npar
        WRITE(530,*) vperp(iv), vpar(ip), f(iv,ip)
      END DO
    END DO
    CLOSE(530)

  END SUBROUTINE write_vdf_snapshot

END MODULE time_comps_mod
