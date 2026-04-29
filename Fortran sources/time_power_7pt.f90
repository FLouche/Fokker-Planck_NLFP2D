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
    

SUBROUTINE time_power_7pt(f, dens, pcoll, pRF, psource, plosses, pcoll_self)

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
  REAL(dp), PARAMETER :: pmass = 1.6726d-27   ! proton mass [kg]

  INTEGER :: iv, ip, imu, ix, ib, row, k

  !================================================================
  ! 0.  Preliminary
  !================================================================
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

    CALL apply_operator(colin20_sp(:,:,ib), colin02_sp(:,:,ib), &
                        colin11_sp(:,:,ib), colin10_sp(:,:,ib), &
                        colin01_sp(:,:,ib), colin00_sp(:,:,ib), &
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
  IF (isc == -1) THEN

    CALL apply_operator(sc20, sc02, sc11, sc10, sc01, sc00, f, Lf)

    DO iv = 1, nperp
      DO ip = 1, npar
        ix = index_mat(iv, ip)
        fint(iv,ip) = ekin(iv,ip) * Lf(ix) * jacob(iv,ip)
      END DO
    END DO

    CALL ncint_2d(fint, pcoll_self)

  ELSE
    pcoll_self = 0.0_dp
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

    taum = 0.0_dp   ! restore

  ELSE
    plosses = 0.0_dp
    psource = 0.0_dp
  END IF


END SUBROUTINE time_power_7pt
