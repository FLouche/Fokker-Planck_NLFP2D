module coulomb_log_mod
  !-----------------------------------------------------------------
  ! Module: coulomb_log_mod
  !
  ! Provides two subroutines:
  !   coulomb_log_ae  : Coulomb log for ion a with thermal electrons
  !   coulomb_log_ab  : Coulomb log for ion a with thermal ion b
  !
  ! Reference: NRL Plasma Formulary (2019), p. 34
  !
  ! Units convention:
  !   ne, na, nb  [m^-3]    number densities (converted to cm^-3 internally)
  !   Ta, Tb, Te  [eV]      temperatures
  !   ma, mb      [amu]     masses
  !   Za, Zb      [-]       charge numbers (integers)
  !-----------------------------------------------------------------
  implicit none
  private
  public :: coulomb_log_ae, coulomb_log_ab

contains

  !================================================================
  subroutine coulomb_log_ae(ne, Za, Te, lnae)
  !----------------------------------------------------------------
  ! Coulomb logarithm for ion species a with thermal electrons
  !
  ! NRL formula (electron-ion):
  !   Te < 10*Za^2 eV  (classical, b_min = distance of closest approach):
  !     lnae = 23 - ln( sqrt(ne) * Za * Te^{-3/2} )
  !
  !   Te > 10*Za^2 eV  (quantum, b_min = de Broglie wavelength):
  !     lnae = 24 - ln( sqrt(ne) / Te )
  !
  ! Inputs:
  !   ne   [m^-3]   electron density
  !   Za   [-]      charge number of ion a
  !   Te   [eV]     electron temperature
  !
  ! Output:
  !   lnae [-]      Coulomb logarithm
  !----------------------------------------------------------------
    double precision, intent(in)  :: ne
    double precision, intent(in)  :: Za
    double precision, intent(in)  :: Te
    double precision, intent(out) :: lnae

    double precision, parameter :: m3_to_cm3 = 1.0d-6   ! n[cm^-3] = n[m^-3] * 1e-6
    double precision :: ne_cgs, Za_r

    ! Sanity checks
    if (ne <= 0.0d0) then
      write(*,*) "ERROR coulomb_log_ae: ne must be > 0"
      stop
    end if
    if (Te <= 0.0d0) then
      write(*,*) "ERROR coulomb_log_ae: Te must be > 0"
      stop
    end if
    if (Za <= 0) then
      write(*,*) "ERROR coulomb_log_ae: Za must be > 0"
      stop
    end if

    Za_r   = real(Za, 8)
    ne_cgs = ne * m3_to_cm3   ! convert to cm^-3 for NRL formula

    if (Te < 10.0d0 * Za_r**2) then
      ! Classical regime
      lnae = 23.0d0 - log( sqrt(ne_cgs) * Za_r * Te**(-1.5d0) )
    else
      ! Quantum (de Broglie) regime
      lnae = 24.0d0 - log( sqrt(ne_cgs) / Te )
    end if

  end subroutine coulomb_log_ae


  !================================================================
  subroutine coulomb_log_ab(Za, ma, Ta, na, Zb, mb, Tb, nb, lnab)
  !----------------------------------------------------------------
  ! Coulomb logarithm for ion species a with thermal ion species b
  !
  ! NRL formula (unlike-ion collisions):
  !   lnab = 23 - ln[ Za^2*Zb^2*(ma+mb) / (ma*Tb + mb*Ta)
  !                   * sqrt( na*Za^2/Ta + nb*Zb^2/Tb ) ]
  !
  ! Inputs:
  !   Za, Zb   [-]      charge numbers
  !   ma, mb   [amu]    masses
  !   Ta, Tb   [eV]     temperatures
  !   na, nb   [m^-3]   densities
  !
  ! Output:
  !   lnab [-]   Coulomb logarithm
  !----------------------------------------------------------------
    double precision, intent(in)  :: Za, Zb
    double precision, intent(in)  :: ma, mb
    double precision, intent(in)  :: Ta, Tb
    double precision, intent(in)  :: na, nb
    double precision, intent(out) :: lnab

    double precision, parameter :: m3_to_cm3 = 1.0d-6   ! n[cm^-3] = n[m^-3] * 1e-6
    double precision :: na_cgs, nb_cgs
    double precision :: Za_r, Zb_r, Za2, Zb2
    double precision :: mass_factor, inner_sum

    ! Sanity checks
    if (Ta <= 0.0d0 .or. Tb <= 0.0d0) then
      write(*,*) "ERROR coulomb_log_ab: temperatures must be > 0"
      stop
    end if
    if (na <= 0.0d0 .or. nb <= 0.0d0) then
      write(*,*) "ERROR coulomb_log_ab: densities must be > 0"
      stop
    end if
    if (Za <= 0 .or. Zb <= 0) then
      write(*,*) "ERROR coulomb_log_ab: charge numbers must be > 0"
      stop
    end if
    if (ma <= 0.0d0 .or. mb <= 0.0d0) then
      write(*,*) "ERROR coulomb_log_ab: masses must be > 0"
      stop
    end if

    Za_r   = real(Za, 8)
    Zb_r   = real(Zb, 8)
    Za2    = Za_r**2
    Zb2    = Zb_r**2
    na_cgs = na * m3_to_cm3   ! convert to cm^-3 for NRL formula
    nb_cgs = nb * m3_to_cm3

    mass_factor = Za2 * Zb2 * (ma + mb) / (ma * Tb + mb * Ta)
    inner_sum   = sqrt( na_cgs * Za2 / Ta  +  nb_cgs * Zb2 / Tb )

    lnab = 23.0d0 - log( mass_factor * inner_sum )

  end subroutine coulomb_log_ab

end module coulomb_log_mod