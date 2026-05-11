!***********************************************************************
!  mod_ss_check — rolling-window steady-state convergence monitor
!
!  Call ss_check once per time step after diagnostics are computed.
!  Returns converged=.TRUE. when the relative change of three quantities
!  over the last n_ss_window steps has fallen below ss_tol.
!
!  Monitored quantities
!  --------------------
!    tk       : total kinetic energy (any consistent unit, > 0)
!    dens     : particle density (> 0)
!    p_net    : net power deposition = SUM(pcoll) + pRF + pcoll_self + ...
!    p_drive  : driving-power scale: max of |individual power terms|,
!               floored at 1.0 so the criterion is always meaningful
!
!  Convergence criteria (all three must be satisfied simultaneously)
!  ----------------------------------------------------------------
!    |tk(now)   - tk(now-W)|   / tk_scale   < ss_tol
!    |dens(now) - dens(now-W)| / dens_scale < ss_tol
!    |p_net(now)- p_net(now-W)|/ p_drive    < ss_tol
!
!  where the scale for energy and density is the arithmetic mean of the
!  current and the window-old value; p_drive is passed by the caller.
!
!  Parameters (from shared_timer, set via namelist)
!  -------------------------------------------------
!    i_ss_check   : 0 = disabled; -1 = enabled (consistent with irf, ifd7 convention)
!    n_ss_window  : window width in time steps (default 50)
!    ss_tol       : relative tolerance (default 1e-3)
!
!  Notes
!  -----
!  For sourceless runs (isource=0), density drifts slowly due to
!  finite-difference truncation; ss_tol may need to be relaxed slightly
!  or n_ss_window widened so the drift rate is captured accurately.
!***********************************************************************

module mod_ss_check

  use shared_timer, only: n_ss_window, ss_tol

  implicit none
  private
  public :: ss_check

  integer,  parameter :: dp = kind(1.0d0)

  ! Module-level circular buffers (allocated/reallocated on first call
  ! or whenever n_ss_window changes between runs).
  integer               :: w_size = 0
  real(dp), allocatable :: hist_ek(:), hist_n(:), hist_pnet(:)

contains

  !--------------------------------------------------------------------
  ! ss_check — call once per time step, after all diagnostics are done.
  !
  !   itime      : step index starting at 1
  !   tk         : total kinetic energy
  !   dens       : particle density
  !   p_net      : net power (sum of all signed power terms)
  !   p_drive    : positive power scale used to normalise the p_net
  !                criterion; call site should pass
  !                max(|pRF|, |pcoll_max|, |pcoll_self|, |psource|, 1.0)
  !   converged  : .TRUE. when all criteria are satisfied
  !--------------------------------------------------------------------
  subroutine ss_check(itime, tk, dens, p_net, p_drive, converged)

    integer,  intent(in)  :: itime
    real(dp), intent(in)  :: tk, dens, p_net, p_drive
    logical,  intent(out) :: converged

    integer  :: ptr
    real(dp) :: ek_old, n_old, pnet_old
    real(dp) :: rel_ek, rel_n, rel_pnet, scale

    ! (Re)allocate if this is the first call or window size changed.
    if (n_ss_window /= w_size) then
      if (allocated(hist_ek)) deallocate(hist_ek, hist_n, hist_pnet)
      allocate(hist_ek(n_ss_window), hist_n(n_ss_window), hist_pnet(n_ss_window))
      hist_ek   = 0.0_dp
      hist_n    = 0.0_dp
      hist_pnet = 0.0_dp
      w_size    = n_ss_window
      write(*,'(A,I5,A,ES8.1)') &
        '  [SS check] window=', w_size, '  tol=', ss_tol
    end if

    converged = .false.

    ! Circular buffer slot for this step (1-based).
    ! The same slot holds the value from exactly w_size steps ago —
    ! read it before overwriting.
    ptr = mod(itime - 1, w_size) + 1

    ek_old   = hist_ek(ptr)
    n_old    = hist_n(ptr)
    pnet_old = hist_pnet(ptr)

    hist_ek(ptr)   = tk
    hist_n(ptr)    = dens
    hist_pnet(ptr) = p_net

    ! Do not evaluate until the buffer has been populated at least once.
    if (itime <= w_size) return

    ! --- Relative change in kinetic energy over the window -----------
    scale  = 0.5_dp * (abs(tk) + abs(ek_old))
    if (scale > 0.0_dp) then
      rel_ek = abs(tk - ek_old) / scale
    else
      rel_ek = 0.0_dp
    end if

    ! --- Relative change in density over the window ------------------
    scale  = 0.5_dp * (abs(dens) + abs(n_old))
    if (scale > 0.0_dp) then
      rel_n = abs(dens - n_old) / scale
    else
      rel_n = 0.0_dp
    end if

    ! --- Relative change in net power --------------------------------
    ! Normalised to the caller-supplied driving-power scale.
    rel_pnet = abs(p_net - pnet_old) / max(p_drive, 1.0_dp)

    ! --- Periodic status print (every w_size steps) ------------------
    if (mod(itime, w_size) == 0) then
      write(*,'(A,I7,3(2X,A,ES9.2),2X,A,ES9.2,A)') &
        '  [SS]', itime, &
        'dE/E=',  rel_ek, &
        'dn/n=',  rel_n, &
        'dP/Pd=', rel_pnet, &
        '[tol=',  ss_tol, ']'
    end if

    converged = (rel_ek < ss_tol) .and. (rel_n < ss_tol) .and. (rel_pnet < ss_tol)

    if (converged) then
      write(*,'(A,I7,3(2X,A,ES9.2))') &
        '  [SS] CONVERGED at step', itime, &
        'dE/E=',  rel_ek, &
        'dn/n=',  rel_n, &
        'dP/Pd=', rel_pnet
    end if

  end subroutine ss_check

end module mod_ss_check
