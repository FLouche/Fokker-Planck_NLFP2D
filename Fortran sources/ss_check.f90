!***********************************************************************
!  mod_ss_check — rolling-window steady-state convergence monitor
!
!  Call ss_check once per time step after diagnostics are computed.
!  Returns converged=.TRUE. when the relative *rate of change* of four
!  quantities over the last n_ss_window steps has fallen below ss_tol.
!
!  Rate-based criterion (dt-independent)
!  -------------------------------------
!  Each relative change is divided by the physical time the window spans
!  (time_now - time_(now-W)), so the criterion measures a per-second
!  relative rate.  This makes ss_tol independent of the time step: a run
!  with dt=1 ms and one with dt=200 ms that are equally close to steady
!  state report the same rate.  The window duration is taken from the
!  actual simulation time, so it is also correct across phase boundaries
!  where dt changes.
!
!  Monitored quantities
!  --------------------
!    tk       : total kinetic energy (any consistent unit, > 0)
!    tkperp   : perpendicular kinetic energy
!               RF drive acts primarily in the perpendicular direction;
!               tkperp converges last and is therefore the most
!               discriminating energy convergence indicator.
!    pRF      : RF power absorbed by the plasma
!               Sensitive to changes in the resonant part of the VDF;
!               replaces the density criterion which drifts numerically
!               in sourceless runs.
!    p_net    : net power deposition = SUM(pcoll) + pRF + pcoll_self + ...
!    p_drive  : driving-power scale: max of |individual power terms|,
!               floored at 1.0 so all power criteria are always finite
!
!  Convergence criteria (all four must be satisfied simultaneously)
!  ----------------------------------------------------------------
!    |tk(now)     - tk(now-W)|     / tk_scale     / dt_win < ss_tol
!    |tkperp(now) - tkperp(now-W)| / tkperp_scale / dt_win < ss_tol
!    |pRF(now)    - pRF(now-W)|    / p_drive      / dt_win < ss_tol
!    |p_net(now)  - p_net(now-W)|  / p_drive      / dt_win < ss_tol
!
!  where the scale for energies is the arithmetic mean of the current
!  and the window-old value; p_drive is passed by the caller; and
!  dt_win = time(now) - time(now-W) is the window's physical duration.
!  When irf/=-1 (no RF), pRF=0 at every step so the pRF criterion is
!  trivially satisfied and the other three drive convergence.
!
!  Parameters (from shared_timer, set via namelist)
!  -------------------------------------------------
!    i_ss_check   : 0 = disabled; -1 = enabled (consistent with irf convention)
!    n_ss_window  : window width in time steps (default 50)
!    ss_tol       : relative tolerance PER SECOND (default 1e-3 /s)
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
  real(dp), allocatable :: hist_ek(:), hist_ekperp(:), hist_prf(:), &
                           hist_pnet(:), hist_time(:)

contains

  !--------------------------------------------------------------------
  ! ss_check — call once per time step, after all diagnostics are done.
  !
  !   itime      : step index starting at 1
  !   time       : current simulation time (s); used to measure the
  !                window's physical duration for the rate criterion
  !   tk         : total kinetic energy
  !   tkperp     : perpendicular kinetic energy
  !   pRF        : RF power absorbed by the plasma (0 when irf/=-1)
  !   p_net      : net power (sum of all signed power terms)
  !   p_drive    : positive power scale used to normalise the pRF and
  !                p_net criteria; call site should pass
  !                max(|pRF|, |pcoll_max|, |pcoll_self|, |psource|, 1.0)
  !   converged  : .TRUE. when all four per-second rates are below ss_tol
  !--------------------------------------------------------------------
  subroutine ss_check(itime, time, tk, tkperp, pRF, p_net, p_drive, converged)

    integer,  intent(in)  :: itime
    real(dp), intent(in)  :: time, tk, tkperp, pRF, p_net, p_drive
    logical,  intent(out) :: converged

    integer  :: ptr
    real(dp) :: ek_old, ekperp_old, prf_old, pnet_old, time_old, dt_win
    real(dp) :: rel_ek, rel_ekperp, rel_prf, rel_pnet, scale
    real(dp) :: rate_ek, rate_ekperp, rate_prf, rate_pnet

    ! (Re)allocate if this is the first call or window size changed.
    if (n_ss_window /= w_size) then
      if (allocated(hist_ek)) &
        deallocate(hist_ek, hist_ekperp, hist_prf, hist_pnet, hist_time)
      allocate(hist_ek(n_ss_window), hist_ekperp(n_ss_window), &
               hist_prf(n_ss_window), hist_pnet(n_ss_window), &
               hist_time(n_ss_window))
      hist_ek     = 0.0_dp
      hist_ekperp = 0.0_dp
      hist_prf    = 0.0_dp
      hist_pnet   = 0.0_dp
      hist_time   = 0.0_dp
      w_size      = n_ss_window
      write(*,'(A,I5,A,ES8.1,A)') &
        '  [SS check] window=', w_size, '  tol=', ss_tol, ' /s (rate-based)'
    end if

    converged = .false.

    ! Circular buffer slot for this step (1-based).
    ! The same slot holds the value from exactly w_size steps ago —
    ! read it before overwriting.
    ptr = mod(itime - 1, w_size) + 1

    ek_old     = hist_ek(ptr)
    ekperp_old = hist_ekperp(ptr)
    prf_old    = hist_prf(ptr)
    pnet_old   = hist_pnet(ptr)
    time_old   = hist_time(ptr)

    hist_ek(ptr)     = tk
    hist_ekperp(ptr) = tkperp
    hist_prf(ptr)    = pRF
    hist_pnet(ptr)   = p_net
    hist_time(ptr)   = time

    ! Do not evaluate until the buffer has been populated at least once.
    if (itime <= w_size) return

    ! Physical duration spanned by the window (correct across phase
    ! boundaries where dt changes). Guard against a non-positive span.
    dt_win = time - time_old
    if (dt_win <= 0.0_dp) return

    ! --- Relative change in total kinetic energy over the window -----
    scale  = 0.5_dp * (abs(tk) + abs(ek_old))
    rel_ek = merge(abs(tk - ek_old) / scale, 0.0_dp, scale > 0.0_dp)

    ! --- Relative change in perpendicular kinetic energy -------------
    scale      = 0.5_dp * (abs(tkperp) + abs(ekperp_old))
    rel_ekperp = merge(abs(tkperp - ekperp_old) / scale, 0.0_dp, scale > 0.0_dp)

    ! --- Relative change in RF power ---------------------------------
    ! Normalised to the caller-supplied driving-power scale.
    rel_prf = abs(pRF - prf_old) / max(p_drive, 1.0_dp)

    ! --- Relative change in net power --------------------------------
    rel_pnet = abs(p_net - pnet_old) / max(p_drive, 1.0_dp)

    ! --- Convert to per-second rates (dt-independent criterion) ------
    rate_ek     = rel_ek     / dt_win
    rate_ekperp = rel_ekperp / dt_win
    rate_prf    = rel_prf    / dt_win
    rate_pnet   = rel_pnet   / dt_win

    ! --- Periodic status print (every w_size steps) ------------------
    if (mod(itime, w_size) == 0) then
      write(*,'(A,I7,4(2X,A,ES9.2),2X,A,ES9.2,A)') &
        '  [SS]', itime, &
        'dE/E/s=',   rate_ek, &
        'dEp/Ep/s=', rate_ekperp, &
        'dPRF/Pd/s=',rate_prf, &
        'dP/Pd/s=',  rate_pnet, &
        '[tol=',   ss_tol, ' /s]'
    end if

    converged = (rate_ek     < ss_tol) .and. &
                (rate_ekperp < ss_tol) .and. &
                (rate_prf    < ss_tol) .and. &
                (rate_pnet   < ss_tol)

    if (converged) then
      write(*,'(A,I7,4(2X,A,ES9.2))') &
        '  [SS] CONVERGED at step', itime, &
        'dE/E/s=',   rate_ek, &
        'dEp/Ep/s=', rate_ekperp, &
        'dPRF/Pd/s=',rate_prf, &
        'dP/Pd/s=',  rate_pnet
    end if

  end subroutine ss_check

end module mod_ss_check
