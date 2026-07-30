!***********************************************************************
!  mod_conv_diag — physically consistent convergence diagnostics for the
!  time-dependent Fokker-Planck solvers.
!
!  Call conv_diag_init once, after make_grid (the weights depend on the
!  grid), then conv_diag_step once per time step after the solve, while
!  f^n is still available.
!
!  WHY THE JACOBIAN MATTERS
!  ------------------------
!  f lives in CYLINDRICAL velocity space, so the phase-space volume
!  element is
!
!      d3v = 2 pi vperp dvperp dvpar
!
!  Every norm and every moment here is integrated against that measure.
!  A plain sum over grid points would over-weight the near-axis region
!  (which holds almost no particles, since its phase-space volume goes to
!  zero) and under-weight the bulk, so a solution could look converged
!  while the physically dominant part of f was still moving.
!
!  QUADRATURE
!  ----------
!  The discretisation is finite difference on NODAL values (Fornberg
!  stencils in fd_stencil_2d act on point values f(i,j)), NOT finite
!  volume, so composite-trapezoidal weights on the non-uniform grids are
!  the consistent choice:
!
!      w_perp(1)     = (vperp(2)     - vperp(1))       / 2
!      w_perp(i)     = (vperp(i+1)   - vperp(i-1))     / 2     1<i<nperp
!      w_perp(nperp) = (vperp(nperp) - vperp(nperp-1)) / 2
!
!  (same on vpar; neither axis is assumed uniform).  The Jacobian is then
!  folded in:
!
!      W_perp(i) = 2 pi vperp(i) w_perp(i)
!      W_par(j)  = w_par(j)
!
!  and the 2D weight is W_perp(i)*W_par(j).  Because the trapezoidal rule
!  is exact for linear integrands and 2 pi vperp is linear in vperp,
!
!      SUM_i W_perp(i) == pi ( vperp(nperp)^2 - vperp(1)^2 )
!
!  holds to machine precision.  conv_diag_init checks exactly that, which
!  is a sharp test of both the Jacobian and the non-uniform weights.
!
!  WINDOWED EVALUATION
!  -------------------
!  The diagnostic is EVALUATED ONCE EVERY W = n_ss_window STEPS, not at
!  every step.  At each evaluation the current state is compared with the
!  state stored at the previous evaluation, and the difference is divided
!  by the PHYSICAL duration the window spans,
!
!      dt_win = time(now) - time(previous evaluation)
!
!  This is the same philosophy as the scalar rate-based criterion this
!  test replaces, and it buys three things.  Taking dt_win from the
!  simulation clock rather than W*dt keeps the rate correct across phase
!  boundaries where the time step changes.  A window wider than one step
!  is far less sensitive to short-term fluctuation: a single-step
!  increment can be small by accident (or oscillate about a mean that is
!  still drifting), whereas a window-long difference cannot.  And because
!  the full-grid norms are only formed on evaluation steps, the cost is
!  amortised over W steps instead of paid every step.
!
!  Only ONE snapshot of f has to be retained (the previous evaluation's
!  state), so the memory overhead is a single nperp*npar array regardless
!  of W.  W = 1 evaluates every step, which is what most unit tests use.
!
!  DIAGNOSTICS
!  -----------
!  1. eps      = ||f^now - f^prev||_2 / ( dt_win ( ||f^now||_2 + eps_abs ) )
!
!     Dividing by dt_win is essential: the scheme is implicit, so the raw
!     increment scales with the step size and an undivided criterion
!     reports false convergence whenever dt is reduced.  As written, eps
!     estimates ||df/dt||/||f|| and carries units of 1/time, so it is
!     dt-independent.  Dividing by ||f^n|| removes the dimensions of f and
!     any arbitrary normalisation of the density.  eps_abs is an absolute
!     floor guarding the f -> 0 case.
!
!     Since eps has units of 1/time it is only meaningful against a
!     physical rate, so eps_ratio = eps/nu_ref is reported too:
!     convergence then means "f changes negligibly on the timescale of
!     the physics" rather than being below an arbitrary number.
!
!  2. eps_tail: same, but with the extra weight
!
!         T(i,j) = ( 1 + (vperp(i)^2 + vpar(j)^2)/vth^2 )^k_tail
!
!     The Jacobian-weighted L2 norm is dominated by the bulk, so a run can
!     look converged there while the high-energy tail still evolves.
!     eps_tail is reported SEPARATELY and never merged into eps: the
!     interesting failure mode is precisely eps small with eps_tail large.
!
!  3. Moment drifts (density, parallel flow, energy), each as
!
!         drift_Q = |Q^now - Q^prev| / ( dt_win ( |Q^now| + eps_abs ) )
!
!     Norm convergence is necessary but not sufficient.  If eps
!     is small but a moment is still drifting, the fault is a flux or
!     boundary condition error at vperp = 0 or at the domain edges, NOT a
!     convergence failure -- look there first.
!
!     u_par uses a physical floor (1e-3 vth) in place of eps_abs, because
!     unlike n and E it is signed and legitimately zero for a vpar-
!     symmetric f; see the comment at its definition in conv_diag_step.
!
!  AXIS BLIND SPOT
!  ---------------
!  If vperp(1) = 0 then W_perp(1) = 0 and the whole i=1 row contributes
!  nothing to any norm or moment above.  That is physically correct (the
!  row has zero phase-space volume) but it means errors ON the axis are
!  invisible to eps, eps_tail and the moments.  axis_rate below is an
!  extra UNWEIGHTED max-norm over that row, reported as a safety check.
!
!  ROUND-OFF
!  ---------
!  All sums use Kahan compensated summation.  On a 200x150 grid a norm
!  accumulates 3e4 terms whose magnitudes span many orders (bulk vs tail),
!  and naive summation loses digits in exactly the small tail terms that
!  eps_tail is meant to resolve; compensated summation recovers them at
!  negligible cost and without needing a wider real kind.
!
!  07/2026: F. Louche
!***********************************************************************

module mod_conv_diag

  use shared_grid,  only: nperp, npar, vperp, vpar
  use shared_timer, only: outfile, n_ss_window, &
                          ss_tol_eps, ss_tol_tail, ss_tol_moment

  implicit none
  private

  public :: conv_diag_t
  public :: conv_diag_init, conv_diag_step, conv_diag_converged, conv_diag_finalize

  integer, parameter :: dp = kind(1.0d0)

  real(dp), parameter :: pi_c = 3.141592653589793238462643d0

  !--- Results of one step -----------------------------------------
  type :: conv_diag_t
    real(dp) :: eps            = 0.0_dp   ! normalised convergence rate   [1/s]
    real(dp) :: eps_ratio      = 0.0_dp   ! eps / nu_ref                  [-]
    real(dp) :: eps_tail       = 0.0_dp   ! tail-weighted rate            [1/s]
    real(dp) :: eps_tail_ratio = 0.0_dp   ! eps_tail / nu_ref             [-]
    real(dp) :: dens           = 0.0_dp   ! density                       [m^-3]
    real(dp) :: upar           = 0.0_dp   ! parallel flow velocity        [m/s]
    real(dp) :: energy         = 0.0_dp   ! kinetic energy density        [keV m^-3]
    real(dp) :: drift_dens     = 0.0_dp   ! relative drift rates          [1/s]
    real(dp) :: drift_upar     = 0.0_dp
    real(dp) :: drift_energy   = 0.0_dp
    real(dp) :: axis_rate      = 0.0_dp   ! unweighted max-norm rate on i=1 [1/s]
    real(dp) :: dt_win         = 0.0_dp   ! physical duration of the window [s]
    logical  :: valid          = .false.  ! .true. only on an evaluation step
                                          ! that has a previous checkpoint
  end type conv_diag_t

  !--- Module state (built once by conv_diag_init) -----------------
  logical               :: initialised = .false.
  real(dp), allocatable :: wperp(:)      ! W_perp(i) = 2 pi vperp(i) w_perp(i)
  real(dp), allocatable :: wpar(:)       ! W_par(j)  = w_par(j)
  real(dp), allocatable :: tailw(:,:)    ! T(i,j)

  !--- Reference scales and knobs ----------------------------------
  real(dp) :: nu_ref_m     = 1.0_dp      ! reference physical rate        [1/s]
  real(dp) :: vth_m        = 1.0_dp      ! reference thermal speed        [m/s]
  integer  :: k_tail_m     = 2           ! tail weight exponent
  real(dp) :: eps_abs_m    = 1.0d-300    ! absolute floor
  real(dp) :: tol_m        = 1.0d-3      ! tolerance on eps_ratio
  real(dp) :: tol_tail_m   = 1.0d-2      ! tolerance on eps_tail_ratio
  real(dp) :: tol_moment_m = 1.0d-3      ! tolerance on the moment drift ratios

  !--- Window checkpoint (state at the previous evaluation) ---------
  ! Because evaluation happens only every W steps, a single snapshot is
  ! enough -- no ring buffer, and the memory cost does not grow with W.
  integer               :: w_size = 0
  logical               :: have_chk = .false.
  real(dp), allocatable :: fchk(:,:)
  real(dp)              :: t_chk = 0.0_dp
  real(dp)              :: dens_chk = 0.0_dp, upar_chk = 0.0_dp, ener_chk = 0.0_dp

  !--- Output files ------------------------------------------------
  ! 519: conv_eps_vs_time.txt -- the epsilon time trace (time, eps, eps_tail),
  !      in the same 'time value ...' layout as the other *_vs_time.txt files
  !      so fp2d_plot.py picks it up automatically.
  ! 518: conv_diag_vs_time.csv -- the full history (every diagnostic), for
  !      post-processing.  Note outfile() only routes .txt to the null device
  !      under notxt, so the .csv is written regardless of notxt.
  integer, parameter :: eps_unit  = 519
  integer, parameter :: hist_unit = 518
  logical            :: eps_open  = .false.
  logical            :: hist_open = .false.

contains

  !--------------------------------------------------------------------
  ! kadd — one Kahan compensated-summation update: s = s + x, with c
  ! carrying the running round-off compensation.  See the ROUND-OFF note
  ! in the module header for why this is used rather than a naive sum.
  !--------------------------------------------------------------------
  pure subroutine kadd(s, c, x)
    real(dp), intent(inout) :: s, c
    real(dp), intent(in)    :: x
    real(dp) :: y, t
    y = x - c
    t = s + y
    c = (t - s) - y
    s = t
  end subroutine kadd

  !--------------------------------------------------------------------
  ! conv_diag_init — build and validate the quadrature weights and the
  ! tail weight.  Must be called AFTER make_grid and before the time
  ! loop; nothing is allocated or recomputed after this point.
  !
  ! The window width and the three tolerances come from the namelist
  ! (n_ss_window, ss_tol_eps, ss_tol_tail, ss_tol_moment in shared_timer);
  ! the optional arguments below override them, which is what the unit
  ! tests use.
  !
  !   vth        : reference thermal speed for the tail weight    [m/s]
  !   nu_ref     : reference physical rate (e.g. a collision       [1/s]
  !                frequency or inverse slowing-down time) used to
  !                render eps dimensionless
  !   k_tail     : optional, tail weight exponent            (default 2)
  !   eps_abs    : optional, absolute floor           (default 1e-300)
  !   tol        : optional, overrides ss_tol_eps
  !   tol_tail   : optional, overrides ss_tol_tail
  !   tol_moment : optional, overrides ss_tol_moment
  !   write_hist : optional, open the history files       (default .false.)
  !--------------------------------------------------------------------
  subroutine conv_diag_init(vth, nu_ref, k_tail, eps_abs, tol, tol_tail, &
                            tol_moment, write_hist)

    real(dp), intent(in)           :: vth, nu_ref
    integer,  intent(in), optional :: k_tail
    real(dp), intent(in), optional :: eps_abs, tol, tol_tail, tol_moment
    logical,  intent(in), optional :: write_hist

    integer  :: i, j
    real(dp) :: sum_wperp, sum_wpar, exact_perp, exact_par
    real(dp) :: err_perp, err_par, scale_p, scale_a, cc, v2

    if (allocated(wperp)) deallocate(wperp)
    if (allocated(wpar))  deallocate(wpar)
    if (allocated(tailw)) deallocate(tailw)
    allocate(wperp(nperp), wpar(npar), tailw(nperp,npar))

    !--- Window checkpoint -----------------------------------------
    w_size = max(n_ss_window, 1)
    if (allocated(fchk)) deallocate(fchk)
    allocate(fchk(nperp,npar))
    fchk     = 0.0_dp
    have_chk = .false.

    !--- Reference scales and knobs --------------------------------
    ! Tolerances default to the namelist values; the optionals override.
    vth_m        = vth
    nu_ref_m     = nu_ref
    tol_m        = ss_tol_eps
    tol_tail_m   = ss_tol_tail
    tol_moment_m = ss_tol_moment
    if (present(k_tail))     k_tail_m     = k_tail
    if (present(eps_abs))    eps_abs_m    = eps_abs
    if (present(tol))        tol_m        = tol
    if (present(tol_tail))   tol_tail_m   = tol_tail
    if (present(tol_moment)) tol_moment_m = tol_moment

    if (vth_m    <= 0.0_dp) vth_m    = 1.0_dp     ! keep T(i,j) finite
    if (nu_ref_m <= 0.0_dp) nu_ref_m = 1.0_dp     ! ratios degrade to raw rates

    !--- Composite-trapezoidal weights x Jacobian ------------------
    if (nperp == 1) then
      wperp(1) = 0.0_dp
    else
      wperp(1)     = 0.5_dp * (vperp(2)     - vperp(1))
      do i = 2, nperp - 1
        wperp(i)   = 0.5_dp * (vperp(i+1)   - vperp(i-1))
      end do
      wperp(nperp) = 0.5_dp * (vperp(nperp) - vperp(nperp-1))
    end if
    do i = 1, nperp
      wperp(i) = 2.0_dp * pi_c * vperp(i) * wperp(i)
    end do

    if (npar == 1) then
      wpar(1) = 0.0_dp
    else
      wpar(1)    = 0.5_dp * (vpar(2)    - vpar(1))
      do j = 2, npar - 1
        wpar(j)  = 0.5_dp * (vpar(j+1)  - vpar(j-1))
      end do
      wpar(npar) = 0.5_dp * (vpar(npar) - vpar(npar-1))
    end if

    !--- Tail weight T(i,j) ----------------------------------------
    do i = 1, nperp
      do j = 1, npar
        v2 = vperp(i)**2 + vpar(j)**2
        tailw(i,j) = (1.0_dp + v2 / vth_m**2)**k_tail_m
      end do
    end do

    !--- Validation: both identities are EXACT, not approximate ----
    sum_wperp = 0.0_dp; cc = 0.0_dp
    do i = 1, nperp
      call kadd(sum_wperp, cc, wperp(i))
    end do
    sum_wpar = 0.0_dp; cc = 0.0_dp
    do j = 1, npar
      call kadd(sum_wpar, cc, wpar(j))
    end do

    exact_perp = pi_c * (vperp(nperp)**2 - vperp(1)**2)
    exact_par  = vpar(npar) - vpar(1)

    scale_p  = max(abs(exact_perp), tiny(1.0_dp))
    scale_a  = max(abs(exact_par),  tiny(1.0_dp))
    err_perp = abs(sum_wperp - exact_perp) / scale_p
    err_par  = abs(sum_wpar  - exact_par ) / scale_a

    write(*,'(A,I5,A,F8.2,A)') '  [conv_diag] evaluated every', w_size, &
      ' steps   (one snapshot of f: ', &
      real(nperp,dp)*real(npar,dp)*8.0_dp/1.048576d6, ' MB)'
    write(*,'(A,3(2X,A,ES9.2))') '  [conv_diag] tolerances:', &
      'eps/nu <', tol_m, 'epsT/nu <', tol_tail_m, 'moments <', tol_moment_m
    write(*,'(A)')        '  [conv_diag] quadrature weights built:'
    write(*,'(A,ES12.5,A,ES12.5,A,ES9.2)') &
      '    sum(W_perp) =', sum_wperp, '   exact =', exact_perp, '   rel.err =', err_perp
    write(*,'(A,ES12.5,A,ES12.5,A,ES9.2)') &
      '    sum(W_par ) =', sum_wpar,  '   exact =', exact_par,  '   rel.err =', err_par
    if (err_perp > 1.0d-10 .or. err_par > 1.0d-10) &
      write(*,'(A)') '    *** WARNING: weight identity violated -- check the grid ***'

    ! Axis blind spot: flag it explicitly at run time as well as in the
    ! module header, since it changes how the numbers should be read.
    if (vperp(1) == 0.0_dp) &
      write(*,'(A)') '    note: vperp(1)=0 -> W_perp(1)=0; the axis row is invisible' // &
                     ' to the weighted norms (see axis_rate).'

    if (present(write_hist)) then
      if (write_hist) then
        open(eps_unit, file=TRIM(outfile('conv_eps_vs_time.txt')), status='unknown')
        write(eps_unit,'(A)') '# time[s]   eps[1/s]   eps_tail[1/s]'
        eps_open = .true.

        open(hist_unit, file=TRIM(outfile('conv_diag_vs_time.csv')), status='unknown')
        write(hist_unit,'(A)') 'time,dt,eps,eps_ratio,eps_tail,eps_tail_ratio,' // &
                               'dens,upar,energy,drift_dens,drift_upar,drift_energy,axis_rate'
        hist_open = .true.
      end if
    end if

    initialised = .true.

  end subroutine conv_diag_init

  !--------------------------------------------------------------------
  ! conv_diag_step — all diagnostics for one step.
  !
  ! Call it every step: the routine itself decides when to evaluate, and
  ! returns immediately (silently, d%valid = .false.) on the W-1 steps
  ! between evaluations, so the full-grid norms are formed only once per
  ! window.
  !
  !   itime  : global step index (1 on the first step)
  !   time   : current simulation time                            [s]
  !   f_new  : f at the current step, on the (nperp,npar) grid
  !   d      : all diagnostics; d%valid is .true. only on an evaluation
  !            step that has a previous checkpoint to compare against and
  !            a positive dt_win
  !
  ! The comparison state is the module's own checkpoint from the previous
  ! evaluation, so no previous iterate needs to be passed in.
  !--------------------------------------------------------------------
  subroutine conv_diag_step(itime, time, f_new, d)

    integer,           intent(in)  :: itime
    real(dp),          intent(in)  :: time
    real(dp),          intent(in)  :: f_new(nperp,npar)
    type(conv_diag_t), intent(out) :: d

    real(dp), parameter :: pmass    = 1.6726d-27   ! proton mass [kg]
    real(dp), parameter :: kev_in_J = 1.60218d-16  ! keV -> J

    integer  :: i, j
    real(dp) :: w, df, fn, v2, dt_win
    real(dp) :: s_dnorm, c_dnorm, s_fnorm, c_fnorm
    real(dp) :: s_dtail, c_dtail, s_ftail, c_ftail
    real(dp) :: s_dens,  c_dens,  s_mom,   c_mom,  s_ener, c_ener
    real(dp) :: nrm_df, nrm_f, nrm_df_t, nrm_f_t
    real(dp) :: dens, upar, ener, amax_df, amax_f
    real(dp) :: dens_old, upar_old, ener_old

    d = conv_diag_t()          ! default-initialised: valid = .false.

    if (.not. initialised) then
      write(*,'(A)') '  [conv_diag] ERROR: conv_diag_step called before conv_diag_init.'
      return
    end if

    ! Evaluate only on window boundaries.  Returning here -- BEFORE the
    ! grid loops -- is what makes the diagnostic cost one full-grid pass
    ! per window rather than one per step.
    if (mod(itime, w_size) /= 0) return

    s_dnorm = 0.0_dp; c_dnorm = 0.0_dp
    s_fnorm = 0.0_dp; c_fnorm = 0.0_dp
    s_dtail = 0.0_dp; c_dtail = 0.0_dp
    s_ftail = 0.0_dp; c_ftail = 0.0_dp
    s_dens  = 0.0_dp; c_dens  = 0.0_dp
    s_mom   = 0.0_dp; c_mom   = 0.0_dp
    s_ener  = 0.0_dp; c_ener  = 0.0_dp
    amax_df = 0.0_dp; amax_f  = 0.0_dp

    do i = 1, nperp
      do j = 1, npar
        w  = wperp(i) * wpar(j)
        fn = f_new(i,j)
        df = fn - fchk(i,j)           ! change since the previous evaluation
        v2 = vperp(i)**2 + vpar(j)**2

        call kadd(s_dnorm, c_dnorm, w * df * df)
        call kadd(s_fnorm, c_fnorm, w * fn * fn)
        call kadd(s_dtail, c_dtail, w * tailw(i,j) * df * df)
        call kadd(s_ftail, c_ftail, w * tailw(i,j) * fn * fn)

        call kadd(s_dens,  c_dens,  w * fn)
        call kadd(s_mom,   c_mom,   w * vpar(j) * fn)
        call kadd(s_ener,  c_ener,  w * v2 * fn)

        ! Unweighted axis-row max-norm: the ONLY diagnostic here that
        ! can see the i=1 row when vperp(1)=0 (W_perp(1)=0 there).
        if (i == 1) then
          amax_df = max(amax_df, abs(df))
          amax_f  = max(amax_f,  abs(fn))
        end if
      end do
    end do

    nrm_df   = sqrt(max(s_dnorm, 0.0_dp))
    nrm_f    = sqrt(max(s_fnorm, 0.0_dp))
    nrm_df_t = sqrt(max(s_dtail, 0.0_dp))
    nrm_f_t  = sqrt(max(s_ftail, 0.0_dp))

    dens = s_dens
    if (abs(dens) > 0.0_dp) then
      upar = s_mom / dens
    else
      upar = 0.0_dp
    end if
    ! Energy density in keV m^-3, using the codebase convention
    ! (0.5 m_p A v^2 converted to keV; cf. time_energy in time_comps_mod).
    ener = 0.5_dp * pmass * aa_ref() * s_ener / kev_in_J

    d%dens   = dens
    d%upar   = upar
    d%energy = ener

    ! Read the previous checkpoint, then refresh it with the current
    ! state.  Order matters: the checkpoint is both the comparison state
    ! and this evaluation's destination.
    dens_old = dens_chk
    upar_old = upar_chk
    ener_old = ener_chk
    dt_win   = time - t_chk

    fchk     = f_new
    t_chk    = time
    dens_chk = dens
    upar_chk = upar
    ener_chk = ener

    d%dt_win = dt_win

    ! The first evaluation only lays down the checkpoint -- there is
    ! nothing to compare against yet.  The window must also span a
    ! positive duration.
    if (.not. have_chk .or. dt_win <= 0.0_dp) then
      have_chk = .true.
      d%valid  = .false.
      call conv_diag_log(time, d)
      return
    end if

    d%eps      = nrm_df   / (dt_win * (nrm_f   + eps_abs_m))
    d%eps_tail = nrm_df_t / (dt_win * (nrm_f_t + eps_abs_m))
    d%eps_ratio      = d%eps      / nu_ref_m
    d%eps_tail_ratio = d%eps_tail / nu_ref_m

    ! Density and energy are positive-definite and O(1) in their own units,
    ! so the bare eps_abs floor is enough for them.
    d%drift_dens   = abs(dens - dens_old) / (dt_win * (abs(dens) + eps_abs_m))
    d%drift_energy = abs(ener - ener_old) / (dt_win * (abs(ener) + eps_abs_m))

    ! u_par is DIFFERENT: it is a signed quantity that legitimately passes
    ! through zero (any distribution symmetric in vpar has u_par = 0 up to
    ! round-off).  Normalising it by |u_par| would then divide a round-off
    ! numerator by a round-off denominator and return a large meaningless
    ! number that never decays -- so conv_diag_converged could never be
    ! satisfied for a symmetric case.  A single global eps_abs cannot floor
    ! both a density of ~1e19 and a velocity of ~0, so u_par gets its own
    ! physical floor: a thousandth of the reference thermal speed.  Below
    ! that, this measures d(u_par/u_floor)/dt; above it, it degrades
    ! smoothly to the same relative rate as the other moments.
    d%drift_upar   = abs(upar - upar_old) &
                     / (dt_win * (abs(upar) + max(eps_abs_m, 1.0d-3 * vth_m)))

    d%axis_rate = amax_df / (dt_win * (amax_f + eps_abs_m))

    d%valid = .true.

    call conv_diag_log(time, d)

  end subroutine conv_diag_step

  !--------------------------------------------------------------------
  ! aa_ref — minority mass number, read from shared_plasma.  Isolated in
  ! a function so the module's only compile-time coupling to the plasma
  ! data is this one line.
  !--------------------------------------------------------------------
  function aa_ref() result(a)
    use shared_plasma, only: aa
    real(dp) :: a
    a = aa
  end function aa_ref

  !--------------------------------------------------------------------
  ! conv_diag_converged — all criteria simultaneously.  Deliberately
  ! separate from conv_diag_step so the caller can log without stopping.
  !--------------------------------------------------------------------
  logical function conv_diag_converged(d) result(ok)

    type(conv_diag_t), intent(in) :: d

    ok = .false.
    if (.not. d%valid) return

    ok =       (d%eps_ratio      <  tol_m)        &
         .and. (d%eps_tail_ratio <  tol_tail_m)   &
         .and. (d%drift_dens     <  tol_moment_m) &
         .and. (d%drift_upar     <  tol_moment_m) &
         .and. (d%drift_energy   <  tol_moment_m)

  end function conv_diag_converged

  !--------------------------------------------------------------------
  ! conv_diag_log — one concise line to stdout, plus the CSV history row
  ! when the history file was requested at init.
  !--------------------------------------------------------------------
  subroutine conv_diag_log(time, d)

    real(dp),          intent(in) :: time
    type(conv_diag_t), intent(in) :: d

    if (d%valid) then
      ! eps and eps_tail are the raw rates [1/s]; the /nu columns are the
      ! same numbers made dimensionless against the reference rate, which
      ! is what the tolerances are actually applied to.
      write(*,'(A,ES11.4,A,ES10.3,A,ES10.3,A,ES9.2,A,ES9.2,A,ES9.2,A,ES9.2,A,ES9.2,A,ES9.2)') &
        '  [conv] t=', time,                    &
        ' s  eps=',         d%eps,              &
        ' /s  epsT=',       d%eps_tail,         &
        ' /s  eps/nu=',     d%eps_ratio,        &
        '  epsT/nu=',       d%eps_tail_ratio,   &
        '  dn=',            d%drift_dens,       &
        '  du=',            d%drift_upar,       &
        '  dE=',            d%drift_energy,     &
        '  axis=',          d%axis_rate
    else
      write(*,'(A,ES11.4,A,I0,A)') &
        '  [conv] t=', time, ' s  (first checkpoint laid down; next ', &
        w_size, '-step window gives the first verdict)'
    end if

    ! The .txt trace carries only the converged-on quantities and starts at
    ! the first VALID step: a leading eps=0 row would otherwise plot as a
    ! spurious dropout on the log axis fp2d_plot.py uses for it.
    if (eps_open .and. d%valid) then
      write(eps_unit,'(3(1X,ES15.7))') time, d%eps, d%eps_tail
    end if

    if (hist_open) then
      write(hist_unit,'(ES15.7,12(",",ES15.7))') &
        time, d%dt_win, d%eps, d%eps_ratio, d%eps_tail, d%eps_tail_ratio, &
        d%dens, d%upar, d%energy, &
        d%drift_dens, d%drift_upar, d%drift_energy, d%axis_rate
    end if

  end subroutine conv_diag_log

  !--------------------------------------------------------------------
  ! conv_diag_finalize — close the history file and release the weights.
  !--------------------------------------------------------------------
  subroutine conv_diag_finalize()

    if (eps_open) then
      close(eps_unit)
      eps_open = .false.
    end if
    if (hist_open) then
      close(hist_unit)
      hist_open = .false.
    end if
    if (allocated(wperp)) deallocate(wperp)
    if (allocated(wpar))  deallocate(wpar)
    if (allocated(tailw)) deallocate(tailw)
    if (allocated(fchk))  deallocate(fchk)
    initialised = .false.
    have_chk    = .false.
    w_size      = 0

  end subroutine conv_diag_finalize

end module mod_conv_diag
