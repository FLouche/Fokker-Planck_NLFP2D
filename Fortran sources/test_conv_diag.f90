!***********************************************************************
!  test_conv_diag — standalone verification of mod_conv_diag.
!
!  NOT part of the FP2D_QLRF_NL.vfproj build: it contains its own PROGRAM
!  unit and would clash with main-FP_Coll_2D.  Build and run it on its own:
!
!    ifx /nologo /module:<objdir> /object:<objdir> ^
!        "shared_data.f90" "conv_diag.f90" "test_conv_diag.f90" ^
!        /exe:<objdir>\test_conv_diag.exe
!
!  It sets up shared_grid itself (nperp, npar, vperp, vpar) so it needs
!  no solver, no PARDISO and no MKL.
!
!  Tests
!  -----
!   1  sum(W_perp) == pi ( vperp(N)^2 - vperp(1)^2 )   exactly
!   2  sum(W_par)  == vpar(N) - vpar(1)                exactly
!   3  density and energy of an analytic Maxwellian
!   4  dt-independence of eps
!   5  tail vs bulk sensitivity of eps_tail relative to eps
!
!  07/2026: F. Louche
!***********************************************************************

program test_conv_diag

  use shared_grid,  only: nperp, npar, vperp, vpar
  use shared_plasma, only: aa
  use shared_timer,  only: n_ss_window
  use mod_conv_diag

  implicit none

  integer,  parameter :: dp = kind(1.0d0)
  real(dp), parameter :: pi_c = 3.141592653589793238462643d0

  integer  :: i, j
  integer  :: nfail
  real(dp) :: vth, nu_ref, vpmax, vpamax
  real(dp) :: sum_wp, sum_wa, exact_wp, exact_wa, err

  real(dp), allocatable :: f0(:,:), f1(:,:), fflat(:)
  type(conv_diag_t) :: d, d_a, d_b, d_tail, d_bulk

  nfail = 0

  !===================================================================
  ! Grid: deliberately NON-UNIFORM in both directions, so the tests
  ! exercise the general weights rather than a uniform special case.
  !===================================================================
  nperp  = 160
  npar   = 121
  vpmax  = 3.0d6
  vpamax = 3.0d6
  aa     = 2.0d0                       ! deuterium, as in the D-beam cases

  allocate(vperp(nperp), vpar(npar))

  ! vperp: quadratic stretching from 0 (mimics ising=+1), vperp(1)=0 so
  ! the axis blind spot is exercised too.
  do i = 1, nperp
    vperp(i) = vpmax * (real(i-1,dp) / real(nperp-1,dp))**2
  end do

  ! vpar: smoothly stretched and asymmetric about 0 -- NOT uniform.
  do j = 1, npar
    vpar(j) = -vpamax + 2.0_dp*vpamax * &
              ( real(j-1,dp)/real(npar-1,dp) )**1.3_dp
  end do

  vth    = 5.0d5
  nu_ref = 1.0d2

  allocate(f0(nperp,npar), f1(nperp,npar), fflat(nperp*npar))

  write(*,'(A)') '======================================================'
  write(*,'(A)') ' test_conv_diag'
  write(*,'(A,I0,A,I0,A)') '   grid ', nperp, ' x ', npar, &
                           '  (vperp quadratic, vpar stretched)'
  write(*,'(A)') '======================================================'

  ! W = 1 makes the test evaluate every step (plain step-to-step),
  ! which is what tests 3-5 assume; test 6 re-inits with a wider window.
  n_ss_window = 1
  call conv_diag_init(vth, nu_ref, k_tail=2, tol=1.0d-3, tol_tail=1.0d-2, &
                      tol_moment=1.0d-3)

  !===================================================================
  ! 1 & 2 — weight-sum identities.
  ! The trapezoidal rule is EXACT for linear integrands and 2 pi vperp is
  ! linear, so these are identities, not approximations.  Tolerance is
  ! pure round-off: a few ULP amplified by the ~1e13 dynamic range of the
  ! accumulation, so 1e-14 relative is the right bar.
  !===================================================================
  sum_wp = 0.0_dp
  do i = 1, nperp
    sum_wp = sum_wp + wperp_of(i)
  end do
  sum_wa = 0.0_dp
  do j = 1, npar
    sum_wa = sum_wa + wpar_of(j)
  end do

  exact_wp = pi_c * (vperp(nperp)**2 - vperp(1)**2)
  exact_wa = vpar(npar) - vpar(1)

  err = abs(sum_wp - exact_wp) / abs(exact_wp)
  call check('1  sum(W_perp) = pi (vperp_N^2 - vperp_1^2)', err, 1.0d-14, nfail)
  write(*,'(A,ES22.15)') '      computed = ', sum_wp
  write(*,'(A,ES22.15)') '      exact    = ', exact_wp

  err = abs(sum_wa - exact_wa) / abs(exact_wa)
  call check('2  sum(W_par)  = vpar_N - vpar_1', err, 1.0d-14, nfail)
  write(*,'(A,ES22.15)') '      computed = ', sum_wa
  write(*,'(A,ES22.15)') '      exact    = ', exact_wa

  !===================================================================
  ! 3 — moment recovery for an analytic Maxwellian
  !
  !   f = n0 / ( (2 pi)^{3/2} vth^3 ) exp( -v^2 / (2 vth^2) )
  !
  ! Exact moments:  density = n0
  !                 energy  = 0.5 m A (3 vth^2) n0   [J m^-3] -> keV
  !
  ! Tolerance: this is a DISCRETISATION test, not a round-off test.  The
  ! trapezoidal rule on a smooth, well-resolved integrand converges as
  ! O(h^2), and the domain truncates the Maxwellian at 6 vth (where the
  ! integrand is ~1e-8 of its peak).  With the quadratic vperp grid the
  ! near-axis spacing is fine but the outer spacing is coarse, so 1e-3
  ! relative is the honest bar; anything much tighter would be testing
  ! the grid, not the weights.
  !===================================================================
  call maxwellian(f0, 1.0d19, vth)
  call flatten(f0, fflat)
  call conv_diag_step(1, 0.0_dp, f0, d)

  err = abs(d%dens - 1.0d19) / 1.0d19
  call check('3a density of an analytic Maxwellian', err, 1.0d-3, nfail)
  write(*,'(A,ES22.15)') '      computed = ', d%dens
  write(*,'(A,ES22.15)') '      exact    = ', 1.0d19

  err = abs(d%energy - maxw_energy(1.0d19, vth)) / maxw_energy(1.0d19, vth)
  call check('3b energy  of an analytic Maxwellian', err, 1.0d-3, nfail)
  write(*,'(A,ES22.15)') '      computed = ', d%energy
  write(*,'(A,ES22.15)') '      exact    = ', maxw_energy(1.0d19, vth)

  !===================================================================
  ! 4 — dt-independence.
  ! Known smooth evolution f(t+dt) = f(t) (1 + alpha dt), so
  !
  !     ||df|| = alpha dt ||f_old||,   ||f_new|| = (1+alpha dt) ||f_old||
  !     =>  eps = alpha / (1 + alpha dt)
  !
  ! eps is therefore dt-independent to O(alpha dt) -- NOT exactly, because
  ! the definition normalises by ||f_new|| (the new iterate) rather than
  ! ||f_old||.  That is the specified form and the right one: it is what
  ! makes eps a backward-difference estimate of ||df/dt||/||f||.  The test
  ! therefore probes the regime the criterion is meant for, alpha dt << 1,
  ! where "rate of change" is meaningful at all; here alpha dt <= 3e-4, so
  ! eps must agree to ~1e-3 relative across FOUR decades of dt.
  !
  ! The contrast printed below is the real point: the UNDIVIDED increment
  ! ||df||/||f_new|| changes by four decades over the same dt range, which
  ! is exactly the false convergence an undivided criterion would report
  ! whenever dt is reduced.
  !===================================================================
  call run_pair(f0, 1.0d-6, vth, d_a)
  call run_pair(f0, 1.0d-2, vth, d_b)

  err = abs(d_a%eps - d_b%eps) / d_a%eps
  call check('4  eps invariant over 4 decades of dt', err, 1.0d-3, nfail)
  write(*,'(A,ES22.15)') '      eps(dt=1e-6) = ', d_a%eps
  write(*,'(A,ES22.15)') '      eps(dt=1e-2) = ', d_b%eps
  write(*,'(A,ES12.5,A,ES12.5,A,ES9.2)')                     &
        '      undivided ||df||/||f||: ', d_a%eps*1.0d-6,    &
        '  vs ', d_b%eps*1.0d-2, '   ratio ',                &
        (d_b%eps*1.0d-2)/(d_a%eps*1.0d-6)

  !===================================================================
  ! 5 — tail vs bulk sensitivity.
  ! A perturbation localised at high v must give eps_tail >> eps; one
  ! localised in the bulk must give eps_tail <~ eps.  This is the failure
  ! mode the tail norm exists to catch, so both orderings are asserted.
  !===================================================================
  call perturb(f0, f1, 4.5_dp, 0.5_dp, vth)      ! shell at v = 4.5 vth
  call flatten(f0, fflat)
  call conv_diag_step(1, 0.0_dp, f0, d)     ! prime
  call conv_diag_step(2, 1.0d-3, f1, d_tail)

  call perturb(f0, f1, 0.5_dp, 0.5_dp, vth)      ! shell at v = 0.5 vth
  call flatten(f0, fflat)
  call conv_diag_step(1, 0.0_dp, f0, d)     ! prime
  call conv_diag_step(2, 1.0d-3, f1, d_bulk)

  write(*,'(A)') ''
  write(*,'(A)') '  5  tail vs bulk perturbation'
  write(*,'(A,ES12.5,A,ES12.5,A,F9.3)') '      tail: eps=', d_tail%eps, &
        '  eps_tail=', d_tail%eps_tail, '   ratio=', d_tail%eps_tail/d_tail%eps
  write(*,'(A,ES12.5,A,ES12.5,A,F9.3)') '      bulk: eps=', d_bulk%eps, &
        '  eps_tail=', d_bulk%eps_tail, '   ratio=', d_bulk%eps_tail/d_bulk%eps

  if (d_tail%eps_tail / d_tail%eps > 5.0_dp) then
    write(*,'(A)') '      PASS  tail perturbation: eps_tail >> eps'
  else
    write(*,'(A)') '      FAIL  tail perturbation did not amplify eps_tail'
    nfail = nfail + 1
  end if

  if (d_bulk%eps_tail / d_bulk%eps < d_tail%eps_tail / d_tail%eps) then
    write(*,'(A)') '      PASS  bulk perturbation gives the opposite ordering'
  else
    write(*,'(A)') '      FAIL  bulk perturbation did not give the opposite ordering'
    nfail = nfail + 1
  end if

  !===================================================================
  ! 6 — windowed evaluation.
  ! Re-init with W = 10 and drive a known steady exponential-in-time
  ! growth f(t) = f0 (1 + alpha t) sampled every dt.  Two assertions:
  !   (a) no verdict is issued while the window is still filling
  !       (itime <= W), so an early accidental small increment can never
  !       stop the run;
  !   (b) once full, the reported rate matches the analytic window rate
  !         ||f(t) - f(t-W dt)|| / (W dt ||f(t)||)
  !         = alpha W dt / (W dt (1 + alpha t)) = alpha / (1 + alpha t),
  !       i.e. the window measures the SAME physical rate as a single
  !       step, which is the point of dividing by the window duration
  !       rather than by dt.
  !===================================================================
  call conv_diag_finalize()
  n_ss_window = 10
  call conv_diag_init(vth, nu_ref, k_tail=2, tol=1.0d-3, tol_tail=1.0d-2, &
                      tol_moment=1.0d-3)
  call maxwellian(f0, 1.0d19, vth)
  call window_run(f0, 10, 1.0d-3, 3.0d-2, nfail)

  !===================================================================
  ! 7 — shape variant under a pure drain.
  ! Reproduces the failure mode found in the sourceless RF run: f decays
  ! as a FIXED shape with falling amplitude, f(t) = phi(v) exp(-t/tau).
  ! The amplitude rate eps must then sit at the drain rate 1/tau and stay
  ! there, while the shape rate eps_shape must be ~0 (machine epsilon),
  ! since normalising removes the amplitude entirely.  Without this the
  ! run could never satisfy a tolerance below 1/tau.
  !===================================================================
  call conv_diag_finalize()
  n_ss_window = 5
  call conv_diag_init(vth, nu_ref, k_tail=2, tol=1.0d-3, tol_tail=1.0d-2, &
                      tol_moment=1.0d-3)
  call maxwellian(f0, 1.0d19, vth)
  call drain_run(f0, 5, 1.0d-2, 1.7393d-3, nfail)

  !===================================================================
  write(*,'(A)') '======================================================'
  if (nfail == 0) then
    write(*,'(A)') ' ALL TESTS PASSED'
  else
    write(*,'(A,I0,A)') ' ', nfail, ' TEST(S) FAILED'
  end if
  write(*,'(A)') '======================================================'

  call conv_diag_finalize()

  if (nfail /= 0) stop 1

contains

  !--- Weight accessors: rebuild the same trapezoid x Jacobian weights
  !    INDEPENDENTLY of the module, so tests 1 and 2 verify the identity
  !    against a separate implementation rather than against themselves.
  function wperp_of(i) result(w)
    integer, intent(in) :: i
    real(dp) :: w
    if (i == 1) then
      w = 0.5_dp * (vperp(2) - vperp(1))
    else if (i == nperp) then
      w = 0.5_dp * (vperp(nperp) - vperp(nperp-1))
    else
      w = 0.5_dp * (vperp(i+1) - vperp(i-1))
    end if
    w = 2.0_dp * pi_c * vperp(i) * w
  end function wperp_of

  function wpar_of(j) result(w)
    integer, intent(in) :: j
    real(dp) :: w
    if (j == 1) then
      w = 0.5_dp * (vpar(2) - vpar(1))
    else if (j == npar) then
      w = 0.5_dp * (vpar(npar) - vpar(npar-1))
    else
      w = 0.5_dp * (vpar(j+1) - vpar(j-1))
    end if
  end function wpar_of

  subroutine maxwellian(f, n0, vt)
    real(dp), intent(out) :: f(nperp,npar)
    real(dp), intent(in)  :: n0, vt
    real(dp) :: nrm, v2
    integer  :: ii, jj
    nrm = n0 / ((2.0_dp*pi_c)**1.5_dp * vt**3)
    do ii = 1, nperp
      do jj = 1, npar
        v2 = vperp(ii)**2 + vpar(jj)**2
        f(ii,jj) = nrm * exp(-v2 / (2.0_dp*vt**2))
      end do
    end do
  end subroutine maxwellian

  ! Exact energy density of the above Maxwellian, in the codebase's keV
  ! convention: E = 0.5 m_p A <v^2> n0 / kev_in_J with <v^2> = 3 vth^2.
  function maxw_energy(n0, vt) result(e)
    real(dp), intent(in) :: n0, vt
    real(dp) :: e
    real(dp), parameter :: pmass = 1.6726d-27, kev_in_J = 1.60218d-16
    e = 0.5_dp * pmass * aa * 3.0_dp * vt**2 * n0 / kev_in_J
  end function maxw_energy

  subroutine flatten(f, fv)
    real(dp), intent(in)  :: f(nperp,npar)
    real(dp), intent(out) :: fv(:)
    integer :: ii, jj
    do ii = 1, nperp
      do jj = 1, npar
        fv((ii-1)*npar + jj) = f(ii,jj)
      end do
    end do
  end subroutine flatten

  ! f_new = f_old (1 + alpha dt): a known smooth evolution whose exact
  ! eps is |alpha| independent of dt.
  subroutine run_pair(fbase, dt, vt, dd)
    real(dp), intent(in)  :: fbase(nperp,npar), dt, vt
    type(conv_diag_t), intent(out) :: dd
    real(dp), parameter :: alpha = 3.0d-2
    real(dp) :: fn(nperp,npar), fv(nperp*npar)
    type(conv_diag_t) :: dprime
    fn = fbase * (1.0_dp + alpha*dt)
    call flatten(fbase, fv)
    call conv_diag_step(1, 0.0_dp, fbase, dprime)   ! prime (fills the W=1 window)
    call conv_diag_step(2, dt,          fn,    dd)
  end subroutine run_pair

  ! Gaussian shell at v = vc*vth, width wid*vth, amplitude 1% of local f.
  subroutine perturb(fin, fout_, vc, wid, vt)
    real(dp), intent(in)  :: fin(nperp,npar), vc, wid, vt
    real(dp), intent(out) :: fout_(nperp,npar)
    integer  :: ii, jj
    real(dp) :: v, g
    do ii = 1, nperp
      do jj = 1, npar
        v = sqrt(vperp(ii)**2 + vpar(jj)**2) / vt
        g = exp(-((v - vc)/wid)**2)
        fout_(ii,jj) = fin(ii,jj) * (1.0_dp + 1.0d-2 * g)
      end do
    end do
  end subroutine perturb

  ! Drive f(t) = fbase (1 + alpha t) for 3W steps and check the window
  ! behaviour: invalid while filling, analytic rate once full.
  subroutine window_run(fbase, w, dt, alpha, nf)
    real(dp), intent(in)    :: fbase(nperp,npar), dt, alpha
    integer,  intent(in)    :: w
    integer,  intent(inout) :: nf
    real(dp) :: fn(nperp,npar)
    type(conv_diag_t) :: dd
    real(dp) :: t, expect, relerr, worst
    integer  :: n, nverdict
    logical  :: early_verdict
    early_verdict = .false.
    worst   = 0.0_dp
    nverdict = 0
    write(*,'(A)') ''
    write(*,'(A,I0,A)') '  6  windowed evaluation (every W = ', w, ' steps)'
    do n = 1, 3*w
      t  = n*dt
      fn = fbase * (1.0_dp + alpha*t)
      call conv_diag_step(n, t, fn, dd)
      ! Nothing may be reported on non-window steps, and the first window
      ! only lays down the checkpoint.
      if (mod(n, w) /= 0 .and. dd%valid) early_verdict = .true.
      if (n <= w .and. dd%valid)         early_verdict = .true.
      if (dd%valid) then
        nverdict = nverdict + 1
        ! f(t) = fbase (1 + alpha t) compared over [t-W dt, t]:
        !   ||df|| / (dt_win ||f(t)||) = alpha / (1 + alpha t)
        expect = alpha / (1.0_dp + alpha*t)
        relerr = abs(dd%eps - expect) / expect
        worst  = max(worst, relerr)
      end if
    end do
    if (.not. early_verdict) then
      write(*,'(A,I0,A)') '      PASS  no verdict off a window boundary or in the first ', &
                          w, ' steps'
    else
      write(*,'(A)') '      FAIL  a verdict was issued off a window boundary'
      nf = nf + 1
    end if
    ! 3W steps, evaluated every W, first one only checkpoints -> 2 verdicts
    if (nverdict == 2) then
      write(*,'(A,I0,A)') '      PASS  exactly ', nverdict, &
                          ' verdicts over 3 windows (first only checkpoints)'
    else
      write(*,'(A,I0,A)') '      FAIL  expected 2 verdicts over 3 windows, got ', &
                          nverdict, ''
      nf = nf + 1
    end if
    if (worst <= 1.0d-3) then
      write(*,'(A,ES9.2,A)') '      PASS  windowed rate matches analytic (worst rel.err ', &
                             worst, ')'
    else
      write(*,'(A,ES9.2,A)') '      FAIL  windowed rate off analytic (worst rel.err ', &
                             worst, ')'
      nf = nf + 1
    end if
  end subroutine window_run

  ! Pure exponential drain at rate gam: f(t) = fbase exp(-gam t).  Shape
  ! is constant, so eps -> gam and eps_shape -> 0.
  subroutine drain_run(fbase, w, dt, gam, nf)
    real(dp), intent(in)    :: fbase(nperp,npar), dt, gam
    integer,  intent(in)    :: w
    integer,  intent(inout) :: nf
    real(dp) :: fn(nperp,npar)
    type(conv_diag_t) :: dd
    real(dp) :: t, worst_amp, worst_shape, expect
    integer  :: n
    worst_amp = 0.0_dp; worst_shape = 0.0_dp
    write(*,'(A)') ''
    write(*,'(A,ES9.2,A)') '  7  pure drain at rate ', gam, ' /s (fixed shape)'
    do n = 1, 6*w
      t  = n*dt
      fn = fbase * exp(-gam*t)
      call conv_diag_step(n, t, fn, dd)
      if (dd%valid) then
        ! amplitude rate: ||df||/(dt_win ||f_now||) with f_now the SMALLER
        ! of the pair, so it reads (exp(gam*dt_win)-1)/dt_win
        expect      = (exp(gam*w*dt) - 1.0_dp) / (w*dt)
        worst_amp   = max(worst_amp, abs(dd%eps - expect)/expect)
        worst_shape = max(worst_shape, dd%eps_shape)
      end if
    end do
    if (worst_amp <= 1.0d-6) then
      write(*,'(A,ES9.2,A)') '      PASS  eps sits at the drain rate (rel.err ', &
                             worst_amp, ')'
    else
      write(*,'(A,ES9.2,A)') '      FAIL  eps not at the drain rate (rel.err ', &
                             worst_amp, ')'
      nf = nf + 1
    end if
    ! eps_shape should be at round-off, i.e. utterly negligible next to
    ! the amplitude rate it has to see through.
    if (worst_shape <= 1.0d-10 * gam) then
      write(*,'(A,ES9.2,A)') '      PASS  eps_shape ~ 0 despite the drain (max ', &
                             worst_shape, ' /s)'
    else
      write(*,'(A,ES9.2,A)') '      FAIL  eps_shape polluted by the drain (max ', &
                             worst_shape, ' /s)'
      nf = nf + 1
    end if
  end subroutine drain_run

  subroutine check(label, err, tol, nf)
    character(len=*), intent(in)    :: label
    real(dp),         intent(in)    :: err, tol
    integer,          intent(inout) :: nf
    write(*,'(A)') ''
    if (err <= tol) then
      write(*,'(A,A,A,ES9.2,A,ES9.2,A)') '  PASS  ', label, &
        '   (rel.err ', err, ' <= ', tol, ')'
    else
      write(*,'(A,A,A,ES9.2,A,ES9.2,A)') '  FAIL  ', label, &
        '   (rel.err ', err, ' >  ', tol, ')'
      nf = nf + 1
    end if
  end subroutine check

end program test_conv_diag
