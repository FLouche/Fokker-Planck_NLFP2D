!***********************************************************************
!*  mod_timefp_ms -- coupled time evolution of all species              *
!*  (stage 3 of docs/Multispecies_design.pdf; namelist i_coupled = -1)  *
!*                                                                     *
!*  Every species s (s = 1 heated, then electrons, then bulk ions) has  *
!*  its own grid (mod_species) and evolves under                        *
!*                                                                     *
!*     df_s/dt = sum_j C(f_s, f_j) + [Q_RF + S - f/tau_p]_(s=1) + L_E f_s *
!*                                                                     *
!*  Coupling (design decision 1): the potentials of every species are   *
!*  lagged by one step (the iteration scheme of Louche 2004, extended   *
!*  to all pairs).  With n_picard > 0, steps 1-5 are repeated with     *
!*  f* = theta f^(n+1,k) + (1-theta) f^n in place of f^n, so that each   *
!*  pair is evaluated at the time level of the theta scheme and the    *
!*  energy exchange is antisymmetric at convergence.  At step n:        *
!*    1. moments, temperatures and Coulomb logarithms from f^n;         *
!*    2. Legendre potentials of every f_j^n on its own grid (lmax_leg   *
!*       for the heated species, lmax_bulk for the bulk species; L >= 1 *
!*       is needed for momentum conservation);                          *
!*    3. for each target s, the coefficients of every pair (s <- j)     *
!*       on the grid of s: pair_coefs (mod_legendre_pot), except the    *
!*       heated species' self-collisions, which keep main_nlterm (FFT   *
!*       kernel, v3.2);                                                 *
!*    4. one theta-scheme step per species, each with its own CSR       *
!*       pattern and PARDISO handle (block-diagonal system, never       *
!*       assembled);                                                    *
!*    5. (after the Picard iterations) the rates of every pair on       *
!*       f^theta = theta f^(n+1) + (1-theta) f^n, the conservation       *
!*       defects (namelist i_ms_cons), and their removal by a moment    *
!*       projection onto the implicit responses to f, div(v f) and      *
!*       df/dv_par: particle numbers, total energy and total momentum   *
!*       are then conserved to round-off.                               *
!*                                                                     *
!*  Quadrature of the coupled loop: q(i,j) = w_perp(i) w_par(j) 2 pi    *
!*  v_perp(i) with trapezoidal weights (the natural weights of the      *
!*  stencil's discrete adjoint; Simpson's alternate).  The              *
!*  axis row (Neumann constraint, sum_k w_k f(k,j) = 0) is eliminated:  *
!*  f(1,j) = sum_{k>=2} a_k f(k,j), so moments and rates are sums over   *
!*  the interior unknowns with folded weights.                          *
!*                                                                     *
!*  Diagnostics (every step, files ms_*_vs_time.txt): moments of every  *
!*  species, the raw particle, power and momentum rates of every pair,  *
!*  the NRL Maxwellian prediction of the power, the Coulomb             *
!*  logarithms, the RF and confinement powers and the corrections.      *
!*                                                                     *
!*  23/11/2026                                                          *
!***********************************************************************

module mod_timefp_ms

  use pardiso_solver
  implicit none
  private
  public :: timefp_ms

  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: pmass = 1.6726d-27      ! proton mass [kg]
  real(dp), parameter :: qe    = 1.60218d-19     ! J per eV

  ! Wall-clock profile of the coupled step (printed at the end of the run).
  integer, parameter :: nt_prof = 10
  real(dp)   :: tprof(nt_prof) = 0.d0
  integer(8) :: tclk(nt_prof) = 0
  character(len=24), parameter :: tname(nt_prof) = [character(len=24) :: &
       'moments', 'Legendre build', 'main_nlterm (a<-a)', 'pair_coefs (cross)', &
       'theta step (solve)', 'projection responses', 'pair rates on f^theta', &
       'RF / L_E rates', 'select_species', 'total']

  ! Per-species workspace: operator, LHS, PARDISO handle, coefficients,
  ! quadrature weights and the axis elimination.
  type ms_ws_t
    integer :: n = 0, nnz = 0
    integer,  allocatable :: ia(:), ja(:), ialhs(:), jalhs(:)
    real(dp), allocatable :: aL(:), alhs(:), rhs(:), x(:), xprev(:), Lf(:), f(:), fn(:), vals(:), colsum(:)
    real(dp), allocatable :: fs(:,:)         ! f* at which the coefficients are evaluated
    real(dp), allocatable, dimension(:,:) :: c20, c02, c11, c10, c01, c00
    real(dp), allocatable :: pc(:,:,:,:)     ! coefficients of every pair (s <- j): (i,j,1:6,j)
    real(dp), allocatable :: g(:,:), fth(:)  ! projection: responses to the unit operators, f^theta
    real(dp), allocatable :: uval(:,:)       ! CSR values of the dilation and shift (constant)
    real(dp) :: cproj(3) = 0.d0              ! projection amplitudes [1/s, 1/s, m/s^2]
    real(dp), allocatable :: q(:,:)          ! trapezoidal weights times 2 pi v_perp
    real(dp) :: ax(7) = 0.d0                 ! f(1,j) = sum_{k=2..7} ax(k) f(k,j)
    type(pardiso_handle_t) :: h
  end type ms_ws_t

contains

subroutine tic(k)
  integer, intent(in) :: k
  call system_clock(tclk(k))
end subroutine tic

subroutine toc(k)
  integer, intent(in) :: k
  integer(8) :: c, r
  call system_clock(c, r)
  tprof(k) = tprof(k) + dble(c - tclk(k)) / dble(r)
end subroutine toc

subroutine sel(s)
  use mod_species, only: select_species
  integer, intent(in) :: s
  call tic(9)
  call select_species(s)
  call toc(9)
end subroutine sel

!=======================================================================
subroutine timefp_ms(fstart_a, fout_a, otime)

  use mod_species
  use mod_legendre_pot
  use mod_fd_stencil_2d
  use mod_ncint
  use shared_grid
  use shared_plasma
  use shared_FPterms
  use shared_timer
  use shared_beam
  use shared_RF
  use shared_multi
  use nlterm
  use assemble_FP_lin, only: add_lossE_terms
  use mod_timefp_7pt_nl, only: phi_kernel_prepare
  use mod_phi_kernel, only: phi_kern, phi_khat

  real(dp), intent(inout) :: fstart_a(:)          ! initial f of the heated species (vector)
  real(dp), intent(out)   :: fout_a(:,:)          ! final f of the heated species
  real(dp), intent(in)    :: otime

  type(ms_ws_t),  allocatable :: ws(:)
  type(legpot_t), allocatable :: pot(:)
  real(dp), allocatable :: dens(:), wkin(:), wperp(:), ppar(:), tev(:), tperp(:), tpar(:)
  real(dp), allocatable :: lnl(:,:), ppair(:,:), rpair(:,:), npair(:,:), pnrl(:,:)
  real(dp), allocatable :: pE(:), dNc(:), dPc(:), dRc(:)
  real(dp), allocatable, dimension(:,:) :: t20, t02, t11, t10, t01, t00, zero2
  real(dp) :: theta, time, dt, taum_save, pRF, mu, kap, wtot, ptot, rdum
  integer  :: s, j, ns, iphase, itime, itime_global, error, k, l, ip, nk, nsys, row, col, ix, info
  integer  :: ipiv(60)
  real(dp) :: amat(60,60), bvec(60), vsx
  real(dp), allocatable :: uu(:,:,:), rr(:,:,:,:), wgt(:,:,:)
  logical,  allocatable :: active(:)
  real(dp) :: pres, pres_s
  logical  :: sc_skip

  ns = nspec
  if (icn == -1) then
    theta = 0.5_dp
  else if (icn == 1) then
    theta = 0.75_dp
  else
    theta = 1.0_dp
  end if

  allocate(ws(ns), pot(ns), dens(ns), wkin(ns), wperp(ns), ppar(ns), tev(ns), tperp(ns), &
           tpar(ns), lnl(ns,ns), ppair(ns,ns), rpair(ns,ns), npair(ns,ns), pnrl(ns,ns), pE(ns), &
           dNc(ns), dPc(ns), dRc(ns), uu(3,3,ns), rr(3,0:3,ns,ns), wgt(3,ns,ns), active(ns))
  if (3 * nspec > 60) stop 'timefp_ms: too many species for the projection system'
  uu = 0.d0;  rr = 0.d0;  wgt = 0.d0;  active = .false.
  dNc = 0.d0;  dPc = 0.d0;  dRc = 0.d0

  !--- initial distribution of the heated species: from main (istart), or a
  !    Maxwellian at ta_init
  call sel(1)
  if (ta_init > 0.d0) then
    call maxwellian_on_grid(spec(1)%aa, spec(1)%dens, ta_init, spec(1)%f)
  else
    spec(1)%f = reshape_vec(fstart_a)
  end if

  !--- per-species workspace, pattern and PARDISO symbolic phase ----------
  taum_save = taum
  taum = 0.d0
  do s = 1, ns
    call sel(s)
    call ws_setup(ws(s), s, ns, theta, timestep(1))
  end do
  taum = taum_save

  !--- phi kernel of the heated species' grid (main_nlterm) -------------
  call sel(1)
  call phi_kernel_prepare

  call open_outputs()
  write(*,'(A,I0,A,I0)') '  coupled loop: conservative construction i_ms_cons = ', i_ms_cons, &
                         ',  Picard iterations n_picard = ', n_picard

  itime_global = 0
  time = otime
  phase_loop: do iphase = 1, 3
    if (ntimes(iphase) == 0) cycle phase_loop
    dt = timestep(iphase)
    timestep_cur = dt
    if (iphase > 1) write(*,'(A,I0,A,ES12.4,A)') '  Phase ', iphase, ': dt = ', dt, ' s'

    time_loop: do itime = 1, ntimes(iphase)
      itime_global = itime_global + 1
      call tic(10)

      ! f* (at which potentials, coefficients and corrections are evaluated):
      ! f^n, then theta f^(n+1,k) + (1-theta) f^n for n_picard iterations.
      do s = 1, ns
        call sel(s)
        ws(s)%fs = spec(s)%f
        ws(s)%fn = vec(spec(s)%f)
      end do
      pres = 0.d0
      picard: do ip = 0, n_picard
        !--- 1. moments and Coulomb logarithms of f^n ----------------------
        call tic(1); call all_moments(); call toc(1)
        do s = 1, ns
          do j = 1, ns
            lnl(s,j) = coulomb_log_pair(s, j, tev(s), tev(j))
          end do
        end do
        if (ip == 0) call write_moments(time)   ! state at t^n (first line: the initial state)

        !--- 2. Legendre potentials of every species -------------------------
        do j = 1, ns
          call sel(j)
          call tic(2)
          call legpot_build(pot(j), ws(j)%fs, merge(lmax_leg, lmax_bulk, spec(j)%kind == 0), nr_leg, nxi_leg)
          call toc(2)
        end do

        !--- 3. coefficients of every pair (s <- j), stored per pair --------
        taum_save = taum
        taum = 0.d0                    ! operators below are collisional only
        do s = 1, ns
          call sel(s)
          allocate(t20(nperp,npar), t02(nperp,npar), t11(nperp,npar), t10(nperp,npar), &
                   t01(nperp,npar), t00(nperp,npar))
          ws(s)%f = vec(ws(s)%fs)
          do j = 1, ns
            if (s == 1 .and. j == 1) then
              ! heated species' self-collisions: v3.2 kernel route (same skip as
              ! timefp_7pt_nl while a beam-filled distribution is still empty)
              sc_skip = (isource == -1 .and. dens(1) < 0.05d0 * npart)
              if (sc_skip) then
                t20 = 0.d0; t02 = 0.d0; t11 = 0.d0; t10 = 0.d0; t01 = 0.d0; t00 = 0.d0
              else
                call tic(3)
                call main_nlterm(ws(1)%f, tev(1) * 1.d-3, t00, t10, t01, t20, t11, t02)
                call toc(3)
              end if
            else
              mu  = spec(s)%aa / spec(j)%aa
              kap = pair_coef_factor(spec(s)%za, spec(s)%aa, spec(j)%za, lnl(s,j))
              call tic(4)
              call pair_coefs(pot(j), mu, kap, t20, t02, t11, t10, t01, t00)
              call toc(4)
            end if
            ws(s)%pc(:,:,1,j) = t20;  ws(s)%pc(:,:,2,j) = t02;  ws(s)%pc(:,:,3,j) = t11
            ws(s)%pc(:,:,4,j) = t10;  ws(s)%pc(:,:,5,j) = t01;  ws(s)%pc(:,:,6,j) = t00
          end do
          ws(s)%c20 = sum(ws(s)%pc(:,:,1,:), dim=3);  ws(s)%c02 = sum(ws(s)%pc(:,:,2,:), dim=3)
          ws(s)%c11 = sum(ws(s)%pc(:,:,3,:), dim=3);  ws(s)%c10 = sum(ws(s)%pc(:,:,4,:), dim=3)
          ws(s)%c01 = sum(ws(s)%pc(:,:,5,:), dim=3);  ws(s)%c00 = sum(ws(s)%pc(:,:,6,:), dim=3)
          ! RF (heated species only)
          if (s == 1 .and. irf == -1) then
            ws(s)%c20 = ws(s)%c20 + rf20;  ws(s)%c02 = ws(s)%c02 + rf02
            ws(s)%c11 = ws(s)%c11 + rf11;  ws(s)%c10 = ws(s)%c10 + rf10
            ws(s)%c01 = ws(s)%c01 + rf01
          end if
          ! energy-confinement loss (every species, only with a source)
          if (tau_E > 0.d0) call add_lossE_terms(ws(s)%c00, ws(s)%c10, ws(s)%c01)
          deallocate(t20, t02, t11, t10, t01, t00)
        end do
        taum = taum_save

        !--- 4. one theta step per species ----------------------------------
        pres = 0.d0
        do s = 1, ns
          call sel(s)
          if (s /= 1) taum = 0.d0   ! the particle loss -f/tau_p acts on the heated species only
          ws(s)%xprev = ws(s)%x
          call tic(5)
          call ws_step(ws(s), s, theta, dt, itime_global)
          call toc(5)
          taum = taum_save
          if (ip > 0) then
            pres_s = maxval(abs(ws(s)%x - ws(s)%xprev)) / max(maxval(abs(ws(s)%x)), tiny(1.d0))
            pres = max(pres, pres_s)
          end if
          ws(s)%fs = theta * reshape_vec(ws(s)%x) + (1.d0 - theta) * spec(s)%f
        end do
      end do picard

      !--- 5. rates of every pair on f^theta, conservation defects, projection
      ! The step applied each pair operator to f_s^theta = theta f^(n+1) +
      ! (1-theta) f^n.  Conservation violations V(s,X), X = N, W, P_par
      ! (i_ms_cons >= 1: N; = 2: also W and P_par):
      !   V(s,N) = sum_j N(s<-j)
      !   V(s,W) = P(s<-s) + sum_{j/=s} w_sj [P(s<-j) + P(j<-s)]
      !   V(s,P) = R(s<-s) + sum_{j/=s} u_sj [R(s<-j) + R(j<-s)]
      ! with w_sj = n_s/(n_s+n_j) (equal temperature corrections) and
      ! u_sj = rho_s/(rho_s+rho_j) (equal velocity corrections; an equal split
      ! would drive the stiff e-i friction).  They are removed IMPLICITLY:
      ! f_s^(n+1) += sum_k c_sk g_sk, g_sk = dt (I - theta dt L_s)^-1 U_k f_s^theta,
      ! U = (1, (1/2) div(v .), d/dv_par), with the amplitudes c chosen so that
      !   sum_k c_sk rate_X(U_k f_s^theta) + V(s,X)(f^theta + theta sum c g) = 0,
      ! i.e. the conditions hold for the corrected state.  Evaluating V on the
      ! uncorrected f^theta (an explicit correction) is unstable: with lagged
      ! potentials the self-pair defect of the electrons is a stiff
      ! nu_ee (u^theta - u*) and its explicit removal amplifies.
      taum_save = taum
      taum = 0.d0
      pnrl = 0.d0;  pE = 0.d0;  pRF = 0.d0
      nk = merge(3, 1, i_ms_cons == 2)
      do s = 1, ns
        call sel(s)
        allocate(zero2(nperp,npar), t00(nperp,npar), t10(nperp,npar), t01(nperp,npar))
        zero2 = 0.d0
        ws(s)%fth = theta * ws(s)%x + (1.d0 - theta) * ws(s)%fn
        ws(s)%f = ws(s)%fth
        active(s) = (i_ms_cons >= 1 .and. dens(s) > 1.d-6 * spec(s)%dens)
        call tic(6)
        if (active(s)) call responses(ws(s), s, theta, dt, nk, uu(:,:,s))
        call toc(6)
        do j = 1, ns
          call tic(7)
          call pair_rates(ws(s), s, ws(s)%pc(:,:,:,j), merge(nk, 0, active(s)), rr(:,0:,j,s))
          call toc(7)
          npair(s,j) = rr(1,0,j,s);  ppair(s,j) = rr(2,0,j,s);  rpair(s,j) = rr(3,0,j,s)
          if (s /= j) pnrl(s,j) = p_nrl(s, j, tev(s), tev(j), lnl(s,j))
        end do
        call tic(8)
        if (s == 1 .and. irf == -1) &
          call op_rates(ws(s), s, rf20, rf02, rf11, rf10, rf01, zero2, rdum, pRF, rdum)
        if (tau_E > 0.d0) then
          t00 = 0.d0; t10 = 0.d0; t01 = 0.d0
          call add_lossE_terms(t00, t10, t01)
          call op_rates(ws(s), s, zero2, zero2, zero2, t10, t01, t00, rdum, pE(s), rdum)
        end if
        call toc(8)
        deallocate(zero2, t00, t10, t01)
      end do
      dNc = 0.d0;  dPc = 0.d0;  dRc = 0.d0
      if (i_ms_cons >= 1) then
        ! weights of each pair defect (X = 1 N: none; 2 W: n; 3 P: rho)
        do s = 1, ns
          do j = 1, ns
            if (j == s) then
              wgt(:,s,j) = 1.d0
            else
              wgt(1,s,j) = 0.d0
              wgt(2,s,j) = dens(s) / (dens(s) + dens(j))
              wgt(3,s,j) = dens(s) * spec(s)%aa / (dens(s) * spec(s)%aa + dens(j) * spec(j)%aa)
            end if
          end do
        end do
        ! linear system, unknowns (s,k), rows (s,X)
        nsys = ns * nk
        amat = 0.d0;  bvec = 0.d0
        do s = 1, ns
          do ix = 1, nk
            row = (s - 1) * nk + ix
            ! V(s,X) on f^theta
            vsx = 0.d0
            do j = 1, ns
              if (ix == 1) then
                vsx = vsx + rr(1,0,j,s)                            ! particles: every pair
              else
                vsx = vsx + wgt(ix,s,j) * rr(ix,0,j,s)
                if (j /= s) vsx = vsx + wgt(ix,s,j) * rr(ix,0,s,j)
              end if
            end do
            if (ix == 1) dNc(s) = -vsx
            if (ix == 2) dPc(s) = -vsx
            if (ix == 3) dRc(s) = -vsx
            bvec(row) = -vsx
            if (.not. active(s)) then
              amat(row,row) = 1.d0;  bvec(row) = 0.d0
              cycle
            end if
            do k = 1, nk
              col = (s - 1) * nk + k
              amat(row,col) = amat(row,col) + uu(ix,k,s)
              do j = 1, ns                                   ! own response through every pair (s<-j)
                if (ix == 1) then
                  amat(row,col) = amat(row,col) + theta * rr(1,k,j,s)
                else
                  amat(row,col) = amat(row,col) + theta * wgt(ix,s,j) * rr(ix,k,j,s)
                end if
              end do
            end do
            if (ix >= 2) then
              do j = 1, ns                                   ! response of partner j through (j<-s)
                if (j == s .or. .not. active(j)) cycle
                do k = 1, nk
                  col = (j - 1) * nk + k
                  amat(row,col) = amat(row,col) + theta * wgt(ix,s,j) * rr(ix,k,s,j)
                end do
              end do
            end if
          end do
        end do
        call dgesv(nsys, 1, amat, size(amat,1), ipiv, bvec, size(bvec), info)
        if (info /= 0) then
          write(*,'(A,I0)') '  *** coupled loop: projection system singular, info = ', info
          stop
        end if
        do s = 1, ns
          ws(s)%cproj = 0.d0
          if (.not. active(s)) cycle
          do k = 1, nk
            ws(s)%cproj(k) = bvec((s - 1) * nk + k)
            ws(s)%x = ws(s)%x + ws(s)%cproj(k) * ws(s)%g(:,k)
          end do
        end do
      end if
      taum = taum_save
      call write_rates(time)
      do s = 1, ns
        call sel(s)
        spec(s)%f = reshape_vec(ws(s)%x)
      end do
      if (n_picard > 0 .and. mod(itime_global, 50) == 0) write(*,'(A,I0,A,ES10.3)') &
        '  step ', itime_global, ': last Picard change of f^(n+1) (relative) = ', pres

      call toc(10)
      time = time + dt
    end do time_loop
  end do phase_loop

  ! final state
  do s = 1, ns
    ws(s)%fs = spec(s)%f
  end do
  call all_moments()
  call write_moments(time)
  call close_outputs()

  !--- profile -------------------------------------------------------------
  if (itime_global > 0) then
    write(*,'(/,A,I0,A)') ' COUPLED LOOP PROFILE (wall clock per step, ', itime_global, ' steps)'
    do k = 1, nt_prof
      write(*,'(4X,A24,F10.4,A,F6.1,A)') tname(k), tprof(k) / itime_global, ' s  (', &
           100.d0 * tprof(k) / max(tprof(nt_prof), tiny(1.d0)), ' %)'
    end do
  end if

  !--- summary -------------------------------------------------------------
  write(*,'(/,A)') ' COUPLED MULTISPECIES RUN (i_coupled = -1): final state'
  write(*,'(A)')   '    s  tag        n [m^-3]       Teff [eV]     T_perp [eV]    T_par [eV]'
  do s = 1, ns
    write(*,'(I5,2X,A4,4ES15.6)') s, spec(s)%tag, dens(s), tev(s), tperp(s), tpar(s)
  end do

  !--- release and write the final distributions ---------------------------
  do s = 1, ns
    call pardiso_solve_finalize(ws(s)%h, ws(s)%ialhs, ws(s)%jalhs, error)
    call legpot_free(pot(s))
  end do
  deallocate(phi_kern)
  if (allocated(phi_khat)) deallocate(phi_khat)

  do s = ns, 1, -1                 ! ends with the heated species selected
    call sel(s)
    open(40, file=trim(outfile('fout.txt')), status='unknown')
    open(42, file=trim(outfile('xout.dat')), status='unknown')
    write(42,*) time, nstep_restart + itime_global
    do k = 1, nperp
      do l = 1, npar
        write(40,*) vperp(k), vpar(l), spec(s)%f(k,l)
        write(42,*) spec(s)%f(k,l)
      end do
    end do
    close(42);  close(40)
  end do
  fout_a = spec(1)%f
  fstart_a = vec(spec(1)%f)

contains

  !---------------------------------------------------------------------
  function reshape_vec(v) result(a)
    real(dp), intent(in) :: v(:)
    real(dp) :: a(nperp,npar)
    integer :: i, jj
    do i = 1, nperp
      do jj = 1, npar
        a(i,jj) = v((i-1)*npar + jj)
      end do
    end do
  end function reshape_vec

  function vec(a) result(v)
    real(dp), intent(in) :: a(:,:)
    real(dp) :: v(nperp*npar)
    integer :: i, jj
    do i = 1, nperp
      do jj = 1, npar
        v((i-1)*npar + jj) = a(i,jj)
      end do
    end do
  end function vec

  ! unit dilation (1/2) div(v f) = (1/2)[3 f + v_perp f_perp + v_par f_par]
  subroutine dilation(a10, a01, a00)
    real(dp), intent(out) :: a10(:,:), a01(:,:), a00(:,:)
    integer :: i, jj
    do jj = 1, npar
      do i = 1, nperp
        a00(i,jj) = 1.5d0;  a10(i,jj) = 0.5d0 * vperp(i);  a01(i,jj) = 0.5d0 * vpar(jj)
      end do
    end do
  end subroutine dilation

  !---------------------------------------------------------------------
  ! Maxwellian at T [eV], normalised to n on the grid now selected.
  ! (drifting along v_par at ua_init [m/s], namelist; 0 by default)
  subroutine maxwellian_on_grid(amass, n0, tev0, fm)
    real(dp), intent(in)  :: amass, n0, tev0
    real(dp), intent(out) :: fm(:,:)
    real(dp) :: vth2, dn
    real(dp), allocatable :: fint(:,:)
    integer :: i, jj
    vth2 = tev0 * qe / (amass * pmass)
    do i = 1, nperp
      do jj = 1, npar
        fm(i,jj) = exp(-(vperp(i)**2 + (vpar(jj) - ua_init)**2) / (2.d0 * vth2))
      end do
    end do
    allocate(fint(nperp,npar))
    fint = fm * jacob
    call ncint_2d(fint, dn)
    fm = fm * n0 / dn
    deallocate(fint)
  end subroutine maxwellian_on_grid

  !---------------------------------------------------------------------
  ! Moments of every species with the weights q of the coupled loop.
  subroutine all_moments()
    integer :: s1
    do s1 = 1, ns
      call sel(s1)
      call moments(ws(s1), s1, ws(s1)%fs, dens(s1), wkin(s1), wperp(s1), ppar(s1))
      if (dens(s1) > 1.d-6 * spec(s1)%dens) then
        tev(s1)   = (2.d0/3.d0) * wkin(s1) / dens(s1) / qe
        tperp(s1) = wperp(s1) / dens(s1) / qe
        tpar(s1)  = 2.d0 * (wkin(s1) - wperp(s1)) / dens(s1) / qe
      else                                      ! empty (beam fill-up): nominal T
        tev(s1) = spec(s1)%t0;  tperp(s1) = tev(s1);  tpar(s1) = tev(s1)
      end if
    end do
  end subroutine all_moments

  !---------------------------------------------------------------------
  ! Coulomb logarithm of the pair (s <- j), NRL formulary; symmetric in
  ! (s, j), so both directions of a pair use the same value.
  real(dp) function coulomb_log_pair(s1, j1, ts, tj)
    use coulomb_log_mod
    integer,  intent(in) :: s1, j1
    real(dp), intent(in) :: ts, tj
    real(dp) :: r
    if (spec(s1)%kind == 1 .and. spec(j1)%kind == 1) then
      r = coulomb_log_ee(spec(s1)%dens, ts)
    else if (spec(s1)%kind == 1) then
      call coulomb_log_ae(spec(s1)%dens, spec(j1)%za, ts, r)
    else if (spec(j1)%kind == 1) then
      call coulomb_log_ae(spec(j1)%dens, spec(s1)%za, tj, r)
    else
      call coulomb_log_ab(spec(s1)%za, spec(s1)%aa, ts, spec(s1)%dens, &
                          spec(j1)%za, spec(j1)%aa, tj, spec(j1)%dens, r)
    end if
    coulomb_log_pair = r
  end function coulomb_log_pair

  ! Electron-electron (NRL formulary 2019, p. 34): n_e [m^-3], T_e [eV]
  real(dp) function coulomb_log_ee(ne_m3, te)
    real(dp), intent(in) :: ne_m3, te
    coulomb_log_ee = 23.5d0 - log(sqrt(ne_m3 * 1.d-6) * te**(-1.25d0)) &
                   - sqrt(1.d-5 + (log(te) - 2.d0)**2 / 16.d0)
  end function coulomb_log_ee

  !---------------------------------------------------------------------
  ! Maxwellian energy exchange (Trubnikov 1965; NRL formulary form):
  !   dT_s/dt = nu (T_j - T_s),
  !   nu = C0 sqrt(m_s m_j) Z_s^2 Z_j^2 n_j lnL / (m_s T_j + m_j T_s)^1.5
  ! m in g, n in cm^-3, T in eV, C0 = (8 sqrt(2 pi)/3) e^4 / eV^1.5 (Gaussian
  ! units) = 1.7544e-19 (the formulary rounds it to 1.8e-19);
  ! P(s<-j) = 3/2 n_s nu (T_j - T_s) [W/m^3].
  real(dp) function p_nrl(s1, j1, ts, tj, lnlsj)
    integer,  intent(in) :: s1, j1
    real(dp), intent(in) :: ts, tj, lnlsj
    real(dp), parameter :: c0 = 1.75436d-19
    real(dp) :: ms, mj, nu
    ms = spec(s1)%aa * pmass * 1.d3;  mj = spec(j1)%aa * pmass * 1.d3
    nu = c0 * sqrt(ms * mj) * spec(s1)%za**2 * spec(j1)%za**2 * spec(j1)%dens * 1.d-6 &
         * lnlsj / (ms * tj + mj * ts)**1.5d0
    p_nrl = 1.5d0 * dens(s1) * nu * (tj - ts) * qe
  end function p_nrl

  !---------------------------------------------------------------------
  subroutine open_outputs()
    integer :: s1, j1
    character(len=600) :: hdr
    open(620, file=trim(outfile('ms_moments_vs_time.txt')),        status='unknown')
    open(621, file=trim(outfile('ms_power_pairs_vs_time.txt')),    status='unknown')
    open(622, file=trim(outfile('ms_power_nrl_vs_time.txt')),      status='unknown')
    open(623, file=trim(outfile('ms_momentum_pairs_vs_time.txt')), status='unknown')
    open(624, file=trim(outfile('ms_density_pairs_vs_time.txt')),  status='unknown')
    open(625, file=trim(outfile('ms_coulomb_log_vs_time.txt')),    status='unknown')
    open(626, file=trim(outfile('ms_power_other_vs_time.txt')),    status='unknown')
    open(627, file=trim(outfile('ms_corrections_vs_time.txt')),    status='unknown')
    hdr = '# time [s]; per species s ('
    do s1 = 1, ns
      hdr = trim(hdr)//' '//trim(stag(s1))
    end do
    write(620,'(A)') trim(hdr)//'): n [m^-3], W [J/m^3], Teff, T_perp, T_par [eV], P_par [kg/m^2/s];'// &
                     ' then W_tot, P_par_tot'
    hdr = '# time [s]; then (s<-j) for s = 1..ns, j = 1..ns:'
    do s1 = 1, ns
      do j1 = 1, ns
        hdr = trim(hdr)//' '//trim(stag(s1))//'<'//trim(stag(j1))
      end do
    end do
    write(621,'(A)') trim(hdr)//'; then the sum  [MW/m^3]  (raw pair rates on f^theta, before the projection)'
    write(622,'(A)') trim(hdr)//'  (Maxwellian exchange at the same T and lnL)  [MW/m^3]'
    write(623,'(A)') trim(hdr)//'; then the sum  [kg m^-2 s^-2]  (raw)'
    write(624,'(A)') trim(hdr)//'  [m^-3 s^-1]'
    write(625,'(A)') trim(hdr)//'  (Coulomb logarithms)'
    write(626,'(A)') '# time [s]; P_RF (heated species); P_E(s), s = 1..ns  [MW/m^3]'
    write(627,'(A)') '# time [s]; per species: defects removed dN [m^-3/s], dP [MW/m^3], dR [kg m^-2 s^-2],'// &
                     ' projection amplitudes c_gamma, c_alpha [1/s], c_beta [m/s^2]'
  end subroutine open_outputs

  subroutine close_outputs()
    integer :: u
    do u = 620, 627
      close(u)
    end do
  end subroutine close_outputs

  character(len=4) function stag(s1)
    integer, intent(in) :: s1
    stag = spec(s1)%tag
    if (s1 == 1) stag = 'a'
  end function stag

  subroutine write_moments(tt)
    real(dp), intent(in) :: tt
    integer :: s1, j1
    wtot = sum(wkin);  ptot = sum(ppar)
    write(620,'(ES16.8,*(ES25.16E3))') tt, (dens(s1), wkin(s1), tev(s1), tperp(s1), tpar(s1), &
                                         ppar(s1), s1 = 1, ns), wtot, ptot
    write(625,'(ES16.8,*(ES14.6))') tt, ((lnl(s1,j1), j1 = 1, ns), s1 = 1, ns)
  end subroutine write_moments

  subroutine write_rates(tt)
    real(dp), intent(in) :: tt
    integer :: s1, j1
    write(621,'(ES16.8,*(ES17.9))') tt, ((ppair(s1,j1)*1.d-6, j1 = 1, ns), s1 = 1, ns), sum(ppair)*1.d-6
    write(622,'(ES16.8,*(ES17.9))') tt, ((pnrl(s1,j1)*1.d-6, j1 = 1, ns), s1 = 1, ns)
    write(623,'(ES16.8,*(ES17.9))') tt, ((rpair(s1,j1), j1 = 1, ns), s1 = 1, ns), sum(rpair)
    write(624,'(ES16.8,*(ES17.9))') tt, ((npair(s1,j1), j1 = 1, ns), s1 = 1, ns)
    write(626,'(ES16.8,*(ES17.9))') tt, pRF*1.d-6, (pE(s1)*1.d-6, s1 = 1, ns)
    write(627,'(ES16.8,*(ES17.9))') tt, (dNc(s1), dPc(s1)*1.d-6, dRc(s1), ws(s1)%cproj, s1 = 1, ns)
  end subroutine write_rates

end subroutine timefp_ms

!=======================================================================
! Workspace of species s (selected): CSR pattern of the full 7-point
! stencil (all six coefficients set to 1, so every entry the step can
! produce is present), LHS copy, quadrature weights, axis elimination,
! PARDISO symbolic phase.  taum must be 0 (set by the caller).
!=======================================================================
subroutine ws_setup(w, s, nsp, theta, dt0)
  use shared_grid
  use mod_fd_stencil_2d
  integer,       intent(in)    :: s, nsp
  type(ms_ws_t), intent(inout) :: w
  real(dp),      intent(in)    :: theta, dt0
  integer  :: col_idx(49), n_entries, i, jj, row, ptr, k, error
  real(dp) :: coeff(49), rhs_ij, wpe(nperp), wpa(npar), w1

  w%n = nperp * npar
  allocate(w%ia(w%n+1), w%ja(49*w%n), w%aL(49*w%n))
  allocate(w%rhs(w%n), w%x(w%n), w%xprev(w%n), w%Lf(w%n), w%f(w%n), w%fn(w%n), w%colsum(w%n))
  allocate(w%fs(nperp,npar), w%g(w%n,3), w%fth(w%n))
  allocate(w%pc(nperp,npar,6,nsp))
  allocate(w%c20(nperp,npar), w%c02(nperp,npar), w%c11(nperp,npar), &
           w%c10(nperp,npar), w%c01(nperp,npar), w%c00(nperp,npar), w%q(nperp,npar))
  w%ia(1) = 1
  ptr = 1
  do i = 1, nperp
    do jj = 1, npar
      row = (i-1)*npar + jj
      call fd_stencil_2d(i, jj, nperp, npar, vperp, dvpar, 1.0_dp, 1.0_dp, 1.0_dp, &
                         1.0_dp, 1.0_dp, 1.0_dp, col_idx, coeff, n_entries, rhs_ij)
      call sort_stencil(col_idx, coeff, n_entries)
      do k = 1, n_entries
        w%ja(ptr) = col_idx(k);  w%aL(ptr) = coeff(k);  ptr = ptr + 1
      end do
      w%ia(row+1) = ptr
    end do
  end do
  w%nnz = w%ia(w%n+1) - 1
  allocate(w%ialhs(w%n+1), w%jalhs(w%nnz), w%alhs(w%nnz), w%vals(w%nnz))
  w%ialhs = w%ia;  w%jalhs = w%ja(1:w%nnz)

  ! trapezoidal weights times 2 pi v_perp
  call trap_weights(vperp, nperp, wpe)
  call trap_weights(vpar,  npar,  wpa)
  do i = 1, nperp
    do jj = 1, npar
      w%q(i,jj) = wpe(i) * wpa(jj) * jacob(i,jj)
    end do
  end do

  ! axis row (Neumann): sum_{k=1..7} c_k f(k,j) = 0, the same for every j
  call fd_stencil_2d(1, 2, nperp, npar, vperp, dvpar, 1.0_dp, 1.0_dp, 1.0_dp, &
                     1.0_dp, 1.0_dp, 1.0_dp, col_idx, coeff, n_entries, rhs_ij)
  w%ax = 0.d0
  w1 = 0.d0
  do k = 1, n_entries
    if ((col_idx(k) - 2) / npar + 1 == 1) w1 = coeff(k)
  end do
  do k = 1, n_entries
    i = (col_idx(k) - 2) / npar + 1
    if (i >= 2 .and. i <= 7) w%ax(i) = -coeff(k) / w1
  end do

  ! constant unit operators of the projection: dilation (1/2) div(v f), shift df/dv_par
  allocate(w%uval(w%nnz,2))
  block
    real(dp), allocatable, dimension(:,:) :: z, u10, u01, u00
    allocate(z(nperp,npar), u10(nperp,npar), u01(nperp,npar), u00(nperp,npar))
    z = 0.d0
    do jj = 1, npar
      do i = 1, nperp
        u00(i,jj) = 1.5d0;  u10(i,jj) = 0.5d0 * vperp(i);  u01(i,jj) = 0.5d0 * vpar(jj)
      end do
    end do
    call fill_vals(w, z, z, z, u10, u01, u00)
    w%uval(:,1) = w%vals(1:w%nnz)
    u01 = 1.d0
    call fill_vals(w, z, z, z, z, u01, z)
    w%uval(:,2) = w%vals(1:w%nnz)
    deallocate(z, u10, u01, u00)
  end block

  call lhs_from_L(w, theta, dt0)
  call pardiso_solve_init(w%h, w%n, w%alhs, w%ialhs, w%jalhs, a_constant=.false., mtype=11, &
                          msglvl=0, error=error, &
                          lu_reuse_L=merge(nint(-log10(lu_reuse_tol)), 0, i_lu_reuse == -1))
  if (error /= 0) then
    write(*,'(A,I0,A,I0)') ' timefp_ms: pardiso_solve_init failed for species ', s, ', error=', error
    stop
  end if
  write(*,'(A,I0,A,I0,A,I0)') '  species ', s, ': ', w%n, ' unknowns, nnz = ', w%nnz
end subroutine ws_setup

subroutine trap_weights(x, n, wt)
  integer,  intent(in)  :: n
  real(dp), intent(in)  :: x(n)
  real(dp), intent(out) :: wt(n)
  integer :: i
  wt(1) = 0.5d0 * (x(2) - x(1))
  wt(n) = 0.5d0 * (x(n) - x(n-1))
  do i = 2, n - 1
    wt(i) = 0.5d0 * (x(i+1) - x(i-1))
  end do
end subroutine trap_weights

!=======================================================================
! Values of the operator with coefficients a** on the pattern of w
! (taum = 0 expected), in w%vals.
!=======================================================================
subroutine fill_vals(w, a20, a02, a11, a10, a01, a00)
  use shared_grid
  use mod_fd_stencil_2d
  type(ms_ws_t), intent(inout) :: w
  real(dp), intent(in), dimension(nperp,npar) :: a20, a02, a11, a10, a01, a00
  integer  :: col_idx(49), n_entries, i, jj, row, k
  real(dp) :: coeff(49), rhs_ij
  logical  :: ok
  ok = .true.
  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) &
  !$OMP   PRIVATE(i, jj, k, row, col_idx, coeff, n_entries, rhs_ij) REDUCTION(.AND.:ok)
  do i = 1, nperp
    do jj = 1, npar
      row = (i-1)*npar + jj
      call fd_stencil_2d(i, jj, nperp, npar, vperp, dvpar, a00(i,jj), a10(i,jj), a01(i,jj), &
                         a20(i,jj), a11(i,jj), a02(i,jj), col_idx, coeff, n_entries, rhs_ij)
      call sort_stencil(col_idx, coeff, n_entries)
      call place_row(w, row, col_idx, coeff, n_entries, w%vals, ok)
    end do
  end do
  !$OMP END PARALLEL DO
  if (.not. ok) then
    write(*,'(A)') '  *** PATTERN MISMATCH in fill_vals (coupled loop)'
    stop
  end if
end subroutine fill_vals

!-----------------------------------------------------------------------
! Folded weights: g(1,j) carried onto the interior unknowns of column j
! through the axis elimination, so that sum_all g f = sum_interior gt f
! for every f obeying the axis constraint (rows i >= 2; gt(1,:) = 0).
subroutine folded_weights(w, g, gt)
  use shared_grid
  type(ms_ws_t), intent(in) :: w
  real(dp), intent(in)  :: g(nperp,npar)
  real(dp), intent(out) :: gt(nperp,npar)
  integer :: k
  gt = g
  gt(1,:) = 0.d0
  do k = 2, 7
    gt(k,:) = gt(k,:) + w%ax(k) * g(1,:)
  end do
end subroutine folded_weights

!=======================================================================
! Rates of an operator applied to f (in w%f) on the interior rows, with
! the folded weights of 1, 1/2 m v^2 and m v_par: dN/dt, dW/dt, dP/dt.
!=======================================================================
subroutine op_rates(w, s, a20, a02, a11, a10, a01, a00, rn, rw, rp)
  use shared_grid
  type(ms_ws_t), intent(inout) :: w
  integer,       intent(in)    :: s
  real(dp), intent(in), dimension(nperp,npar) :: a20, a02, a11, a10, a01, a00
  real(dp), intent(out) :: rn, rw, rp
  call fill_vals(w, a20, a02, a11, a10, a01, a00)
  call csr_matvec(w%n, w%ia, w%ja, w%vals, w%f, w%Lf)
  call lf_rates(w, s, rn, rw, rp)
end subroutine op_rates

!-----------------------------------------------------------------------
! dN/dt, dW/dt, dP/dt of the vector w%Lf (= L f), interior rows, folded
! weights.
subroutine lf_rates(w, s, rn, rw, rp)
  use shared_grid
  use mod_species, only: spec
  type(ms_ws_t), intent(inout) :: w
  integer,       intent(in)    :: s
  real(dp), intent(out) :: rn, rw, rp
  real(dp), allocatable :: g(:,:), gt(:,:), lf2(:,:)
  real(dp) :: m
  integer :: i, jj, row

  allocate(g(nperp,npar), gt(nperp,npar), lf2(nperp,npar))
  lf2 = 0.d0
  do i = 2, nperp - 1
    do jj = 2, npar - 1
      row = (i-1)*npar + jj
      lf2(i,jj) = w%Lf(row)
    end do
  end do
  m = spec(s)%aa * 1.6726d-27
  call folded_weights(w, w%q, gt)
  rn = sum(gt * lf2)
  do i = 1, nperp
    do jj = 1, npar
      g(i,jj) = w%q(i,jj) * 0.5d0 * m * (vperp(i)**2 + vpar(jj)**2)
    end do
  end do
  call folded_weights(w, g, gt)
  rw = sum(gt * lf2)
  do i = 1, nperp
    do jj = 1, npar
      g(i,jj) = w%q(i,jj) * m * vpar(jj)
    end do
  end do
  call folded_weights(w, g, gt)
  rp = sum(gt * lf2)
  deallocate(g, gt, lf2)
end subroutine lf_rates

!=======================================================================
! Implicit responses of species s to the unit operators of the projection:
!   g_k = dt (I - theta dt L_s)^-1 U_k f^theta,  U = (1, (1/2) div(v .), d/dv_par)
! (boundary rows homogeneous, so f + c g still obeys the constraints), and
! the rates (N, W, P_par) of U_k f^theta in uk(:,k).  Uses w%fth and the
! factorised LHS of the step.
!=======================================================================
subroutine responses(w, s, theta, dt, nk, uk)
  use shared_grid
  type(ms_ws_t), intent(inout) :: w
  integer,       intent(in)    :: s, nk
  real(dp),      intent(in)    :: theta, dt
  real(dp),      intent(out)   :: uk(3,3)
  integer :: k, i, jj, error
  real(dp) :: mx(3), lr(3), dum
  real(dp), allocatable :: g2(:,:)

  allocate(g2(nperp,npar))
  uk = 0.d0
  do k = 1, nk
    ! U_1 = identity (interior rows); U_2, U_3 from the CSR values cached at setup
    if (k == 1) then
      w%Lf = w%fth
    else
      call csr_matvec(w%n, w%ia, w%ja, w%uval(:,k-1), w%fth, w%Lf)
    end if
    w%rhs = dt * w%Lf
    do i = 1, nperp
      do jj = 1, npar
        if (i == 1 .or. i == nperp .or. jj == 1 .or. jj == npar) w%rhs((i-1)*npar + jj) = 0.0_dp
      end do
    end do
    ! one back-substitution with the current factors (approximate when LU
    ! reuse has kept older factors): the projection only needs a good
    ! response direction, since it uses what g_k actually does (below)
    call pardiso_solve_fixed(w%h, w%alhs, w%ialhs, w%jalhs, w%rhs, w%g(:,k), error=error)
    ! constraints exactly: Dirichlet walls 0, axis row from the interior
    do jj = 1, npar
      w%g((nperp-1)*npar + jj, k) = 0.d0
      if (jj == 1 .or. jj == npar) then
        do i = 1, nperp
          w%g((i-1)*npar + jj, k) = 0.d0
        end do
      else
        w%g(jj, k) = sum(w%ax(2:7) * w%g((/ (((i-1)*npar + jj), i = 2, 7) /), k))
      end if
    end do
    ! effective unit rate: moment change of g_k per dt, minus the part the
    ! step operator produces on it, = X(U_k f^theta) for an exact solve
    do i = 1, nperp
      do jj = 1, npar
        g2(i,jj) = w%g((i-1)*npar + jj, k)
      end do
    end do
    call moments(w, s, g2, mx(1), mx(2), dum, mx(3))
    call csr_matvec(w%n, w%ia, w%ja, w%aL, w%g(:,k), w%Lf)
    call lf_rates(w, s, lr(1), lr(2), lr(3))
    uk(:,k) = mx / dt - theta * lr
  end do
  deallocate(g2)
end subroutine responses

!-----------------------------------------------------------------------
! Rates (N, W, P_par) of the pair operator with coefficients pc(:,:,1:6)
! applied to f^theta (r(:,0)) and to the responses g_k (r(:,k), k <= nk).
subroutine pair_rates(w, s, pc, nk, r)
  use shared_grid
  type(ms_ws_t), intent(inout) :: w
  integer,       intent(in)    :: s, nk
  real(dp),      intent(in)    :: pc(:,:,:)
  real(dp),      intent(out)   :: r(3,0:3)
  integer :: k
  r = 0.d0
  call fill_vals(w, pc(:,:,1), pc(:,:,2), pc(:,:,3), pc(:,:,4), pc(:,:,5), pc(:,:,6))
  call csr_matvec(w%n, w%ia, w%ja, w%vals, w%fth, w%Lf)
  call lf_rates(w, s, r(1,0), r(2,0), r(3,0))
  do k = 1, nk
    call csr_matvec(w%n, w%ia, w%ja, w%vals, w%g(:,k), w%Lf)
    call lf_rates(w, s, r(1,k), r(2,k), r(3,k))
  end do
end subroutine pair_rates

!-----------------------------------------------------------------------
! Density, kinetic energy, its perpendicular part [J/m^3] and parallel
! momentum [kg m^-2 s^-1] with the weights q (all nodes; f obeys the
! axis constraint, so this equals the folded interior sum).
subroutine moments(w, s, fa, n0, w0, wp0, p0)
  use shared_grid
  use mod_species, only: spec
  type(ms_ws_t), intent(in) :: w
  integer,  intent(in)  :: s
  real(dp), intent(in)  :: fa(:,:)
  real(dp), intent(out) :: n0, w0, wp0, p0
  real(dp) :: m
  integer :: i, jj
  m = spec(s)%aa * 1.6726d-27
  n0 = 0.d0;  w0 = 0.d0;  wp0 = 0.d0;  p0 = 0.d0
  do jj = 1, npar
    do i = 1, nperp
      n0  = n0  + w%q(i,jj) * fa(i,jj)
      w0  = w0  + w%q(i,jj) * fa(i,jj) * 0.5d0 * m * (vperp(i)**2 + vpar(jj)**2)
      wp0 = wp0 + w%q(i,jj) * fa(i,jj) * 0.5d0 * m * vperp(i)**2
      p0  = p0  + w%q(i,jj) * fa(i,jj) * m * vpar(jj)
    end do
  end do
end subroutine moments

!=======================================================================
! One theta-scheme step of species s (selected), coefficients in w%c**:
!   (I - theta dt L) f^{n+1} = (I + (1-theta) dt L) f^n + dt S
! with the boundary rows kept as constraints (as in timefp_7pt_nl).
!=======================================================================
subroutine ws_step(w, s, theta, dt, istep)
  use shared_grid
  use shared_beam, only: source_v
  use mod_fd_stencil_2d
  type(ms_ws_t), intent(inout) :: w
  integer,       intent(in)    :: s, istep
  real(dp),      intent(in)    :: theta, dt
  integer  :: col_idx(49), n_entries, i, jj, row, k, error
  real(dp) :: coeff(49), rhs_ij
  logical  :: pattern_ok

  pattern_ok = .true.
  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) &
  !$OMP   PRIVATE(i, jj, k, row, col_idx, coeff, n_entries, rhs_ij) REDUCTION(.AND.:pattern_ok)
  do i = 1, nperp
    do jj = 1, npar
      row = (i-1)*npar + jj
      call fd_stencil_2d(i, jj, nperp, npar, vperp, dvpar, &
                         w%c00(i,jj), w%c10(i,jj), w%c01(i,jj), &
                         w%c20(i,jj), w%c11(i,jj), w%c02(i,jj), &
                         col_idx, coeff, n_entries, rhs_ij)
      call sort_stencil(col_idx, coeff, n_entries)
      call place_row(w, row, col_idx, coeff, n_entries, w%aL, pattern_ok)
    end do
  end do
  !$OMP END PARALLEL DO
  if (.not. pattern_ok) then
    write(*,'(A,I0,A,I0)') '  *** PATTERN MISMATCH (coupled loop), species ', s, ', step ', istep
    stop
  end if

  call lhs_from_L(w, theta, dt)

  call csr_matvec(w%n, w%ia, w%ja, w%aL, w%fn, w%Lf)
  w%rhs = w%fn + (1.0_dp - theta) * dt * w%Lf
  if (s == 1) w%rhs = w%rhs + dt * source_v
  do i = 1, nperp
    do jj = 1, npar
      if (i == 1 .or. i == nperp .or. jj == 1 .or. jj == npar) w%rhs((i-1)*npar + jj) = 0.0_dp
    end do
  end do

  call pardiso_solve_step(w%h, w%alhs, w%ialhs, w%jalhs, w%rhs, w%x, a_changed=.true., error=error)
  if (error /= 0) then
    write(*,'(A,I0,A,I0)') ' timefp_ms: pardiso_solve_step failed for species ', s, ', error=', error
    stop
  end if
end subroutine ws_step

!-----------------------------------------------------------------------
! LHS = I - theta dt L, the boundary rows (i = 1 Neumann, walls Dirichlet)
! kept as the L stencil (constraints).
subroutine lhs_from_L(w, theta, dt)
  use shared_grid
  type(ms_ws_t), intent(inout) :: w
  real(dp),      intent(in)    :: theta, dt
  integer :: row, ptr, i, jj
  !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(row, ptr, i, jj)
  do row = 1, w%n
    i  = (row - 1) / npar + 1
    jj = row - (i - 1) * npar
    if (i == 1 .or. i == nperp .or. jj == 1 .or. jj == npar) then
      do ptr = w%ia(row), w%ia(row+1) - 1
        w%alhs(ptr) = w%aL(ptr)
      end do
    else
      do ptr = w%ia(row), w%ia(row+1) - 1
        w%alhs(ptr) = -theta * dt * w%aL(ptr)
        if (w%ja(ptr) == row) w%alhs(ptr) = w%alhs(ptr) + 1.0_dp
      end do
    end if
  end do
  !$OMP END PARALLEL DO
end subroutine lhs_from_L

!-----------------------------------------------------------------------
! Put the sorted stencil entries of a row at their positions in the fixed
! pattern (entries absent from the stencil, e.g. a coefficient that is
! exactly zero, are set to 0).  ok = .false. if an entry is not in the
! pattern.
subroutine place_row(w, row, col_idx, coeff, n_entries, vals, ok)
  type(ms_ws_t), intent(in)    :: w
  integer,       intent(in)    :: row, n_entries, col_idx(:)
  real(dp),      intent(in)    :: coeff(:)
  real(dp),      intent(inout) :: vals(:)
  logical,       intent(inout) :: ok
  integer :: ptr, k
  ptr = w%ia(row)
  vals(w%ia(row):w%ia(row+1)-1) = 0.0_dp
  do k = 1, n_entries
    do while (ptr < w%ia(row+1))
      if (w%ja(ptr) == col_idx(k)) exit
      ptr = ptr + 1
    end do
    if (ptr >= w%ia(row+1)) then
      ok = .false.
      return
    end if
    vals(ptr) = coeff(k)
  end do
end subroutine place_row

!-----------------------------------------------------------------------
! 3x3 linear solve (Gaussian elimination, partial pivoting).
subroutine solve3(a_in, b_in, x)
  real(dp), intent(in)  :: a_in(3,3), b_in(3)
  real(dp), intent(out) :: x(3)
  real(dp) :: a(3,3), b(3), t(3), tb, fct
  integer :: k, i, p
  a = a_in;  b = b_in
  do k = 1, 3
    p = k - 1 + maxloc(abs(a(k:3,k)), 1)
    if (p /= k) then
      t = a(k,:);  a(k,:) = a(p,:);  a(p,:) = t
      tb = b(k);   b(k) = b(p);      b(p) = tb
    end if
    do i = k + 1, 3
      fct = a(i,k) / a(k,k)
      a(i,k:3) = a(i,k:3) - fct * a(k,k:3)
      b(i) = b(i) - fct * b(k)
    end do
  end do
  do k = 3, 1, -1
    x(k) = (b(k) - sum(a(k,k+1:3) * x(k+1:3))) / a(k,k)
  end do
end subroutine solve3

subroutine csr_matvec(n, ia, ja, aa, x, y)
  integer,  intent(in)  :: n, ia(n+1), ja(:)
  real(dp), intent(in)  :: aa(:), x(n)
  real(dp), intent(out) :: y(n)
  integer :: row, ptr
  !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(row, ptr)
  do row = 1, n
    y(row) = 0.0_dp
    do ptr = ia(row), ia(row+1) - 1
      y(row) = y(row) + aa(ptr) * x(ja(ptr))
    end do
  end do
  !$OMP END PARALLEL DO
end subroutine csr_matvec

subroutine sort_stencil(idx, val, n)
  integer,  intent(in)    :: n
  integer,  intent(inout) :: idx(:)
  real(dp), intent(inout) :: val(:)
  integer  :: p, q, tmp_i
  real(dp) :: tmp_v
  do p = 2, n
    tmp_i = idx(p);  tmp_v = val(p)
    q = p - 1
    do while (q >= 1)
      if (idx(q) <= tmp_i) exit
      idx(q+1) = idx(q);  val(q+1) = val(q)
      q = q - 1
    end do
    idx(q+1) = tmp_i;  val(q+1) = tmp_v
  end do
end subroutine sort_stencil

end module mod_timefp_ms
