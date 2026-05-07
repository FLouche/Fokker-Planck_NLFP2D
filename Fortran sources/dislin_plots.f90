!   dislin_plots.f90
!
!   End-of-run and time-trace visualisation via DISLIN 11.5 (real64)
!   PNG output; files read from outfile() names already written to disk.
!
!   plot_endof_run()    -- always called; reads fout*, Ekin_perp* files
!   plot_time_traces()  -- called when ntimes/=0 and iplot_traces==-1

module mod_dislin_plots

  use dislin
  use shared_grid    ! nperp, npar, vperp, vpar
  use shared_timer   ! outfile, isc, iplot_traces
  use shared_plasma  ! nbulk
  use shared_beam    ! isource
  use shared_rf      ! irf

  implicit none
  private
  public :: plot_endof_run, plot_time_traces

  integer, parameter :: dp = kind(1.0d0)

contains

  !-----------------------------------------------------------------
  ! Count non-empty lines in a text file; 0 if file absent/unreadable
  !-----------------------------------------------------------------
  integer function count_lines(fname)
    character(len=*), intent(in) :: fname
    integer :: io, lun
    logical :: exists
    count_lines = 0
    inquire(file=trim(fname), exist=exists)
    if (.not. exists) return
    open(newunit=lun, file=trim(fname), status='old', action='read', iostat=io)
    if (io /= 0) return
    do
      read(lun, *, iostat=io)
      if (io /= 0) exit
      count_lines = count_lines + 1
    end do
    close(lun)
  end function count_lines

  !-----------------------------------------------------------------
  ! 2-column file -> single 1D PNG plot
  !-----------------------------------------------------------------
  subroutine plot_1d(fname, xlabel, ylabel, title_str, pngname)
    character(len=*), intent(in) :: fname, xlabel, ylabel, title_str, pngname
    real(dp), allocatable :: xd(:), yd(:)
    real(dp) :: xmin, xmax, ymin, ymax, dx, dy
    integer :: n, i, io, lun

    n = count_lines(fname)
    if (n < 2) return

    allocate(xd(n), yd(n))
    open(newunit=lun, file=trim(fname), status='old', action='read')
    do i = 1, n
      read(lun, *, iostat=io) xd(i), yd(i)
      if (io /= 0) then; n = i - 1; exit; end if
    end do
    close(lun)
    if (n < 2) then; deallocate(xd, yd); return; end if

    xmin = minval(xd(1:n)); xmax = maxval(xd(1:n))
    ymin = minval(yd(1:n)); ymax = maxval(yd(1:n))
    dx = (xmax - xmin) / 5.0_dp; if (dx == 0.0_dp) dx = 1.0_dp
    dy = (ymax - ymin) / 5.0_dp; if (dy == 0.0_dp) dy = 1.0_dp

    call metafl('PNG')
    call setfil(trim(pngname))
    call disini()
    call pagera()
    call hwfont()
    call axspos(450, 1800)
    call axslen(2200, 1200)
    call name(trim(xlabel), 'X')
    call name(trim(ylabel), 'Y')
    call titlin(trim(title_str), 1)
    call graf(xmin, xmax, xmin, dx, ymin, ymax, ymin, dy)
    call title()
    call curve(xd(1:n), yd(1:n), n)
    call disfin()

    deallocate(xd, yd)
  end subroutine plot_1d

  !-----------------------------------------------------------------
  ! 3-column file (outer loop=vperp, inner loop=vpar) ->
  ! filled contour PNG with v_par on x-axis, v_perp on y-axis
  !-----------------------------------------------------------------
  subroutine plot_2d(fname, title_str, pngname)
    character(len=*), intent(in) :: fname, title_str, pngname
    real(dp), allocatable :: zmat(:,:)
    real(dp) :: zlev(20)
    real(dp) :: xmin, xmax, ymin, ymax, dx, dy
    real(dp) :: zmin, zmax, dz, dum1, dum2, val
    integer :: i, j, io, nc, lun
    logical :: exists

    inquire(file=trim(fname), exist=exists)
    if (.not. exists) return

    allocate(zmat(npar, nperp))   ! ZMAT(x_idx, y_idx) = ZMAT(vpar_j, vperp_i)
    open(newunit=lun, file=trim(fname), status='old', action='read')
    do i = 1, nperp
      do j = 1, npar
        read(lun, *, iostat=io) dum1, dum2, val
        if (io /= 0) then
          close(lun); deallocate(zmat); return
        end if
        zmat(j, i) = val
      end do
    end do
    close(lun)

    xmin = vpar(1);  xmax = vpar(npar)
    ymin = vperp(1); ymax = vperp(nperp)
    dx = (xmax - xmin) / 5.0_dp; if (dx == 0.0_dp) dx = 1.0_dp
    dy = (ymax - ymin) / 5.0_dp; if (dy == 0.0_dp) dy = 1.0_dp

    zmin = minval(zmat); zmax = maxval(zmat)
    if (zmin == zmax) zmax = zmin + 1.0_dp
    nc = 20
    dz = (zmax - zmin) / real(nc, dp)
    do i = 1, nc
      zlev(i) = zmin + (i - 0.5_dp) * dz
    end do

    call metafl('PNG')
    call setfil(trim(pngname))
    call disini()
    call pagera()
    call hwfont()
    call setvlt('RAIN')
    call zscale(zmin, zmax)
    call axspos(450, 1800)
    call axslen(2200, 1200)
    call name('v_par (v_th)', 'X')
    call name('v_perp (v_th)', 'Y')
    call titlin(trim(title_str), 1)
    call graf(xmin, xmax, xmin, dx, ymin, ymax, ymin, dy)
    call title()
    call conshd(vpar, npar, vperp, nperp, zmat, zlev, nc)
    call disfin()

    deallocate(zmat)
  end subroutine plot_2d

  !-----------------------------------------------------------------
  ! 5 end-of-run plots (all files must exist on disk)
  !-----------------------------------------------------------------
  subroutine plot_endof_run()
    character(len=256) :: fn, pn

    fn = trim(outfile('fout_at_vpar0.txt'))
    pn = trim(outfile('fout_at_vpar0.png'))
    call plot_1d(fn, 'v_perp (v_th)', 'f', 'VDF at v_par = 0', pn)

    fn = trim(outfile('fout_at_vperp0.txt'))
    pn = trim(outfile('fout_at_vperp0.png'))
    call plot_1d(fn, 'v_par (v_th)', 'f', 'VDF at v_perp = 0', pn)

    fn = trim(outfile('Ekin_perp_at_vpar0.txt'))
    pn = trim(outfile('Ekin_perp_at_vpar0.png'))
    call plot_1d(fn, 'v_perp (v_th)', 'E_kin_perp (keV)', &
                 'Perp. kinetic energy at v_par = 0', pn)

    fn = trim(outfile('fout.txt'))
    pn = trim(outfile('fout.png'))
    call plot_2d(fn, '2D VDF f(v_perp, v_par)', pn)

    fn = trim(outfile('Ekin_perp.txt'))
    pn = trim(outfile('Ekin_perp.png'))
    call plot_2d(fn, 'Perp. kinetic energy (keV)', pn)

  end subroutine plot_endof_run

  !-----------------------------------------------------------------
  ! 3 time-trace plots: temperature, collisional power, power balance
  !-----------------------------------------------------------------
  subroutine plot_time_traces()
    real(dp), allocatable :: t(:), y1(:), y2(:)
    real(dp), allocatable :: pcoll_e(:), pcoll_ions(:,:), pcoll_self(:)
    real(dp), allocatable :: pcoll_tot(:), pRF(:), psource(:), plosses(:)
    real(dp) :: xmin, xmax, ymin, ymax, dx, dy, tmp, tmp2
    character(len=256) :: fn, pn, dynfn
    character(len=512) :: legbuf
    character(len=8)   :: ion_label
    character(len=2)   :: ibstr
    integer :: nt, n, ib, io, i, lun, n_leg

    !-- Plot 1: kinetic energy vs time --------------------------------
    fn = trim(outfile('energy_vs_time.txt'))
    nt = count_lines(fn)
    if (nt >= 2) then
      allocate(t(nt), y1(nt), y2(nt))
      open(newunit=lun, file=trim(fn), status='old', action='read')
      do i = 1, nt
        read(lun, *, iostat=io) t(i), y1(i), y2(i)
        if (io /= 0) then; nt = i - 1; exit; end if
      end do
      close(lun)

      if (nt >= 2) then
        xmin = t(1); xmax = t(nt)
        ymin = min(minval(y1(1:nt)), minval(y2(1:nt)))
        ymax = max(maxval(y1(1:nt)), maxval(y2(1:nt)))
        dx = (xmax - xmin) / 5.0_dp; if (dx == 0.0_dp) dx = 1.0_dp
        dy = (ymax - ymin) / 5.0_dp; if (dy == 0.0_dp) dy = 1.0_dp

        pn = trim(outfile('temperature_vs_time.png'))
        call metafl('PNG')
        call setfil(trim(pn))
        call disini()
        call pagera()
        call hwfont()
        call axspos(450, 1800)
        call axslen(2200, 1200)
        call name('Time (s)', 'X')
        call name('Energy (keV)', 'Y')
        call titlin('Kinetic energy vs time', 1)
        call legini(legbuf, 2, 8)
        call leglin(legbuf, 'Temp', 1)
        call leglin(legbuf, 'Tperp', 2)
        call graf(xmin, xmax, xmin, dx, ymin, ymax, ymin, dy)
        call title()
        call color('RED')
        call curve(t(1:nt), y1(1:nt), nt)
        call color('BLUE')
        call curve(t(1:nt), y2(1:nt), nt)
        call color('FORE')
        call legend(legbuf, 1)
        call disfin()
      end if
      deallocate(t, y1, y2)
    end if

    !-- Plot 2: collisional power density vs time ---------------------
    fn = trim(outfile('power_coll_e_vs_time.txt'))
    nt = count_lines(fn)
    if (nt >= 2) then
      allocate(t(nt), pcoll_e(nt))
      if (nbulk > 1) allocate(pcoll_ions(nt, nbulk-1))
      if (isc /= 0) allocate(pcoll_self(nt))

      open(newunit=lun, file=trim(fn), status='old', action='read')
      do i = 1, nt
        read(lun, *, iostat=io) t(i), pcoll_e(i)
        if (io /= 0) then; nt = i - 1; exit; end if
      end do
      close(lun)

      if (nbulk > 1) then
        do ib = 1, nbulk - 1
          write(ibstr, '(i2)') ib
          dynfn = 'power_coll_ion'//trim(adjustl(ibstr))//'_vs_time.txt'
          fn = trim(outfile(trim(dynfn)))
          n = min(count_lines(fn), nt)
          if (n > 0) then
            open(newunit=lun, file=trim(fn), status='old', action='read')
            do i = 1, n
              read(lun, *, iostat=io) tmp, pcoll_ions(i, ib)
              if (io /= 0) then; n = i - 1; exit; end if
            end do
            close(lun)
          end if
          if (n < nt) pcoll_ions(n+1:nt, ib) = 0.0_dp
        end do
      end if

      if (isc /= 0) then
        fn = trim(outfile('power_coll_self_vs_time.txt'))
        n = min(count_lines(fn), nt)
        if (n > 0) then
          open(newunit=lun, file=trim(fn), status='old', action='read')
          do i = 1, n
            read(lun, *, iostat=io) tmp, pcoll_self(i)
            if (io /= 0) then; n = i - 1; exit; end if
          end do
          close(lun)
        end if
        if (n < nt) pcoll_self(n+1:nt) = 0.0_dp
      end if

      if (nt >= 2) then
        xmin = t(1); xmax = t(nt)
        ymin = minval(pcoll_e(1:nt)); ymax = maxval(pcoll_e(1:nt))
        if (nbulk > 1) then
          ymin = min(ymin, minval(pcoll_ions(1:nt,:)))
          ymax = max(ymax, maxval(pcoll_ions(1:nt,:)))
        end if
        if (isc /= 0) then
          ymin = min(ymin, minval(pcoll_self(1:nt)))
          ymax = max(ymax, maxval(pcoll_self(1:nt)))
        end if
        dx = (xmax - xmin) / 5.0_dp; if (dx == 0.0_dp) dx = 1.0_dp
        dy = (ymax - ymin) / 5.0_dp; if (dy == 0.0_dp) dy = 1.0_dp

        n_leg = nbulk + merge(1, 0, isc /= 0)

        pn = trim(outfile('power_coll_vs_time.png'))
        call metafl('PNG')
        call setfil(trim(pn))
        call disini()
        call pagera()
        call hwfont()
        call axspos(450, 1800)
        call axslen(2200, 1200)
        call name('Time (s)', 'X')
        call name('Power density (MW/m3)', 'Y')
        call titlin('Collisional power density vs time', 1)
        call legini(legbuf, n_leg, 8)
        call leglin(legbuf, 'e-', 1)
        do ib = 1, nbulk - 1
          write(ion_label, '(a,i0)') 'ion ', ib
          call leglin(legbuf, trim(ion_label), 1 + ib)
        end do
        if (isc /= 0) call leglin(legbuf, 'self', n_leg)
        call graf(xmin, xmax, xmin, dx, ymin, ymax, ymin, dy)
        call title()
        call color('RED')
        call curve(t(1:nt), pcoll_e(1:nt), nt)
        if (nbulk > 1) then
          do ib = 1, nbulk - 1
            select case (ib)
              case (1); call color('BLUE')
              case (2); call color('GREEN')
              case (3); call color('CYAN')
              case default; call color('MAGEN')
            end select
            call curve(t(1:nt), pcoll_ions(1:nt, ib), nt)
          end do
        end if
        if (isc /= 0) then
          call color('YELLO')
          call curve(t(1:nt), pcoll_self(1:nt), nt)
        end if
        call color('FORE')
        call legend(legbuf, 1)
        call disfin()
      end if

      deallocate(t, pcoll_e)
      if (allocated(pcoll_ions)) deallocate(pcoll_ions)
      if (allocated(pcoll_self)) deallocate(pcoll_self)
    end if

    !-- Plot 3: power balance vs time ---------------------------------
    fn = trim(outfile('power_coll_tot_vs_time.txt'))
    nt = count_lines(fn)
    if (nt >= 2) then
      allocate(t(nt), pcoll_tot(nt))
      if (irf    == -1) allocate(pRF(nt))
      if (isource == -1) allocate(psource(nt), plosses(nt))

      open(newunit=lun, file=trim(fn), status='old', action='read')
      do i = 1, nt
        read(lun, *, iostat=io) t(i), pcoll_tot(i)
        if (io /= 0) then; nt = i - 1; exit; end if
      end do
      close(lun)

      if (irf == -1) then
        fn = trim(outfile('power_RF_vs_time.txt'))
        n = min(count_lines(fn), nt)
        if (n > 0) then
          open(newunit=lun, file=trim(fn), status='old', action='read')
          do i = 1, n
            read(lun, *, iostat=io) tmp, pRF(i)
            if (io /= 0) then; n = i - 1; exit; end if
          end do
          close(lun)
        end if
        if (n < nt) pRF(n+1:nt) = 0.0_dp
      end if

      if (isource == -1) then
        fn = trim(outfile('power_NBI_vs_time.txt'))
        n = min(count_lines(fn), nt)
        if (n > 0) then
          open(newunit=lun, file=trim(fn), status='old', action='read')
          do i = 1, n
            read(lun, *, iostat=io) tmp, psource(i), plosses(i)
            if (io /= 0) then; n = i - 1; exit; end if
          end do
          close(lun)
        end if
        if (n < nt) then
          psource(n+1:nt) = 0.0_dp; plosses(n+1:nt) = 0.0_dp
        end if
      end if

      if (nt >= 2) then
        xmin = t(1); xmax = t(nt)
        ymin = minval(pcoll_tot(1:nt)); ymax = maxval(pcoll_tot(1:nt))
        if (irf == -1) then
          ymin = min(ymin, minval(pRF(1:nt)))
          ymax = max(ymax, maxval(pRF(1:nt)))
        end if
        if (isource == -1) then
          ymin = min(ymin, minval(psource(1:nt)))
          ymax = max(ymax, maxval(psource(1:nt)))
        end if
        dx = (xmax - xmin) / 5.0_dp; if (dx == 0.0_dp) dx = 1.0_dp
        dy = (ymax - ymin) / 5.0_dp; if (dy == 0.0_dp) dy = 1.0_dp

        n_leg = 1 + merge(1, 0, irf == -1) + merge(1, 0, isource == -1)

        pn = trim(outfile('power_balance_vs_time.png'))
        call metafl('PNG')
        call setfil(trim(pn))
        call disini()
        call pagera()
        call hwfont()
        call axspos(450, 1800)
        call axslen(2200, 1200)
        call name('Time (s)', 'X')
        call name('Power density (MW/m3)', 'Y')
        call titlin('Power balance vs time', 1)
        call legini(legbuf, n_leg, 8)
        call leglin(legbuf, 'coll', 1)
        if (irf == -1 .and. isource == -1) then
          call leglin(legbuf, 'RF',   2)
          call leglin(legbuf, 'beam', 3)
        else if (irf == -1) then
          call leglin(legbuf, 'RF',   2)
        else if (isource == -1) then
          call leglin(legbuf, 'beam', 2)
        end if
        call graf(xmin, xmax, xmin, dx, ymin, ymax, ymin, dy)
        call title()
        call color('RED')
        call curve(t(1:nt), pcoll_tot(1:nt), nt)
        if (irf == -1) then
          call color('BLUE')
          call curve(t(1:nt), pRF(1:nt), nt)
        end if
        if (isource == -1) then
          call color('GREEN')
          call curve(t(1:nt), psource(1:nt), nt)
        end if
        call color('FORE')
        call legend(legbuf, 1)
        call disfin()
      end if

      deallocate(t, pcoll_tot)
      if (allocated(pRF))     deallocate(pRF)
      if (allocated(psource)) deallocate(psource)
      if (allocated(plosses)) deallocate(plosses)
    end if

  end subroutine plot_time_traces

end module mod_dislin_plots
