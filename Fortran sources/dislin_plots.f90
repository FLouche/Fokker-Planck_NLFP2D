! dislin_plots.f90 – no-op stub.
! All plotting is now handled by fp2d_plot.py (Python wrapper).
! The module is kept so the Visual Studio project compiles without modification.
module mod_dislin_plots
  implicit none
  private
  public :: plot_endof_run, plot_time_traces
contains
  subroutine plot_endof_run()
  end subroutine plot_endof_run
  subroutine plot_time_traces()
  end subroutine plot_time_traces
end module mod_dislin_plots
