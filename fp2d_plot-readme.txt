================================================================================
  fp2d_plot.py  --  Run & plot wrapper for the FP2D_QLRF_NL Fokker-Planck solver
  F. Louche -- LPP-ERM/KMS
  June 2026
================================================================================

OVERVIEW
--------
fp2d_plot.py is a command-line tool that can run the FP2D_QLRF_NL solver and
immediately visualise its output, or plot the output files produced by a
previous run.  It reads the Fortran text output files directly and produces
matplotlib figures (interactive and/or saved as PNG).

Three sub-commands are available:

  run      Execute the solver, then plot all outputs.
  plot     Plot the output files in an existing directory.
  compare  Overlay the same-type outputs from several simulation cases on
           shared axes, for direct visual comparison.


DEPENDENCIES
------------
  Python >= 3.8
  numpy, matplotlib, scipy  (install with: pip install numpy matplotlib scipy)


================================================================================
SUBCOMMAND: run
================================================================================
Execute the solver and plot all output files produced by the run.

SYNTAX
------
  python fp2d_plot.py run <exe> <input.dat> [options]

POSITIONAL ARGUMENTS
  <exe>         Path to the FP2D_QLRF_NL executable.
  <input.dat>   Fortran namelist input file.

SPECIFIC OPTIONS
  --outdir DIR  Working directory for the solver (default: same folder as
                <input.dat>).  Output files are written here, and plots read
                them from here.
  --out FILE    Redirect the solver's standard output (Fortran WRITE(*,*))
                to FILE instead of the terminal.

AUTO-DETECTION FROM THE NAMELIST
  The following namelist parameters are read automatically from <input.dat>
  and influence which plots are produced:

  casename      Appended to every plot title and to output filenames
                (e.g. "JET-beam7-TD0-NLSC").  Can be overridden with
                --casename.
  ntimes        If ntimes(1) = 0 the run is steady-state; all time-trace
                plots are skipped automatically (equivalent to --steady-state).
  isc           If isc = 0 there are no self-collisions; self-collision power
                plots are suppressed automatically (equivalent to --no-sc).
  iplot_pow     If iplot_pow = 0 the solver writes no power-vs-time files;
                the power and SC-power-split figures are suppressed
                automatically (equivalent to --no-pow).
  iplot_mom     If iplot_mom = 0 the solver writes no momentum-vs-time files;
                momentum breakdown and balance figures are suppressed
                automatically (equivalent to --no-mom).

EXAMPLES
  # Run and display interactively
  python fp2d_plot.py run FP2D_QLRF_NL.exe jet_beam_case7.txt --show

  # Run, save PNGs, logarithmic scale for distribution functions
  python fp2d_plot.py run FP2D_QLRF_NL.exe jet_beam_case7.txt ^
      --outdir x64/Release --save plots/beam7 --log

  # Redirect solver output to a log file
  python fp2d_plot.py run FP2D_QLRF_NL.exe jet_beam_case7.txt ^
      --out run_beam7.log --save plots/beam7


================================================================================
SUBCOMMAND: plot
================================================================================
Plot the output files already present in a directory (no solver is launched).

SYNTAX
------
  python fp2d_plot.py plot <outdir> [options]

POSITIONAL ARGUMENTS
  <outdir>  Directory containing the solver output .txt files.

CHOOSING THE CASE(S)
  --cases takes one or more casenames and draws a SEPARATE set of figures for
  each; use the compare subcommand if you want them overlaid on shared axes
  instead.  With --show, the windows of every case open together.  --casename
  is accepted as an alias, so older command lines keep working.

  A casename you type also selects WHICH FILES ARE READ, not merely how they
  are labelled.  So

      --cases ITER_EDA-RF-NLSC_2 --files fout

  plots fout-ITER_EDA-RF-NLSC_2.txt alone, and not the fout of every other case
  in the directory.  Note that a longer name is a different case: the command
  above does NOT pick up fout-ITER_EDA-RF-NLSC_2-Grid_1.txt.

  Omit --cases and the casename is auto-detected from the filenames.  An
  auto-detected name only labels the plots, so  --files fout  with no --cases
  still plots the fout of every case in a multi-case directory.

EXAMPLES
  # Plot everything, display interactively
  python fp2d_plot.py plot x64/Release --show

  # Save PNGs for a named case
  python fp2d_plot.py plot x64/Release ^
      --cases JET-beam7-TD0-NLSC --save plots/beam7

  # Plot only two specific files
  python fp2d_plot.py plot x64/Release ^
      --files fout.txt energy_vs_time-JET-beam7-TD0-NLSC.txt --show

  # The same quantity for two cases, one set of windows each
  python fp2d_plot.py plot . ^
      --cases ITER_EDA-RF-NLSC_2 ITER_EDA-RF-NLSC_2-Grid_1 ^
      --files fout --show

  # Skip momentum plots (e.g. they were not written)
  python fp2d_plot.py plot x64/Release --no-mom --show

  # Zoom every plot to the low-velocity region 0 - 5e6 m/s
  python fp2d_plot.py plot x64/Release --xrange 0:5e6 --show


================================================================================
SUBCOMMAND: compare
================================================================================
Overlay the same-type output files from multiple simulation cases on shared
axes, making it easy to compare the effect of self-collisions (Lin vs NLSC),
different injection angles, different timesteps, etc.

SYNTAX
------
  python fp2d_plot.py compare <outdir> --cases CASE1 CASE2 [...] [options]

POSITIONAL ARGUMENTS
  <outdir>  Directory containing the output .txt files for all listed cases.

REQUIRED OPTION
  --cases CASE1 CASE2 [CASE3 ...]
            List of casenames to overlay.  The script looks for files of the
            form  <stem>-<casename>.txt  in <outdir>.

BEHAVIOR
  For each known file type (time traces, 1D profiles) that is present for at
  least one of the listed cases, one figure is produced with all available
  cases on the same axes.

  - Each case is drawn in a distinct colour.
  - For multi-column files (e.g. energy_vs_time which contains E_total and
    E_perp), the columns within the same case keep the same colour but use
    different linestyles (solid / dashed / dotted).
  - 2D distribution files (fout, Ekin, beam) are skipped -- superposing
    contour maps on the same axes is not meaningful.
  - Output files are named  <stem>-compare.png.

--files OPTION IN COMPARE MODE
  The --files option accepts any of:
    (a) a bare stem key    e.g.  anisotropy_vs_time
    (b) a full filename    e.g.  anisotropy_vs_time-JET-beam7-TD0-Lin.txt
    (c) a glob             e.g.  "power_*"      (quote it, or the shell expands
                                                 it against the local folder)
  Form (a) is shorter and recommended when --cases is already given.

  An entry that matches nothing produces a warning naming it, instead of
  quietly plotting less than you asked for.

FILES WHOSE NAME CONTAINS A SPACE
  The per-background-species outputs carry the species index in the stem, with
  a space before it:

      power_coll_ion 1_vs_time-<casename>.txt
      momentum_coll_ion 1_vs_time-<casename>.txt

  All of these spellings work, in both plot and compare mode, quoted or not:

      --files power_coll_ion 1_vs_time
      --files "power_coll_ion 1_vs_time"
      --files power_coll_ion1_vs_time
      --files power_coll_ion_1_vs_time

  (Matching ignores spaces, underscores and letter case, and entries that a
  shell split at the space are rejoined automatically.)

EXAMPLES
  # Compare Lin and NLSC for all available file types
  python fp2d_plot.py compare x64/Release ^
      --cases JET-beam7-TD0-Lin JET-beam7-TD0-NLSC --show

  # Compare only anisotropy and effective temperature
  python fp2d_plot.py compare x64/Release ^
      --cases JET-beam7-TD0-Lin JET-beam7-TD0-NLSC ^
      --files anisotropy_vs_time Teff_vs_time --show

  # Save comparison PNGs to a sub-folder
  python fp2d_plot.py compare x64/Release ^
      --cases JET-beam7-TD0-Lin JET-beam7-TD0-NLSC ^
      --save plots/Lin_vs_NLSC

  # Compare 1D profiles, zoomed to the perpendicular tail
  python fp2d_plot.py compare x64/Release ^
      --cases ITER-RF1-TD1-NLMax1 ITER-RF1-TD1-NLSC ^
      --files sc_Fpe_at_vpar0 --xrange 0:5e6 --show


================================================================================
COMMON OPTIONS  (all three sub-commands)
================================================================================
  --show            Open interactive matplotlib windows.  This is the default
                    when --save is not given.

  --save DIR        Save each figure as a PNG file in DIR (created if absent).
                    Can be combined with --show.

  --log             Use logarithmic y-scale for 1D profiles, or logarithmic
                    colour scale for 2D contour maps (distribution functions).

  --xrange xmin:xmax
                    (plot and compare only)  Zoom the x-axis of every figure
                    to the interval [xmin, xmax], and rescale the y-axis to the
                    data within that window.  Applies to whatever lies on the
                    x-axis of each plot: v_par for 2D maps, v_perp or v_par for
                    1D profiles, and time for *_vs_time traces.  Bounds are in
                    the file's own x-axis units: m/s for velocities (the grid
                    spans ~0 to a few 1e7 m/s, so use e.g. 0:5e6, NOT 0:1) and
                    seconds for time traces.  Accepts a ':' or ',' separator
                    and requires xmin < xmax.  If the window contains no data a
                    warning is printed (the plot would otherwise be blank).
                    Example: --xrange 0:5e6

  --cases CASE ...  (plot and compare)  Case(s) to work on.  In plot mode each
                    case gets its own set of figures and the name also selects
                    which files are read; in compare mode the cases are
                    overlaid on shared axes.  Omitted in plot mode, the case is
                    auto-detected from the filenames and then only labels the
                    plots.  --casename is accepted as an alias in plot mode.

  --casename STR    (run only)  Case label for the output filenames and plot
                    titles.  Read from the namelist when omitted.

  --files F ...     Restrict plotting to the listed files.  A bare stem key
                    (--files anisotropy_vs_time), a full filename or a quoted
                    glob are all accepted, in both plot and compare mode.

  --steady-state    (run and plot only)  Skip all time-trace plots.  Set
                    automatically in run mode when ntimes(1) = 0.

  --no-sc           Suppress self-collision power plots.  Set automatically
                    in run mode when isc = 0.

  --no-pow          Suppress power-vs-time and SC-power-split plots.  Set
                    automatically in run mode when iplot_pow = 0.

  --no-mom          Suppress momentum breakdown and balance plots.  Set
                    automatically in run mode when iplot_mom = 0.

  --3d              (run and plot only)  Add 3D surface plots alongside each
                    2D contour map (fout, Ekin, beam, ...).


================================================================================
OUTPUT FILES AND THEIR PLOTS
================================================================================

INDIVIDUAL FIGURES
(one figure per file; produced by run, plot, and compare)

  1D velocity-space profiles
  --------------------------
  fout_at_vpar0.txt           VDF f(v_perp) at v_par = 0
  fout_at_vperp0.txt          VDF f(v_par)  at v_perp = 0
  fout_at_vperpmax.txt        VDF f(v_par)  at v_perp = v_perp,max
  fout_at_vparmax.txt         VDF f(v_perp) at v_par  = v_par,max
  Ekin_perp_at_vpar0.txt      Perpendicular kinetic energy E_perp(v_perp)
                                at v_par = 0
  fstix_at_vpar0.txt          Stix Maxwellian at v_par = 0
  fstix_at_vperp0.txt         Stix Maxwellian at v_perp = 0
  fsc_maxw_at_vpar0.txt       SC Maxwellian background at v_par = 0, evaluated
                                at the final T_eff  (written when isc = 2)

  Self-collision diagnostics  (v_par = 0 slices; SC_diagnostics build)
  --------------------------
  sc_Dpepe_at_vpar0.txt       SC perpendicular diffusion  D_perp,perp(v_perp)
  sc_Fpe_at_vpar0.txt         SC perpendicular friction   F_perp(v_perp)
  sc_power_density_at_vpar0.txt  SC power density  dP_SC/d^3v = 1/2 m v^2 C_SC[f]
                                at v_par = 0 (>0 source, <0 sink)

  2D distribution maps
  --------------------
  fout.txt                    Full 2D VDF f(v_perp, v_par)  -- filled contour
  Ekin.txt                    Total kinetic energy map (keV)
  Ekin_perp.txt               Perpendicular kinetic energy map (keV)
  Ekin_par.txt                Parallel kinetic energy map (keV)
  beam.txt                    Beam source S(v_perp, v_par)
  fsc_maxw.txt                SC Maxwellian background f_M(v_perp, v_par) at the
                                final T_eff  (written when isc = 2)
  sc_power_density.txt        SC power density map dP_SC/d^3v (signed; drawn with
                                a diverging colour scale: red > 0 source,
                                blue < 0 sink)  -- SC_diagnostics build

  Time traces -- single curve
  ----------------------------
  anisotropy_vs_time.txt      Perpendicular anisotropy (%):
                                0%  = all energy in parallel motion
                                50% = Maxwellian (isotropic)
                               100% = all energy in perpendicular motion
  Teff_vs_time.txt            Effective temperature T_eff (keV)

  Time traces -- multiple curves
  --------------------------------
  energy_vs_time.txt          Two curves: E_total and E_perp (keV) vs time
  power_coll_self_vs_time.txt Three curves: P_SC total, P_SC_perp, P_SC_par
                                (MW/m^3) -- self-collision power split
  momentum_coll_e_vs_time.txt Two curves: perp. and par. electron collisional
                                momentum transfer rate (N/m^3) vs time
  momentum_coll_self_vs_time.txt  Two curves: perp. and par. self-collision
                                   momentum transfer rate (N/m^3) vs time
  momentum_RF_vs_time.txt     Two curves: perp. and par. RF momentum transfer
  momentum_coll_tot_vs_time.txt   Two curves: total perp. and par. collisional
                                   momentum transfer rate vs time


COMPOSITE FIGURES
(multi-curve or multi-panel; produced by run and plot only, not by compare)

  power_vs_time               Two-panel figure:
                                Top:    collisional power breakdown by species
                                        (electrons, bulk ions, self-collisions)
                                Bottom: full power balance
                                        (collisions + RF + beam source
                                         + particle losses + net sum)
                              Controlled by iplot_pow / --no-pow.

  power_sc_split_vs_time      Self-collision power split:
                                P_SC_perp and P_SC_par vs time, with the
                                near-zero total as a reference line.
                                Quantifies the rate of energy transfer from
                                perpendicular to parallel motion by
                                self-collisions.
                              Only produced when isc != 0 and iplot_pow = -1.

  coulomb_log_all_vs_time     Coulomb logarithm vs time for each background
                                species and (if isc != 0) for self-collisions.

  momentum_breakdown_vs_time  Grid of panels: one column per operator
                                (electrons, ions, RF, beam, self-collisions),
                                perpendicular (top row) and parallel (bottom
                                row) momentum transfer rate vs time.
                              Controlled by iplot_mom / --no-mom.

  momentum_balance_vs_time    Momentum balance for both perp. and par.
                                directions: collisions + RF + beam source
                                + losses + net sum vs time.
                              Controlled by iplot_mom / --no-mom.


================================================================================
CASENAME CONVENTION
================================================================================
When the namelist contains a non-empty casename (e.g. "JET-beam7-TD0-NLSC"),
the solver appends it to every output filename:

  density_vs_time-JET-beam7-TD0-NLSC.txt
  anisotropy_vs_time-JET-beam7-TD0-NLSC.txt
  fout-JET-beam7-TD0-NLSC.txt
  ...

fp2d_plot.py detects the casename automatically in run mode by reading the
namelist.  In plot and compare mode, name the case(s) with --cases; in plot
mode it is auto-detected from the filenames in the directory if you do not.

The distinction matters in a directory holding several cases.  A casename you
give with --cases selects which files are read; an auto-detected one only
labels the plots.  So

  plot . --cases NLSC_2 --files fout   plots fout-NLSC_2.txt only
  plot . --files fout                  plots the fout of every case present

Two further points on these names.  A longer name is a different case:
--cases NLSC_2 does not match fout-NLSC_2-Grid_1.txt.  And the per-species
files carry a space in the stem ("power_coll_ion 1_vs_time-<casename>.txt");
--files accepts that name with the space, with an underscore, or with nothing
in its place, quoted or not.


================================================================================
NAMELIST PARAMETERS THAT AFFECT PLOTTING
================================================================================
  iplot_pow = -1   (default)  Power-vs-time files are written; power figures
                               are produced.
  iplot_pow =  0              Power files are NOT written; power figures are
                               suppressed automatically in run mode.

  iplot_mom = -1              Momentum-vs-time files are written; momentum
                               figures are produced.
  iplot_mom =  0   (default)  Momentum files are NOT written; momentum figures
                               are suppressed automatically in run mode.

  isc = -1 or != 0            Self-collision output files are written; SC
                               plots are produced.
  isc =  0                    No self-collision files; SC plots suppressed.

  ntimes(1) = 0               Steady-state run; all time traces suppressed.


================================================================================
END OF MANUAL
================================================================================
