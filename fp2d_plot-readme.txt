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

  For --movie only:
  Pillow           Writes GIF movies.  Installed with matplotlib.
  ffmpeg           Needed for AVI and MP4 movies.  Either an ffmpeg on the
                   PATH, or the copy bundled with the imageio-ffmpeg package:
                       pip install imageio-ffmpeg
                   fp2d_plot.py finds that copy by itself.  Without any ffmpeg,
                   movies fall back to GIF.


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
  n_snap        With --movie only.  If n_snap = 0 the solver would write no
                snapshots, so the command stops with an error BEFORE the solver
                is launched, rather than after a run that cannot give a movie.

EXAMPLES
  # Run and display interactively
  python fp2d_plot.py run FP2D_QLRF_NL.exe jet_beam_case7.txt --show

  # Run, save PNGs, logarithmic scale for distribution functions
  python fp2d_plot.py run FP2D_QLRF_NL.exe jet_beam_case7.txt ^
      --outdir x64/Release --save plots/beam7 --logf

  # Redirect solver output to a log file
  python fp2d_plot.py run FP2D_QLRF_NL.exe jet_beam_case7.txt ^
      --out run_beam7.log --save plots/beam7

  # Run, then make a movie of the snapshots (needs n_snap > 0 in the namelist)
  python fp2d_plot.py run FP2D_QLRF_NL.exe jet_beam_case7.txt ^
      --outdir x64/Release --movie --logf --save movies


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

  # Log scales: the VDF on a log colour scale over a log v_perp axis
  python fp2d_plot.py plot x64/Release --files fout --logf --logy --show

  # A time trace over several decades, both axes logarithmic
  python fp2d_plot.py plot x64/Release --files energy_vs_time --logx --logy --show

  # Fix both axes: the first second, energies up to 400 keV
  python fp2d_plot.py plot x64/Release --files energy_vs_time ^
      --xrange 0:1 --yrange 0:400 --show

  # The VDF tail, six decades below the peak
  python fp2d_plot.py plot x64/Release --files fout_at_vpar0 ^
      --logy --yrange 1e-6:10 --show

  # Skip momentum plots (e.g. they were not written)
  python fp2d_plot.py plot x64/Release --no-mom --show

  # Zoom every plot to the low-velocity region 0 - 5e6 m/s
  python fp2d_plot.py plot x64/Release --xrange 0:5e6 --show

  # Movie (AVI) of f and Ekin from the snapshots of one case, f on a log scale
  python fp2d_plot.py plot x64/Release --cases ITER_EDA-RF-NLSC_2 ^
      --movie --logf

  # The same as an MP4 at 8 frames/s, zoomed, with the namelist named
  # explicitly (it is not in the output folder)
  python fp2d_plot.py plot x64/Release --cases ITER_EDA-RF-NLSC_2 ^
      --movie --logf --movie-format mp4 --fps 8 --xrange=-5e6:5e6 ^
      --namelist inputs/ITER_EDA_N=2.txt --save movies


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

  # Overlay the VDF tails on a logarithmic y-axis
  python fp2d_plot.py compare x64/Release ^
      --cases JET-beam7-TD0-Lin JET-beam7-TD0-NLSC ^
      --files fout_at_vpar0 --logy --show

  # Compare 1D profiles, zoomed to the perpendicular tail
  python fp2d_plot.py compare x64/Release ^
      --cases ITER-RF1-TD1-NLMax1 ITER-RF1-TD1-NLSC ^
      --files sc_Fpe_at_vpar0 --xrange 0:5e6 --show


================================================================================
MOVIES FROM THE VDF SNAPSHOTS  (--movie;  run and plot)
================================================================================
With the namelist parameter n_snap = N > 0, the solver saves the distribution
function every N time steps (both time-dependent solvers):

  vdf_snap_<step>-<casename>.txt      (vdf_snap_<step>.txt without casename)

<step> is the 6-digit total step count.  Each file has a "# time = <t>" header
line followed by three columns v_perp, v_par, f, as in fout.txt.  A restart
(iold = -1) continues the numbering of the run it restarts from, so its
snapshots do not overwrite the earlier ones.

--movie turns these snapshots into one animation per case, INSTEAD of the
usual figures.  Any of the movie-only options (--movie-scale, --movie-format,
--fps) implies --movie, and a message says so: they mean nothing without a
movie, so  plot . --cases X --movie-scale frame  makes the movie rather than
silently drawing the ordinary figures.

WHAT IS SHOWN
  Two panels per frame, side by side:
    left    f(v_perp, v_par)
    right   the kinetic-energy density (keV), the quantity of Ekin.txt
  The frame title gives the case, the time, the step and the frame number, and
  each panel title the maximum of that panel in that frame.  One frame per
  snapshot, in step order.

  COLOUR SCALES (--movie-scale)
    fixed   (default)  One scale for the whole movie, from the largest value
            of any frame, so frames can be compared directly.  A frame of
            much smaller amplitude than the largest then looks nearly empty.
    frame   Each frame scaled to its own maximum: 0 .. max for a linear f and
            for the energy, eight decades below max with --logf.  Every frame
            shows its structure whatever its amplitude; the colour bars follow
            the frame, so read the amplitude from the bar or the panel title.

  The snapshots store f only.  The energy panel is computed from each
  snapshot exactly as the solver computes Ekin.txt (it reproduces the solver's
  own Ekin.txt to round-off); this needs the ion mass number aa, read from the
  namelist.

  The snapshot f is the raw solution: the end-of-run renormalisation to npart
  that sourceless runs apply to fout.txt is not applied to it.

THE NAMELIST, AND THE n_snap CHECK
  The movie needs the namelist of the run, for n_snap and aa.
    run    The namelist is <input.dat>.  If n_snap = 0 there, the command
           stops BEFORE launching the solver.
    plot   The namelist is given with --namelist FILE; if omitted, the
           output folder is searched for a namelist (a file starting with
           &INPUT) whose casename is the case being animated.
  The movie of a case is skipped, with a message saying why, when no namelist
  is found, when n_snap = 0, when no snapshot files exist for the case, or
  when there is only one.

OPTIONS ACTING ON THE MOVIE
  --logf                f on a logarithmic colour scale covering eight decades
                        below its maximum; anything smaller is left blank,
                        rather than painting the round-off of the far tail in
                        the lowest colour.  The energy panel stays linear.
  --xrange / --yrange   Zoom both panels (v_par and v_perp, in m/s).
  --fps N               Frames per second (default 5).
  --movie-scale S       fixed (default) or frame.  See COLOUR SCALES above.
  --movie-format F      avi (default), mp4 or gif.  See below.
  --save DIR            Folder for the movie (default: the output folder).

FORMATS
  avi    MPEG-4 video, tagged XviD; plays in VLC and in the Windows media
         players.  Needs ffmpeg.
  mp4    H.264 video; the smallest file.  Needs ffmpeg.
  gif    Animated GIF, written with Pillow; no ffmpeg needed.
  ffmpeg is taken from the PATH, or else from the imageio-ffmpeg package
  (pip install imageio-ffmpeg).  Without either, avi and mp4 fall back to gif
  and a message says so.  The movie is written as  movie-<casename>.<format>.

CHOOSING n_snap
  At 150 x 150 each snapshot is about 0.7 MB, and each one is a frame.  For a
  run of 1000 steps, n_snap = 25 gives 40 frames (28 MB of snapshots), or 8 s
  of movie at the default 5 frames/s.  The step size changes between phases
  (timestep(1..3)), so equal step intervals are not equal time intervals; the
  time in each frame title is always the true time.

EXAMPLES
  # AVI of one case, f on a log scale
  python fp2d_plot.py plot x64/Release --cases ITER_EDA-RF-NLSC_2 ^
      --movie --logf

  # GIF, 2 frames/s, namelist given explicitly
  python fp2d_plot.py plot x64/Release --cases ITER_EDA-RF-NLSC_2 ^
      --movie --movie-format gif --fps 2 --namelist inputs/ITER_EDA_N=2.txt

  # Each frame on its own colour scale (e.g. a beam filling up from zero)
  python fp2d_plot.py plot x64/Release --cases JET-beam7-TD0-NLSC ^
      --movie --movie-scale frame


================================================================================
COMMON OPTIONS  (all three sub-commands)
================================================================================
  --show            Open interactive matplotlib windows.  This is the default
                    when --save is not given.

  --save DIR        Save each figure as a PNG file in DIR (created if absent).
                    Can be combined with --show.

  --logf            Logarithmic scale for the plotted QUANTITY: the colour
                    scale of a 2D contour map, the z of a 3D surface, the y of
                    a 1D profile.  --log is an accepted alias.

  --logx            Logarithmic x-axis.

  --logy            Logarithmic y-axis.  On a 1D profile the y-axis IS the
                    function, so --logy and --logf coincide there.

                    The three combine freely, e.g.
                      --files fout --logy --logf
                    gives a 2D map with a log v_perp axis and a log colour
                    scale.

                    An axis that genuinely CROSSES zero is left linear and a
                    note is printed: v_par is signed, and so is a power that
                    changes sign, so a log scale would discard half the data
                    with nothing on the figure to say so.

                    An axis that merely STARTS at zero (v_perp, or a time trace
                    at t = 0), or that dips below it only by round-off (the far
                    tail of a VDF, a few 1e-9 of the peak), is drawn
                    logarithmically; the points that cannot be shown are
                    counted in a note.  The threshold between the two is a
                    negative excursion of 1e-6 of the positive range.

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
                    A NEGATIVE lower bound (e.g. a v_par window) must be
                    attached with '=', or it is read as an option name:
                      --xrange=-5e6:5e6      (works)
                      --xrange -5e6:5e6      (error: expected one argument)
                    The same holds for --yrange.

  --yrange ymin:ymax
                    (plot and compare only)  The same for the y-axis: bounds in
                    the file's own units, a ':' or ',' separator, ymin < ymax,
                    and a warning when the window holds no data.
                    Example: --yrange 1e-6:10

                    Note the interaction with --xrange.  Given alone, --xrange
                    rescales the y-axis to the data inside the x-window, which
                    is usually what you want when zooming a trace.  Giving
                    --yrange suppresses that rescale and uses the bounds you
                    asked for, so the two can be combined to fix both axes:
                      --xrange 0:1 --yrange 0:400

                    It is applied after --logy, so it also sets the limits of a
                    logarithmic axis:
                      --files fout_at_vpar0 --logy --yrange 1e-6:10
                    plots the VDF tail down to six decades below the peak.

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

  --movie           (run and plot only)  Instead of the usual figures, make a
                    movie of f and the kinetic-energy density from the VDF
                    snapshots.  Needs n_snap > 0 in the namelist.  See MOVIES
                    FROM THE VDF SNAPSHOTS above.

  --fps N           (run and plot only)  Frames per second of the movie
                    (default 5).  Implies --movie.

  --movie-scale S   (run and plot only)  Colour scales of the movie: fixed
                    over the whole movie (default) or scaled to each frame's
                    own maximum (frame).  Implies --movie.

  --movie-format F  (run and plot only)  avi (default), mp4 or gif.  avi and
                    mp4 need ffmpeg (on the PATH, or via pip install
                    imageio-ffmpeg) and fall back to gif without it.  Implies
                    --movie.

  --namelist FILE   (plot only)  Namelist of the run, read by --movie for
                    n_snap and aa.  Default: the namelist in the output folder
                    whose casename matches the case.


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

  VDF snapshots  (not plotted one by one)
  -------------
  vdf_snap_<step>.txt         f every n_snap steps.  Skipped by the normal
                                plotting, which would otherwise draw every
                                snapshot as a separate map; animated by
                                --movie instead.

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

  n_snap =  0      (default)  No VDF snapshots; --movie has nothing to show
                               (run --movie stops before the solver starts).
  n_snap =  N > 0             f is saved every N steps to vdf_snap_<step>.txt,
                               which --movie animates.

  aa                          Ion mass number; read by --movie to compute the
                               kinetic-energy panel.


================================================================================
END OF MANUAL
================================================================================
