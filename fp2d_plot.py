#!/usr/bin/env python3
"""
fp2d_plot.py – Run & plot wrapper for the FP2D_QLRF_NL Fokker-Planck solver.

Subcommands
-----------
  run  <exe> <input.dat>   Execute the solver then plot all outputs.
  plot <outdir>            Plot from an existing output directory.

Options (both subcommands)
--------------------------
  --save DIR      Write PNG files to DIR (created if absent).
  --show          Open interactive matplotlib windows
                  (default when --save is not given).
  --logf          Logarithmic scale for the plotted QUANTITY: the colour scale
                  of a 2D map, the z of a 3D surface, the y of a 1D profile.
                  (--log is an accepted alias.)
  --logx          Logarithmic x-axis.
  --logy          Logarithmic y-axis (on a 1D profile, same as --logf).
                  An axis that genuinely crosses zero -- v_par, or a power that
                  changes sign -- is left linear, with a note saying so.  One
                  that starts at zero (v_perp, t=0) or dips below it only by
                  round-off (a VDF tail) is drawn logarithmically, and the
                  points not shown are counted in a note.
  --3d            Add 3D surface plots for 2D distribution files.
  --files F ...   Plot only these filenames (basenames, e.g. fout.txt).
  --movie         Instead of the usual plots, animate f and the kinetic-energy
                  density over the VDF snapshots vdf_snap_<step>.txt.  Checks
                  n_snap > 0 in the namelist first (for 'run', before the
                  solver starts).  One frame per snapshot, colour scales fixed
                  over the movie; --logf puts f on a log scale, --xrange and
                  --yrange zoom both panels.  Written as movie-<case>.<fmt> to
                  --save DIR or the output directory.
  --fps N         Frames per second of the movie (default 5).
                  --fps, --movie-scale and --movie-format each imply --movie.
  --movie-scale {fixed,frame}
                  Colour scales of the movie: fixed over the whole movie, so
                  frames compare (default), or scaled to each frame's own
                  maximum.  Each panel title gives the frame's maximum.
  --movie-format {avi,mp4,gif}
                  Movie file format (default avi).  avi (MPEG-4) and mp4
                  (H.264) need ffmpeg: one on the PATH, or the binary bundled
                  with  pip install imageio-ffmpeg.  Without it the movie
                  falls back to an animated GIF.

plot / compare only
-------------------
  --cases CASE ... Cases to plot (plot) or overlay (compare).  In plot mode
                   each case gets its own set of figures, and the case also
                   selects which files are read, not merely how they are
                   labelled; --casename is accepted as an alias.  Omit it and
                   the case is auto-detected from the filenames.
  --xrange xmin:xmax   Zoom the x-axis of every plot to [xmin, xmax]
                       (e.g. --xrange 0:5e6). Accepts ':' or ',' separators.
  --yrange ymin:ymax   The same for the y-axis. Given alone, --xrange rescales
                       y to the data in the window; --yrange suppresses that
                       and uses the bounds asked for.
  --namelist FILE      (plot only) Namelist of the run, read by --movie for
                       n_snap and aa.  Default: the namelist in the output
                       directory whose casename matches the case.

run only
--------
  --casename STR  Case label for the output filenames and plot titles
                  (read from the namelist if omitted).

run-only options
----------------
  --outdir DIR    Solver working directory (default: folder of <input.dat>).
  --out FILE      Redirect solver stdout (Fortran write(*,*)) to FILE.

Examples
--------
  # Run then display interactively
  python fp2d_plot.py run FP2D_QLRF_NL.exe inputs/jet_RF_case1.dat --show

  # Run, save PNGs, log scale for f
  python fp2d_plot.py run FP2D_QLRF_NL.exe inputs/jet_RF_case1.dat \\
      --outdir x64/Release --save plots/jet_run1 --log

  # Plot an existing output directory
  python fp2d_plot.py plot x64/Release --save plots/jet_run1

  # Plot specific files only
  python fp2d_plot.py plot x64/Release --files fout.txt energy_vs_time.txt --show

  # Same files for two cases, one set of windows each
  python fp2d_plot.py plot . --cases RF-NLSC_2 RF-NLSC_2-Grid_1 --files fout --show

  # Movie of f (log scale) and Ekin from the snapshots of one case
  python fp2d_plot.py plot x64/Release --cases RF-NLSC_2 --movie --logf --fps 8
"""

import argparse
import fnmatch
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 – registers 3d projection
from scipy.interpolate import RectBivariateSpline

# ---------------------------------------------------------------------------
# Per-file metadata: ptype in {"1d", "2d", "ts", "ts2"}
#   1d  – (x, y) line plot
#   2d  – (vperp, vpar, z) filled-contour plot
#   ts  – (time, y) time-series line plot
#   ts2 – (time, y1, y2, ...) multi-curve time series
# ---------------------------------------------------------------------------
FILE_META = {
    # 1-D slices ---------------------------------------------------------------
    "fout_at_vpar0":      {"ptype": "1d",  "xlabel": "v⊥ (v_th)", "ylabel": "f",
                           "title": "VDF at v∥ = 0",              "sci_y": True},
    "fout_at_vperp0":     {"ptype": "1d",  "xlabel": "v∥ (v_th)", "ylabel": "f",
                           "title": "VDF at v⊥ = 0",             "sci_y": True},
    "fout_at_vperpmax":   {"ptype": "1d",  "xlabel": "v∥ (v_th)", "ylabel": "f",
                           "title": "VDF at v⊥ = v⊥,max",   "sci_y": True},
    "fout_at_vparmax":    {"ptype": "1d",  "xlabel": "v⊥ (v_th)", "ylabel": "f",
                           "title": "VDF at v∥ = v∥,max",   "sci_y": True},
    "Ekin_perp_at_vpar0": {"ptype": "1d",  "xlabel": "v⊥ (v_th)",
                           "ylabel": "E_kin,⊥ (keV)",
                           "title": "Perp. kinetic energy at v∥ = 0"},
    "fstix_at_vpar0":     {"ptype": "1d",  "xlabel": "v⊥ (v_th)", "ylabel": "f_Stix",
                           "title": "Stix Maxwellian at v∥ = 0", "sci_y": True},
    "fstix_at_vperp0":    {"ptype": "1d",  "xlabel": "v∥ (v_th)", "ylabel": "f_Stix",
                           "title": "Stix Maxwellian at v⊥ = 0", "sci_y": True},
    "fstix_at_vperpmax":  {"ptype": "1d",  "xlabel": "v∥ (v_th)", "ylabel": "f_Stix",
                           "title": "Stix Maxwellian at v⊥ = v⊥,max", "sci_y": True},
    "fstix_at_vparmax":   {"ptype": "1d",  "xlabel": "v⊥ (v_th)", "ylabel": "f_Stix",
                           "title": "Stix Maxwellian at v∥ = v∥,max", "sci_y": True},
    # 2-D grids ----------------------------------------------------------------
    "fout":               {"ptype": "2d",  "title": "2D VDF  f(v⊥, v∥)",
                           "sci_z": True},
    "Ekin":               {"ptype": "2d",  "title": "Kinetic energy (keV)"},
    "Ekin_perp":          {"ptype": "2d",  "title": "Perp. kinetic energy (keV)"},
    "Ekin_par":           {"ptype": "2d",  "title": "Par. kinetic energy (keV)"},
    "beam":               {"ptype": "2d",  "title": "Beam source  S(v⊥, v∥)",
                           "sci_z": True},
    "fsc_maxw":           {"ptype": "2d",  "title": "SC Maxwellian  f_M(v⊥, v∥)  at final T_eff",
                           "sci_z": True},
    # 1-D slices of SC Maxwellian ------------------------------------------
    "fsc_maxw_at_vpar0":  {"ptype": "1d",  "xlabel": "v⊥ (v_th)", "ylabel": "f_M",
                           "title": "SC Maxwellian at v∥ = 0  (final T_eff)", "sci_y": True},
    # SC coefficient diagnostics (DISABLED) --------------------------------
    # The solver no longer writes sc_*_at_vpar0.txt: the two calls that
    # produced them are commented out in TimeFP_7pt.f90 (isc=1,2,3) and
    # TimeFP_7pt_NL.f90 (isc=-1).  These entries are commented out with
    # them so that an old run directory that still holds the files does not
    # plot them either.  Uncomment both sides to get the diagnostic back;
    # the plot functions themselves need no change.
    # "sc_Dpepe_at_vpar0":  {"ptype": "1d",  "xlabel": "v⊥ (v_th)",
    #                         "ylabel": "D⊥⊥ (m² s⁻³)",
    #                         "title":  "SC diffusion D⊥⊥ at v∥ = 0", "sci_y": True},
    # "sc_Dpapa_at_vpar0":  {"ptype": "1d",  "xlabel": "v⊥ (v_th)",
    #                         "ylabel": "D∥∥ (m² s⁻³)",
    #                         "title":  "SC diffusion D∥∥ at v∥ = 0", "sci_y": True},
    # "sc_Dpepa_at_vpar0":  {"ptype": "1d",  "xlabel": "v⊥ (v_th)",
    #                         "ylabel": "D⊥∥ (m² s⁻³)",
    #                         "title":  "SC cross diffusion D⊥∥ at v∥ = 0", "sci_y": True},
    # "sc_Fpe_at_vpar0":    {"ptype": "1d",  "xlabel": "v⊥ (v_th)",
    #                         "ylabel": "F⊥ (m s⁻²)",
    #                         "title":  "SC friction F⊥ at v∥ = 0",   "sci_y": True},
    # "sc_Fpa_at_vpar0":    {"ptype": "1d",  "xlabel": "v⊥ (v_th)",
    #                         "ylabel": "F∥ (m s⁻²)",
    #                         "title":  "SC friction F∥ at v∥ = 0",   "sci_y": True},
    # "sc_psi_at_vpar0":    {"ptype": "1d",  "xlabel": "v⊥ (v_th)",
    #                         "ylabel": "ψ",
    #                         "title":  "Rosenbluth potential ψ at v∥ = 0", "sci_y": True},
    # "sc_phi_at_vpar0":    {"ptype": "1d",  "xlabel": "v⊥ (v_th)",
    #                         "ylabel": "φ",
    #                         "title":  "Rosenbluth potential φ at v∥ = 0", "sci_y": True},
    # Simple time series -------------------------------------------------------
    "density_vs_time":         {"ptype": "ts",  "ylabel": "Density (m⁻³)",
                                "title": "Particle density vs time"},
    "anisotropy_vs_time":      {"ptype": "ts",  "ylabel": "Perp. anisotropy (%, 50=Maxwell)",
                                "title": "Anisotropy factor vs time"},
    "power_coll_tot_vs_time":  {"ptype": "ts",
                                "ylabel": "Power density (MW·m⁻³)",
                                "title": "Total collisional power density"},
    "power_coll_e_vs_time":    {"ptype": "ts",
                                "ylabel": "Power density (MW·m⁻³)",
                                "title": "Electron collisional power density"},
    "power_RF_vs_time":        {"ptype": "ts",
                                "ylabel": "Power density (MW·m⁻³)",
                                "title": "RF power density"},
    "power_coll_self_vs_time": {"ptype": "ts2",
                                "ylabel": "Power density (MW·m⁻³)",
                                "title": "Self-collision power split (⊥ / ∥)",
                                "labels": ["total", "⊥", "∥"],
                                "skip_cols": [1]},
    # Two-curve time series ----------------------------------------------------
    "energy_vs_time":               {"ptype": "ts2", "ylabel": "Energy (keV)",
                                     "title": "Kinetic energy vs time",
                                     "labels": ["E_total", "E_⊥"]},
    # Effective temperature
    "Teff_vs_time":                 {"ptype": "ts",  "ylabel": "T_eff (keV)",
                                     "title":  "Effective temperature vs time"},
    # Density-characteristic (log-slope) temperature
    "Tn_vs_time":                   {"ptype": "ts",  "ylabel": "T_n (keV)",
                                     "title":  "Density-characteristic temperature vs time"},
    # RF tail formation time, tau_RF = npart*Teff/P_RF (time_power_7pt): the
    # energy stored in the tail divided by the rate the RF supplies it.  Written
    # only when irf = -1; zero-filled otherwise, hence absent for non-RF runs.
    # Kept in FILE_META so 'compare' can overlay them across cases; both are in
    # _SKIP_STEMS because the composite plot_timescales draws them together.
    "tau_rf_vs_time":               {"ptype": "ts",  "ylabel": "τ_RF (s)",
                                     "title":  "RF tail formation time vs time"},
    # Effective collisional times at the running Teff, npart*Teff/|P_coll|
    # (time_power_7pt).  Log y: they span decades as the tail heats, and at
    # steady state 1/tau_RF = 1/tau_ii + 1/tau_ie.
    "tau_coll_vs_time":             {"ptype": "ts2", "ylabel": "τ (s)",
                                     "title":  "Effective collision times vs time",
                                     "labels": ["τ_ii  (to background ions)",
                                                "τ_ie  (to electrons)"],
                                     "logy": True},
    # Convergence rate epsilon = ||f^n - f^(n-1)|| / (dt ||f^n||), Jacobian-weighted
    # (mod_conv_diag).  Log y: epsilon decays over orders of magnitude as the run
    # converges.  eps_tail uses the extra tail weight and is plotted alongside --
    # the diagnostic failure mode is eps small while eps_tail is not.
    "conv_eps_vs_time":             {"ptype": "ts2", "ylabel": "ε  (1/s)",
                                     "title":  "Convergence rate ε vs time",
                                     "labels": ["ε  (amplitude, Jacobian-weighted)",
                                                "ε_tail  (amplitude, tail-weighted)",
                                                "ε_shape  (f normalised)",
                                                "ε_tail,shape  (f normalised, tail)"],
                                     "logy": True},
    # Coulomb logarithms (data files behind the coulomb_log_all_vs_time plot).
    # Listed here so 'compare' can overlay them; still skipped in plot mode
    # (the composite plot_coulomb_log handles them) via _SKIP_STEMS.
    "coulomb_log_self_vs_time":     {"ptype": "ts",  "ylabel": "ln Λ (self)",
                                     "title":  "Self-collision Coulomb logarithm vs time"},
    "coulomb_log_vs_time":          {"ptype": "ts2", "ylabel": "ln Λ",
                                     "title":  "Background-ion Coulomb logarithm vs time",
                                     "labels": ["ion 1", "ion 2", "ion 3"]},
    # Momentum transfer rate (⊥ and ∥ per file) --------------------------------
    "momentum_coll_tot_vs_time":    {"ptype": "ts2",
                                     "ylabel": "Momentum transfer rate (N·m⁻³)",
                                     "title": "Total collisional momentum transfer",
                                     "labels": ["⊥", "∥"]},
    "momentum_coll_e_vs_time":      {"ptype": "ts2",
                                     "ylabel": "Momentum transfer rate (N·m⁻³)",
                                     "title": "Electron collisional momentum transfer",
                                     "labels": ["⊥", "∥"]},
    "momentum_RF_vs_time":          {"ptype": "ts2",
                                     "ylabel": "Momentum transfer rate (N·m⁻³)",
                                     "title": "RF momentum transfer",
                                     "labels": ["⊥", "∥"]},
    "momentum_coll_self_vs_time":   {"ptype": "ts2",
                                     "ylabel": "Momentum transfer rate (N·m⁻³)",
                                     "title": "Self-collision momentum transfer",
                                     "labels": ["⊥", "∥"]},
}

# Files to skip (unusual format or not useful for plotting)
_SKIP_STEMS = {"RF_dirac", "fstix",
               "density_vs_time",
               "power_RF_vs_time",
               "power_coll_e_vs_time", "power_coll_ion", "power_coll_self_vs_time",
               "power_NBI_vs_time", "power_coll_tot_vs_time",
               "momentum_coll_e_vs_time", "momentum_coll_ion",
               "momentum_RF_vs_time", "momentum_coll_self_vs_time",
               "momentum_NBI_vs_time", "momentum_coll_tot_vs_time",
               "coulomb_log_vs_time", "coulomb_log_self_vs_time",
               "tau_coll_vs_time", "tau_rf_vs_time",
               # VDF snapshots (n_snap > 0): one 2-D file every n_snap steps.
               # Auto-detection would draw each of them as a map; they are
               # meant for --movie, which animates them instead.
               "vdf_snap"}

# SC coefficient / Rosenbluth-potential diagnostics, switched off.  The solver
# no longer writes them (the two calls are commented out in TimeFP_7pt.f90 and
# TimeFP_7pt_NL.f90), but older run directories still hold the files, and an
# unknown stem would otherwise be auto-detected and plotted with default
# labels.  To re-enable the diagnostic: uncomment the calls in the solver, the
# FILE_META entries above, and this update().
_SKIP_STEMS.update({"sc_Dpepe_at_vpar0", "sc_Dpapa_at_vpar0", "sc_Dpepa_at_vpar0",
                    "sc_Fpe_at_vpar0", "sc_Fpa_at_vpar0",
                    "sc_psi_at_vpar0", "sc_phi_at_vpar0",
                    "sc_phi_raw_at_vpar0", "sc_Fpe_raw_at_vpar0",
                    "sc_Fpa_raw_at_vpar0",
                    # No writer left in the sources for these two; they survive
                    # only in run directories from before that writer was removed.
                    "sc_power_density", "sc_power_density_at_vpar0"})

_PALETTE = ["#8B1A1A", "#1A1A8B", "#1A8B1A", "#8B8B1A", "#8B1A8B", "#1A8B8B"]

# Axis zoom for the 'plot' and 'compare' subcommands; (lo, hi) or None.
# Set in main() from --xrange / --yrange and applied to every figure in
# _finish().  An explicit _YRANGE also suppresses the automatic y-rescale that
# _XRANGE would otherwise perform, so the two can be combined.
_XRANGE = None
_YRANGE = None

# Logarithmic AXES, set in main() from --logx/--logy and applied to every
# figure in _finish().  The third log option, --logf, scales the plotted
# QUANTITY (the colour scale of a 2D map, the z of a surface, the y of a 1D
# profile); that cannot be a post-hoc axis change, so it stays the `log`
# argument threaded through the individual plot functions.
_LOGX = False
_LOGY = False


def _parse_range(s: str, opt: str = "--xrange"):
    """Parse a 'lo:hi' (or 'lo,hi') bound pair, as given to --xrange/--yrange."""
    txt = s.strip().lstrip("[").rstrip("]")
    sep = ":" if ":" in txt else ","
    parts = txt.split(sep)
    a, b = ("xmin", "xmax") if opt == "--xrange" else ("ymin", "ymax")
    if len(parts) != 2:
        sys.exit(f"Error: {opt} expects '{a}:{b}', got '{s}'")
    try:
        lo, hi = float(parts[0]), float(parts[1])
    except ValueError:
        sys.exit(f"Error: {opt} bounds must be numbers, got '{s}'")
    if hi <= lo:
        sys.exit(f"Error: {opt} requires {a} < {b}, got '{s}'")
    return (lo, hi)


def _parse_xrange(s: str):
    """Kept for callers and tests that use the original name."""
    return _parse_range(s, "--xrange")


def _get_meta(stem: str) -> dict:
    """Return FILE_META entry for stem, falling back to prefix match for case-named files."""
    if stem in FILE_META:
        return FILE_META[stem]
    # Case-named files: e.g. "fout_at_vpar0-JET-1HarmH_1-TD1-Lin"
    for key in sorted(FILE_META, key=len, reverse=True):
        if stem.startswith(key + "-"):
            return dict(FILE_META[key])
    return {}


# ---------------------------------------------------------------------------
# Data loading & type detection
# ---------------------------------------------------------------------------

def _load(path: Path):
    try:
        data = np.loadtxt(path)
    except Exception:
        return None
    if data.ndim == 1:
        data = data.reshape(-1, 1)
    return data


def _load_ncol(path: Path, ncols: int):
    """Load a Fortran list-directed file whose records may wrap across lines."""
    try:
        tokens = path.read_text().split()
        arr = np.array(tokens, dtype=float)
    except Exception:
        return None
    if arr.size == 0 or arr.size % ncols != 0:
        return None
    return arr.reshape(-1, ncols)


def _load_auto(path: Path):
    """Load a text file; fall back to token-based reshape for line-wrapped records."""
    data = _load(path)
    if data is not None:
        return data
    try:
        tokens = path.read_text().split()
        arr = np.array(tokens, dtype=float)
        if arr.size == 0:
            return None
        for ncols in (4, 3, 2):
            if arr.size % ncols == 0:
                return arr.reshape(-1, ncols)
    except Exception:
        pass
    return None


def _detect_type(path: Path, data: np.ndarray) -> str:
    """Infer plot type from filename and column count."""
    stem  = path.stem
    ncols = data.shape[1]
    if "_vs_time" in stem:
        return "ts2" if ncols >= 3 else "ts"
    if ncols == 2:
        return "1d"
    if ncols >= 3:
        # 2-D grid: first column (vperp) contains many repeated values
        n_unique = len(np.unique(data[:, 0]))
        return "2d" if n_unique < len(data) * 0.9 else "ts"
    return "unknown"


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def _autoscale_y_to_xrange(ax, xrange) -> None:
    """Rescale an axes' y-limits to the line data lying within *xrange*.

    Only line plots are affected: contour maps, colorbars and 3D surfaces have
    no Line2D data and are left untouched. Reference lines (axhline/axvline,
    exactly 2 points) are ignored so they cannot corrupt the range.
    """
    xmin, xmax = xrange
    is_log = ax.get_yscale() == "log"
    ys = []
    for line in ax.get_lines():
        xd = np.asarray(line.get_xdata(), dtype=float)
        yd = np.asarray(line.get_ydata(), dtype=float)
        if xd.size <= 2:                       # skip axhline/axvline references
            continue
        m = (xd >= xmin) & (xd <= xmax) & np.isfinite(yd)
        if is_log:
            m &= yd > 0
        if np.any(m):
            ys.append(yd[m])
    if not ys:
        return
    yall = np.concatenate(ys)
    ylo, yhi = float(yall.min()), float(yall.max())
    if is_log:
        lo, hi = np.log10(ylo), np.log10(yhi)
        pad = 0.05 * (hi - lo) if hi > lo else 0.1
        ax.set_ylim(10.0**(lo - pad), 10.0**(hi + pad))
    else:
        pad = 0.05 * (yhi - ylo) if yhi > ylo else (abs(yhi) * 0.05 or 1.0)
        ax.set_ylim(ylo - pad, yhi + pad)


# A negative excursion smaller than this fraction of the positive range is
# round-off, not signal: the far tail of a VDF dips a few 1e-9 of its peak
# below zero, and refusing a log axis over that would block the commonest use
# of --logy there.  Anything larger is real structure that a log scale would
# hide, so the axis is left linear instead.
_LOG_NEG_TOL = 1e-6


def _positive_span(ax, which: str):
    """(smallest positive, largest, n_dropped) of the data on one axis of *ax*.

    Line data is preferred, ignoring the 2-point reference lines drawn by
    axhline/axvline.  Where there are no lines -- a contour map -- the axes'
    own data limits are used instead, so a log scale still works on the
    v_perp axis of a 2D plot.  Returns None when a log scale would hide real
    data, which is the signal to leave the axis linear.
    """
    vals = []
    for line in ax.get_lines():
        d = np.asarray(line.get_xdata() if which == "x" else line.get_ydata(),
                       dtype=float)
        if d.size <= 2:
            continue
        d = d[np.isfinite(d)]
        if d.size:
            vals.append(d)
    if vals:
        allv = np.concatenate(vals)
        lo, hi = float(allv.min()), float(allv.max())
        if hi <= 0 or _hides_data(lo, hi):
            return None
        pos = allv[allv > 0]
        if not pos.size:
            return None
        return float(pos.min()), hi, int((allv <= 0).sum())

    box = ax.dataLim
    lo, hi = (box.x0, box.x1) if which == "x" else (box.y0, box.y1)
    if not np.isfinite([lo, hi]).all() or hi <= 0 or _hides_data(lo, hi):
        return None
    return (lo if lo > 0 else hi * 1e-4), float(hi), 0


def _hides_data(lo: float, hi: float) -> bool:
    """True if a log scale on this range would conceal real data.

    An axis that merely STARTS at zero (v_perp, t=0 on a time trace) is fine:
    one end point is lost.  So is one whose negative excursion is round-off,
    like the far tail of a VDF.  One that is genuinely SIGNED (v_par, or a
    power that changes sign) is not -- half the data would vanish with nothing
    on the figure to say so, which is worse than no log scale at all.
    """
    return lo < -_LOG_NEG_TOL * abs(hi)


def _apply_log_axes(fig, stem: str) -> None:
    """Apply --logx / --logy to the data axes of a figure.

    A log axis needs positive data, and several of these axes legitimately
    reach or cross zero: v_par is signed, and every time trace starts at t=0.
    Where a log scale would hide real data the axis is left linear and the
    reason is printed -- a half-empty plot with no explanation is worse than a
    linear one.  Otherwise the axis is set to log and its lower limit pinned to
    the smallest positive sample, which stops matplotlib from padding down to
    an arbitrary decade; any non-positive points dropped are reported.
    """
    if not (_LOGX or _LOGY):
        return
    for ax in fig.axes:
        if ax.get_label() == "<colorbar>":
            continue
        for which, want in (("x", _LOGX), ("y", _LOGY)):
            if not want:
                continue
            span = _positive_span(ax, which)
            if span is None:
                # e.g. --logx on a v_par axis, which crosses zero
                print(f"    (--log{which}: the {which}-axis of {stem} crosses "
                      f"zero, left linear)")
                continue
            lo, hi, dropped = span
            getattr(ax, f"set_{which}scale")("log")
            lim = getattr(ax, f"get_{which}lim")()
            if lim[0] <= 0:
                getattr(ax, f"set_{which}lim")(lo * 0.9, max(hi * 1.1, lim[1]))
            if dropped:
                print(f"    (--log{which}: {dropped} non-positive point(s) "
                      f"not shown on the {which}-axis of {stem})")


def _warn_if_range_empty(fig, stem: str, rng, which: str) -> None:
    """Warn when --xrange/--yrange selects no line data (otherwise the plot is
    blank with no hint why — usually a units mismatch, e.g. 0:1 on an m/s
    axis)."""
    lo, hi = rng
    vals, in_window = [], False
    for ax in fig.axes:
        for line in ax.get_lines():
            d = np.asarray(line.get_xdata() if which == "x" else line.get_ydata(),
                           dtype=float)
            if d.size <= 2:                    # skip reference lines
                continue
            vals.append(d)
            if np.any((d >= lo) & (d <= hi)):
                in_window = True
    if vals and not in_window:
        allv = np.concatenate(vals)
        print(f"    WARNING: --{which}range [{lo:g}, {hi:g}] selects no data for "
              f"'{stem}' ({which} spans [{allv.min():g}, {allv.max():g}]); "
              f"plot is empty.")


def _finish(fig, stem: str, save_dir, show: bool) -> None:
    if _XRANGE is not None and fig.axes:
        # Zoom the data axes' x-range. fig.axes[0] is always the main plot:
        # colorbars are appended afterwards, and every multi-panel figure is
        # built with sharex=True, so a single set_xlim propagates to all panels.
        fig.axes[0].set_xlim(_XRANGE)
        # Rescale y to the data now visible in the x-window (line plots only),
        # unless the user has asked for a y-window of their own.
        if _YRANGE is None:
            for ax in fig.axes:
                _autoscale_y_to_xrange(ax, _XRANGE)
        _warn_if_range_empty(fig, stem, _XRANGE, "x")
    _apply_log_axes(fig, stem)
    if _YRANGE is not None:
        # After _apply_log_axes, which may have set its own limits.  Applied to
        # every data axes rather than just the first: multi-panel figures share
        # x but not y, so each panel needs it.
        for ax in fig.axes:
            if ax.get_label() != "<colorbar>":
                ax.set_ylim(_YRANGE)
        _warn_if_range_empty(fig, stem, _YRANGE, "y")
    if save_dir is not None:
        out = Path(save_dir).resolve() / f"{stem}.png"
        fig.savefig(out, dpi=150, bbox_inches="tight")
        print(f"    -> {out}")
    if not show:
        plt.close(fig)


def _title(meta: dict, fallback: str, casename: str) -> str:
    base = meta.get("title", fallback)
    return f"{base}  [{casename}]" if casename else base


def _outfile(outdir: Path, stem: str, casename: str) -> Path:
    """Construct the output filename the Fortran code would produce for this stem."""
    name = f"{stem}-{casename}.txt" if casename else f"{stem}.txt"
    return outdir / name


def _matches_casename(path: Path, casename: str) -> bool:
    """True when path belongs to the current job (identified by casename)."""
    stem = path.stem
    if casename:
        return stem.endswith(f"-{casename}")
    # Empty casename: accept only files with no '-' suffix (base filenames only).
    # BASE stems in FILE_META use underscores, never hyphens, so any hyphen means
    # the file belongs to a named case that isn't the current one.
    return "-" not in stem


def _canon_stem(s: str) -> str:
    """Canonical form of a stem for --files matching.

    The per-background-species files carry a space in the stem, e.g.
    'power_coll_ion 1_vs_time-<case>.txt' (the Fortran builds the name from the
    species index).  A shell splits an unquoted argument at that space, so every
    plausible spelling -- with the space, with an underscore, or with nothing at
    all -- has to match the same file.  Case is ignored too.
    """
    if s.endswith(".txt"):
        s = s[:-4]
    return s.replace(" ", "").replace("_", "").lower()


def _rejoin_split_stems(pats, keys) -> list:
    """Undo the shell's splitting of a stem that contains a space.

    '--files power_coll_ion 1_vs_time' reaches argparse as two entries.  Where
    two adjacent entries rejoin into a stem that actually exists, treat them as
    the single name the user meant; everything else passes through untouched.
    """
    canon = {_canon_stem(k) for k in keys}

    def known(p: str) -> bool:
        # keys may be bare stems or whole filenames carrying a '-<case>' suffix
        c = _canon_stem(p)
        return c in canon or any(k.startswith(c + "-") for k in canon)

    out, i = [], 0
    while i < len(pats):
        if i + 1 < len(pats) and not known(pats[i]) \
                and known(pats[i] + pats[i + 1]):
            out.append(f"{pats[i]} {pats[i + 1]}")
            i += 2
        else:
            out.append(pats[i])
            i += 1
    return out


def _restrict_match(name: str, pat: str) -> bool:
    """True if output filename *name* matches a --files entry *pat*.

    *pat* may be (a) a glob (contains * ? [), used with fnmatch; (b) a full
    filename ending in .txt, matched exactly; or (c) a bare stem key such as
    'Teff_vs_time', which matches both the base file 'Teff_vs_time.txt' and any
    case-named 'Teff_vs_time-<case>.txt' (so plot mode now behaves like compare).
    Stems also match canonically, so spaces and underscores need not line up.
    """
    if any(c in pat for c in "*?["):
        return fnmatch.fnmatch(name, pat)
    if pat.endswith(".txt"):
        return name == pat
    stem = name[:-4] if name.endswith(".txt") else name
    if stem == pat or stem.startswith(pat + "-"):
        return True
    cpat, cstem = _canon_stem(pat), _canon_stem(stem)
    return cstem == cpat or cstem.startswith(cpat + "-")


# ---------------------------------------------------------------------------
# Individual plot functions
# ---------------------------------------------------------------------------

def plot_1d(data, meta, stem, save_dir, show, log, casename):
    x, y = data[:, 0], data[:, 1]
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(x, y, color=_PALETTE[0], linewidth=1.5)
    ax.set_xlabel(meta.get("xlabel", "x"))
    ax.set_ylabel(meta.get("ylabel", "y"))
    ax.set_title(_title(meta, stem, casename))
    ax.grid(True, alpha=0.3)
    if x.min() >= 0:
        ax.set_xlim(left=0)
    if log and np.any(y > 0):
        ax.set_yscale("log")
    elif meta.get("sci_y"):
        ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    fig.tight_layout()
    _finish(fig, stem, save_dir, show)


def plot_2d(data, meta, stem, save_dir, show, log, casename):
    vperp_u = np.unique(data[:, 0])
    vpar_u  = np.unique(data[:, 1])
    nperp, npar = len(vperp_u), len(vpar_u)
    try:
        Z = data[:, 2].reshape(nperp, npar)   # Z[i_vperp, j_vpar]
    except ValueError:
        print(f"    Cannot reshape {len(data[:, 2])} rows → {nperp}×{npar}; skipping")
        return

    zmin, zmax = float(Z.min()), float(Z.max())
    if zmin == zmax:
        zmax = zmin + 1.0

    diverging = meta.get("diverging", False)
    if diverging:
        # Signed field: symmetric range about 0 with a diverging colormap so
        # sources (>0) and sinks (<0) are immediately distinguishable.
        zabs   = max(abs(zmin), abs(zmax)) or 1.0
        norm   = None
        levels = np.linspace(-zabs, zabs, 21)
        cmap   = "RdBu_r"
    elif log and np.any(Z > 0):
        pos    = Z[Z > 0]
        norm   = mcolors.LogNorm(vmin=float(pos.min()), vmax=zmax)
        levels = np.logspace(np.log10(float(pos.min())), np.log10(zmax), 20)
        cmap   = "rainbow"
    else:
        norm   = None
        levels = np.linspace(zmin, zmax, 20)
        cmap   = "rainbow"

    fig, ax = plt.subplots(figsize=(8, 6))
    cf   = ax.contourf(vpar_u, vperp_u, Z, levels=levels, norm=norm, cmap=cmap)
    cbar = fig.colorbar(cf, ax=ax)
    if meta.get("sci_z") and norm is None:
        cbar.formatter.set_powerlimits((0, 0))
        cbar.update_ticks()
    ax.set_xlabel("v∥ (v_th)")
    ax.set_ylabel("v⊥ (v_th)")
    ax.set_title(_title(meta, stem, casename))
    ax.grid(True, alpha=0.2)
    fig.tight_layout()
    _finish(fig, stem, save_dir, show)


def plot_3d(data, meta, stem, save_dir, show, log, casename):
    vperp_u = np.unique(data[:, 0])
    vpar_u  = np.unique(data[:, 1])
    nperp, npar = len(vperp_u), len(vpar_u)
    try:
        Z = data[:, 2].reshape(nperp, npar)
    except ValueError:
        print(f"    Cannot reshape {len(data[:, 2])} rows → {nperp}×{npar}; skipping 3D")
        return

    # Upsample to a finer grid for a smoother surface (cap at 300 per axis)
    nperp_f = min(3 * nperp, 300)
    npar_f  = min(3 * npar,  300)
    vperp_f = np.linspace(vperp_u[0], vperp_u[-1], nperp_f)
    vpar_f  = np.linspace(vpar_u[0],  vpar_u[-1],  npar_f)

    if log and np.any(Z > 0):
        Z_pos  = np.where(Z > 0, Z, np.nanmin(Z[Z > 0]))
        spline = RectBivariateSpline(vperp_u, vpar_u, np.log10(Z_pos))
        Z_plot = spline(vperp_f, vpar_f)
        zlabel = "log₁₀(z)"
    else:
        spline = RectBivariateSpline(vperp_u, vpar_u, Z)
        Z_plot = spline(vperp_f, vpar_f)
        zlabel = "z"

    X, Y = np.meshgrid(vpar_f, vperp_f)

    fig = plt.figure(figsize=(9, 6))
    ax  = fig.add_subplot(111, projection="3d")
    surf = ax.plot_surface(X, Y, Z_plot, cmap="rainbow",
                           linewidth=0, antialiased=True)
    fig.colorbar(surf, ax=ax, shrink=0.5, aspect=10)
    ax.set_xlabel("v∥")
    ax.set_ylabel("v⊥")
    ax.set_zlabel(zlabel)
    ax.set_title(_title(meta, stem, casename))
    fig.tight_layout()
    _finish(fig, stem + "_3d", save_dir, show)


def plot_ts(data, meta, stem, save_dir, show, casename):
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(data[:, 0], data[:, 1], color=_PALETTE[0], linewidth=1.5)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel(meta.get("ylabel", ""))
    ax.set_title(_title(meta, stem, casename))
    ax.set_xlim(left=0)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    _finish(fig, stem, save_dir, show)


def plot_ts2(data, meta, stem, save_dir, show, casename):
    labels = meta.get("labels", [])
    fig, ax = plt.subplots(figsize=(8, 5))
    for i in range(1, data.shape[1]):
        lbl = labels[i - 1] if i - 1 < len(labels) else f"col {i}"
        ax.plot(data[:, 0], data[:, i], color=_PALETTE[(i - 1) % len(_PALETTE)],
                linewidth=1.5, label=lbl)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel(meta.get("ylabel", ""))
    ax.set_title(_title(meta, stem, casename))
    ax.set_xlim(left=0)
    ax.legend()
    ax.grid(True, alpha=0.3)
    # Opt-in log y-axis (e.g. the convergence rate, which decays over orders
    # of magnitude).  Guarded: log scale needs at least one positive sample.
    if meta.get("logy") and np.any(data[:, 1:] > 0):
        ax.set_yscale("log")
    fig.tight_layout()
    _finish(fig, stem, save_dir, show)


# ---------------------------------------------------------------------------
# Composite plots (multi-curve, mirroring dislin's power plots)
# ---------------------------------------------------------------------------

def _ts_col(path: Path, col: int = 1):
    """Load a time-series file and return (time_array, col_array) or (None, None)."""
    data = _load(path)
    if data is None or data.shape[1] <= col:
        return None, None
    return data[:, 0], data[:, col]


def _ts_col2(path: Path):
    """Load a 3-column time-series; return (t, col1, col2) or (None, None, None)."""
    data = _load(path)
    if data is None or data.shape[1] < 3:
        return None, None, None
    return data[:, 0], data[:, 1], data[:, 2]


def plot_power_coll(outdir: Path, save_dir, show: bool, casename: str,
                    show_sc: bool = True) -> None:
    """Collisional power breakdown: electrons + bulk ions + self-collisions."""
    base = _outfile(outdir, "power_coll_e_vs_time", casename)
    if not base.exists():
        return
    fig, ax = plt.subplots(figsize=(8, 5))
    plotted = False

    t, y = _ts_col(base)
    if t is not None:
        ax.plot(t, y, color=_PALETTE[0], linewidth=1.5, label="e⁻")
        plotted = True

    for ib in range(1, 10):
        # Fortran writes ibstr with i2 format: " 1", " 2", … (leading space)
        f = _outfile(outdir, f"power_coll_ion {ib}_vs_time", casename)
        if not f.exists():
            f = _outfile(outdir, f"power_coll_ion{ib}_vs_time", casename)
        if not f.exists():
            break
        t, y = _ts_col(f)
        if t is not None:
            ax.plot(t, y, color=_PALETTE[ib % len(_PALETTE)], linewidth=1.5, label=f"ion {ib}")
            plotted = True

    if show_sc:
        f = _outfile(outdir, "power_coll_self_vs_time", casename)
        if f.exists():
            _sc = _load_ncol(f, 4)
            if _sc is not None:
                ax.plot(_sc[:, 0], _sc[:, 1], color=_PALETTE[4],
                        linewidth=1.5, linestyle="--", label="self")
                plotted = True

    if not plotted:
        plt.close(fig)
        return
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Power density (MW·m⁻³)")
    ax.set_title(_title({}, "Collisional power density vs time", casename))
    ax.set_xlim(left=0)
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    stem_out = f"power_coll_vs_time-{casename}" if casename else "power_coll_vs_time"
    _finish(fig, stem_out, save_dir, show)


def plot_power_balance(outdir: Path, save_dir, show: bool, casename: str) -> None:
    """Power balance: total collisions + RF + NBI + net sum."""
    base = _outfile(outdir, "power_coll_tot_vs_time", casename)
    if not base.exists():
        return
    fig, ax = plt.subplots(figsize=(8, 5))
    plotted = False

    # Accumulate net sum on the common time grid
    t_ref, y_sum = None, None

    def _add_to_sum(t, y):
        nonlocal t_ref, y_sum
        if t is None:
            return
        if t_ref is None:
            t_ref = t
            y_sum = y.copy()
        elif len(t) == len(t_ref):
            y_sum += y

    t, y = _ts_col(base)
    if t is not None:
        ax.plot(t, y, color=_PALETTE[0], linewidth=1.5, label="collisional")
        _add_to_sum(t, y)
        plotted = True

    f = _outfile(outdir, "power_RF_vs_time", casename)
    if f.exists():
        t, y = _ts_col(f)
        if t is not None:
            ax.plot(t, y, color=_PALETTE[1], linewidth=1.5, label="RF")
            _add_to_sum(t, y)
            plotted = True

    f = _outfile(outdir, "power_NBI_vs_time", casename)
    if f.exists():
        data = _load(f)
        if data is not None and data.shape[1] >= 3:
            ax.plot(data[:, 0], data[:, 1], color=_PALETTE[2], linewidth=1.5, label="NBI source")
            ax.plot(data[:, 0], data[:, 2], color=_PALETTE[2], linewidth=1.5,
                    linestyle="--", label="NBI losses")
            _add_to_sum(data[:, 0], data[:, 1])
            _add_to_sum(data[:, 0], data[:, 2])
            plotted = True

    if not plotted:
        plt.close(fig)
        return

    if t_ref is not None and y_sum is not None:
        ax.plot(t_ref, y_sum, color="black", linewidth=2.0, linestyle="--", label="net")

    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Power density (MW·m⁻³)")
    ax.set_title(_title({}, "Power balance vs time", casename))
    ax.set_xlim(left=0)
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    stem_out = f"power_balance_vs_time-{casename}" if casename else "power_balance_vs_time"
    _finish(fig, stem_out, save_dir, show)


def plot_power_combined(outdir: Path, save_dir, show: bool, casename: str,
                        show_sc: bool = True) -> None:
    """Two-panel: collisional power breakdown (top) and power balance (bottom)."""
    base_coll = _outfile(outdir, "power_coll_e_vs_time", casename)
    base_tot  = _outfile(outdir, "power_coll_tot_vs_time", casename)
    if not base_coll.exists() and not base_tot.exists():
        return
    fig, (ax_top, ax_bot) = plt.subplots(2, 1, figsize=(8, 8), sharex=True)

    # --- Top: collisional breakdown ---
    plotted_top = False
    if base_coll.exists():
        t, y = _ts_col(base_coll)
        if t is not None:
            ax_top.plot(t, y, color=_PALETTE[0], linewidth=1.5, label="e⁻")
            plotted_top = True

    for ib in range(1, 10):
        f = _outfile(outdir, f"power_coll_ion {ib}_vs_time", casename)
        if not f.exists():
            f = _outfile(outdir, f"power_coll_ion{ib}_vs_time", casename)
        if not f.exists():
            break
        t, y = _ts_col(f)
        if t is not None:
            ax_top.plot(t, y, color=_PALETTE[ib % len(_PALETTE)], linewidth=1.5, label=f"ion {ib}")
            plotted_top = True

    if show_sc:
        f = _outfile(outdir, "power_coll_self_vs_time", casename)
        if f.exists():
            _sc = _load_ncol(f, 4)
            if _sc is not None and _sc.shape[1] >= 4:
                ax_top.plot(_sc[:, 0], _sc[:, 2], color=_PALETTE[4],
                            linewidth=1.5, linestyle="--", label="SC ⊥")
                ax_top.plot(_sc[:, 0], _sc[:, 3], color=_PALETTE[5],
                            linewidth=1.5, linestyle=":",  label="SC ∥")
                plotted_top = True

    ax_top.set_ylabel("Power density (MW·m⁻³)")
    ax_top.set_title(_title({}, "Power vs time", casename))
    if plotted_top:
        ax_top.legend()
    ax_top.grid(True, alpha=0.3)
    ax_top.set_xlim(left=0)

    # --- Bottom: power balance ---
    plotted_bot = False
    t_ref, y_sum = None, None

    def _add(t, y):
        nonlocal t_ref, y_sum
        if t is None:
            return
        if t_ref is None:
            t_ref = t; y_sum = y.copy()
        elif len(t) == len(t_ref):
            y_sum += y

    if base_tot.exists():
        t, y = _ts_col(base_tot)
        if t is not None:
            ax_bot.plot(t, y, color=_PALETTE[0], linewidth=1.5, label="collisional")
            _add(t, y); plotted_bot = True

    f = _outfile(outdir, "power_RF_vs_time", casename)
    if f.exists():
        t, y = _ts_col(f)
        if t is not None:
            ax_bot.plot(t, y, color=_PALETTE[1], linewidth=1.5, label="RF")
            _add(t, y); plotted_bot = True

    f = _outfile(outdir, "power_NBI_vs_time", casename)
    if f.exists():
        data = _load(f)
        if data is not None and data.shape[1] >= 3:
            ax_bot.plot(data[:, 0], data[:, 1], color=_PALETTE[2], linewidth=1.5, label="NBI source")
            ax_bot.plot(data[:, 0], data[:, 2], color=_PALETTE[2], linewidth=1.5,
                        linestyle="--", label="NBI losses")
            _add(data[:, 0], data[:, 1]); _add(data[:, 0], data[:, 2])
            plotted_bot = True

    if t_ref is not None:
        ax_bot.plot(t_ref, y_sum, color="black", linewidth=2.0, linestyle="--", label="net")

    ax_bot.set_xlabel("Time (s)")
    ax_bot.set_ylabel("Power density (MW·m⁻³)")
    if plotted_bot:
        ax_bot.legend()
    ax_bot.grid(True, alpha=0.3)
    ax_bot.set_xlim(left=0)

    fig.tight_layout()
    stem_out = f"power_vs_time-{casename}" if casename else "power_vs_time"
    _finish(fig, stem_out, save_dir, show)


def plot_coulomb_log(outdir: Path, save_dir, show: bool, casename: str,
                     show_sc: bool = True) -> None:
    """Coulomb logarithm vs time: background ions + self-collision on one axes."""
    base    = _outfile(outdir, "coulomb_log_vs_time", casename)
    sc_file = _outfile(outdir, "coulomb_log_self_vs_time", casename)
    if not base.exists() and not sc_file.exists():
        return
    fig, ax = plt.subplots(figsize=(8, 5))
    plotted = False

    if base.exists():
        data = _load(base)
        if data is not None and data.shape[1] >= 2:
            for i in range(1, data.shape[1]):
                ax.plot(data[:, 0], data[:, i],
                        color=_PALETTE[(i - 1) % len(_PALETTE)],
                        linewidth=1.5, label=f"ion {i}")
                plotted = True

    if show_sc and sc_file.exists():
        t, y = _ts_col(sc_file)
        if t is not None:
            ax.plot(t, y, color=_PALETTE[4], linewidth=1.5,
                    linestyle="--", label="self")
            plotted = True

    if not plotted:
        plt.close(fig)
        return
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Coulomb logarithm")
    ax.set_title(_title({}, "Coulomb logarithm vs time", casename))
    ax.set_xlim(left=0)
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    stem_out = f"coulomb_log_all_vs_time-{casename}" if casename else "coulomb_log_all_vs_time"
    _finish(fig, stem_out, save_dir, show)


def plot_timescales(outdir: Path, save_dir, show: bool, casename: str) -> None:
    """Characteristic timescales vs time, all on one axes.

    tau_ii and tau_ie (tau_coll_vs_time) are the effective collisional times
    npart*Teff/|P_coll| for the background-ion and electron channels; tau_RF
    (tau_rf_vs_time) is npart*Teff/P_RF and exists only when irf = -1, so it is
    added only if its file is present.

    Log y: the three span decades while the tail forms.  At steady state
    1/tau_RF = 1/tau_ii + 1/tau_ie, so tau_RF sits below both collisional
    curves; where they cross tells which channel governs the balance.
    """
    coll = _outfile(outdir, "tau_coll_vs_time", casename)
    rf   = _outfile(outdir, "tau_rf_vs_time", casename)
    if not coll.exists() and not rf.exists():
        return

    fig, ax = plt.subplots(figsize=(8, 5))
    plotted = False

    if coll.exists():
        data = _load(coll)
        if data is not None and data.shape[1] >= 3:
            ax.plot(data[:, 0], data[:, 1], color=_PALETTE[0], linewidth=1.5,
                    label="τ_ii  (to background ions)")
            ax.plot(data[:, 0], data[:, 2], color=_PALETTE[1], linewidth=1.5,
                    label="τ_ie  (to electrons)")
            plotted = True

    if rf.exists():
        t, y = _ts_col(rf)
        if t is not None:
            ax.plot(t, y, color=_PALETTE[2], linewidth=1.8, linestyle="--",
                    label="τ_RF  (tail formation)")
            plotted = True

    if not plotted:
        plt.close(fig)
        return
    ax.set_yscale("log")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("τ (s)")
    ax.set_title(_title({}, "Characteristic timescales vs time", casename))
    ax.set_xlim(left=0)
    ax.legend()
    ax.grid(True, alpha=0.3, which="both")
    fig.tight_layout()
    stem_out = f"timescales_vs_time-{casename}" if casename else "timescales_vs_time"
    _finish(fig, stem_out, save_dir, show)


def plot_sc_power_split(outdir: Path, save_dir, show: bool, casename: str) -> None:
    """Self-collision power split: perpendicular and parallel components vs time.

    Reads power_coll_self_vs_time (4-column format: time, total, perp, par).
    Plots P_SC_⊥ and P_SC_∥; also draws the total as a thin reference so
    readers can verify P_SC_⊥ + P_SC_∥ ≈ 0 (energy conservation).
    Silently skips if the file is absent or has fewer than 4 columns (old format).
    """
    f = _outfile(outdir, "power_coll_self_vs_time", casename)
    if not f.exists():
        return
    data = _load_ncol(f, 4)   # list-directed WRITE wraps lines; tokenise first
    if data is None:
        return
    t = data[:, 0]
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(t, data[:, 2], color=_PALETTE[0], linewidth=1.5, label="P_SC ⊥")
    ax.plot(t, data[:, 3], color=_PALETTE[1], linewidth=1.5, label="P_SC ∥")
    ax.plot(t, data[:, 1], color="grey",       linewidth=0.8,
            linestyle="--", label="total (≈0)")
    ax.axhline(0, color="black", linewidth=0.5, linestyle=":")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Power density (MW·m⁻³)")
    ax.set_title(_title({}, "Self-collision power split (⊥ / ∥)", casename))
    ax.set_xlim(left=0)
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    stem_out = f"power_sc_split_vs_time-{casename}" if casename else "power_sc_split_vs_time"
    _finish(fig, stem_out, save_dir, show)


def plot_momentum_coll(outdir: Path, save_dir, show: bool, casename: str,
                       show_sc: bool = True) -> None:
    """Collisional momentum breakdown: electrons + bulk ions + self-collisions (⊥ and ∥)."""
    base = _outfile(outdir, "momentum_coll_e_vs_time", casename)
    if not base.exists():
        return
    fig, (ax_perp, ax_par) = plt.subplots(2, 1, figsize=(8, 8), sharex=True)
    plotted = False

    t, yp, yl = _ts_col2(base)
    if t is not None:
        ax_perp.plot(t, yp, color=_PALETTE[0], linewidth=1.5, label="e⁻")
        ax_par.plot(t, yl, color=_PALETTE[0], linewidth=1.5, label="e⁻")
        plotted = True

    for ib in range(1, 10):
        f = _outfile(outdir, f"momentum_coll_ion {ib}_vs_time", casename)
        if not f.exists():
            f = _outfile(outdir, f"momentum_coll_ion{ib}_vs_time", casename)
        if not f.exists():
            break
        t, yp, yl = _ts_col2(f)
        if t is not None:
            ax_perp.plot(t, yp, color=_PALETTE[ib % len(_PALETTE)], linewidth=1.5, label=f"ion {ib}")
            ax_par.plot(t, yl, color=_PALETTE[ib % len(_PALETTE)], linewidth=1.5, label=f"ion {ib}")
            plotted = True

    if show_sc:
        f = _outfile(outdir, "momentum_coll_self_vs_time", casename)
        if f.exists():
            t, yp, yl = _ts_col2(f)
            if t is not None:
                ax_perp.plot(t, yp, color=_PALETTE[4], linewidth=1.5, linestyle="--", label="self")
                ax_par.plot(t, yl, color=_PALETTE[4], linewidth=1.5, linestyle="--", label="self")
                plotted = True

    if not plotted:
        plt.close(fig)
        return
    ax_perp.set_ylabel("⊥ (N·m⁻³)")
    ax_par.set_ylabel("∥ (N·m⁻³)")
    ax_par.set_xlabel("Time (s)")
    ax_perp.set_title(_title({}, "Collisional momentum transfer vs time", casename))
    for ax in (ax_perp, ax_par):
        ax.set_xlim(left=0)
        ax.legend()
        ax.grid(True, alpha=0.3)
    fig.tight_layout()
    stem_out = f"momentum_coll_vs_time-{casename}" if casename else "momentum_coll_vs_time"
    _finish(fig, stem_out, save_dir, show)


def plot_momentum_breakdown(outdir: Path, save_dir, show: bool, casename: str,
                            show_sc: bool = True) -> None:
    """Grid (2 rows × N cols): each operator's ⊥ (top) and ∥ (bottom) momentum transfer."""
    # Collect (label, t, yp, yl) tuples in physical order
    contributions = []

    f = _outfile(outdir, "momentum_coll_e_vs_time", casename)
    if f.exists():
        t, yp, yl = _ts_col2(f)
        if t is not None:
            contributions.append(("e⁻", t, yp, yl))

    for ib in range(1, 10):
        f = _outfile(outdir, f"momentum_coll_ion {ib}_vs_time", casename)
        if not f.exists():
            f = _outfile(outdir, f"momentum_coll_ion{ib}_vs_time", casename)
        if not f.exists():
            break
        t, yp, yl = _ts_col2(f)
        if t is not None:
            contributions.append((f"ion {ib}", t, yp, yl))

    f = _outfile(outdir, "momentum_RF_vs_time", casename)
    if f.exists():
        t, yp, yl = _ts_col2(f)
        if t is not None:
            contributions.append(("RF", t, yp, yl))

    # NBI: 4-data-column file (msrc_perp, msrc_par, mloss_perp, mloss_par)
    f = _outfile(outdir, "momentum_NBI_vs_time", casename)
    if f.exists():
        data = _load_ncol(f, 5)
        if data is not None:
            contributions.append(("NBI src",  data[:, 0], data[:, 1], data[:, 2]))
            contributions.append(("NBI loss", data[:, 0], data[:, 3], data[:, 4]))

    if show_sc:
        f = _outfile(outdir, "momentum_coll_self_vs_time", casename)
        if f.exists():
            t, yp, yl = _ts_col2(f)
            if t is not None:
                contributions.append(("self", t, yp, yl))

    if not contributions:
        return

    n = len(contributions)
    fig, axes = plt.subplots(2, n, figsize=(4 * n, 6), sharex=True)
    if n == 1:
        axes = axes.reshape(2, 1)

    for col, (label, t, yp, yl) in enumerate(contributions):
        color = _PALETTE[col % len(_PALETTE)]
        axes[0, col].plot(t, yp, color=color, linewidth=1.5)
        axes[1, col].plot(t, yl, color=color, linewidth=1.5)
        axes[0, col].set_title(label)
        axes[1, col].set_xlabel("Time (s)")
        for ax in (axes[0, col], axes[1, col]):
            ax.grid(True, alpha=0.3)
            ax.set_xlim(left=0)
            ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

    axes[0, 0].set_ylabel("⊥ (N·m⁻³)")
    axes[1, 0].set_ylabel("∥ (N·m⁻³)")
    fig.suptitle(_title({}, "Momentum transfer by operator", casename))
    fig.tight_layout()
    stem_out = f"momentum_breakdown_vs_time-{casename}" if casename else "momentum_breakdown_vs_time"
    _finish(fig, stem_out, save_dir, show)


def plot_momentum_balance(outdir: Path, save_dir, show: bool, casename: str) -> None:
    """Momentum balance: total collisions + RF + NBI + net sum (⊥ and ∥)."""
    base = _outfile(outdir, "momentum_coll_tot_vs_time", casename)
    if not base.exists():
        return
    fig, (ax_perp, ax_par) = plt.subplots(2, 1, figsize=(8, 8), sharex=True)
    plotted = False

    t_ref, yp_sum, yl_sum = None, None, None

    def _add_to_sum(t, yp, yl):
        nonlocal t_ref, yp_sum, yl_sum
        if t is None:
            return
        if t_ref is None:
            t_ref = t
            yp_sum = yp.copy()
            yl_sum = yl.copy()
        elif len(t) == len(t_ref):
            yp_sum += yp
            yl_sum += yl

    t, yp, yl = _ts_col2(base)
    if t is not None:
        ax_perp.plot(t, yp, color=_PALETTE[0], linewidth=1.5, label="collisional")
        ax_par.plot(t, yl, color=_PALETTE[0], linewidth=1.5, label="collisional")
        _add_to_sum(t, yp, yl)
        plotted = True

    f = _outfile(outdir, "momentum_RF_vs_time", casename)
    if f.exists():
        t, yp, yl = _ts_col2(f)
        if t is not None:
            ax_perp.plot(t, yp, color=_PALETTE[1], linewidth=1.5, label="RF")
            ax_par.plot(t, yl, color=_PALETTE[1], linewidth=1.5, label="RF")
            _add_to_sum(t, yp, yl)
            plotted = True

    f = _outfile(outdir, "momentum_NBI_vs_time", casename)
    if f.exists():
        data = _load_ncol(f, 5)
        if data is not None:
            # columns: t, msrc_perp, msrc_par, mloss_perp, mloss_par
            ax_perp.plot(data[:, 0], data[:, 1], color=_PALETTE[2], linewidth=1.5, label="NBI source")
            ax_perp.plot(data[:, 0], data[:, 3], color=_PALETTE[2], linewidth=1.5,
                         linestyle="--", label="NBI losses")
            ax_par.plot(data[:, 0], data[:, 2], color=_PALETTE[2], linewidth=1.5, label="NBI source")
            ax_par.plot(data[:, 0], data[:, 4], color=_PALETTE[2], linewidth=1.5,
                        linestyle="--", label="NBI losses")
            _add_to_sum(data[:, 0], data[:, 1] + data[:, 3], data[:, 2] + data[:, 4])
            plotted = True

    if not plotted:
        plt.close(fig)
        return

    if t_ref is not None:
        ax_perp.plot(t_ref, yp_sum, color="black", linewidth=2.0, linestyle="--", label="net")
        ax_par.plot(t_ref, yl_sum, color="black", linewidth=2.0, linestyle="--", label="net")

    ax_perp.set_ylabel("⊥ (N·m⁻³)")
    ax_par.set_ylabel("∥ (N·m⁻³)")
    ax_par.set_xlabel("Time (s)")
    ax_perp.set_title(_title({}, "Momentum balance vs time", casename))
    for ax in (ax_perp, ax_par):
        ax.set_xlim(left=0)
        ax.legend()
        ax.grid(True, alpha=0.3)
    fig.tight_layout()
    stem_out = f"momentum_balance_vs_time-{casename}" if casename else "momentum_balance_vs_time"
    _finish(fig, stem_out, save_dir, show)


def plot_fout_vs_maxw_at_vpar0(outdir: Path, save_dir, show: bool, log: bool,
                               casename: str) -> None:
    """Overlay VDF and SC Maxwellian at v∥ = 0 on the same axes (isc=2 runs)."""
    f_maxw = _outfile(outdir, "fsc_maxw_at_vpar0", casename)
    if not f_maxw.exists():
        return
    fig, ax = plt.subplots(figsize=(8, 5))
    plotted = False

    f_fout = _outfile(outdir, "fout_at_vpar0", casename)
    if f_fout.exists():
        d = _load(f_fout)
        if d is not None and d.shape[1] >= 2:
            ax.plot(d[:, 0], d[:, 1], color=_PALETTE[0], linewidth=1.5, label="f  (VDF)")
            plotted = True

    d = _load(f_maxw)
    if d is not None and d.shape[1] >= 2:
        ax.plot(d[:, 0], d[:, 1], color=_PALETTE[1], linewidth=1.5,
                linestyle="--", label="f_M  (SC Maxwellian at T_eff)")
        plotted = True

    if not plotted:
        plt.close(fig)
        return

    ax.set_xlabel("v⊥ (v_th)")
    ax.set_ylabel("f")
    ax.set_title(_title({}, "VDF vs SC Maxwellian at v∥ = 0", casename))
    if log:
        all_y = np.concatenate([l.get_ydata() for l in ax.lines])
        if np.any(all_y > 0):
            ax.set_yscale("log")
    else:
        ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    if ax.get_xlim()[0] < 0:
        ax.set_xlim(left=0)
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    stem_out = f"fout_vs_maxw_at_vpar0-{casename}" if casename else "fout_vs_maxw_at_vpar0"
    _finish(fig, stem_out, save_dir, show)


# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------

def plot_directory(outdir: Path, save_dir, show: bool, log: bool,
                   casename: str, restrict=None, steady_state: bool = False,
                   show_sc: bool = True, show_pow: bool = True,
                   show_mom: bool = True, plot3d: bool = False,
                   strict_case: bool = False, defer_show: bool = False) -> None:
    if save_dir is not None:
        Path(save_dir).mkdir(parents=True, exist_ok=True)

    txt_files = sorted(outdir.glob("*.txt"))
    if restrict:
        # Each entry may be a glob, a full filename, or a bare stem key
        # (e.g. 'Teff_vs_time' matches 'Teff_vs_time-<case>.txt').  Stems that
        # contain a space, such as 'power_coll_ion 1_vs_time', arrive split
        # unless the user quoted them, so rejoin those first.
        restrict = _rejoin_split_stems(restrict, [f.stem for f in txt_files])
        txt_files = [f for f in txt_files
                     if any(_restrict_match(f.name, pat) for pat in restrict)]
        # --files selects WHICH quantities; an explicit --casename still selects
        # WHICH case.  Without this, '--files fout --casename X' plotted the
        # fout of every case in the directory.  An auto-detected casename does
        # not filter, so a multi-case folder can still be plotted across cases.
        if strict_case and casename:
            txt_files = [f for f in txt_files if _matches_casename(f, casename)]
    else:
        txt_files = [f for f in txt_files if _matches_casename(f, casename)]

    if not txt_files:
        # Nothing matched -- help the user rather than silently printing "Done.".
        detected = _detect_casename(outdir)
        hint = (f"\n  Hint: this folder contains case '{detected}' -- "
                f"add  --casename {detected}") if (detected and detected != casename) else ""
        if restrict:
            print(f"  No output files matched the --files pattern(s).{hint}")
        elif casename:
            print(f"  No output files matched casename '{casename}'.{hint}")
        else:
            print(f"  No case-named output files matched (empty casename).{hint}")
        return

    for path in txt_files:
        stem = path.stem
        # _SKIP_STEMS lists files a composite function draws instead, so the
        # per-file loop does not duplicate them.  Those composites only run
        # when the plot set is unrestricted, so a file named explicitly through
        # --files must be drawn here or it is drawn by nobody.
        if restrict is None and any(stem.startswith(s) for s in _SKIP_STEMS):
            continue

        data = _load(path)
        if data is None or data.shape[0] < 2:
            continue

        meta  = _get_meta(stem)
        ptype = meta.get("ptype") or _detect_type(path, data)

        if steady_state and ptype in ("ts", "ts2"):
            continue

        if not show_sc and stem.startswith("power_coll_self_vs_time"):
            continue

        print(f"  [{ptype:3s}]  {path.name}")
        # Use each file's own casename for the title, so a directory holding
        # several cases (e.g. plotted via --files) labels every plot with the
        # case it actually belongs to instead of one auto-detected casename.
        file_case = _case_from_stem(stem, casename)
        try:
            if ptype == "1d":
                plot_1d(data, meta, stem, save_dir, show, log, file_case)
            elif ptype == "2d":
                plot_2d(data, meta, stem, save_dir, show, log, file_case)
                if plot3d:
                    plot_3d(data, meta, stem, save_dir, show, log, file_case)
            elif ptype == "ts":
                plot_ts(data, meta, stem, save_dir, show, file_case)
            elif ptype == "ts2":
                plot_ts2(data, meta, stem, save_dir, show, file_case)
            else:
                print("    (skipped — unknown type)")
        except Exception as exc:
            print(f"    Error: {exc}")

    if restrict is None and not steady_state:
        if show_pow:
            print("  [cmp]  power_vs_time")
            plot_power_combined(outdir, save_dir, show, casename, show_sc=show_sc)
            if show_sc:
                print("  [cmp]  power_sc_split_vs_time")
                plot_sc_power_split(outdir, save_dir, show, casename)
        print("  [cmp]  coulomb_log_all_vs_time")
        plot_coulomb_log(outdir, save_dir, show, casename, show_sc=show_sc)
        print("  [cmp]  timescales_vs_time")
        plot_timescales(outdir, save_dir, show, casename)
        if show_mom:
            print("  [cmp]  momentum_breakdown_vs_time")
            plot_momentum_breakdown(outdir, save_dir, show, casename, show_sc=show_sc)
            print("  [cmp]  momentum_balance_vs_time")
            plot_momentum_balance(outdir, save_dir, show, casename)

    if restrict is None:
        f_check = _outfile(outdir, "fsc_maxw_at_vpar0", casename)
        if f_check.exists():
            print("  [cmp]  fout_vs_maxw_at_vpar0")
            plot_fout_vs_maxw_at_vpar0(outdir, save_dir, show, log, casename)

    # One blocking call opens every window at once.  When several cases are
    # being plotted the caller defers it to the end, so the windows of all of
    # them appear together instead of case n+1 waiting on case n being closed.
    if show and not defer_show:
        plt.show()


# ---------------------------------------------------------------------------
# Movie from the VDF snapshots (--movie)
# ---------------------------------------------------------------------------
# The solver writes f every n_snap steps to vdf_snap_<step>[-<case>].txt
# (write_vdf_snapshot, time_comps_mod.f90): three columns vperp, vpar, f in the
# layout of fout.txt, with a '# time = <t>' header line.  A snapshot holds f
# only; the kinetic-energy density is rebuilt here exactly as analysis.f90
# builds Ekin.txt, which needs the ion mass number aa from the namelist.

_SNAP_RE = re.compile(r"^vdf_snap_(\d{6})(?:-(.*))?\.txt$")

# Movie-only options: (command-line name, args attribute) and their defaults.
# Any of them given on its own implies --movie (see main).
_MOVIE_OPTS     = (("--fps", "fps"), ("--movie-format", "movie_format"),
                   ("--movie-scale", "movie_scale"))
_MOVIE_DEFAULTS = {"fps": 5.0, "movie_format": "avi", "movie_scale": "fixed"}
_PMASS   = 1.6726e-27     # proton mass (kg), as in analysis.f90
_KEV_J   = 1.60218e-16    # keV in J, as in analysis.f90


def _read_nml_number(input_file: Path, key: str):
    """Value of a numeric namelist entry (Fortran 'd' exponents allowed), or None."""
    try:
        text = input_file.read_text(errors="replace")
    except Exception:
        return None
    m = re.search(rf"\b{key}\s*=\s*([+-]?[0-9.]+(?:[dDeE][+-]?\d+)?)", text)
    if not m:
        return None
    return float(m.group(1).replace("d", "e").replace("D", "e"))


def _find_namelist(outdir: Path, casename: str):
    """A namelist in *outdir* whose casename is *casename*, or None.

    Only the head of each candidate is read, since a run folder can hold
    thousands of output .txt files and a namelist starts with '&INPUT'.
    """
    for pat in ("*.txt", "*.dat", "*.nml", "*.in"):
        for p in sorted(outdir.glob(pat)):
            try:
                with open(p, errors="replace") as fh:
                    head = fh.read(4096)
            except Exception:
                continue
            if "&input" not in head.lower():
                continue
            if _read_casename_from_namelist(p) == casename:
                return p
    return None


def _snapshot_files(outdir: Path, casename: str) -> list:
    """(step, path) of the snapshots of *casename*, sorted by step."""
    out = []
    for p in outdir.glob("vdf_snap_*.txt"):
        m = _SNAP_RE.match(p.name)
        if m and (m.group(2) or "") == casename:
            out.append((int(m.group(1)), p))
    return sorted(out)


def _simpson_weights(x: np.ndarray) -> np.ndarray:
    """Composite Simpson weights on an arbitrary grid -- ncint.f90 simpson_weights."""
    n = len(x); w = np.zeros(n); i = 0
    while i + 2 <= n - 1:
        h0, h1 = x[i+1] - x[i], x[i+2] - x[i+1]; hs = h0 + h1
        w[i]   += hs / 6 * (2 - h1 / h0)
        w[i+1] += hs**3 / (6 * h0 * h1)
        w[i+2] += hs / 6 * (2 - h0 / h1)
        i += 2
    if i < n - 1:
        h0 = x[i+1] - x[i]; w[i] += h0 / 2; w[i+1] += h0 / 2
    return w


def _read_snapshot(path: Path):
    """(time, vperp, vpar, f[i_vperp, j_vpar]) of one snapshot file."""
    with open(path) as fh:
        head = fh.readline()
    m = re.search(r"time\s*=\s*([-+0-9.EeDd]+)", head)
    time = float(m.group(1).replace("D", "E").replace("d", "e")) if m else float("nan")
    data = np.loadtxt(path, comments="#")
    vperp, vpar = np.unique(data[:, 0]), np.unique(data[:, 1])
    return time, vperp, vpar, data[:, 2].reshape(len(vperp), len(vpar))


def _ekin_map(f, vperp, vpar, aa: float) -> np.ndarray:
    """Kinetic-energy density (keV) as analysis.f90 writes Ekin.txt."""
    jac  = 2 * np.pi * vperp[:, None] * np.ones_like(f)
    mod0 = (np.outer(_simpson_weights(vperp), _simpson_weights(vpar)) * f * jac).sum()
    v2   = vperp[:, None]**2 + vpar[None, :]**2
    return f * v2 * jac / mod0 * 0.5 * _PMASS * aa / _KEV_J


def _ffmpeg_available() -> bool:
    """True if matplotlib can reach an ffmpeg binary.

    An ffmpeg on the PATH is used as is.  Otherwise the binary bundled with the
    imageio-ffmpeg package (pip install imageio-ffmpeg) is handed to matplotlib,
    which is what makes AVI and MP4 work on a machine without ffmpeg installed.
    """
    from matplotlib import animation
    if animation.writers.is_available("ffmpeg"):
        return True
    try:
        import imageio_ffmpeg
        plt.rcParams["animation.ffmpeg_path"] = imageio_ffmpeg.get_ffmpeg_exe()
    except Exception:
        return False
    return animation.writers.is_available("ffmpeg")


def make_movie(outdir: Path, casename: str, namelist, save_dir, log: bool,
               fps: float = 5.0, fmt: str = "avi", scale: str = "fixed") -> None:
    """Animate f and the kinetic-energy density over the snapshots of one case.

    Checks n_snap in the namelist first: with n_snap = 0 the solver wrote no
    snapshots and there is nothing to animate.  One frame per snapshot.  With
    *scale* 'fixed' the colour scales are fixed over the whole movie, so frames
    compare; with 'frame' each frame is scaled to its own maximum.  *fmt* is avi
    (MPEG-4 Part 2, plays in VLC and Windows Media Player), mp4 (H.264) or gif
    (Pillow); avi and mp4 need ffmpeg and fall back to gif without it.
    """
    from matplotlib import animation

    label = casename or "(no casename)"
    if namelist is None:
        namelist = _find_namelist(outdir, casename)
    if namelist is None or not Path(namelist).exists():
        print(f"  [movie] {label}: no namelist found for this case in {outdir};"
              f" pass it with --namelist FILE.  Movie skipped.")
        return
    namelist = Path(namelist)
    n_snap = _read_nml_number(namelist, "n_snap")
    if not n_snap:
        print(f"  [movie] {label}: n_snap = 0 in {namelist.name} -- the solver wrote no "
              f"snapshots.  Set n_snap > 0 and rerun.  Movie skipped.")
        return
    aa = _read_nml_number(namelist, "aa")
    if aa is None:
        print(f"  [movie] {label}: no 'aa' in {namelist.name}; Ekin needs the ion mass "
              f"number.  Movie skipped.")
        return
    snaps = _snapshot_files(outdir, casename)
    if not snaps:
        print(f"  [movie] {label}: n_snap = {int(n_snap)} but no vdf_snap_*.txt files for "
              f"this case in {outdir}.  Movie skipped.")
        return
    if len(snaps) < 2:
        print(f"  [movie] {label}: only one snapshot; a movie needs at least two.  Skipped.")
        return

    print(f"  [movie] {label}: {len(snaps)} snapshots (n_snap = {int(n_snap)}, "
          f"steps {snaps[0][0]}..{snaps[-1][0]}), aa = {aa:g}")
    frames = []
    for step, path in snaps:
        t, vperp, vpar, f = _read_snapshot(path)
        frames.append((step, t, f, _ekin_map(f, vperp, vpar, aa)))

    # Colour scales.  'fixed': one scale for the whole movie, from the largest
    # value of any frame, so frames compare directly.  'frame': each frame is
    # scaled to its own maximum, so its structure stays visible whatever its
    # amplitude; the colour bars and the panel titles follow the frame.
    per_frame = (scale == "frame")
    fmax_all = max(fr[2].max() for fr in frames)
    fmin_all = min(fr[2].min() for fr in frames)
    emax_all = max(fr[3].max() for fr in frames)

    def f_limits(f):
        hi = f.max() if per_frame else fmax_all
        hi = hi if hi > 0 else 1.0
        if log:
            # Eight decades.  Anything below the floor -- the far tail, and the
            # round-off of either sign around it -- is left blank rather than
            # drawn in the lowest colour, which would paint that noise as a
            # checkerboard.
            return hi * 1e-8, hi
        return min(0.0, f.min() if per_frame else fmin_all), hi

    def e_limit(e):
        hi = e.max() if per_frame else emax_all
        return hi if hi > 0 else 1.0

    def prep(z, lo):
        return np.where(z > lo, z, np.nan) if log else z

    fig, (axf, axe) = plt.subplots(1, 2, figsize=(13, 5.2))
    step0, t0, f0, e0 = frames[0]
    lo0, hi0 = f_limits(f0)
    fnorm = (mcolors.LogNorm if log else mcolors.Normalize)(vmin=lo0, vmax=hi0)
    enorm = mcolors.Normalize(vmin=0.0, vmax=e_limit(e0))
    mf = axf.pcolormesh(vpar, vperp, prep(f0, lo0), norm=fnorm, cmap="rainbow", shading="nearest")
    me = axe.pcolormesh(vpar, vperp, e0, norm=enorm, cmap="rainbow", shading="nearest")
    tail = "  (scaled to each frame)" if per_frame else ""
    fig.colorbar(mf, ax=axf, label="f" + ("  (log)" if log else "") + tail)
    fig.colorbar(me, ax=axe, label="kinetic-energy density (keV)" + tail)
    for ax in (axf, axe):
        ax.set_xlabel("v∥ (m/s)"); ax.set_ylabel("v⊥ (m/s)")
        if _XRANGE is not None: ax.set_xlim(_XRANGE)
        if _YRANGE is not None: ax.set_ylim(_YRANGE)
    sup = fig.suptitle("")
    fig.tight_layout(rect=(0, 0, 1, 0.94))

    def draw(k):
        step, t, f, e = frames[k]
        lo, hi = f_limits(f)
        mf.set_clim(lo, hi)
        mf.set_array(prep(f, lo).ravel())
        me.set_clim(0.0, e_limit(e))
        me.set_array(e.ravel())
        # The frame's own maximum in the titles: with --movie-scale frame it is
        # the top of the colour bar, with fixed scales it shows the amplitude.
        axf.set_title(f"VDF  f(v⊥, v∥)      max = {f.max():.3g}")
        axe.set_title(f"Kinetic energy      max = {e.max():.3g} keV")
        sup.set_text(f"{label}    t = {t:.4g} s    (step {step}, frame {k + 1}/{len(frames)})")
        return mf, me, sup

    anim = animation.FuncAnimation(fig, draw, frames=len(frames), blit=False)
    dest = Path(save_dir) if save_dir is not None else outdir
    dest.mkdir(parents=True, exist_ok=True)
    stem = f"movie-{casename}" if casename else "movie"
    if fmt in ("avi", "mp4") and not _ffmpeg_available():
        print(f"    no ffmpeg found (pip install imageio-ffmpeg provides one): "
              f"writing a GIF instead of {fmt.upper()}")
        fmt = "gif"
    if fmt == "gif":
        out = (dest / f"{stem}.gif").resolve()
        anim.save(out, writer=animation.PillowWriter(fps=fps), dpi=90)
    else:
        # yuv420p for player compatibility; it needs even frame dimensions,
        # which the scale filter guarantees whatever the figure size and dpi.
        codec = "mpeg4" if fmt == "avi" else "libx264"
        extra = ["-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2", "-pix_fmt", "yuv420p"]
        if fmt == "avi":
            # MPEG-4 quality scale (1 best .. 31), and the XVID fourcc: ffmpeg
            # labels the stream FMP4 by default, which some Windows players
            # do not recognise, whereas XVID-tagged MPEG-4 plays natively.
            extra += ["-q:v", "3", "-vtag", "xvid"]
        out = (dest / f"{stem}.{fmt}").resolve()
        anim.save(out, writer=animation.FFMpegWriter(fps=fps, codec=codec,
                                                     extra_args=extra), dpi=120)
    plt.close(fig)
    print(f"    -> {out}")


# ---------------------------------------------------------------------------
# Casename helpers
# ---------------------------------------------------------------------------

def _read_casename_from_namelist(input_file: Path) -> str:
    """Return the casename value from a Fortran namelist, or '' if absent/unset."""
    try:
        text = input_file.read_text(errors="replace")
        m = re.search(r'\bcasename\s*=\s*["\']([^"\']*)["\']', text, re.IGNORECASE)
        if m:
            return m.group(1).strip()
    except Exception:
        pass
    return ""


def _read_ntimes_from_namelist(input_file: Path) -> int:
    """Return ntimes from a Fortran namelist, or -1 if not found."""
    try:
        text = input_file.read_text(errors="replace")
        m = re.search(r'\bntimes\s*=\s*([+-]?\d+)', text, re.IGNORECASE)
        if m:
            return int(m.group(1))
    except Exception:
        pass
    return -1


def _read_isc_from_namelist(input_file: Path) -> int:
    """Return isc from a Fortran namelist, or 0 if not found."""
    try:
        text = input_file.read_text(errors="replace")
        m = re.search(r'\bisc\s*=\s*([+-]?\d+)', text, re.IGNORECASE)
        if m:
            return int(m.group(1))
    except Exception:
        pass
    return 0


def _read_iplot_pow_from_namelist(input_file: Path) -> int:
    """Return iplot_pow from a Fortran namelist; default -1 (enabled)."""
    try:
        text = input_file.read_text(errors="replace")
        m = re.search(r'\biplot_pow\s*=\s*([+-]?\d+)', text, re.IGNORECASE)
        if m:
            return int(m.group(1))
    except Exception:
        pass
    return -1


def _read_iplot_mom_from_namelist(input_file: Path) -> int:
    """Return iplot_mom from a Fortran namelist; default 0 (disabled)."""
    try:
        text = input_file.read_text(errors="replace")
        m = re.search(r'\biplot_mom\s*=\s*([+-]?\d+)', text, re.IGNORECASE)
        if m:
            return int(m.group(1))
    except Exception:
        pass
    return 0


def _case_from_stem(stem: str, default: str = "") -> str:
    """Casename suffix of a single output-file stem, or *default* if none.

    'Teff_vs_time-JET-...-Grid3-test3' -> 'JET-...-Grid3-test3'.
    Lets each file carry its own casename when several cases are plotted
    together (e.g. via --files), so titles match the file they describe.
    """
    for key in sorted(FILE_META, key=len, reverse=True):
        if stem.startswith(key + "-"):
            return stem[len(key) + 1:]
    return default


def _detect_casename(outdir: Path, names=None) -> str:
    """Infer casename from output files by matching known FILE_META stems.

    If *names* is given (a list of bare filenames), search those instead of
    scanning the whole directory — used when --files restricts the plot set.
    """
    if names is not None:
        stems = sorted(Path(n).stem for n in names)
    else:
        stems = (p.stem for p in sorted(outdir.glob("*.txt")))
    for stem in stems:
        case = _case_from_stem(stem)
        if case:
            return case
    return ""


# ---------------------------------------------------------------------------
# Solver runner
# ---------------------------------------------------------------------------

def run_solver(exe: Path, input_file: Path, run_dir: Path,
               out_file: Path = None) -> int:
    print(f"Running : {exe}")
    print(f"Input   : {input_file}")
    print(f"Workdir : {run_dir}")
    if out_file is not None:
        print(f"Stdout  : {out_file}")
    with open(input_file) as fin:
        if out_file is not None:
            with open(out_file, "w") as fout:
                result = subprocess.run([str(exe)], stdin=fin, stdout=fout,
                                        cwd=run_dir)
        else:
            result = subprocess.run([str(exe)], stdin=fin, cwd=run_dir)
    return result.returncode


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _add_common(p: argparse.ArgumentParser, multi: bool = False) -> None:
    p.add_argument("--save",        type=Path, default=None, dest="save_dir", metavar="DIR",
                   help="save PNG files to DIR")
    p.add_argument("--show",        action="store_true",
                   help="display plots interactively")
    p.add_argument("--logf", "--log", action="store_true", dest="log",
                   help="logarithmic scale for the plotted QUANTITY: the "
                        "colour scale of a 2D map, the z of a 3D surface, the "
                        "y of a 1D profile.  --log is an accepted alias")
    p.add_argument("--logx",        action="store_true",
                   help="logarithmic x-axis (left linear where the data reach "
                        "zero, e.g. v_par or a trace starting at t=0)")
    p.add_argument("--logy",        action="store_true",
                   help="logarithmic y-axis; on a 1D profile this is the same "
                        "as --logf")
    if multi:
        # --cases matches the spelling used by the compare subcommand.
        # --casename is kept as an alias: it appears in existing scripts.
        p.add_argument("--cases", "--casename", nargs="+", default=None,
                       dest="cases", metavar="CASE",
                       help="case(s) to plot, one set of figures each; selects "
                            "the files and labels them (auto-detected if "
                            "omitted).  --casename is an accepted alias")
    else:
        p.add_argument("--casename", default="", metavar="STR",
                       help="case label for the output filenames and titles "
                            "(read from the namelist if omitted)")
    p.add_argument("--files",       nargs="+", default=None, metavar="F",
                   help="plot only these file types (stem key, glob or full "
                        "filename); still restricted to the case(s) when given")
    p.add_argument("--steady-state", action="store_true", dest="steady_state",
                   help="skip time-trace plots (for ntimes=0 runs)")
    p.add_argument("--no-sc",        action="store_true", dest="no_sc",
                   help="suppress self-collision power plots (for isc=0 runs)")
    p.add_argument("--no-pow",       action="store_true", dest="no_pow",
                   help="suppress power vs time plots (for iplot_pow=0 runs)")
    p.add_argument("--no-mom",       action="store_true", dest="no_mom",
                   help="suppress momentum vs time plots (for iplot_mom=0 runs)")
    p.add_argument("--3d",           action="store_true", dest="plot3d",
                   help="add 3D surface plots for 2D distribution files")
    p.add_argument("--movie",        action="store_true",
                   help="make a movie of f and the kinetic-energy density from "
                        "the VDF snapshots (vdf_snap_*.txt) instead of the usual "
                        "plots; needs n_snap > 0 in the namelist")
    # The three options below only make sense for a movie, so giving any of
    # them implies --movie (main).  Their defaults are therefore None here, to
    # tell an option the user typed from an untouched one, and are filled in
    # from _MOVIE_DEFAULTS afterwards.
    p.add_argument("--fps",          type=float, default=None, metavar="N",
                   help="frames per second of the movie (default 5); "
                        "implies --movie")
    p.add_argument("--movie-format", choices=("avi", "mp4", "gif"), default=None,
                   dest="movie_format",
                   help="file format of the movie (default avi; avi and mp4 "
                        "need ffmpeg, e.g. from pip install imageio-ffmpeg, and "
                        "fall back to gif without it); implies --movie")
    p.add_argument("--movie-scale", choices=("fixed", "frame"), default=None,
                   dest="movie_scale",
                   help="colour scales of the movie: 'fixed' over the whole "
                        "movie, from its largest value, so frames compare "
                        "(default); 'frame' scaled to each frame's own maximum; "
                        "implies --movie")
    if multi:
        p.add_argument("--namelist", type=Path, default=None, metavar="FILE",
                       help="namelist of the run, read by --movie for n_snap and "
                            "aa (default: the namelist in DIR whose casename "
                            "matches the case)")


_COMPARE_ALIASES = {
    # power_sc_split_vs_time is a plot-only name; data lives in power_coll_self_vs_time
    "power_sc_split_vs_time": "power_coll_self_vs_time",
    # coulomb_log_all_vs_time is a plot-only composite; for compare, default to
    # the self-collision log (the one that differs between isc models). Use
    # coulomb_log_vs_time explicitly for the background-ion logs.
    "coulomb_log_all_vs_time": "coulomb_log_self_vs_time",
}


def _species_file_meta(outdir: Path, cases: list) -> dict:
    """Metadata for the per-background-species files, discovered from disk.

    The Fortran writes one power and one momentum file per bulk species, with
    the species index built into the stem: 'power_coll_ion 1_vs_time'.  They are
    therefore not static FILE_META keys, and in plot mode a composite function
    draws them.  Compare mode has no such composite, so without this they were
    invisible to it however they were spelled on the command line.

    Both the spaced and the unspaced spelling are probed, mirroring
    plot_power_coll, and only stems that exist are registered.
    """
    meta = {}
    for kind, ylabel, what in (("power",    "Power density (MW·m⁻³)",
                                "collisional power density"),
                               ("momentum", "Momentum transfer rate (N·m⁻³)",
                                "momentum transfer")):
        for ib in range(1, 10):
            for stem in (f"{kind}_coll_ion {ib}_vs_time",
                         f"{kind}_coll_ion{ib}_vs_time"):
                if any(_outfile(outdir, stem, c).exists() for c in cases):
                    meta[stem] = {"ptype": "ts", "ylabel": ylabel,
                                  "title": f"Ion species {ib} {what}"}
                    break
    return meta


def compare_directory(outdir: Path, cases: list, save_dir, show: bool,
                      log: bool, restrict=None) -> None:
    """Overlay same-type output files from multiple casenames on shared axes.

    For each known stem that has at least one matching file, one figure is
    produced with one curve per case (ts/ts2) or one profile per case (1d).
    2D distribution files are skipped — they are not meaningful to superpose.
    """
    if save_dir is not None:
        Path(save_dir).mkdir(parents=True, exist_ok=True)

    # FILE_META plus the per-species files, which carry their index in the stem
    all_meta = {**FILE_META, **_species_file_meta(outdir, cases)}

    if restrict:
        restrict = [_COMPARE_ALIASES.get(r, r) for r in restrict]
        restrict = _rejoin_split_stems(restrict, all_meta)

    # Collect (case, data) pairs for every known stem key
    stem_cases: dict = {}
    matched_names: set = set()
    for case in cases:
        for key in all_meta:
            path = _outfile(outdir, key, case)
            if not path.exists():
                continue
            # same matching as plot mode: globs, full filenames, or bare stems
            if restrict and not any(_restrict_match(path.name, pat)
                                    for pat in restrict):
                continue
            data = _load_auto(path)
            if data is None or data.shape[0] < 2:
                continue
            matched_names.add(path.name)
            stem_cases.setdefault(key, []).append((case, data))

    # Say so when an entry matched nothing, rather than silently plotting less
    # than was asked for.
    if restrict:
        missed = [r for r in restrict
                  if not any(_restrict_match(n, r) for n in matched_names)]
        if missed:
            known = ", ".join(sorted(k for k in all_meta
                                     if k not in _SKIP_STEMS)[:6])
            print(f"  Warning: no files matched {', '.join(repr(m) for m in missed)}"
                  f"\n           for case(s) {', '.join(cases)}."
                  f"\n           Known stems include: {known}, ...")

    if not stem_cases:
        print("  (no matching files found for the given casenames)")
        return

    _lstyles = ["-", "--", ":", "-."]
    _markers = ["o", "s", "^", "D", "v", "*", "P", "X"]

    for key in sorted(stem_cases):
        entries = stem_cases[key]
        meta   = all_meta[key]
        ptype  = meta.get("ptype", "")
        if ptype not in ("ts", "ts2", "1d"):
            continue

        labels    = meta.get("labels", [])
        skip_cols = set(meta.get("skip_cols", []))
        print(f"  [cmp]  {key}")
        fig, ax = plt.subplots(figsize=(8, 5))
        plotted = False

        for ci, (case, data) in enumerate(entries):
            color  = _PALETTE[ci % len(_PALETTE)]
            lstyle = _lstyles[ci % len(_lstyles)]   # distinguish each case by line style
            marker = _markers[ci % len(_markers)]   # ... and by marker shape
            # ~12 markers spread along the curve (not one per data point)
            mevery = max(1, int(round(data.shape[0] / 12.0)))
            if ptype in ("ts", "ts2"):
                ncols = data.shape[1]
                if ncols == 2:
                    ax.plot(data[:, 0], data[:, 1], color=color, linestyle=lstyle,
                            marker=marker, markevery=mevery, markersize=5,
                            linewidth=1.5, label=case)
                else:
                    # Multi-column file: columns within a case differ by line
                    # style; cases are told apart by colour + marker shape.
                    for j in range(1, ncols):
                        if j in skip_cols:
                            continue
                        col_lbl = labels[j - 1] if j - 1 < len(labels) else f"col{j}"
                        ax.plot(data[:, 0], data[:, j], color=color,
                                linestyle=_lstyles[(j - 1) % len(_lstyles)],
                                marker=marker, markevery=mevery, markersize=5,
                                linewidth=1.5, label=f"{case}  [{col_lbl}]")
                ax.set_xlabel("Time (s)")
                ax.set_xlim(left=0)
            elif ptype == "1d":
                ax.plot(data[:, 0], data[:, 1], color=color, linestyle=lstyle,
                        marker=marker, markevery=mevery, markersize=5,
                        linewidth=1.5, label=case)
                ax.set_xlabel(meta.get("xlabel", ""))
                if log:
                    ax.set_yscale("log")
            plotted = True

        if not plotted:
            plt.close(fig)
            continue

        ax.set_ylabel(meta.get("ylabel", ""))
        ax.set_title(_title(meta, key, ""))
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        _finish(fig, f"{key}-compare", save_dir, show)

    if show:
        plt.show()


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = p.add_subparsers(dest="command", required=True)

    run_p = sub.add_parser("run", help="execute solver then plot outputs")
    run_p.add_argument("exe",        type=Path, help="path to FP2D_QLRF_NL.exe")
    run_p.add_argument("input_file", type=Path, help="Fortran namelist input (.dat)")
    run_p.add_argument("--outdir",   type=Path, default=None, metavar="DIR",
                       help="solver working directory (default: folder of input_file)")
    run_p.add_argument("--out",      type=Path, default=None, metavar="FILE",
                       help="redirect solver stdout (Fortran write(*,*)) to FILE")
    _add_common(run_p)

    plot_p = sub.add_parser("plot", help="plot from an existing output directory")
    plot_p.add_argument("outdir", type=Path, help="directory containing .txt output files")
    _add_common(plot_p, multi=True)
    plot_p.add_argument("--xrange", default=None, metavar="xmin:xmax",
                        help="zoom the x-axis of every plot to [xmin, xmax]")
    plot_p.add_argument("--yrange", default=None, metavar="ymin:ymax",
                        help="zoom the y-axis of every plot to [ymin, ymax] "
                             "(suppresses the automatic y-rescale that "
                             "--xrange performs)")

    cmp_p = sub.add_parser("compare",
                            help="overlay same-type outputs from multiple cases")
    cmp_p.add_argument("outdir", type=Path,
                       help="directory containing the output .txt files")
    cmp_p.add_argument("--cases", nargs="+", required=True, metavar="CASE",
                       help="casenames to overlay (e.g. --cases Lin NLSC)")
    cmp_p.add_argument("--save",  type=Path, default=None, dest="save_dir",
                       metavar="DIR", help="save PNG files to DIR")
    cmp_p.add_argument("--show",  action="store_true",
                       help="display plots interactively")
    cmp_p.add_argument("--logf", "--log", action="store_true", dest="log",
                       help="logarithmic y-scale for 1D profile plots "
                            "(--log is an accepted alias)")
    cmp_p.add_argument("--logx",  action="store_true",
                       help="logarithmic x-axis")
    cmp_p.add_argument("--logy",  action="store_true",
                       help="logarithmic y-axis")
    cmp_p.add_argument("--files", nargs="+", default=None, metavar="F",
                       help="restrict to these file types (stem key, e.g. anisotropy_vs_time)"
                            " or full filenames")
    cmp_p.add_argument("--xrange", default=None, metavar="xmin:xmax",
                       help="zoom the x-axis of every plot to [xmin, xmax]")
    cmp_p.add_argument("--yrange", default=None, metavar="ymin:ymax",
                       help="zoom the y-axis of every plot to [ymin, ymax]")

    return p


def main(argv=None):
    global _XRANGE, _YRANGE, _LOGX, _LOGY
    args = build_parser().parse_args(argv)

    # Default to interactive display when no save directory is given
    if not args.show and args.save_dir is None:
        args.show = True

    # axis zoom (plot / compare only; run defines neither)
    if getattr(args, "xrange", None):
        _XRANGE = _parse_range(args.xrange, "--xrange")
    if getattr(args, "yrange", None):
        _YRANGE = _parse_range(args.yrange, "--yrange")

    # logarithmic axes; --logf travels separately as the `log` argument
    _LOGX = bool(getattr(args, "logx", False))
    _LOGY = bool(getattr(args, "logy", False))

    # A movie option on its own means a movie: without this, e.g.
    # "plot . --movie-scale frame" silently produced the ordinary figures.
    if hasattr(args, "movie"):
        given = [opt for opt, dest in _MOVIE_OPTS if getattr(args, dest) is not None]
        if given and not args.movie:
            args.movie = True
            print(f"{', '.join(given)} given: making a movie (--movie implied).")
        for dest, default in _MOVIE_DEFAULTS.items():
            if getattr(args, dest) is None:
                setattr(args, dest, default)

    # ----------------------------------------------------------------
    # compare command — handled entirely here, then return
    # ----------------------------------------------------------------
    if args.command == "compare":
        outdir = args.outdir.resolve()
        if not outdir.is_dir():
            sys.exit(f"Error: output directory not found: {outdir}")
        if not args.show:
            plt.switch_backend("Agg")
        print(f"\nComparing {len(args.cases)} case(s) in: {outdir}")
        for c in args.cases:
            print(f"  • {c}")
        compare_directory(outdir, args.cases,
                          save_dir=args.save_dir,
                          show=args.show,
                          log=args.log,
                          restrict=args.files)
        print("Done.")
        return

    if args.command == "run":
        exe        = args.exe.resolve()
        input_file = args.input_file.resolve()
        run_dir    = (args.outdir or input_file.parent).resolve()
        if not exe.exists():
            sys.exit(f"Error: executable not found: {exe}")
        if not input_file.exists():
            sys.exit(f"Error: input file not found: {input_file}")
        if not args.casename:
            args.casename = _read_casename_from_namelist(input_file)
            if args.casename:
                print(f"Casename (from namelist): {args.casename}")
        if not args.steady_state:
            ntimes = _read_ntimes_from_namelist(input_file)
            if ntimes == 0:
                args.steady_state = True
                print("Steady-state run (ntimes=0): time traces will be skipped.")
        if not args.no_sc:
            isc = _read_isc_from_namelist(input_file)
            if isc == 0:
                args.no_sc = True
                print("No self-collisions (isc=0): self-collision power plots will be skipped.")
        if not args.no_pow:
            if _read_iplot_pow_from_namelist(input_file) == 0:
                args.no_pow = True
                print("iplot_pow=0: power vs time plots will be skipped.")
        if not args.no_mom:
            if _read_iplot_mom_from_namelist(input_file) == 0:
                args.no_mom = True
                print("iplot_mom=0: momentum vs time plots will be skipped.")
        if args.movie and not _read_nml_number(input_file, "n_snap"):
            # Checked before the solve: without snapshots the run would be wasted.
            sys.exit("Error: --movie needs VDF snapshots, but n_snap = 0 in "
                     f"{input_file.name}.  Set n_snap > 0 in the namelist.")
        rc = run_solver(exe, input_file, run_dir, out_file=args.out)
        if rc != 0:
            sys.exit(f"Solver exited with code {rc} — aborting (no plots written).")
        outdir = run_dir
    else:
        outdir = args.outdir.resolve()

    if not outdir.is_dir():
        sys.exit(f"Error: output directory not found: {outdir}")

    # 'plot' takes one or more --cases; 'run' carries the single --casename it
    # read from the namelist.  Cases the user typed restrict which files are
    # plotted; one we guessed only labels them (see plot_directory).
    cases = list(getattr(args, "cases", None) or [])
    if not cases and getattr(args, "casename", ""):
        cases = [args.casename]
    strict_case = bool(cases)
    if not cases:
        # Auto-detect the casename from the output filenames. With --files, search
        # only those names; otherwise scan the whole directory, so a plain
        # `plot <dir>` works for a single-case folder without requiring --cases.
        detected = _detect_casename(outdir, names=args.files)
        if detected:
            print(f"Casename (auto-detected): {detected}")
        cases = [detected]          # may be "" — the no-casename file set

    if args.movie:
        # Movie only: the snapshots replace the usual figures.  For 'run' the
        # namelist is the input file; for 'plot' it is --namelist or found in
        # outdir by casename (make_movie).
        plt.switch_backend("Agg")
        namelist = input_file if args.command == "run" else args.namelist
        print(f"\nMaking movie(s) from the snapshots in: {outdir}")
        for case in cases:
            make_movie(outdir, case, namelist, save_dir=args.save_dir,
                       log=args.log, fps=args.fps, fmt=args.movie_format,
                       scale=args.movie_scale)
        print("Done.")
        return

    if not args.show:
        plt.switch_backend("Agg")

    print(f"\nPlotting output files in: {outdir}")
    for n, case in enumerate(cases):
        if len(cases) > 1:
            print(f"\n=== case {n + 1}/{len(cases)}: {case} ===")
        plot_directory(
            outdir,
            save_dir=args.save_dir,
            show=args.show,
            log=args.log,
            casename=case,
            restrict=args.files,
            strict_case=strict_case,
            steady_state=args.steady_state,
            show_sc=not args.no_sc,
            show_pow=not args.no_pow,
            show_mom=not args.no_mom,
            plot3d=args.plot3d,
            # hold the blocking show() until every case has been drawn
            defer_show=(n < len(cases) - 1),
        )
    print("Done.")


if __name__ == "__main__":
    main()
