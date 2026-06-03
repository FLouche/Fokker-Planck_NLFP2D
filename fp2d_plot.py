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
  --log           Logarithmic colour/y-scale for distribution functions.
  --3d            Add 3D surface plots for 2D distribution files.
  --casename STR  Case label appended to every plot title.
  --files F ...   Plot only these filenames (basenames, e.g. fout.txt).

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
    # Simple time series -------------------------------------------------------
    "density_vs_time":         {"ptype": "ts",  "ylabel": "Density (m⁻³)",
                                "title": "Particle density vs time"},
    "anisotropy_vs_time":      {"ptype": "ts",  "ylabel": "Anisotropy (%)",
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
                                "title": "Self-collision power density",
                                "labels": ["total", "⊥", "∥"]},
    # Two-curve time series ----------------------------------------------------
    "energy_vs_time":               {"ptype": "ts2", "ylabel": "Energy (keV)",
                                     "title": "Kinetic energy vs time",
                                     "labels": ["E_total", "E_⊥"]},
    # Effective temperature
    "Teff_vs_time":                 {"ptype": "ts",  "ylabel": "T_eff (keV)",
                                     "title":  "Effective temperature vs time"},
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
               "coulomb_log_vs_time", "coulomb_log_self_vs_time"}

_PALETTE = ["#8B1A1A", "#1A1A8B", "#1A8B1A", "#8B8B1A", "#8B1A8B", "#1A8B8B"]


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

def _finish(fig, stem: str, save_dir, show: bool) -> None:
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

    if log and np.any(Z > 0):
        pos    = Z[Z > 0]
        norm   = mcolors.LogNorm(vmin=float(pos.min()), vmax=zmax)
        levels = np.logspace(np.log10(float(pos.min())), np.log10(zmax), 20)
    else:
        norm   = None
        levels = np.linspace(zmin, zmax, 20)

    fig, ax = plt.subplots(figsize=(8, 6))
    cf   = ax.contourf(vpar_u, vperp_u, Z, levels=levels, norm=norm, cmap="rainbow")
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
            t, y = _ts_col(f)
            if t is not None:
                ax.plot(t, y, color=_PALETTE[4], linewidth=1.5, linestyle="--", label="self")
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
            t, y = _ts_col(f)
            if t is not None:
                ax_top.plot(t, y, color=_PALETTE[4], linewidth=1.5, linestyle="--", label="self")
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


# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------

def plot_directory(outdir: Path, save_dir, show: bool, log: bool,
                   casename: str, restrict=None, steady_state: bool = False,
                   show_sc: bool = True, plot3d: bool = False) -> None:
    if save_dir is not None:
        Path(save_dir).mkdir(parents=True, exist_ok=True)

    txt_files = sorted(outdir.glob("*.txt"))
    if restrict:
        # Each entry in restrict may be an exact filename or a glob pattern
        txt_files = [f for f in txt_files
                     if any(fnmatch.fnmatch(f.name, pat) for pat in restrict)]
    else:
        txt_files = [f for f in txt_files if _matches_casename(f, casename)]

    for path in txt_files:
        stem = path.stem
        if any(stem.startswith(s) for s in _SKIP_STEMS):
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
        try:
            if ptype == "1d":
                plot_1d(data, meta, stem, save_dir, show, log, casename)
            elif ptype == "2d":
                plot_2d(data, meta, stem, save_dir, show, log, casename)
                if plot3d:
                    plot_3d(data, meta, stem, save_dir, show, log, casename)
            elif ptype == "ts":
                plot_ts(data, meta, stem, save_dir, show, casename)
            elif ptype == "ts2":
                plot_ts2(data, meta, stem, save_dir, show, casename)
            else:
                print("    (skipped — unknown type)")
        except Exception as exc:
            print(f"    Error: {exc}")

    if restrict is None and not steady_state:
        print("  [cmp]  power_vs_time")
        plot_power_combined(outdir, save_dir, show, casename, show_sc=show_sc)
        print("  [cmp]  coulomb_log_all_vs_time")
        plot_coulomb_log(outdir, save_dir, show, casename, show_sc=show_sc)
        print("  [cmp]  momentum_breakdown_vs_time")
        plot_momentum_breakdown(outdir, save_dir, show, casename, show_sc=show_sc)
        print("  [cmp]  momentum_balance_vs_time")
        plot_momentum_balance(outdir, save_dir, show, casename)

    if show:
        plt.show()  # single blocking call — all windows open simultaneously


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
        for key in sorted(FILE_META, key=len, reverse=True):
            if stem.startswith(key + "-"):
                return stem[len(key) + 1:]
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

def _add_common(p: argparse.ArgumentParser) -> None:
    p.add_argument("--save",        type=Path, default=None, dest="save_dir", metavar="DIR",
                   help="save PNG files to DIR")
    p.add_argument("--show",        action="store_true",
                   help="display plots interactively")
    p.add_argument("--log",         action="store_true",
                   help="logarithmic scale for distribution functions")
    p.add_argument("--casename",    default="", metavar="STR",
                   help="case label for plot titles")
    p.add_argument("--files",       nargs="+", default=None, metavar="F",
                   help="plot only these filenames (basenames)")
    p.add_argument("--steady-state", action="store_true", dest="steady_state",
                   help="skip time-trace plots (for ntimes=0 runs)")
    p.add_argument("--no-sc",        action="store_true", dest="no_sc",
                   help="suppress self-collision power plots (for isc=0 runs)")
    p.add_argument("--3d",           action="store_true", dest="plot3d",
                   help="add 3D surface plots for 2D distribution files")


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
    _add_common(plot_p)

    return p


def main(argv=None):
    args = build_parser().parse_args(argv)

    # Default to interactive display when no save directory is given
    if not args.show and args.save_dir is None:
        args.show = True

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
        rc = run_solver(exe, input_file, run_dir, out_file=args.out)
        if rc != 0:
            print(f"Warning: solver exited with code {rc}", file=sys.stderr)
        outdir = run_dir
    else:
        outdir = args.outdir.resolve()

    if not outdir.is_dir():
        sys.exit(f"Error: output directory not found: {outdir}")

    if not args.casename and args.files:
        args.casename = _detect_casename(outdir, names=args.files)
        if args.casename:
            print(f"Casename (auto-detected): {args.casename}")

    if not args.show:
        plt.switch_backend("Agg")

    print(f"\nPlotting output files in: {outdir}")
    plot_directory(
        outdir,
        save_dir=args.save_dir,
        show=args.show,
        log=args.log,
        casename=args.casename,
        restrict=args.files,
        steady_state=args.steady_state,
        show_sc=not args.no_sc,
        plot3d=args.plot3d,
    )
    print("Done.")


if __name__ == "__main__":
    main()
