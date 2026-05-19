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
  --casename STR  Case label appended to every plot title.
  --files F ...   Plot only these filenames (basenames, e.g. fout.txt).

run-only options
----------------
  --outdir DIR    Solver working directory (default: folder of <input.dat>).

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
    "power_coll_self_vs_time": {"ptype": "ts",
                                "ylabel": "Power density (MW·m⁻³)",
                                "title": "Self-collision power density"},
    # Two-curve time series ----------------------------------------------------
    "energy_vs_time":          {"ptype": "ts2", "ylabel": "Energy (keV)",
                                "title": "Kinetic energy vs time",
                                "labels": ["E_total", "E_⊥"]},
}

# Files to skip (unusual format or not useful for plotting)
_SKIP_STEMS = {"RF_dirac", "fstix"}

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
        out = Path(save_dir) / f"{stem}.png"
        fig.savefig(out, dpi=150, bbox_inches="tight")
        print(f"    -> {out.name}")
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


# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------

def plot_directory(outdir: Path, save_dir, show: bool, log: bool,
                   casename: str, restrict=None, steady_state: bool = False,
                   show_sc: bool = True) -> None:
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
            elif ptype == "ts":
                plot_ts(data, meta, stem, save_dir, show, casename)
            elif ptype == "ts2":
                plot_ts2(data, meta, stem, save_dir, show, casename)
            else:
                print("    (skipped — unknown type)")
        except Exception as exc:
            print(f"    Error: {exc}")

    if restrict is None and not steady_state:
        print("  [cmp]  power_coll_vs_time")
        plot_power_coll(outdir, save_dir, show, casename, show_sc=show_sc)
        print("  [cmp]  power_balance_vs_time")
        plot_power_balance(outdir, save_dir, show, casename)

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


def _detect_casename(outdir: Path) -> str:
    """Infer casename from output files by matching known FILE_META stems."""
    for path in sorted(outdir.glob("*.txt")):
        stem = path.stem
        for key in sorted(FILE_META, key=len, reverse=True):
            if stem.startswith(key + "-"):
                return stem[len(key) + 1:]
    return ""


# ---------------------------------------------------------------------------
# Solver runner
# ---------------------------------------------------------------------------

def run_solver(exe: Path, input_file: Path, run_dir: Path) -> int:
    print(f"Running : {exe}")
    print(f"Input   : {input_file}")
    print(f"Workdir : {run_dir}")
    with open(input_file) as fin:
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
        rc = run_solver(exe, input_file, run_dir)
        if rc != 0:
            print(f"Warning: solver exited with code {rc}", file=sys.stderr)
        outdir = run_dir
    else:
        outdir = args.outdir.resolve()

    if not outdir.is_dir():
        sys.exit(f"Error: output directory not found: {outdir}")

    if not args.casename:
        args.casename = _detect_casename(outdir)
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
    )
    print("Done.")


if __name__ == "__main__":
    main()
