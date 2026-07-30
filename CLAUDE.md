# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Git Workflow

After every meaningful unit of work — a bug fix, a new feature, a refactor, or any change that compiles and leaves the code in a coherent state — commit the changed source files and push to GitHub:

```
git add <changed files>
git commit -m "<concise description of what changed and why>"
git push
```

Commit messages should state **what** changed and **why** (the physics or numerical motivation), not just repeat the diff. Examples of good messages:
- `Fix Neumann BC weight normalisation in fd_stencil_2d at i=1`
- `Add self-collision term to timefp_7pt time loop`
- `Increase nnz_max estimate to accommodate denser stencil near boundaries`

Never batch unrelated changes into a single commit. Never commit build outputs, `.obj`/`.mod`/`.exe` files, or simulation result files — those are covered by `.gitignore`.

## Build and Run

**Build:** Open `FP2D_QLRF_NL/FP2D_QLRF_NL.sln` in Visual Studio and build the `Debug|x64` configuration. This uses the Intel `ifx` compiler with Intel MKL (sequential). The executable lands at `FP2D_QLRF_NL/x64/Debug/FP2D_QLRF_NL.exe`.

**Run:** The program reads a Fortran namelist (`&INPUT`) from stdin and writes output files to the working directory (i.e. `x64/Debug/` when launched from Visual Studio, or wherever the shell is):
```
FP2D_QLRF_NL.exe < <case>.dat
```
No example namelists are checked into the repo; supply your own `&INPUT` file (see the `namelist /INPUT/` declaration in `main-FP_Coll_2D.f90` for the full parameter list).

There is no test suite; correctness is verified by inspecting the output files (`density_vs_time.txt`, `energy_vs_time.txt`, `fout.txt`, etc.) and checking power/density balance diagnostics printed to stdout.

## Architecture

### The PDE
The code solves the 2D Fokker-Planck equation for the ion velocity distribution function f(v⊥, v∥) in cylindrical velocity space:

```
A·f + B·∂f/∂v⊥ + C·∂f/∂v∥ + D·∂²f/∂v⊥² + E·∂²f/∂v⊥∂v∥ + F·∂²f/∂v∥² = source
```

The six coefficient arrays (`all00`=A, `all10`=B, `all01`=C, `all20`=D, `all11`=E, `all02`=F) are assembled in `main-FP_Coll_2D.f90` as sums of three independently computed contributions, all stored in `shared_FPterms`:
- **Coulomb collisions** (`colin**`): computed by `cblin` per background species, summed over `nbulk` species
- **Quasi-linear RF** (`rf**`): computed by `qlrfterm` using complex Bessel functions; active when `irf=-1`
- **Self-collisions** (`sc**`): computed by `self_coll_max`; active when `isc≠0`

### Grid and Indexing
Grid is `nperp × npar` in (v⊥, v∥). The global 1D index used in all sparse matrix operations is:
```
ix = (i-1)*npar + j      ! i=v⊥ index (1..nperp), j=v∥ index (1..npar)
```
`index_mat(i,j)` and `index_mat_inv(ix,i,j)` in `func_index` (`functions.f90`) perform the mapping.

The v∥ grid is always uniform. The v⊥ grid has three modes controlled by `ising` (`shared_grid`):
- `ising=0`: uniform
- `ising=-1`: two-domain (fine for v⊥ < `vbound`, coarser above)
- `ising=+1`: quadratic spacing (more resolution near v⊥=0)

Boundary conditions are encoded in `fd_stencil_2d`:
- `i=1` (v⊥=0): Neumann, df/dv⊥=0 (axis symmetry)
- `i=nperp`, `j=1`, `j=npar`: Dirichlet, f=0

### Domain Sizing

The velocity domain must contain not just the beam injection point but the full slowing-down tail. For nearly-parallel beams the tail extends in v∥ well beyond the injection point due to pitch-angle scattering during deceleration.

**Confirmed case (2026-06-08):** 120 keV D beam at 10°, xH=50% (cases #3/#5). Beam injection at v∥ ≈ 3.34 Mm/s is within vpar_max=4 Mm/s, but the slowing-down tail reaches the Dirichlet boundary, forcing a strong gradient there. This corrupts the self-collision flux integral, causing apparent SC energy non-conservation (|Psc_perp + Psc_par| / |Psc_par| ~ 80%). Extending to vpar_max=5 Mm/s gives |Psc_perp + Psc_par| < 0.2% of the component magnitude — the SC operator is conservative when the domain is adequate.

**Rule of thumb:** set vpar_max ≥ v_beam × cos(θ) + 3 × max(beam_dvpar, v_th_bulk), and similarly for vperp_max. Always verify SC energy conservation via Psc_perp + Psc_par ≈ 0 in power_coll_self_vs_time.txt when running near-parallel or near-perpendicular beam cases.

### Solver Path
All solvers use the 7-point Fornberg stencil (`fd_stencil_2d`). The solver selected by `isc`:

| `isc` | Time-dependent solver |
|-------|-----------------------|
| `-1`  | `timefp_7pt_nl` (nonlinear self-collisions) |
| `≠-1` | `timefp_7pt` (linear or Maxwellian-SC) |

The steady-state solver is always `FP_steady_state` ← `fd_stencil_2d` (in `mod_linear`).

### Sparse Matrix and PARDISO
All sparse matrices are stored in 1-based CSR format. The `pardiso_solver` module (`pardiso_solver (2).f90`) wraps Intel MKL PARDISO with four entry points that exploit phased factorisation:
- `pardiso_solve_steady` — phases 11+22+33 then release (used by steady-state solver)
- `pardiso_solve_init` — phases 11+22 once (called once before the time loop in `timefp_7pt`)
- `pardiso_solve_step` — phase 33 only (called each time step; the LHS matrix is constant)
- `pardiso_solve_finalize` — releases PARDISO memory

Matrix type `mtype=11` (real non-symmetric general) is used throughout.

### Row Assembly: fd_stencil_2d
`fd_stencil_2d` (`fd_stencil_2d.f90`) is the central stencil routine. For each grid point (i,j) it returns at most 49 sparse entries (7×7 stencil). The Fornberg weights come from `fornberg_weights` in `derivatives_2d.f90`, which handles non-uniform v⊥ spacing exactly.

Both `solve_fp_pardiso` (steady-state) and `timefp_7pt` (time-dependent) call `fd_stencil_2d` in a two-pass loop: first pass counts non-zeros to build `ia`, second pass fills `ja` and `aa`. Entries in each row must be sorted by column index before passing to PARDISO — `sort_stencil` (insertion sort, defined inside `timefp_7pt`) does this.

### Time-Dependent Solver (timefp_7pt)
Implements Crank-Nicolson (`icn=-1`, θ=0.5) or fully implicit (`icn≠-1`, θ=1.0). Each step solves:
```
(I - θ·dt·L)·f^{n+1} = (I + (1-θ)·dt·L)·f^n + dt·S
```
The LHS `(I - θ·dt·L)` is factorised once at the start. Each step only applies `L` to `f^n` via a hand-written sparse mat-vec (`sparse_matvec_csr`, inside `timefp_7pt`), builds the RHS, and calls `pardiso_solve_step` (phase 33 only).

A particle-conservation rescaling (`fout = fout * npart / dens_tmp`) is applied each step when `isource=0` to correct numerical drift.

### Namelist Parameters (key physics flags)
| Parameter | Values | Effect |
|-----------|--------|--------|
| `ntimes` | 0 / >0 | Steady-state / time-dependent |
| `icn` | -1 / else | Crank-Nicolson / fully implicit |
| `isc` | 0 / 1 / -1 | No self-coll / Maxwellian approx / neglected in TD |
| `irf` | -1 / else | Include QL-RF term / no RF |
| `isource` | -1 / 0 | NBI beam source / no source |
| `iold` | -1 / 0 | Restart from `xout.dat` / fresh start |
| `ising` | 0 / -1 / +1 | Uniform / two-domain / quadratic v⊥ grid |

### Output Files
All output is written to the run directory. Key files:
- `fout.txt` / `xout.dat` — final VDF (text / binary restart)
- `fstix.dat` — Stix reference Maxwellian
- `density_vs_time.txt`, `energy_vs_time.txt` — time traces
- `power_coll_*_vs_time.txt`, `power_RF_vs_time.txt` — power balance
- `fout_at_vpar0.txt`, `fout_at_vperp0.txt`, etc. — 1D cuts of f at boundaries

### Argument Ordering Note
The two main solver calls use different argument orders than each other. In `main`:
- `solve_fp_pardiso(A, B, C, D, E, F, …)` → `(all00, all10, all01, all20, all11, all02, …)`
- `timefp_7pt(all00, all10, all01, all11, all20, all02, …)` — note E and D are swapped relative to `solve_fp_pardiso`'s `(A,B,C,D,E,F)`

## Directory Layout

All `.f90` sources live in `Fortran sources/` at the **project root** — not inside the
nested `FP2D_QLRF_NL/` Visual Studio project folder (which holds the `.sln`/`.vfproj` and
the `x64/Debug`, `x64/Release` build outputs). Run outputs land in whichever `x64/*`
directory the exe is launched from.

## Fortran Internals

### Shared variable semantics

| Variable | Module | Meaning |
|----------|--------|---------|
| `npart` | `shared_plasma` | Particle density from namelist (m⁻³); used as SC-operator density and renorm target |
| `dens_tmp` | local in solvers | Actual density of `fout` at current step (diverges from `npart` during NBI fill) |
| `vteff` | `shared_plasma` | Thermal velocity at Stix temperature (fixed; used by `isc=1`) |
| `teff` | local in `timefp_7pt` | Effective temperature (keV) updated unconditionally every step via `time_energy(fout,…,teff=teff)` |
| `ta_eV` | local in `timefp_7pt` | `teff * 1.0d3` — Teff in eV, recomputed each step for isc=2 |
| `vteff_t` | local in `timefp_7pt` | Thermal velocity used for SC this step: `vteff` (isc=1) or `9.79d3*sqrt(ta_eV/aa)` (isc=2) |
| `jmid` | `shared_grid` | v∥ grid index for v∥ = 0 (used for all `_at_vpar0` output slices) |
| `gamma0` | PARAMETER in solvers | `2.390775d-1` — collision frequency pre-factor |
| `pi15` | PARAMETER in `timefp_7pt` | `5.5683279968` = π^(3/2), used for SC Maxwellian normalisation |

### Maxwellian background formula (isc=2, isc=3)

The SC operator uses a Maxwellian background with density `npart` and thermal velocity `vteff_t`:
```
f_M(v⊥, v∥) = npart / ((2π)^{3/2} · vth³) · exp(−(v⊥² + v∥²) / (2·vth²))
              where  vth = 9.79×10³ · sqrt(T[eV] / aa)  m/s  ( = sqrt(T/m) )
```
Here `vth` is the per-degree-of-freedom rms speed `sqrt(T/m)` (NOT the most-probable
speed `sqrt(2T/m)`), so the normalisation constant is `(2π)^{3/2}` and the exponent
carries a factor `1/(2·vth²)`. The temperature `T` is `Teff` for `isc=2` and `Tn` for
`isc=3`. See `TimeFP_7pt.f90` (`twopi15` parameter and the `fM_ij` assembly, ~line 631).
Written to `fsc_maxw.txt` and `fsc_maxw_at_vpar0.txt` at the end of `timefp_7pt` when
`isc==2` or `isc==3`.

### Output file-unit registry (`timefp_7pt`)

New output files must use units not in this table:

| Unit(s) | File(s) | Condition |
|---------|---------|-----------|
| 40 | `fout.txt`, 1D slices (main) | always |
| 45–47 | density / energy / anisotropy vs time | always |
| 470–479 | power_coll_tot / per-species | `iplot_pow=-1` |
| 480 | power_RF_vs_time | RF |
| 490 | power_NBI_vs_time | NBI |
| 500 | power_coll_self_vs_time | `isc≠0` |
| 505 | coulomb_log_vs_time | `nbulk>1` |
| 506 | coulomb_log_self_vs_time | `isc=1,2` |
| 507 | Teff_vs_time | always |
| 508–509 | fsc_maxw / fsc_maxw_at_vpar0 | `isc=2` |
| 514 | Tn_vs_time | always |
| 570–579 | momentum_coll_tot / per-species | `iplot_mom=-1` |
| 580 | momentum_RF_vs_time | RF |
| 590 | momentum_NBI_vs_time | NBI |
| 600 | momentum_coll_self_vs_time | `isc≠0` |

### Output filename convention
All output files go through `outfile(name)` (defined in `shared_timer`):
- `casename = ''` → writes `name` unchanged
- `casename = 'foo'` → writes `stem-foo.ext` (e.g. `fout-foo.txt`)

## fp2d_plot.py Structure

The script is at the project root. Key extension points:

- **`FILE_META`** — dict keyed by file stem; `ptype` in `{"1d","2d","ts","ts2"}` controls which plot function is used. Add new entries here to make a new output file auto-discovered. `_SKIP_STEMS` lists stems that are handled by composite functions instead of the per-file loop.
- **Per-file plot functions** — `plot_1d`, `plot_2d`, `plot_3d`, `plot_ts`, `plot_ts2` (called by the main loop in `plot_directory`)
- **Composite plot functions** — `plot_power_combined`, `plot_power_coll`, `plot_sc_power_split`, `plot_coulomb_log`, `plot_momentum_coll`, `plot_momentum_breakdown`, `plot_momentum_balance`, `plot_fout_vs_maxw_at_vpar0` — called at the end of `plot_directory` when `restrict is None`
- **`plot_directory`** — main orchestrator: iterates txt files, dispatches to per-file functions, then calls composite functions
- **`_outfile(outdir, stem, casename)`** — constructs the expected filename for a given stem+casename (mirrors the Fortran `outfile()` function)
- **Namelist readers** — `_read_ntimes_from_namelist`, `_read_isc_from_namelist`, `_read_iplot_pow_from_namelist`, `_read_iplot_mom_from_namelist` used by `run` subcommand to set plotting flags before the solver runs
