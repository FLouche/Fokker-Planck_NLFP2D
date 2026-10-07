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

Never batch unrelated changes into a single commit. Never commit build outputs, `.obj`/`.mod`/`.exe` files, or simulation result files — those are covered by `.gitignore`. The one exception is `Benchmark runs/` (see *Directory Layout*): reference results committed on purpose; add to it only when asked.

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
- `pardiso_solve_init` — phases 11+22 once (called once before the time loop in `timefp_7pt`; the symbolic analysis is reused for the whole run because the sparsity pattern never changes)
- `pardiso_solve_step` — phase 33, preceded by phase 22 when `a_changed=.TRUE.`. Both `timefp_7pt` and `timefp_7pt_nl` always pass `.TRUE.`, so every time step does a numerical refactorisation
- `pardiso_solve_finalize` — releases PARDISO memory

Matrix type `mtype=11` (real non-symmetric general) is used throughout.

### Row Assembly: fd_stencil_2d
`fd_stencil_2d` (`fd_stencil_2d.f90`) is the central stencil routine. For each grid point (i,j) it returns at most 49 sparse entries (7×7 stencil). The Fornberg weights come from `fornberg_weights` in `derivatives_2d.f90`, which handles non-uniform v⊥ spacing exactly.

Both `solve_fp_pardiso` (steady-state) and `timefp_7pt` (time-dependent) call `fd_stencil_2d` in a two-pass loop: first pass counts non-zeros to build `ia`, second pass fills `ja` and `aa`. Entries in each row must be sorted by column index before passing to PARDISO — `sort_stencil` (insertion sort, defined inside `timefp_7pt`) does this.

### Time-Dependent Solver (timefp_7pt)
Implements the θ-scheme: Crank-Nicolson (`icn=-1`, θ=0.5), intermediate (`icn=1`, θ=0.75) or fully implicit (any other `icn`, θ=1.0), in both `timefp_7pt` and `timefp_7pt_nl`. Note `icn=1` is **not** fully implicit: its stiff components flip sign each step (amplification → −(1−θ)/θ = −1/3), which is what makes per-step diagnostics such as `density_terms_vs_time` scatter while the density trace stays smooth. Each step solves:
```
(I - θ·dt·L)·f^{n+1} = (I + (1-θ)·dt·L)·f^n + dt·S
```
**The operator is rebuilt and refactorised every step; it is not constant.** Each step recomputes `Teff` (and `Tn`) from `f^n`, updates the Coulomb logarithms, reassembles the coefficient arrays (`assemble_FP_terms`, plus the `sc**` self-collision terms when `isc=1,2,3`), refills the values of `aa_L` and of the LHS `(I - θ·dt·L)` on the fixed CSR pattern, applies `L` to `f^n` via a hand-written sparse mat-vec (`sparse_matvec_csr`, inside `timefp_7pt`) to build the RHS, and calls `pardiso_solve_step` with `a_changed=.TRUE.` (phases 22+33). This holds for every case, including RF-off/SC-off ones. The only exception is the NBI fill-up (`isource=-1`, `iold=0`, density below 5% of `npart`): there the Coulomb-log update is skipped and the precomputed `all**_lin` arrays are used, but the matrix is still refilled and refactorised.

Consequence (stencil study, 2026-09-16): self-collisions (`isc=2`) add no matrix work. They cost only the `self_coll_max` evaluation and six extra N² arrays: +3–13% run time and about +1 MB at 161×161 for the JET RF case.

**No particle-conservation rescaling is applied.** A `fout = fout * npart / dens_tmp`
renormalisation exists but is commented out (`fstart = x_vec!*npart/dens_tmp`,
`TimeFP_7pt.f90`), and with `isource=0` the loss term is off too (`taum = 0`,
`consts.f90`). Nothing therefore counteracts particles absorbed at the Dirichlet
boundaries, so a sourceless RF run drains: once the VDF shape has converged, `f` decays
as a fixed shape with slowly falling amplitude. Verified 2026-07-30 on a 60×41 `isc=2`
RF case — density fell at a steady 1.74×10⁻³ s⁻¹ with the shape frozen. This puts a
**floor** under the amplitude-based convergence rate `eps`; use `i_conv_shape=-1`
(shape criterion) for such cases, or set `ss_tol_eps` above the drain rate.

### Namelist Parameters (key physics flags)
| Parameter | Values | Effect |
|-----------|--------|--------|
| `ntimes` | 0 / >0 | Steady-state / time-dependent |
| `icn` | -1 / 1 / else | Crank-Nicolson (θ=0.5) / θ=0.75 / fully implicit (θ=1) |
| `isc` | 0 / 1 / -1 | No self-coll / Maxwellian approx / neglected in TD |
| `irf` | -1 / else | Include QL-RF term / no RF |
| `isource` | -1 / 0 | NBI beam source / no source |
| `iold` | -1 / 0 | Restart from `xout.dat` / fresh start |
| `ising` | 0 / -1 / +1 | Uniform / two-domain / quadratic v⊥ grid |
| `n_snap` | 0 / N>0 | No snapshots (default) / write f every N steps to `vdf_snap_<step>.txt` (both time solvers) |

### Output Files
All output is written to the run directory. Key files:
- `fout.txt` / `xout.dat` — final VDF (text / restart file). The first line of `xout.dat` is `time  nstep` (total steps, including those of any run it was restarted from); `main` also accepts older files whose first line holds only the time, and then counts from 0.
- `vdf_snap_<step>.txt` — f every `n_snap` steps (`write_vdf_snapshot`, `time_comps_mod`), same layout as `fout.txt` with a `# time =` header. `<step>` is the total step count, so a restart (`iold=-1`) continues the previous run's numbering and cadence instead of overwriting its snapshots. Raw solution: the end-of-run renormalisation of sourceless runs is not applied. Snapshots are not restart files.
- `fstix.dat` — Stix reference Maxwellian
- `density_vs_time.txt`, `energy_vs_time.txt` — time traces
- `power_coll_*_vs_time.txt`, `power_RF_vs_time.txt` — power balance
- `fout_at_vpar0.txt`, `fout_at_vperp0.txt`, etc. — 1D cuts of f at boundaries

The full list is in **Diagnostics** below.

### Diagnostics

Everything the code writes, by category. All files go through `outfile()` (so they
take the `-<casename>` suffix) unless marked otherwise. Unit numbers are in the
file-unit registry under *Fortran Internals*. Time traces are written by both
`timefp_7pt` and `timefp_7pt_nl`.

**Distribution function** (`main-FP_Coll_2D.f90`, `consts.f90`, solvers)

| File | Content | Condition |
|------|---------|-----------|
| `fout.txt`, `xout.dat` | final f (text / restart) | always |
| `fout_at_{vpar0,vparmax,vperp0,vperpmax}.txt` | 1D cuts of f | always |
| `fstix.dat`, `fstix_at_{vpar0,vparmax,vperp0,vperpmax}.txt` | Stix reference Maxwellian and its cuts | always |
| `vdf_snap_<step>.txt` | f every `n_snap` steps | `n_snap>0` |
| `fsc_maxw.txt`, `fsc_maxw_at_vpar0.txt` | Maxwellian background of the SC operator | `isc=2,3` |
| `beam.txt` | NBI source | NBI |

**Energy maps** (`analysis.f90`, end of run): `Ekin.txt`, `Ekin_perp.txt`,
`Ekin_par.txt`, `Ekin_perp_at_vpar0.txt`.

**Time traces** (`*_vs_time.txt`)

| Group | Files | Condition |
|-------|-------|-----------|
| moments | `density`, `energy`, `anisotropy`, `Teff`, `Tn` | always |
| collision times | `tau_coll` (effective τ_ii, τ_ie) | always |
| power balance | `power_coll_tot`, `power_coll_e`, `power_coll_ion<k>` (one per species), `power_coll_self` (`isc≠0`), `power_RF` (RF), `power_NBI` (NBI) | `iplot_pow=-1` |
| density balance | `density_terms`: dn/dt split into total, coll per species, SC, RF, source, losses (m⁻³/s) | `iplot_pow=-1` |
| momentum balance | `momentum_coll_tot`, `momentum_coll_e`, `momentum_coll_ion<k>`, `momentum_coll_self`, `momentum_RF`, `momentum_NBI` (same conditions as power) | `iplot_mom=-1` |
| Coulomb logs | `coulomb_log` (`nbulk>1`), `coulomb_log_self` (`isc=1,2`) | — |
| RF | `tau_rf` (RF tail formation time) | RF |
| convergence | `conv_eps_vs_time.txt` (time, eps, eps_tail), `conv_diag_vs_time.csv` (full history; `conv_diag.f90`) | `i_ss_check=-1` |

The `power_coll_ion<k>` / `momentum_coll_ion<k>` names are built from the species
index and may carry a space (see *`--files` matching* below).

**End-of-run balance checks** (`analysis.f90`, stdout only): power, particle
density and momentum balance (`test_power_balance_7pt`,
`test_density_balance_7pt`, `test_momentum_balance_7pt`), always run.

**Cache, not a diagnostic:** `phi_kern-<case>.dat`, the compressed φ-kernel of
`timefp_7pt_nl` (`mod_phi_kernel`), reused when the grid matches.

**Not routed through `outfile()`** (no casename, overwritten by every run):
`qlrfterm.f90` writes `RF_dirac.txt` (`notxt=0`).

**Removed 2026-10-05 (branch `last-dev`, after `1845cf1`):** the namelist
parameter `idiag` and everything it switched on — `sc_density_map.txt`,
`coef_map.txt`, `sc_*_at_vpar0.txt` (`sc_components_diag`,
`sc_components_maxw_diag`) — and the QL-tensor test output of `qlrfterm.f90`.
Commit `1845cf1` holds the last version with them. A namelist that still sets
`idiag` now stops at read time with `forrtl: severe (19)`; delete the entry.

### Argument Ordering Note
The two main solver calls use different argument orders than each other. In `main`:
- `solve_fp_pardiso(A, B, C, D, E, F, …)` → `(all00, all10, all01, all20, all11, all02, …)`
- `timefp_7pt(all00, all10, all01, all11, all20, all02, …)` — note E and D are swapped relative to `solve_fp_pardiso`'s `(A,B,C,D,E,F)`

## Directory Layout

All `.f90` sources live in `Fortran sources/` at the **project root** — not inside the
nested `FP2D_QLRF_NL/` Visual Studio project folder (which holds the `.sln`/`.vfproj` and
the `x64/Debug`, `x64/Release` build outputs). Run outputs land in whichever `x64/*`
directory the exe is launched from.

### `Benchmark runs/` — reference results under version control

Added 2026-10-07 (`a737981`): a committed copy of the benchmark cases from
`FP2D_QLRF_NL/x64/Release/Benchmark/` (which is git-ignored), with the paths
below `Benchmark/` kept — namelists, outputs, figures and reports:

| Folder | Content |
|--------|---------|
| `ITER-RF/ITER-RF-N=2/nominal/nominal_#2--REF` | ITER RF, N=2 nominal reference case |
| `ITER-RF/ITER-RF-N=1/nominal` | ITER RF, N=1 nominal cases |
| `JET-Beam` | JET NBI benchmark cases |
| `JET-RF` | JET RF benchmark cases (without `RF-Case2/Test isc`) |

Use them to check a solver change against the reference outputs (rerun the
namelist, compare with `fp2d_plot.py compare`). Points to know:

- **Missing on purpose:** the `sum_phi-JET-beam{1,2,3}-NLSC.dat` φ-kernel caches
  of `JET-Beam/Case*/Restart files` (4.2 GB each, over GitHub's 100 MB file
  limit — Git LFS stops at 2 GB too). The solver rebuilds the kernel when
  absent, so the runs stay reproducible. Never try to add files > 100 MB.
- **`results/*.out` are solver stdout logs**, force-added (`git add -f`): they
  match the LaTeX `*.out` ignore pattern only by accident. LaTeX build leftovers
  next to the reports stay ignored.
- The folder is ~2.9 GB; adding more results grows the clone size for good.

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
| 515 | tau_rf_vs_time (RF tail formation time) | RF |
| 516 | tau_coll_vs_time (effective tau_ii, tau_ie) | always |
| 518 | conv_diag_vs_time.csv (full convergence history) | `i_ss_check=-1` |
| 519 | conv_eps_vs_time (time, eps, eps_tail) | `i_ss_check=-1` |
| 520 | density_terms_vs_time (time, total, coll(1:nbulk), SC, RF, source, losses; dn/dt per operator term, m⁻³/s) | `iplot_pow=-1` |
| 530 | vdf_snap_<step> (f every `n_snap` steps; opened and closed per snapshot) | `n_snap>0` |
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
- **`plot_directory`** — main orchestrator: iterates txt files, dispatches to per-file functions, then calls composite functions. Called once per case by `main`, with `defer_show=True` on all but the last so a single blocking `plt.show()` opens every case's windows together
- **`_outfile(outdir, stem, casename)`** — constructs the expected filename for a given stem+casename (mirrors the Fortran `outfile()` function)
- **Namelist readers** — `_read_ntimes_from_namelist`, `_read_isc_from_namelist`, `_read_iplot_pow_from_namelist`, `_read_iplot_mom_from_namelist` used by `run` subcommand to set plotting flags before the solver runs

### Case selection (`--cases`)

Both `plot` and `compare` take `--cases CASE [CASE ...]`; `plot` also accepts
`--casename` as an alias. In `plot` each case gets its own set of figures;
in `compare` the cases are overlaid on shared axes. `run` keeps a **singular**
`--casename` — one run produces one case, and the label is read from the
namelist when omitted.

A case the user *types* also **selects which files are read**, not merely how
they are labelled: `--cases X --files fout` plots `fout-X.txt` alone, not every
`fout-*.txt` in the directory (`strict_case=True`). A case that was
*auto-detected* only labels, so `--files` with no `--cases` still plots across
every case in a multi-case folder.

### `--files` matching

`--files` entries may be a glob, a full filename, or a bare stem; `plot` and
`compare` use the same matcher (`_restrict_match`). Three traps are handled
explicitly, all from the per-species outputs whose stem carries a space
(`power_coll_ion 1_vs_time`, built in Fortran from the species index):

- `_canon_stem` ignores spaces, underscores and case, so `power_coll_ion 1_vs_time`,
  `power_coll_ion1_vs_time` and `power_coll_ion_1_vs_time` all match.
- `_rejoin_split_stems` puts back together adjacent argv entries that a shell
  split at that space, when they name a file that exists.
- `_species_file_meta` discovers these files on disk for `compare`, which
  iterates known stems rather than the directory — they are not `FILE_META`
  keys, so without it they were invisible to it.

An explicitly named file is drawn even when its stem is in `_SKIP_STEMS`: the
composite that would otherwise draw it only runs when the plot set is
unrestricted. A `--files` entry that matches nothing prints a warning rather
than silently plotting less than was asked for.

### Logarithmic scales

Three independent options, available to `plot` and `compare`:

| Option | Scales | Mechanism |
|--------|--------|-----------|
| `--logf` (alias `--log`) | the plotted **quantity** — colour scale of a 2D map, z of a 3D surface, y of a 1D profile | the `log` argument threaded through each plot function |
| `--logx` | the x-axis | `_LOGX` global, applied in `_finish` |
| `--logy` | the y-axis | `_LOGY` global, applied in `_finish` |

`--logf` has to be built in at plot time (it changes contour levels and the
colour norm), so it stays a parameter. The two axis scales are a post-hoc
property of the axes, so they follow the `_XRANGE` pattern: a module global set
in `main`, applied to every figure in `_finish` via `_apply_log_axes`. Adding a
plot function therefore needs no work to support them.

**The zero-crossing rule** (`_hides_data`, used by `_positive_span`): an axis
whose negative excursion exceeds `_LOG_NEG_TOL = 1e-6` of its positive range is
left linear and prints why — `v_par` is signed, and so is a power that changes
sign, and a log scale there drops half the data with no visible indication.
Anything smaller is round-off: the far tail of a VDF dips a few 1e-9 of its
peak below zero, and refusing a log axis over that would block the commonest
use of `--logy`. In that case the axis is drawn logarithmically, the lower
limit is pinned to the smallest positive sample so matplotlib does not pad down
to an arbitrary decade, and the dropped points are counted in a note.

`_positive_span` finds the range from line data where there is any, and from
`ax.dataLim` otherwise — which is what makes `--logy` work on the `v_perp` axis
of a 2D contour map.

### Axis zoom (`--xrange`, `--yrange`)

Both are `plot`/`compare` only, parsed by `_parse_range` and held in the
`_XRANGE` / `_YRANGE` globals, applied in `_finish` like the log scales.
`_warn_if_range_empty` reports a window that selects no data, since the plot
would otherwise be blank with no hint why (usually a units mismatch — the
velocity grid spans ~0 to a few 10⁷ m/s, so `0:1` selects nothing).

Two ordering constraints, both load-bearing:

- `--xrange` alone rescales y to the data inside the window
  (`_autoscale_y_to_xrange`), which is what you want when zooming a trace. That
  rescale is **skipped** when `_YRANGE` is set, or it would immediately undo
  the bounds the user asked for.
- `_YRANGE` is applied **after** `_apply_log_axes`, which sets limits of its
  own when it pins a log axis to the smallest positive sample. So
  `--logy --yrange 1e-6:10` gives the requested decades.

`--xrange` sets `fig.axes[0]` only (multi-panel figures are built `sharex=True`);
`--yrange` sets every non-colorbar axes, since y is not shared.

### Movie from the snapshots (`--movie`)

`run` and `plot` accept `--movie` (plus `--fps`, `--movie-scale`, `--movie-format`, and
`--namelist FILE` for `plot`). It replaces the usual figures with one animation
per case (`make_movie`), a 2×2 array — f and E_kin on top, E_⊥ and E_∥ below
(½m v⊥² f and ½m v∥² f, which add up to E_kin; the solver's `Ekin_par.txt` is
m v∥² f, twice the movie's E_∥) — one
frame per `vdf_snap_<step>.txt`. `--movie-scale fixed` (default) keeps the
colour scales fixed over the whole movie; `frame` rescales each frame to its own
maximum (`set_clim` per frame, so the colour bars follow), and each panel title
shows the frame's maximum in both modes. `--tstop T` ends the movie at time T
(s): snapshots are filtered on the `# time =` header (`_snapshot_time`) before
any is loaded, so later ones cost nothing, and fixed scales come from the kept
frames only. `--fps`, `--movie-scale`, `--movie-format`, `--tstop`, `--logE`,
`--logEperp` and `--logEpar` each **imply `--movie`** (`_MOVIE_OPTS`): their
argparse defaults are `None` (also for the three `store_true` log flags) so `main` can tell a typed option from an untouched one,
and `_MOVIE_DEFAULTS` fills them in afterwards. Without this, `plot . --cases X
--movie-scale frame` silently drew the ordinary figures. Rules:

- **`n_snap` is checked first**, from the namelist. With `run` this happens
  before the solver starts, so a missing `n_snap` does not cost a run. With
  `plot` the namelist is `--namelist`, or the one in the output folder whose
  `casename` matches (`_find_namelist`, which reads only file heads looking for
  `&INPUT`). No namelist, or `n_snap = 0`, skips the movie with a message.
- **Ekin is rebuilt from f** (`_ekin_map`), exactly as `analysis.f90` builds
  `Ekin.txt` (Simpson-weighted density, `aa` from the namelist); it matches the
  solver's `Ekin.txt` to round-off.
- `--logf` gives f eight decades and leaves everything below the floor blank,
  rather than drawing the round-off tail in the lowest colour. `--logE`,
  `--logEperp` and `--logEpar` do the same, independently, for the three
  energy panels (`elog` in `make_movie`; one `limits` function serves all four).
- `--movie-format {avi,mp4,gif}`, default `avi`: MPEG-4 Part 2 tagged `xvid`
  (plays in VLC and Windows' own players), H.264, or an animated GIF (Pillow).
  avi/mp4 need ffmpeg. `_ffmpeg_available` uses one on the PATH, or else the
  binary bundled with `pip install imageio-ffmpeg` (installed on this machine,
  2026-10-01). Without either it falls back to GIF with a message. Written as
  `movie-<case>.<fmt>` to `--save DIR` or the output folder.
- `vdf_snap` is in `_SKIP_STEMS`: without it, a plain `plot` would draw every
  snapshot as a separate 2D map.
