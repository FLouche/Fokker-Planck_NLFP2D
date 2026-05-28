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

**Run:** The program reads a Fortran namelist from stdin and writes output files to the working directory (i.e. `x64/Debug/` when launched from Visual Studio, or wherever the shell is):
```
FP2D_QLRF_NL.exe < inputs/jet_RF_case1.dat
```

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
