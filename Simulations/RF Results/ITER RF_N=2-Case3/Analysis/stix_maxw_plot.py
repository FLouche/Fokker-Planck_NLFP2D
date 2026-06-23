#!/usr/bin/env python3
"""
stix_maxw_plot.py
=================
Plot the isotropic Maxwellian background f_M(v_perp, v_par; T_eff) used by
the isc=2 self-collision operator in FP2D_QLRF_NL, for several T_eff values,
and overlay the actual steady-state VDF from the NLSC (isc=-1) run.

Formula (matches the background used by cblin / self_coll_max, i.e. the
standard Maxwellian at temperature T with v_th = sqrt(T/m)):

    v_th = 9.79e3 * sqrt(T_eff[eV] / A)          [m/s]   (= sqrt(T/m))
    f_M  = n / ((2*pi)^{3/2} * v_th^3)
           * exp(-(v_perp^2 + v_par^2) / (2 * v_th^2))

Case: ITER-EDA RF case 3 (pure ICRF, fundamental N=1, no beam source)
    Minority : T (tritium)   A = 3, Z = 1
    Majority : D (deuterium) A = 2, Z = 1
    n_e = 1.5e20 m^-3,  x_T = 0.50  =>  n_T = 7.5e19 m^-3
    T_bulk = 12 keV  (NLSC steady state: T_eff = 28.1 keV, v_perp box = 35 Mm/s)

Output: stix_maxw_comparison.png (saved alongside this script)
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from pathlib import Path

# ---------------------------------------------------------------------------
# Physics parameters
# ---------------------------------------------------------------------------
A_MIN = 3.0       # tritium atomic mass
N_MIN = 7.5e19    # minority density [m^-3]
C_VTH = 9.79e3    # sqrt(e / m_amu) [m/s]  — matches Fortran convention

# T_eff values for the Maxwellian backgrounds [keV]
TEFF_KEV = [12, 20, 25, 50]

# Path to the NLSC steady-state VDF (fout file, 3 columns: vperp, vpar, f)
FOUT_PATH = Path(
    r'E:\Fokker-Planck\FPCode-CoulColl-2D\FP2D_QLRF_NL'
    r'\FP2D_QLRF_NL\x64\Release\fout-ITER-RF3-TD1-NLSC.txt'
)
FOUT_TEFF_KEV = 28.1    # T_eff of the NLSC steady state [keV]

# Velocity domain for the plots [Mm/s]
VPAR_MAX  = 4.0    # ±Mm/s  (sim range is ±4 Mm/s)
VPERP_MAX = 10.0   # Mm/s   (covers the NLSC perp tail down to ~1e-6 of peak)
NV        = 600

# ---------------------------------------------------------------------------
# Colour palette: cool blue -> warm red as T_eff increases
# ---------------------------------------------------------------------------
CMAP   = plt.cm.plasma
COLORS = CMAP(np.linspace(0.12, 0.88, len(TEFF_KEV)))

# Style for the NLSC VDF overlay
NLSC_COLOR = '#111111'    # near-black
NLSC_LW    = 2.5

# Linestyles for the three Maxwellian iso-contour levels
LEVELS_EXP = [1, 4, 9]   # f = f_peak * exp(-n)
LSTYLES    = ['-', '--', ':']

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def vth(Teff_keV: float) -> float:
    """Thermal velocity [Mm/s] for minority T at T_eff."""
    return C_VTH * np.sqrt(Teff_keV * 1e3 / A_MIN) / 1e6


def fM_peak(Teff_keV: float) -> float:
    """Peak value of f_M at v=0  [m^-3 (m/s)^-3]."""
    vt_SI = vth(Teff_keV) * 1e6
    return N_MIN / ((2.0 * np.pi)**1.5 * vt_SI**3)


def fM_2d(Teff_keV: float,
          VPar_Mms: np.ndarray,
          VPerp_Mms: np.ndarray) -> np.ndarray:
    """2D isotropic Maxwellian on a velocity mesh [m^-3 (m/s)^-3]."""
    vt_SI  = vth(Teff_keV) * 1e6
    v2_SI2 = (VPar_Mms**2 + VPerp_Mms**2) * 1e12   # (Mm/s)^2 -> (m/s)^2
    return fM_peak(Teff_keV) * np.exp(-v2_SI2 / (2.0 * vt_SI**2))


# ---------------------------------------------------------------------------
# Load NLSC VDF
# ---------------------------------------------------------------------------
_d         = np.loadtxt(FOUT_PATH)
NPERP      = len(np.unique(_d[:, 0]))   # 150
NPAR       = len(np.unique(_d[:, 1]))   # 150
FOUT_VPERP = _d[::NPAR,  0] / 1e6      # unique v_perp [Mm/s]
FOUT_VPAR  = _d[:NPAR,   1] / 1e6      # unique v_par  [Mm/s]
FOUT_Z     = _d[:, 2].reshape(NPERP, NPAR)   # f [m^-3 (m/s)^-3], shape (Nperp, Npar)

# 1D cut at v_par closest to 0
JMID       = int(np.argmin(np.abs(FOUT_VPAR)))
FOUT_F1D   = FOUT_Z[:, JMID]           # f(v_perp, v_par≈0)

# ---------------------------------------------------------------------------
# Maxwellian velocity grids
# ---------------------------------------------------------------------------
v_perp_Mms = np.linspace(0,          VPERP_MAX, NV)
v_par_Mms  = np.linspace(-VPAR_MAX,  VPAR_MAX,  NV)
VPar_mesh, VPerp_mesh = np.meshgrid(v_par_Mms, v_perp_Mms)

# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------
fig, (ax2d, ax1d) = plt.subplots(1, 2, figsize=(14, 6.5))
fig.suptitle(
    r'SC Maxwellian backgrounds $f_M$ vs NLSC steady-state VDF'
    ' — ITER-EDA RF case 3 (fundamental, N=1)\n'
    r'Minority T (A=3), $n_T = 7.5\times10^{19}$ m$^{-3}$, '
    r'$T_\mathrm{bulk} = 12$ keV',
    fontsize=12
)

# ── Left panel: 2D iso-contour lines ──────────────────────────────────────

# Maxwellian contours
for i, T in enumerate(TEFF_KEV):
    F   = fM_2d(T, VPar_mesh, VPerp_mesh)
    fpk = fM_peak(T)
    for exp_val, ls in zip(LEVELS_EXP, LSTYLES):
        ax2d.contour(VPar_mesh, VPerp_mesh, F,
                     levels=[fpk * np.exp(-exp_val)],
                     colors=[COLORS[i]], linewidths=1.8, linestyles=[ls])

# NLSC VDF contours
FOUT_VPar_mesh, FOUT_VPerp_mesh = np.meshgrid(FOUT_VPAR, FOUT_VPERP)
fout_pk = FOUT_Z.max()
for exp_val, ls in zip(LEVELS_EXP, LSTYLES):
    ax2d.contour(FOUT_VPar_mesh, FOUT_VPerp_mesh, FOUT_Z,
                 levels=[fout_pk * np.exp(-exp_val)],
                 colors=[NLSC_COLOR], linewidths=NLSC_LW, linestyles=[ls])

# Legends
color_handles = [
    mlines.Line2D([], [], color=COLORS[i], lw=2.5,
                  label=rf'$f_M$: $T_\mathrm{{eff}}$ = {T} keV')
    for i, T in enumerate(TEFF_KEV)
] + [
    mlines.Line2D([], [], color=NLSC_COLOR, lw=NLSC_LW,
                  label=rf'NLSC VDF ($T_\mathrm{{eff}}$ = {FOUT_TEFF_KEV:.0f} keV)')
]
ls_labels = [
    rf'$f = f_\mathrm{{pk}}\,e^{{-{n}}}$  ($v = {int(n**0.5)}\,v_\mathrm{{th}}$)'
    for n in LEVELS_EXP
]
ls_handles = [
    mlines.Line2D([], [], color='k', lw=1.6, ls=ls, label=lbl)
    for ls, lbl in zip(LSTYLES, ls_labels)
]

leg_color = ax2d.legend(handles=color_handles, loc='upper right',
                        fontsize=8.5, framealpha=0.85)
ax2d.add_artist(leg_color)
ax2d.legend(handles=ls_handles, loc='lower right', fontsize=8.5, framealpha=0.85)

ax2d.set_xlabel(r'$v_\parallel$ (Mm s$^{-1}$)', fontsize=12)
ax2d.set_ylabel(r'$v_\perp$ (Mm s$^{-1}$)',     fontsize=12)
ax2d.set_title(r'2D iso-contours in $(v_\parallel,\,v_\perp)$ space', fontsize=11)
ax2d.set_xlim(-VPAR_MAX, VPAR_MAX)
ax2d.set_ylim(0,          VPERP_MAX)
ax2d.axhline(0, color='k', lw=0.6)
ax2d.axvline(0, color='k', lw=0.6, ls=':')
ax2d.set_aspect('equal')
ax2d.grid(True, alpha=0.20)

# ── Right panel: 1D cut at v_par ≈ 0 ──────────────────────────────────────

# Maxwellian curves
for i, T in enumerate(TEFF_KEV):
    vt_SI  = vth(T) * 1e6
    fpk    = fM_peak(T)
    f1d    = fpk * np.exp(-(v_perp_Mms * 1e6)**2 / (2.0 * vt_SI**2))
    vt_Mms = vth(T)
    ax1d.plot(v_perp_Mms, f1d, color=COLORS[i], lw=2.2,
              label=rf'$f_M$: $T_\mathrm{{eff}}$ = {T} keV'
                    rf'  ($v_\mathrm{{th}}$ = {vt_Mms:.2f} Mm/s)')
    ax1d.axvline(vt_Mms, color=COLORS[i], lw=0.9, ls='--', alpha=0.55)

# NLSC VDF curve (limit to plot domain)
mask = FOUT_VPERP <= VPERP_MAX
ax1d.plot(FOUT_VPERP[mask], FOUT_F1D[mask],
          color=NLSC_COLOR, lw=NLSC_LW,
          label=rf'NLSC VDF ($T_\mathrm{{eff}}$ = {FOUT_TEFF_KEV:.0f} keV, $v_\parallel \approx 0$)')

ax1d.set_xlabel(r'$v_\perp$ (Mm s$^{-1}$)', fontsize=12)
ax1d.set_ylabel(r'$f(v_\perp,\,0)$   (m$^{-3}$ (m s$^{-1})^{-3}$)', fontsize=11)
ax1d.set_title(r'1D cut at $v_\parallel \approx 0$  (log scale)', fontsize=11)
ax1d.set_yscale('log')
ax1d.set_xlim(0, VPERP_MAX)
ax1d.legend(fontsize=8.5, loc='lower left')
ax1d.grid(True, alpha=0.25, which='both')
ax1d.text(0.97, 0.55, r'Dashed verticals: $v_\perp = v_\mathrm{th}$',
          transform=ax1d.transAxes, ha='right', va='bottom',
          fontsize=8, color='grey')

plt.tight_layout()

out_path = Path(__file__).with_name('stix_maxw_comparison.png')
fig.savefig(out_path, dpi=150, bbox_inches='tight')
print(f'Saved -> {out_path}')

# ---------------------------------------------------------------------------
# Figure 2: NLSC VDF vs Stix Maxwellian at T_eff = 12 keV only
# ---------------------------------------------------------------------------
T12       = 12.0                          # keV
vt12_SI   = vth(T12) * 1e6               # m/s
fpk12     = fM_peak(T12)
COL_12    = COLORS[0]                     # same colour as in figure 1

fig2, (ax2d2, ax1d2) = plt.subplots(1, 2, figsize=(14, 6.5))
fig2.suptitle(
    r'NLSC steady-state VDF vs SC Maxwellian at $T_\mathrm{eff}$ = 12 keV'
    ' — ITER-EDA RF case 3 (fundamental, N=1)\n'
    r'Minority T (A=3), $n_T = 7.5\times10^{19}$ m$^{-3}$',
    fontsize=12
)

# ── Figure 2 — Left panel: 2D iso-contours ────────────────────────────────

# Maxwellian at 12 keV
F12 = fM_2d(T12, VPar_mesh, VPerp_mesh)
for exp_val, ls in zip(LEVELS_EXP, LSTYLES):
    ax2d2.contour(VPar_mesh, VPerp_mesh, F12,
                  levels=[fpk12 * np.exp(-exp_val)],
                  colors=[COL_12], linewidths=1.8, linestyles=[ls])

# NLSC VDF
fout_pk = FOUT_Z.max()
for exp_val, ls in zip(LEVELS_EXP, LSTYLES):
    ax2d2.contour(FOUT_VPar_mesh, FOUT_VPerp_mesh, FOUT_Z,
                  levels=[fout_pk * np.exp(-exp_val)],
                  colors=[NLSC_COLOR], linewidths=NLSC_LW, linestyles=[ls])

# Legends
color_handles2 = [
    mlines.Line2D([], [], color=COL_12,    lw=2.5,
                  label=rf'$f_M$: $T_\mathrm{{eff}}$ = {T12:.0f} keV'
                        rf'  ($v_\mathrm{{th}}$ = {vth(T12):.2f} Mm/s)'),
    mlines.Line2D([], [], color=NLSC_COLOR, lw=NLSC_LW,
                  label=rf'NLSC VDF ($T_\mathrm{{eff}}$ = {FOUT_TEFF_KEV:.0f} keV)'),
]
ls_handles2 = [
    mlines.Line2D([], [], color='k', lw=1.6, ls=ls, label=lbl)
    for ls, lbl in zip(LSTYLES, ls_labels)
]
leg2_color = ax2d2.legend(handles=color_handles2, loc='upper right',
                          fontsize=9, framealpha=0.85)
ax2d2.add_artist(leg2_color)
ax2d2.legend(handles=ls_handles2, loc='lower right', fontsize=8.5, framealpha=0.85)

ax2d2.set_xlabel(r'$v_\parallel$ (Mm s$^{-1}$)', fontsize=12)
ax2d2.set_ylabel(r'$v_\perp$ (Mm s$^{-1}$)',     fontsize=12)
ax2d2.set_title(r'2D iso-contours in $(v_\parallel,\,v_\perp)$ space', fontsize=11)
ax2d2.set_xlim(-VPAR_MAX, VPAR_MAX)
ax2d2.set_ylim(0,          VPERP_MAX)
ax2d2.axhline(0, color='k', lw=0.6)
ax2d2.axvline(0, color='k', lw=0.6, ls=':')
ax2d2.set_aspect('equal')
ax2d2.grid(True, alpha=0.20)

# ── Figure 2 — Right panel: 1D cut at v_par ≈ 0 ──────────────────────────

# Maxwellian at 12 keV
f1d_12 = fpk12 * np.exp(-(v_perp_Mms * 1e6)**2 / (2.0 * vt12_SI**2))
ax1d2.plot(v_perp_Mms, f1d_12, color=COL_12, lw=2.2,
           label=rf'$f_M$: $T_\mathrm{{eff}}$ = {T12:.0f} keV'
                 rf'  ($v_\mathrm{{th}}$ = {vth(T12):.2f} Mm/s)')
ax1d2.axvline(vth(T12), color=COL_12, lw=0.9, ls='--', alpha=0.55)

# NLSC VDF
mask = FOUT_VPERP <= VPERP_MAX
ax1d2.plot(FOUT_VPERP[mask], FOUT_F1D[mask],
           color=NLSC_COLOR, lw=NLSC_LW,
           label=rf'NLSC VDF ($T_\mathrm{{eff}}$ = {FOUT_TEFF_KEV:.0f} keV, $v_\parallel \approx 0$)')

ax1d2.set_xlabel(r'$v_\perp$ (Mm s$^{-1}$)', fontsize=12)
ax1d2.set_ylabel(r'$f(v_\perp,\,0)$   (m$^{-3}$ (m s$^{-1})^{-3}$)', fontsize=11)
ax1d2.set_title(r'1D cut at $v_\parallel \approx 0$  (log scale)', fontsize=11)
ax1d2.set_yscale('log')
ax1d2.set_xlim(0, VPERP_MAX)
ax1d2.legend(fontsize=9, loc='lower left')
ax1d2.grid(True, alpha=0.25, which='both')
ax1d2.text(0.97, 0.55, r'Dashed vertical: $v_\perp = v_\mathrm{th}$',
           transform=ax1d2.transAxes, ha='right', va='bottom',
           fontsize=8, color='grey')

fig2.tight_layout()

out_path2 = Path(__file__).with_name('stix_maxw_vs_nlsc.png')
fig2.savefig(out_path2, dpi=150, bbox_inches='tight')
print(f'Saved -> {out_path2}')

# ---------------------------------------------------------------------------
# Figure 3: zoom of figure 1's right panel (1D cut at v_par=0) to 0 - 3 Mm/s,
#           with the log y-axis adapted to the data inside that x-window.
# ---------------------------------------------------------------------------
XZOOM_MAX = 3.0   # Mm/s

fig3, ax3 = plt.subplots(figsize=(8, 6))
xz = v_perp_Mms <= XZOOM_MAX
ymins, ymaxs = [], []

for i, T in enumerate(TEFF_KEV):
    vt_SI  = vth(T) * 1e6
    fpk    = fM_peak(T)
    f1d    = fpk * np.exp(-(v_perp_Mms * 1e6)**2 / (2.0 * vt_SI**2))
    vt_Mms = vth(T)
    ax3.plot(v_perp_Mms, f1d, color=COLORS[i], lw=2.2,
             label=rf'$f_M$: $T_\mathrm{{eff}}$ = {T} keV'
                   rf'  ($v_\mathrm{{th}}$ = {vt_Mms:.2f} Mm/s)')
    if vt_Mms <= XZOOM_MAX:
        ax3.axvline(vt_Mms, color=COLORS[i], lw=0.9, ls='--', alpha=0.55)
    yz = f1d[xz]
    ymins.append(yz[yz > 0].min()); ymaxs.append(yz.max())

mask = FOUT_VPERP <= XZOOM_MAX
ax3.plot(FOUT_VPERP[mask], FOUT_F1D[mask],
         color=NLSC_COLOR, lw=NLSC_LW,
         label=rf'NLSC VDF ($T_\mathrm{{eff}}$ = {FOUT_TEFF_KEV:.0f} keV, $v_\parallel \approx 0$)')
yz = FOUT_F1D[mask]
if np.any(yz > 0):
    ymins.append(yz[yz > 0].min()); ymaxs.append(yz.max())

# Adapt the (log) y-range to the data within 0 - XZOOM_MAX, with a small margin
lo, hi = np.log10(min(ymins)), np.log10(max(ymaxs))
pad = 0.05 * (hi - lo)
ax3.set_ylim(10.0**(lo - pad), 10.0**(hi + pad))

ax3.set_xlabel(r'$v_\perp$ (Mm s$^{-1}$)', fontsize=12)
ax3.set_ylabel(r'$f(v_\perp,\,0)$   (m$^{-3}$ (m s$^{-1})^{-3}$)', fontsize=11)
ax3.set_title(rf'1D cut at $v_\parallel \approx 0$  '
              rf'(log scale, zoom $0$--${XZOOM_MAX:.0f}$ Mm s$^{{-1}}$)', fontsize=11)
ax3.set_yscale('log')
ax3.set_xlim(0, XZOOM_MAX)
ax3.legend(fontsize=8.5, loc='lower left')
ax3.grid(True, alpha=0.25, which='both')
ax3.text(0.97, 0.55, r'Dashed verticals: $v_\perp = v_\mathrm{th}$',
         transform=ax3.transAxes, ha='right', va='bottom',
         fontsize=8, color='grey')

fig3.tight_layout()
out_path3 = Path(__file__).with_name('stix_maxw_at_vpar0_zoom.png')
fig3.savefig(out_path3, dpi=150, bbox_inches='tight')
print(f'Saved -> {out_path3}')

plt.show()
