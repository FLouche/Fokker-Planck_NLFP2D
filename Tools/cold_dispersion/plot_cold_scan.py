"""Plot radial scans written by cold_dispersion.exe (cold_dispersion-<case>.txt).

Usage:  python plot_cold_scan.py scan1.txt [scan2.txt ...] [-o out.png]
One column per scan; rows: fast-wave n_perp^2, k_perp, |E-/E+|, each with the
minority in the dielectric and without it (bulk only). Vertical lines: minority
cyclotron harmonics (solid) and the position of B = b0 (dashed); horizontal
dotted lines: the namelist kperp and |emin/eplus|.
"""
import argparse, re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

C = 2.99792458e8
ALL, BULK = '#2a6fdb', '#e8590c'      # with minority / bulk only
INK2, INK3 = '#555555', '#999999'


def load(fn):
    head = ''.join(l for l in open(fn) if l.startswith('#'))
    g = lambda k: float(re.search(k + r'\s*=\s*([-\d.Ee+]+)', head)[1])
    meta = dict(b0=g('b0'), r0=g('r0'), f=g('f') * 1e6, kpar=g('kpar'), kperp=g('namelist kperp'),
                pol=g(r'\|emin/eplus\|'), nharm=int(g('nharm')),
                case=re.search(r'case\s+(\S+)', head)[1])
    return meta, np.loadtxt(fn)


def panel(axc, meta, t):
    R, wOm = t[:, 0], t[:, 2]
    w = 2 * np.pi * meta['f']
    x_all, x_bulk = t[:, 3], t[:, 10]
    k_all = np.where(x_all > 0, w / C * np.sqrt(np.abs(x_all)), np.nan)
    k_bulk = np.where(x_bulk > 0, w / C * np.sqrt(np.abs(x_bulk)), np.nan)
    ax1, ax2, ax3 = axc
    ax1.plot(R, x_all, color=ALL, lw=2, label='with minority')
    ax1.plot(R, x_bulk, color=BULK, lw=2, label='bulk only')
    ax1.set_yscale('symlog', linthresh=100)
    ax1.axhline(0, color=INK3, lw=0.8)
    ax1.set_ylim(-1e5, 1e5)
    ax2.plot(R, k_all, color=ALL, lw=2)
    ax2.plot(R, k_bulk, color=BULK, lw=2)
    ax2.axhline(meta['kperp'], color=INK2, ls=':', lw=1.5)
    ax2.text(R[-1], meta['kperp'], f"namelist kperp = {meta['kperp']:g}", ha='right', va='bottom',
             fontsize=8, color=INK2)
    ax2.set_ylim(0, 120)
    ax3.plot(R, np.abs(t[:, 9]), color=ALL, lw=2)
    ax3.plot(R, np.abs(t[:, 11]), color=BULK, lw=2)
    if meta['pol'] > 0:
        ax3.axhline(meta['pol'], color=INK2, ls=':', lw=1.5)
        ax3.text(R[-1], meta['pol'], f"namelist |E-/E+| = {meta['pol']:g}", ha='right', va='bottom',
                 fontsize=8, color=INK2)
    ax3.set_yscale('log'); ax3.set_ylim(0.1, 100)
    # minority harmonics: omega / Omega_min = N
    for N in range(1, 5):
        i = np.where(np.diff(np.sign(wOm - N)))[0]
        for j in i:
            for ax in axc:
                ax.axvline(R[j], color=INK3, lw=1)
            ax1.text(R[j], 6e4, f' {N}$\\Omega_H$', fontsize=8, color=INK2, va='top')
    for ax in axc:
        ax.axvline(meta['r0'], color=INK2, ls='--', lw=1)
        ax.grid(alpha=.25)
    ax1.text(meta['r0'], -6e4, ' B = b0', fontsize=8, color=INK2)
    ax1.set_title(f"{meta['case']}\nB0 = {meta['b0']:g} T at R0 = {meta['r0']:g} m, "
                  f"f = {meta['f'] / 1e6:g} MHz, k$_\\parallel$ = {meta['kpar']:g} m$^{{-1}}$", fontsize=9)
    ax3.set_xlabel('R (m)')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('scans', nargs='+')
    ap.add_argument('-o', '--out', default='cold_scan.png')
    a = ap.parse_args()
    n = len(a.scans)
    fig, axs = plt.subplots(3, n, figsize=(5.5 * n, 9), sharex='col', squeeze=False)
    for c, fn in enumerate(a.scans):
        panel(axs[:, c], *load(fn))
    axs[0, 0].set_ylabel(r'Re $n_\perp^2$ fast wave (symlog)')
    axs[1, 0].set_ylabel(r'$k_\perp$ fast wave (m$^{-1}$), propagating only')
    axs[2, 0].set_ylabel(r'$|E_-/E_+|$ fast wave')
    axs[0, 0].legend(loc='lower left', fontsize=8)
    fig.tight_layout()
    fig.savefig(a.out, dpi=150)
    print('wrote', a.out)


if __name__ == '__main__':
    main()
