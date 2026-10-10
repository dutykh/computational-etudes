#!/usr/bin/env python3
"""
arctan_tan_sweep.py

Chapter 19: Coordinate Transformations
Computational Etude 19.7: A localised periodic pulse in the right coordinate.

We approximate the same sharply-localised 2-pi-periodic pulse as in
Etude 19.1,

    f(y) = exp( -kappa * (1 - cos y) ),    y in [-pi, pi],

this time sweeping the map parameter ell over a range of values for each
N, so that the student can SEE the arctan/tan map's parameter landscape.
The map used is the 2-pi-periodic variant of Boyd's Eq. (16.27)
(pedagogically identical to his formula but with period 2 pi in y):

    y = 2 * arctan( ell * tan(x / 2) ),    x, y in [-pi, pi].

For each (N, ell) pair we compute the ell-infinity error of the trigonometric
interpolant against the analytic target and record the result.  The
resulting error(N, ell) map is displayed as a heat-map in (N, ell) space
together with slices at fixed N.  The etude's take-away is that ell is
the degree of freedom the method gives you -- and that the range of
GOOD ell is broad, not sharp, so parameter tuning need not be painful.

Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
Part of "Computational Etudes: A Spectral Approach"
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from map_common import (CORAL, NAVY, ORANGE, PURPLE, TEAL, output_dir_for,
                        save_fig, setup_matplotlib)

SCRIPT_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = output_dir_for(SCRIPT_DIR)

KAPPA = 80.0


def target(y):
    return np.exp(-KAPPA * (1.0 - np.cos(y)))


def arctan_tan_map(x, ell):
    return 2.0 * np.arctan(ell * np.tan(x / 2.0))


def arctan_tan_inverse(y, ell):
    return 2.0 * np.arctan(np.tan(y / 2.0) / ell)


def mapped_interp(x_nodes, f_nodes, y_eval, ell):
    """Evaluate the mapped Fourier interpolant at arbitrary y_eval."""
    N = len(x_nodes)
    coeffs = np.fft.fft(f_nodes) / N
    k = np.fft.fftfreq(N, d=1.0 / N)
    x_eval = arctan_tan_inverse(y_eval, ell)
    phase = np.exp(1j * k[:, None] * (x_eval[None, :] + np.pi))
    return (coeffs[:, None] * phase).sum(axis=0).real


def error_for(N, ell, y_eval, truth):
    x = -np.pi + 2.0 * np.pi * np.arange(N) / N
    y = arctan_tan_map(x, ell)
    fv = target(y)
    approx = mapped_interp(x, fv, y_eval, ell)
    return np.max(np.abs(approx - truth))


def make_figure():
    setup_matplotlib()

    Ns = np.array([12, 16, 24, 32, 48, 64, 96])
    Ls = np.array([0.08, 0.10, 0.15, 0.20, 0.30, 0.45, 0.60, 0.80, 1.0, 1.5])
    y_eval = np.linspace(-np.pi, np.pi, 4097)
    y_eval[0] += 1e-9
    y_eval[-1] -= 1e-9
    truth = target(y_eval)

    E = np.zeros((len(Ns), len(Ls)))
    for i, N in enumerate(Ns):
        for j, ell in enumerate(Ls):
            E[i, j] = error_for(N, ell, y_eval, truth)

    fig, axes = plt.subplots(2, 2, figsize=(11.0, 7.6))

    # (a) three mapped grids at N = 32 for different ell
    ax = axes[0, 0]
    y_line = np.linspace(-np.pi, np.pi, 401)
    ax.plot(y_line, target(y_line), color=NAVY, lw=1.4, label=r"$f(y)$")
    Ngrid = 32
    for (ell, c, off) in [(0.1, CORAL, -0.08),
                          (0.3, TEAL, -0.18),
                          (1.0, ORANGE, -0.28)]:
        x = -np.pi + 2.0 * np.pi * np.arange(Ngrid) / Ngrid
        y = arctan_tan_map(x, ell)
        ax.scatter(y, np.full_like(y, off), s=22, color=c, marker="o",
                   facecolors="none", label=fr"$\ell = {ell}$")
    ax.set_xlim(-np.pi, np.pi)
    ax.set_ylim(-0.35, 1.18)
    ax.set_xlabel(r"physical $y$")
    ax.set_title(fr"(a) Mapped grids in physical $y$, $N = {Ngrid}$")
    ax.legend(loc="upper right", frameon=False, fontsize=9)

    # (b) NEW: pulse f(y) vs broadened image f(y(x)) at the chosen ell* = 0.3
    ax = axes[0, 1]
    ELL_STAR = 0.3
    x_line = np.linspace(-np.pi + 1e-9, np.pi - 1e-9, 401)
    f_of_x = target(arctan_tan_map(x_line, ELL_STAR))
    ax.plot(y_line, target(y_line), color="0.65", lw=1.0, ls="--",
            label=r"original $f(y)$")
    ax.plot(x_line, f_of_x, color=TEAL, lw=1.6,
            label=r"$\tilde f(x) = f(y(x))$, $\ell = 0.3$")
    # uniform x-grid below (this becomes the mapped y-grid in panel (a))
    x_grid = -np.pi + 2.0 * np.pi * np.arange(Ngrid) / Ngrid
    ax.scatter(x_grid, np.full_like(x_grid, -0.18),
               s=22, color=TEAL, marker="o", facecolors="none",
               label="uniform $x$-grid")
    ax.set_xlim(-np.pi, np.pi)
    ax.set_ylim(-0.35, 1.18)
    ax.set_xlabel(r"computational coordinate $x$")
    ax.set_title(r"(b) Pulse in computational $x$ at $\ell^* = 0.3$")
    ax.legend(loc="upper right", frameon=False, fontsize=9)

    # (c) heat-map of error(N, ell)
    ax = axes[1, 0]
    im = ax.pcolormesh(Ls, Ns, np.log10(E + 1e-16),
                       shading="auto", cmap="magma_r")
    plt.colorbar(im, ax=ax, label=r"$\log_{10} \|f - f_N\|_\infty$")
    ax.set_xscale("log")
    ax.set_xlabel(r"map parameter $\ell$")
    ax.set_ylabel(r"$N$")
    ax.set_title("(c) Error landscape")

    # (d) slices at fixed N
    ax = axes[1, 1]
    colours = [CORAL, ORANGE, TEAL, NAVY, PURPLE]
    slice_Ns = [16, 24, 32, 48, 64]
    for N_sl, c in zip(slice_Ns, colours):
        i = int(np.where(Ns == N_sl)[0][0])
        ax.semilogy(Ls, E[i, :] + 1e-16, "-o", color=c, ms=4,
                    label=fr"$N = {N_sl}$")
    ax.set_xscale("log")
    ax.set_xlabel(r"$\ell$")
    ax.set_ylabel(r"$\|f - f_N\|_\infty$")
    ax.set_title(r"(d) Slices at fixed $N$")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(frameon=False, fontsize=9)

    fig.tight_layout()
    save_fig(fig, OUTPUT_DIR, "arctan_tan_sweep")
    plt.close(fig)

    print(f"[19.7] saved figure to {OUTPUT_DIR / 'arctan_tan_sweep.pdf'}")
    print("  best ell at each N:")
    for i, N in enumerate(Ns):
        j_best = np.argmin(E[i, :])
        print(f"    N={N:3d}  ell*={Ls[j_best]:.3f}  err={E[i, j_best]:.2e}")


if __name__ == "__main__":
    make_figure()
