#!/usr/bin/env python3
"""
arctan_tan_sweep.py

Chapter 19: Coordinate Transformations
Computational Etude 19.7: A localised periodic pulse in the right coordinate.

We approximate the same sharply-localised 2-pi-periodic pulse as in
Etude 19.1,

    f(y) = exp( -kappa * (1 - cos y) ),    y in [-pi, pi],

this time sweeping the map parameter L over a range of values for each
N, so that the student can SEE the arctan/tan map's parameter landscape.
The map used is the 2-pi-periodic variant of Boyd's Eq. (16.27)
(pedagogically identical to his formula but with period 2 pi in y):

    y = 2 * arctan( L * tan(x / 2) ),    x, y in [-pi, pi].

For each (N, L) pair we compute the L-infinity error of the trigonometric
interpolant against the analytic target and record the result.  The
resulting error(N, L) map is displayed as a heat-map in (N, L) space
together with slices at fixed N.  The etude's take-away is that L is
the degree of freedom the method gives you -- and that the range of
GOOD L is broad, not sharp, so parameter tuning need not be painful.

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


def arctan_tan_map(x, L):
    return 2.0 * np.arctan(L * np.tan(x / 2.0))


def arctan_tan_inverse(y, L):
    return 2.0 * np.arctan(np.tan(y / 2.0) / L)


def mapped_interp(x_nodes, f_nodes, y_eval, L):
    """Evaluate the mapped Fourier interpolant at arbitrary y_eval."""
    N = len(x_nodes)
    coeffs = np.fft.fft(f_nodes) / N
    k = np.fft.fftfreq(N, d=1.0 / N)
    x_eval = arctan_tan_inverse(y_eval, L)
    phase = np.exp(1j * k[:, None] * (x_eval[None, :] + np.pi))
    return (coeffs[:, None] * phase).sum(axis=0).real


def error_for(N, L, y_eval, truth):
    x = -np.pi + 2.0 * np.pi * np.arange(N) / N
    y = arctan_tan_map(x, L)
    fv = target(y)
    approx = mapped_interp(x, fv, y_eval, L)
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
        for j, L in enumerate(Ls):
            E[i, j] = error_for(N, L, y_eval, truth)

    fig, axes = plt.subplots(1, 3, figsize=(13.2, 3.8))

    # (a) three mapped grids at N = 32 for different L
    ax = axes[0]
    y_line = np.linspace(-np.pi, np.pi, 401)
    ax.plot(y_line, target(y_line), color=NAVY, lw=1.1, label=r"$f(y)$")
    Ngrid = 32
    for (L, c, off) in [(0.1, CORAL, -0.08),
                        (0.3, TEAL, -0.18),
                        (1.0, ORANGE, -0.28)]:
        x = -np.pi + 2.0 * np.pi * np.arange(Ngrid) / Ngrid
        y = arctan_tan_map(x, L)
        ax.scatter(y, np.full_like(y, off), s=22, color=c, marker="o",
                   facecolors="none", label=fr"$L={L}$")
    ax.set_xlim(-np.pi, np.pi)
    ax.set_ylim(-0.35, 1.1)
    ax.set_xlabel(r"$y$")
    ax.set_title(f"(a) Mapped grids, $N={Ngrid}$")
    ax.legend(frameon=False, fontsize=9)

    # (b) heat-map of error(N, L)
    ax = axes[1]
    im = ax.pcolormesh(Ls, Ns, np.log10(E + 1e-16), shading="auto", cmap="magma_r")
    plt.colorbar(im, ax=ax, label=r"$\log_{10} \|f - f_N\|_\infty$")
    ax.set_xscale("log")
    ax.set_xlabel(r"map parameter $L$")
    ax.set_ylabel(r"$N$")
    ax.set_title("(b) Error landscape")

    # (c) slices at fixed N
    ax = axes[2]
    colours = [CORAL, ORANGE, TEAL, NAVY, PURPLE]
    slice_Ns = [16, 24, 32, 48, 64]
    for N_sl, c in zip(slice_Ns, colours):
        i = int(np.where(Ns == N_sl)[0][0])
        ax.semilogy(Ls, E[i, :] + 1e-16, "-o", color=c, ms=4,
                    label=fr"$N={N_sl}$")
    ax.set_xscale("log")
    ax.set_xlabel(r"$L$")
    ax.set_ylabel(r"$\|f - f_N\|_\infty$")
    ax.set_title("(c) Slices at fixed $N$")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(frameon=False, fontsize=9)

    fig.tight_layout()
    save_fig(fig, OUTPUT_DIR, "arctan_tan_sweep")
    plt.close(fig)

    print(f"[19.7] saved figure to {OUTPUT_DIR / 'arctan_tan_sweep.pdf'}")
    print("  best L at each N:")
    for i, N in enumerate(Ns):
        j_best = np.argmin(E[i, :])
        print(f"    N={N:3d}  L*={Ls[j_best]:.3f}  err={E[i, j_best]:.2e}")


if __name__ == "__main__":
    make_figure()
