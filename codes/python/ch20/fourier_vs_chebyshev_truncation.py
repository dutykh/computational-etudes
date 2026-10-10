#!/usr/bin/env python3
"""
fourier_vs_chebyshev_truncation.py

Chapter 20: Spectral Methods on Unbounded Intervals.
Computational Etude 20.2: Fourier beats Chebyshev at the far boundary.

For domain truncation on an unbounded interval, Boyd (2000) argues that
an ordinary Fourier basis on [-L, L] is usually more efficient than the
Chebyshev basis:  Fourier has a factor of pi/2 higher density of
interior points than Chebyshev (which clusters near +/-L, i.e. where
the decaying solution is tiny).  The Gibbs oscillation at y = +/- L is
bounded by 0.09 |u(+/-L)| / N, smaller than the domain-truncation
error itself, so it is harmless.

This etude compares Fourier- and Chebyshev-based domain truncation of
f(y) = sech(y) on [-L, L] for a range of N at the optimal
representative L = 10.  We report max-norm error on the window
[-L + 1, L - 1] (to avoid counting the Gibbs overshoot in the Fourier
case, which is in any case smaller than the domain-truncation error).

Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR.parent / "ch07"))
from cheb_matrix import cheb_matrix  # noqa: E402

from unbounded_common import (CORAL, NAVY, TEAL, cheb_eval, dct1_coeffs,
                              output_dir_for, save_fig, setup_matplotlib)

OUTPUT_DIR = output_dir_for(SCRIPT_DIR)


def target(y):
    return 1.0 / np.cosh(y)


def chebyshev_error(N, L):
    _, x = cheb_matrix(N)
    y_nodes = L * x
    f_nodes = target(y_nodes)
    a = dct1_coeffs(f_nodes)
    y_fine = np.linspace(-L + 1.0, L - 1.0, 4001)
    approx = cheb_eval(a, y_fine / L, N)
    return np.max(np.abs(approx - target(y_fine)))


def fourier_error(N, L):
    """Fourier interpolation on [-L, L] with (2 * N)-point uniform grid.

    Pedagogically, we want to compare at MATCHED degrees of freedom.
    Use M = 2 * N points (so the trigonometric polynomial has the same
    dimension as a Chebyshev expansion of degree N).  The grid is
    y_j = -L + 2 * L * j / M  for j = 0, ..., M - 1.
    """
    M = 2 * N
    y_nodes = -L + 2.0 * L * np.arange(M) / M
    f_nodes = target(y_nodes)
    coeffs = np.fft.fft(f_nodes) / M
    k = np.fft.fftfreq(M, d=1.0 / M)         # wave numbers
    y_fine = np.linspace(-L + 1.0, L - 1.0, 4001)
    # re-parametrise: y -> t = pi * (y + L) / L, so that the Fourier
    # basis exp(i * k * t) matches the DFT of samples on a uniform grid.
    t_fine = np.pi * (y_fine + L) / L
    phase = np.exp(1j * k[:, None] * t_fine[None, :])
    approx = (coeffs[:, None] * phase).sum(axis=0).real
    return np.max(np.abs(approx - target(y_fine)))


def make_figure():
    setup_matplotlib()

    L = 10.0
    Ns = np.array([8, 12, 16, 24, 32, 48, 64, 96, 128])
    err_cheb = np.array([chebyshev_error(N, L) for N in Ns])
    err_four = np.array([fourier_error(N, L) for N in Ns])

    # grid densities for visual panel
    Nshow = 32
    _, x_GL = cheb_matrix(Nshow)
    y_GL = L * x_GL
    y_F = -L + 2.0 * L * np.arange(2 * Nshow) / (2 * Nshow)

    fig, axes = plt.subplots(1, 2, figsize=(10.0, 3.8))

    ax = axes[0]
    ax.scatter(y_F, np.full_like(y_F, -0.08), s=14, color=CORAL, marker="x",
               label=f"Fourier, $2N = {2*Nshow}$")
    ax.scatter(y_GL, np.full_like(y_GL, -0.18), s=20, color=TEAL,
               marker="o", facecolors="none", label=f"Chebyshev, $N = {Nshow}$")
    y_plot = np.linspace(-L, L, 401)
    ax.plot(y_plot, target(y_plot), color=NAVY, lw=1.1)
    ax.set_xlim(-L, L)
    ax.set_ylim(-0.25, 1.1)
    ax.set_xlabel(r"$y$")
    ax.set_title(fr"(a) Grid densities on $[-L,L]$, $L={L:.0f}$")
    ax.legend(frameon=False, fontsize=9)

    ax = axes[1]
    ax.semilogy(Ns, err_cheb, "-o", color=TEAL, lw=1.1, mfc="none",
                label="Chebyshev")
    ax.semilogy(Ns, err_four, "-s", color=CORAL, lw=1.1,
                label="Fourier")
    ax.set_xlabel(r"$N$")
    ax.set_ylabel(r"$\|f - f_N\|_\infty$ on $[-L+1, L-1]$")
    ax.set_title(r"(b) Interior error at fixed $L = 10$")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(frameon=False, fontsize=9)

    fig.tight_layout()
    save_fig(fig, OUTPUT_DIR, "fourier_vs_chebyshev_truncation")
    plt.close(fig)

    print(f"[20.2] saved figure to {OUTPUT_DIR / 'fourier_vs_chebyshev_truncation.pdf'}")
    print("  N  Cheb           Four")
    for N, c, f in zip(Ns, err_cheb, err_four):
        print(f"  {N:3d}  {c:.3e}    {f:.3e}")


if __name__ == "__main__":
    make_figure()
