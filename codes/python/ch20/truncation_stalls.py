#!/usr/bin/env python3
"""
truncation_stalls.py

Chapter 20: Spectral Methods on Unbounded Intervals.
Computational Etude 20.1: When spectral convergence stalls.

On an unbounded interval the user must choose TWO parameters, not one:
the truncation order N and a scale parameter, here the half-width L of
the truncation interval.  This etude exhibits the resulting trap.

We approximate  f(y) = sech(y)  on y in (-infty, +infty) by Chebyshev
collocation on [-L, L] with Dirichlet truncation.  Three sweeps:

  (A) fix L = 6, vary N:        error saturates at exp(-L) ~ 2.5e-3;
  (B) fix N = 32, vary L:       error rises with L for small L, then
                                drops with L (Chebyshev poles move
                                towards the real axis for too-large L);
  (C) grow L and N together:    error descends subgeometrically
                                (Boyd's Rule of Thumb 14).

The etude is the chapter's opening shock: ``more modes alone do not
mean better accuracy'' on an unbounded domain.

Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
Part of "Computational Etudes: A Spectral Approach"
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR.parent / "ch07"))
from cheb_matrix import cheb_matrix  # noqa: E402

from unbounded_common import (CORAL, NAVY, TEAL, output_dir_for,
                              save_fig, setup_matplotlib)

OUTPUT_DIR = output_dir_for(SCRIPT_DIR)


def target(y):
    return 1.0 / np.cosh(y)


def chebyshev_approx_error(N, L):
    """Chebyshev-GL interpolant of sech(y) on [-L, L]; max-norm error on
    a fine reference mesh on (-infty, +infty).  The truncation is enforced
    by treating the interpolant as identically zero outside [-L, L]; this
    is the right measure of the TOTAL error of a Dirichlet-truncated
    method.
    """
    _, x = cheb_matrix(N)
    y_nodes = L * x
    f_nodes = target(y_nodes)
    # interpolate on a dense y-mesh
    y_fine = np.linspace(-20.0, 20.0, 4001)
    truth = target(y_fine)
    # Chebyshev coefficients on [-L, L] (via Discrete Cosine -- local copy)
    from unbounded_common import dct1_coeffs, cheb_eval
    coeffs = dct1_coeffs(f_nodes)
    in_window = np.abs(y_fine) <= L
    x_in = y_fine[in_window] / L
    approx = np.zeros_like(y_fine)
    approx[in_window] = cheb_eval(coeffs, x_in, N)
    return np.max(np.abs(approx - truth))


def make_figure():
    setup_matplotlib()

    # (A) fix L, vary N
    L_A = 6.0
    Ns_A = np.array([8, 12, 16, 24, 32, 48, 64, 96, 128])
    err_A = np.array([chebyshev_approx_error(N, L_A) for N in Ns_A])

    # (B) fix N, vary L
    N_B = 32
    Ls_B = np.array([2, 3, 4, 5, 6, 8, 10, 12, 16, 20])
    err_B = np.array([chebyshev_approx_error(N_B, L) for L in Ls_B])

    # (C) increase both
    pairs = [(8, 3), (12, 4), (16, 5), (24, 6), (32, 8), (48, 10), (64, 12), (96, 14)]
    err_C = np.array([chebyshev_approx_error(N, L) for N, L in pairs])
    Ns_C = np.array([p[0] for p in pairs])

    fig, axes = plt.subplots(1, 3, figsize=(13.4, 3.8))

    ax = axes[0]
    ax.semilogy(Ns_A, err_A, "-o", color=CORAL, lw=1.1, label=f"error at $L={L_A:.0f}$")
    ax.axhline(np.exp(-L_A), color=NAVY, ls=":", lw=1.0,
               label=fr"$e^{{-L}} \approx {np.exp(-L_A):.1e}$")
    ax.set_xlabel(r"$N$")
    ax.set_ylabel(r"$\|f - f_N\|_\infty$")
    ax.set_title(r"(a) Fix $L$, vary $N$: plateau")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(frameon=False, fontsize=9)

    ax = axes[1]
    ax.semilogy(Ls_B, err_B, "-s", color=TEAL, lw=1.1, mfc="none")
    ax.set_xlabel(r"$L$")
    ax.set_ylabel(r"$\|f - f_N\|_\infty$")
    ax.set_title(rf"(b) Fix $N={N_B}$, vary $L$: sweet spot")
    ax.grid(True, which="both", alpha=0.3)

    ax = axes[2]
    ax.semilogy(Ns_C, err_C, "-^", color=NAVY, lw=1.1)
    ax.set_xlabel(r"$N$  ($L$ growing with $N$)")
    ax.set_ylabel(r"$\|f - f_N\|_\infty$")
    ax.set_title("(c) Grow both: subgeometric descent")
    ax.grid(True, which="both", alpha=0.3)

    fig.tight_layout()
    save_fig(fig, OUTPUT_DIR, "truncation_stalls")
    plt.close(fig)

    print(f"[20.1] saved figure to {OUTPUT_DIR / 'truncation_stalls.pdf'}")
    print("  (A) fix L=6, vary N:", ["%.2e" % e for e in err_A])
    print("  (C) grow both:      ", ["%.2e" % e for e in err_C])


if __name__ == "__main__":
    make_figure()
