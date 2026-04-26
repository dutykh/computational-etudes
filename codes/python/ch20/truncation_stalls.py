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

    # (D) NEW: full 2D error surface E(N, L) on a coarse grid, so the
    # three sweeps can be overlaid as cuts through one common surface.
    Ns_grid = np.array([8, 12, 16, 24, 32, 48, 64, 96, 128])
    Ls_grid = np.array([2, 3, 4, 5, 6, 8, 10, 12, 16, 20])
    err_grid = np.empty((len(Ls_grid), len(Ns_grid)))
    for i, L in enumerate(Ls_grid):
        for j, N in enumerate(Ns_grid):
            err_grid[i, j] = chebyshev_approx_error(N, L)

    fig, axes = plt.subplots(2, 2, figsize=(11.0, 8.0))

    # ---- (a) Fix L, vary N: plateau ---------------------------------
    ax = axes[0, 0]
    ax.semilogy(Ns_A, err_A, "-o", color=CORAL, lw=1.1, label=f"error at $L={L_A:.0f}$")
    ax.axhline(np.exp(-L_A), color=NAVY, ls=":", lw=1.0,
               label=fr"$e^{{-L}} \approx {np.exp(-L_A):.1e}$")
    ax.set_xlabel(r"$N$")
    ax.set_ylabel(r"$\|f - f_N\|_\infty$")
    ax.set_title(r"(a) Fix $L$, vary $N$: plateau")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(frameon=False, fontsize=9)

    # ---- (b) Fix N, vary L: sweet spot ------------------------------
    ax = axes[0, 1]
    ax.semilogy(Ls_B, err_B, "-s", color=TEAL, lw=1.1, mfc="none")
    ax.set_xlabel(r"$L$")
    ax.set_ylabel(r"$\|f - f_N\|_\infty$")
    ax.set_title(rf"(b) Fix $N={N_B}$, vary $L$: sweet spot")
    ax.grid(True, which="both", alpha=0.3)

    # ---- (c) Grow both: subgeometric descent ------------------------
    ax = axes[1, 0]
    ax.semilogy(Ns_C, err_C, "-^", color=NAVY, lw=1.1)
    ax.set_xlabel(r"$N$  ($L$ growing with $N$)")
    ax.set_ylabel(r"$\|f - f_N\|_\infty$")
    ax.set_title("(c) Grow both: subgeometric descent")
    ax.grid(True, which="both", alpha=0.3)

    # ---- (d) NEW: full error surface E(N, L) with three sweeps -----
    # Three slices through one surface: (a) is a horizontal cut at L=6;
    # (b) is a vertical cut at N=32; (c) is a diagonal trace through
    # the (N_k, L_k) pairs.
    ax = axes[1, 1]
    Lgrid_mesh, Ngrid_mesh = np.meshgrid(Ls_grid, Ns_grid, indexing="ij")
    log_err = np.log10(err_grid + 1e-18)
    im = ax.pcolormesh(Ngrid_mesh, Lgrid_mesh, log_err,
                       cmap="viridis", shading="auto",
                       vmin=log_err.min(), vmax=log_err.max())
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label(r"$\log_{10}\|f - f_N\|_\infty$", fontsize=8)
    cb.ax.tick_params(labelsize=8)

    # Overlay the three sweeps as polylines
    ax.plot(Ns_A, np.full_like(Ns_A, L_A), "-o", color=CORAL, lw=1.4,
            mfc="white", ms=4, label="(a) $L = 6$, vary $N$")
    ax.plot(np.full_like(Ls_B, N_B), Ls_B, "-s", color=TEAL, lw=1.4,
            mfc="white", ms=4, label=fr"(b) $N = {N_B}$, vary $L$")
    Ls_C = np.array([p[1] for p in pairs])
    ax.plot(Ns_C, Ls_C, "-^", color="white", lw=1.6, mfc="white", ms=5,
            label="(c) grow both jointly")

    ax.set_xlabel(r"$N$")
    ax.set_ylabel(r"$L$")
    ax.set_title("(d) the error surface $E(N, L)$ and its three slices",
                 fontsize=10)
    ax.legend(frameon=True, fontsize=8, loc="upper right",
              facecolor="white", edgecolor="none", framealpha=0.85)

    fig.tight_layout()
    save_fig(fig, OUTPUT_DIR, "truncation_stalls")
    plt.close(fig)

    print(f"[20.1] saved figure to {OUTPUT_DIR / 'truncation_stalls.pdf'}")
    print("  (A) fix L=6, vary N:", ["%.2e" % e for e in err_A])
    print("  (C) grow both:      ", ["%.2e" % e for e in err_C])


if __name__ == "__main__":
    make_figure()
