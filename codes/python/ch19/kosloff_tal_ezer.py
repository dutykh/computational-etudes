#!/usr/bin/env python3
"""
kosloff_tal_ezer.py

Chapter 19: Coordinate Transformations
Computational Etude 19.8: Accuracy versus timestep.

The Kosloff-Tal-Ezer (KTE) arcsine map (Boyd Eq. (16.33))

    y = arcsin((1 - beta) x) / arcsin(1 - beta),     x in [-1, 1],

distributes the Chebyshev grid points more uniformly in y and dramatically
reduces the stiffness of the first-derivative matrix, which tightens the
explicit-time-step restriction of spectral methods.

But there is a price.  Boyd's Theorem 31 (BRANCH POINTS OF THE KTE MAP):
the same beta scaling that cures stiffness also drives the branch points
of the arcsin map toward the interval [-1, 1], which means the asymptotic
rate of geometric convergence

    rate ~ exp(-sqrt(2 C))       with beta = C / N^2 + O(N^-3)

DOES NOT DECREASE with N.  Spectral accuracy has been traded for a
bigger timestep.

This etude reproduces the two halves of Boyd's narrative:

 (i)  Minimum grid spacing versus N under beta = C / N^2 and beta = fixed.
 (ii) Chebyshev coefficient magnitudes for the linear function f(y) = y
      as a function of N, for beta_optimum = 1 - cos(1/N)  (aggressive) and
      beta_conservative = 1 - cos(1/2)  (Hesthaven-Dinesen-Lynov 1999).
      The aggressive choice shows the Mead-Renaut asymptote |a_N| ~ 0.488 / N^2
      (Boyd Fig. 16.3 reproduced), while the conservative choice preserves
      geometric coefficient decay at the cost of only doubling the timestep.

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

from map_common import CORAL, NAVY, TEAL, output_dir_for, save_fig, setup_matplotlib

OUTPUT_DIR = output_dir_for(SCRIPT_DIR)


def kte_grid(N, beta):
    """Return the physical y-grid of the KTE mapping."""
    D, xi = cheb_matrix(N)
    denom = np.arcsin(1.0 - beta)
    y = np.arcsin((1.0 - beta) * xi) / denom
    return y, xi, D


def kte_derivative_matrix(N, beta):
    """First-derivative matrix in the physical y-coordinate."""
    y, xi, D = kte_grid(N, beta)
    denom = np.arcsin(1.0 - beta)
    # dy / dxi = (1 - beta) / ( sqrt(1 - (1-beta)^2 xi^2) * arcsin(1 - beta) )
    fp = (1.0 - beta) / (np.sqrt(1.0 - (1.0 - beta) ** 2 * xi ** 2) * denom)
    Dy = np.diag(1.0 / fp) @ D
    return y, xi, D, Dy


def chebyshev_coeffs(v):
    """DCT-I coefficients on GL nodes."""
    N = len(v) - 1
    V = np.concatenate([v, v[N - 1:0:-1]])
    A = np.real(np.fft.fft(V)) / N
    A[0] *= 0.5
    A[N] *= 0.5
    return A[:N + 1]


def min_grid_spacing(y):
    ys = np.sort(y)
    return np.min(np.diff(ys))


def make_figure():
    setup_matplotlib()

    Ns = np.arange(8, 97, 2)

    # (i) min grid spacing under different beta laws
    # Also record stiffness = max |eigenvalue| of the first-derivative matrix
    min_dy_std, min_dy_opt, min_dy_cons = [], [], []
    stiff_std, stiff_opt, stiff_cons = [], [], []
    for N in Ns:
        # standard Chebyshev (beta = 0)
        y, xi, D = kte_grid(N, 0.0)
        min_dy_std.append(min_grid_spacing(xi))
        stiff_std.append(np.max(np.abs(np.linalg.eigvals(D))))
        # Mead-Renaut "optimum" beta ~ 1/(2 N^2)
        beta_opt = 1.0 - np.cos(1.0 / N)
        y, _, _, Dy = kte_derivative_matrix(N, beta_opt)
        min_dy_opt.append(min_grid_spacing(y))
        stiff_opt.append(np.max(np.abs(np.linalg.eigvals(Dy))))
        # Hesthaven conservative beta = 1 - cos(1/2), N-independent
        beta_cons = 1.0 - np.cos(0.5)
        y, _, _, Dy = kte_derivative_matrix(N, beta_cons)
        min_dy_cons.append(min_grid_spacing(y))
        stiff_cons.append(np.max(np.abs(np.linalg.eigvals(Dy))))

    min_dy_std = np.array(min_dy_std)
    min_dy_opt = np.array(min_dy_opt)
    min_dy_cons = np.array(min_dy_cons)
    stiff_std = np.array(stiff_std)
    stiff_opt = np.array(stiff_opt)
    stiff_cons = np.array(stiff_cons)

    # (ii) Chebyshev coefficient |a_N| of f(y) = y under aggressive KTE map
    Ns_coeff = np.arange(3, 64, 2)           # only odd N where a_N != 0
    aN_aggressive, aN_conservative = [], []
    for N in Ns_coeff:
        # sample f(y) = y on the KTE grid, then compute its Chebyshev coefficients
        # (the map parameter is different for each N in the aggressive run)
        beta_agg = 1.0 - np.cos(1.0 / N)
        y_agg, _, _ = kte_grid(N, beta_agg)
        f_agg = y_agg
        a_agg = chebyshev_coeffs(f_agg)
        aN_aggressive.append(np.abs(a_agg[-1]))

        beta_con = 1.0 - np.cos(0.5)
        y_con, _, _ = kte_grid(N, beta_con)
        f_con = y_con
        a_con = chebyshev_coeffs(f_con)
        aN_conservative.append(np.abs(a_con[-1]))

    aN_aggressive = np.array(aN_aggressive)
    aN_conservative = np.array(aN_conservative)

    fig, axes = plt.subplots(2, 2, figsize=(11.0, 7.6))

    # (a) NEW: KTE map shape y(x) for the three beta regimes at a fixed N
    ax = axes[0, 0]
    Ndemo = 32
    _, xi_demo = cheb_matrix(Ndemo)
    x_line = np.linspace(-1, 1, 401)
    # standard Chebyshev (beta = 0): y = x identity
    ax.plot(x_line, x_line, color=NAVY, lw=1.4,
            label=r"standard ($\beta = 0$): identity")
    # aggressive at this N
    beta_a = 1.0 - np.cos(1.0 / Ndemo)
    y_a = np.arcsin((1 - beta_a) * x_line) / np.arcsin(1 - beta_a)
    ax.plot(x_line, y_a, color=CORAL, lw=1.6,
            label=fr"aggressive $\beta \sim 1/N^2$ ($N = {Ndemo}$)")
    # conservative beta = 1 - cos(1/2)
    beta_c = 1.0 - np.cos(0.5)
    y_c = np.arcsin((1 - beta_c) * x_line) / np.arcsin(1 - beta_c)
    ax.plot(x_line, y_c, color=TEAL, lw=1.6,
            label=r"conservative $\beta = 1 - \cos(1/2)$")
    # CGL nodes ticks at the bottom
    ax.scatter(xi_demo, np.full_like(xi_demo, -1.12),
               marker="|", color=NAVY, s=70, lw=1.2)
    ax.text(-1.0, -1.30, "CGL nodes $\\xi_j$", fontsize=8, color=NAVY)
    ax.axhline(0, color="0.85", lw=0.4)
    ax.axvline(0, color="0.85", lw=0.4)
    ax.set_xlim(-1.10, 1.10)
    ax.set_ylim(-1.40, 1.10)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(r"computational coordinate $\xi$")
    ax.set_ylabel(r"physical coordinate $y$")
    ax.set_title(r"(a) KTE map $y = \arcsin((1 - \beta)\xi) / \arcsin(1 - \beta)$")
    ax.legend(loc="upper left", frameon=False, fontsize=8)

    # (b) Minimum spacing vs N
    ax = axes[0, 1]
    ax.loglog(Ns, min_dy_std, "-o", color=NAVY, ms=3,
              label="standard Chebyshev")
    ax.loglog(Ns, min_dy_opt, "-s", color=CORAL, ms=3, mfc="none",
              label=r"KTE $\beta \sim 1/N^2$")
    ax.loglog(Ns, min_dy_cons, "-^", color=TEAL, ms=3, mfc="none",
              label=r"KTE $\beta = 1 - \cos(1/2)$")
    ax.loglog(Ns, 2.0 / Ns, ":", color="gray", lw=0.8, label=r"$2/N$")
    ax.set_xlabel(r"$N$")
    ax.set_ylabel(r"$\min_j |y_{j+1} - y_j|$")
    ax.set_title(r"(b) Minimum spacing vs $N$")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(frameon=False, fontsize=8)

    # (c) Spectral radius of first-derivative matrix
    ax = axes[1, 0]
    ax.loglog(Ns, stiff_std, "-o", color=NAVY, ms=3,
              label="standard Chebyshev")
    ax.loglog(Ns, stiff_opt, "-s", color=CORAL, ms=3, mfc="none",
              label=r"KTE $\beta \sim 1/N^2$")
    ax.loglog(Ns, stiff_cons, "-^", color=TEAL, ms=3, mfc="none",
              label=r"KTE $\beta = 1 - \cos(1/2)$")
    ax.set_xlabel(r"$N$")
    ax.set_ylabel(r"$\rho(D)$")
    ax.set_title(r"(c) Stiffness $\rho(D)$ vs $N$")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(frameon=False, fontsize=8)

    # (d) |a_N| for f(y) = y -- accuracy destroyed by aggressive scaling
    ax = axes[1, 1]
    ax.loglog(Ns_coeff, aN_aggressive + 1e-18, "o", color=CORAL, ms=4,
              label=r"aggressive $\beta = 1 - \cos(1/N)$")
    ax.loglog(Ns_coeff, aN_conservative + 1e-18, "s", color=TEAL, ms=4,
              mfc="none", label=r"conservative $\beta = 1 - \cos(1/2)$")
    ax.loglog(Ns_coeff, 0.488 / Ns_coeff ** 2, ":", color=NAVY, lw=0.9,
              label=r"$0.488 \, / \, N^2$ (Boyd)")
    ax.set_xlabel(r"$N$")
    ax.set_ylabel(r"$|a_N|$ of $f(y) = y$")
    ax.set_title("(d) Accuracy destroyed by aggressive scaling")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(frameon=False, fontsize=8)

    fig.tight_layout()
    save_fig(fig, OUTPUT_DIR, "kosloff_tal_ezer")
    plt.close(fig)

    print(f"[19.8] saved figure to {OUTPUT_DIR / 'kosloff_tal_ezer.pdf'}")
    for N, a, c in zip(Ns_coeff[::4], aN_aggressive[::4], aN_conservative[::4]):
        print(f"  N={N:3d}  |a_N| aggressive={a:.3e}  conservative={c:.3e}")


if __name__ == "__main__":
    make_figure()
