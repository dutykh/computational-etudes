#!/usr/bin/env python3
"""
semi_infinite_compare.py

Chapter 19: Coordinate Transformations
Computational Etude 19.4: Infinity without truncation.

We solve the semi-infinite boundary value problem

    u''(y) - u(y) = 0,    u(0) = 1,    u(y) -> 0 as y -> infinity,

whose exact solution is u_ex(y) = exp(-y), in two ways:

 (I)  Truncation to a finite interval [0, L] with u(L) = 0; Chebyshev
      collocation on [0, L]. This incurs BOTH a spectral-series error
      and a domain-truncation error (Boyd Sec 16.4). Here `L` is a
      truncation half-width.

 (II) Algebraic map   y = ell (1 + x) / (1 - x),   x in (-1, 1],   which
      takes y in [0, infinity) onto x in [-1, 1] and allows Chebyshev
      collocation on the full semi-infinite domain without any truncation.
      The map parameter `ell` (Greek script ell, distinct from the
      truncation L above) controls where the grid clusters.

The etude compares maximum-norm errors across a range of N, sweeping the
truncation length L for method (I) and the map parameter ell for method
(II), producing the "map versus truncation" picture that drives Boyd's
Principle on unbounded domains.

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

from map_common import (CORAL, GOLD, NAVY, ORANGE, PURPLE, TEAL,
                        output_dir_for, save_fig, setup_matplotlib)

OUTPUT_DIR = output_dir_for(SCRIPT_DIR)


def exact_u(y):
    return np.exp(-y)


# ---------------------------------------------------------------- (I) truncation
def solve_truncation(N, L):
    """Chebyshev on [0, L] with u(0) = 1 and u(L) = 0 imposed.
    Here `L` is the truncation half-width."""
    D, x = cheb_matrix(N)        # x in [-1, 1]
    # affine map x in [-1, 1]  ->  y in [0, L]
    y = 0.5 * L * (x + 1.0)
    # chain rule: d/dy = (2 / L) d/dx
    Dy = (2.0 / L) * D
    Dy2 = Dy @ Dy
    # interior block: y in (0, L)
    A = Dy2[1:N, 1:N] - np.eye(N - 1)
    # boundary contributions: u(y=L) at index 0, u(y=0) at index N
    # (because x_0 = 1 -> y = L, x_N = -1 -> y = 0 in this affine map)
    rhs = -Dy2[1:N, 0] * 0.0 - Dy2[1:N, N] * 1.0
    u_int = np.linalg.solve(A, rhs)
    u = np.zeros(N + 1)
    u[0] = 0.0       # y = L
    u[N] = 1.0       # y = 0
    u[1:N] = u_int
    return y, u


# ---------------------------------------------------------- (II) algebraic map
def solve_algebraic_map(N, ell):
    """Algebraic map y = ell (1 + x) / (1 - x). Interior grid only (x in [-1, 1)).
    The x = -1 endpoint maps to y = 0 (where u = 1); the x = 1 endpoint maps
    to y = +infty (where u -> 0); we impose u = 0 there by bordering.
    The map parameter `ell` controls grid clustering near the origin.
    """
    D, x = cheb_matrix(N)
    with np.errstate(divide="ignore", invalid="ignore"):
        # chain rule: dy/dx = 2 ell / (1 - x)^2,  d2y/dx2 = 4 ell / (1 - x)^3
        fp = 2.0 * ell / (1.0 - x) ** 2
        fpp = 4.0 * ell / (1.0 - x) ** 3
        Dx2 = D @ D
        # The x = 1 row is going to be overwritten by the Dirichlet BC u=0
        # anyway, so NaN entries at that row are harmless.
        Dy = np.diag(1.0 / fp) @ D
        Dy2 = np.diag(1.0 / fp ** 2) @ Dx2 - np.diag(fpp / fp ** 3) @ D
        y = ell * (1.0 + x) / (1.0 - x)
    # Impose u(x=-1) = 1 and u(x=1) = 0.  Interior block:
    A = Dy2[1:N, 1:N] - np.eye(N - 1)
    rhs = -Dy2[1:N, 0] * 0.0 - Dy2[1:N, N] * 1.0
    u_int = np.linalg.solve(A, rhs)
    u = np.zeros(N + 1)
    u[0] = 0.0        # x = 1  -> y = +infty
    u[N] = 1.0        # x = -1 -> y = 0
    u[1:N] = u_int
    return y, u


def max_error(y, u, clip_infty=1e6):
    mask = np.isfinite(y) & (y < clip_infty)
    return np.max(np.abs(u[mask] - exact_u(y[mask])))


def make_figure():
    setup_matplotlib()

    Ns = np.array([12, 16, 20, 24, 32, 40, 48, 64])
    L_trunc_list = [10, 20, 40]              # truncation half-widths
    ell_map_list = [1, 2, 4, 8]              # map parameters

    err_trunc = {L: [] for L in L_trunc_list}
    err_map = {ell: [] for ell in ell_map_list}

    for N in Ns:
        for L in L_trunc_list:
            y, u = solve_truncation(N, L)
            err_trunc[L].append(max_error(y, u))
        for ell in ell_map_list:
            y, u = solve_algebraic_map(N, ell)
            err_map[ell].append(max_error(y, u))

    fig, axes = plt.subplots(1, 3, figsize=(13.2, 3.6))

    # (a) sample solution and grids at N = 24
    ax = axes[0]
    y_t, u_t = solve_truncation(24, 20)
    y_m, u_m = solve_algebraic_map(24, 2.0)
    yplot = np.linspace(0, 12, 401)
    ax.plot(yplot, exact_u(yplot), color=NAVY, lw=1.2, label="exact")
    ax.plot(y_t, u_t, "o", color=CORAL, ms=4,
            label=r"truncation, $L=20$")
    # clip algebraic map for visibility
    y_m_clip = np.clip(y_m, 0, 12)
    ax.plot(y_m_clip, u_m, "s", color=TEAL, ms=4, mfc="none",
            label=r"algebraic, $\ell=2$")
    ax.set_xlim(0, 12)
    ax.set_ylim(-0.05, 1.1)
    ax.set_xlabel(r"$y$")
    ax.set_ylabel(r"$u(y)$")
    ax.set_title(r"(a) Solution at $N=24$")
    ax.legend(frameon=False, fontsize=9)
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    colours_trunc = [CORAL, ORANGE, GOLD]
    for L, c in zip(L_trunc_list, colours_trunc):
        ax.semilogy(Ns, err_trunc[L], "-o", color=c, label=fr"trunc. $L={L}$")
    ax.set_xlabel(r"$N$")
    ax.set_ylabel(r"$\|u - u_N\|_\infty$")
    ax.set_title("(b) Truncation: error vs $N$")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(frameon=False, fontsize=9)

    ax = axes[2]
    colours_map = [CORAL, TEAL, PURPLE, NAVY]
    for ell, c in zip(ell_map_list, colours_map):
        ax.semilogy(Ns, err_map[ell], "-s", color=c, mfc="none",
                    label=fr"algebraic $\ell={ell}$")
    ax.set_xlabel(r"$N$")
    ax.set_ylabel(r"$\|u - u_N\|_\infty$")
    ax.set_title(r"(c) Algebraic map: error vs $N$")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(frameon=False, fontsize=9)

    fig.tight_layout()
    save_fig(fig, OUTPUT_DIR, "semi_infinite_compare")
    plt.close(fig)

    print(f"[19.4] saved figure to {OUTPUT_DIR / 'semi_infinite_compare.pdf'}")
    for ell in ell_map_list:
        print(f"  algebraic map ell={ell}: errors = "
              + ", ".join(f"{e:.2e}" for e in err_map[ell]))


if __name__ == "__main__":
    make_figure()
