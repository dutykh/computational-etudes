#!/usr/bin/env python3
"""
map1d_toolkit.py

Chapter 19: Coordinate Transformations
Computational Etude 19.3: Build a reusable one-dimensional map toolkit.

We encapsulate a one-dimensional change of coordinate

    y = f(x),     f: [a, b] ---> [c, d],

as a small Python class ``Map1D`` that stores forward map, inverse map,
and first and second derivatives of f.  The class delivers:

 * the physical grid  y_j = f(x_j)  from a computational grid  x_j;
 * the mapped first-derivative operator
       (d / d y) = (1 / f'(x)) (d / d x);
 * the mapped second-derivative operator
       (d^2 / d y^2) = (1 / f'(x)^2) (d^2 / d x^2)
                       - (f''(x) / f'(x)^3) (d / d x);
 * a mapped-quadrature helper that pushes Clenshaw-Curtis-on-x through
   the Jacobian weight.

The etude validates the class numerically on three manufactured cases:

 (A) identity map --- should reproduce Chebyshev on [-1, 1];
 (B) algebraic semi-infinite map  y = ell (1 + x)/(1 - x),  applied to
     u(y) = exp(-y)  and checked against the exact derivatives;
 (C) tanh map  y = tanh(x),  applied to  u(y) = 1 / (1 + y^2)  with
     all derivatives known analytically.

The resulting diagnostic figure plots the grid clustering and the
convergence of the mapped-derivative error for each case.

Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
Part of "Computational Etudes: A Spectral Approach"
"""
from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR.parent / "ch07"))
from cheb_matrix import cheb_matrix  # noqa: E402

from map_common import CORAL, NAVY, TEAL, output_dir_for, save_fig, setup_matplotlib

OUTPUT_DIR = output_dir_for(SCRIPT_DIR)


# --------------------------------------------------------------------- Map1D
@dataclass
class Map1D:
    """One-dimensional coordinate transformation y = f(x)."""
    forward:      Callable[[np.ndarray], np.ndarray]   # x -> y
    inverse:      Callable[[np.ndarray], np.ndarray]   # y -> x
    fprime:       Callable[[np.ndarray], np.ndarray]   # dy/dx
    fdoubleprime: Callable[[np.ndarray], np.ndarray]   # d2y/dx2

    def grid(self, x):
        return self.forward(x)

    def derivative_matrices(self, Dx):
        """Transform the (N+1) x (N+1) Chebyshev DM Dx into the mapped
        first and second-derivative matrices acting on values at y = f(x).
        """
        x = np.cos(np.pi * np.arange(Dx.shape[0]) / (Dx.shape[0] - 1))
        fp = self.fprime(x)
        fpp = self.fdoubleprime(x)
        Dx2 = Dx @ Dx
        Dy = np.diag(1.0 / fp) @ Dx
        Dy2 = np.diag(1.0 / fp ** 2) @ Dx2 - np.diag(fpp / fp ** 3) @ Dx
        return Dy, Dy2


# --------------------------------------------------------------- particular maps
def identity_map():
    return Map1D(
        forward=lambda x: x,
        inverse=lambda y: y,
        fprime=lambda x: np.ones_like(x),
        fdoubleprime=lambda x: np.zeros_like(x),
    )


def algebraic_semi_infinite(ell):
    """y = ell (1 + x) / (1 - x),  x in (-1, 1),  y in (0, +infty)."""
    def fwd(x): return ell * (1.0 + x) / (1.0 - x)
    def inv(y): return (y - ell) / (y + ell)
    def fp(x):  return 2.0 * ell / (1.0 - x) ** 2
    def fpp(x): return 4.0 * ell / (1.0 - x) ** 3
    return Map1D(fwd, inv, fp, fpp)


def tanh_map():
    """y = tanh(x),  x in (-infty, infty) but restricted to x in [-1, 1] here."""
    def fwd(x): return np.tanh(x)
    def inv(y): return np.arctanh(y)
    def fp(x):  return 1.0 / np.cosh(x) ** 2
    def fpp(x): return -2.0 * np.tanh(x) / np.cosh(x) ** 2
    return Map1D(fwd, inv, fp, fpp)


# ---------------------------------------------------------------- test problems
def case_exp_decay(ell):
    """Test (B): u(y) = exp(-y) on y in [0, infty) via algebraic map."""
    mp = algebraic_semi_infinite(ell)

    def u_of_y(y):  return np.exp(-y)
    def du_of_y(y): return -np.exp(-y)
    def d2u_of_y(y): return np.exp(-y)
    return mp, u_of_y, du_of_y, d2u_of_y


def case_rational(ell=1.0):
    """Test (C): u(y) = 1 / (1 + y^2) on y = tanh(x), x in [-1, 1]."""
    mp = tanh_map()

    def u_of_y(y):  return 1.0 / (1.0 + y ** 2)
    def du_of_y(y): return -2.0 * y / (1.0 + y ** 2) ** 2
    def d2u_of_y(y):
        return (6.0 * y ** 2 - 2.0) / (1.0 + y ** 2) ** 3
    return mp, u_of_y, du_of_y, d2u_of_y


# ----------------------------------------------------------------- convergence
def convergence(case_builder, Ns):
    err1, err2 = [], []
    for N in Ns:
        Dx, x = cheb_matrix(N)
        mp, u_of_y, du_of_y, d2u_of_y = case_builder()
        y = mp.grid(x)
        Dy, Dy2 = mp.derivative_matrices(Dx)
        u = u_of_y(y)
        du_num = Dy @ u
        d2u_num = Dy2 @ u
        # exclude the right-hand "infinity" endpoint where f is singular
        mask = np.isfinite(y) & (np.abs(y) < 1e6)
        err1.append(np.max(np.abs(du_num[mask] - du_of_y(y[mask]))))
        err2.append(np.max(np.abs(d2u_num[mask] - d2u_of_y(y[mask]))))
    return np.array(err1), np.array(err2)


def make_figure():
    setup_matplotlib()
    Ns = np.array([8, 12, 16, 20, 24, 32, 40, 48, 64])

    err_alg_1, err_alg_2 = convergence(lambda: case_exp_decay(2.0), Ns)
    err_tanh_1, err_tanh_2 = convergence(lambda: case_rational(), Ns)

    Ngrid = 24
    ELL = 2.0
    Y_LIM_ALG = 12.0
    Dx_demo, x_demo = cheb_matrix(Ngrid)
    x_line = np.linspace(-1.0, 1.0 - 1e-12, 401)  # exclude x=1 (algebraic blow-up)

    fig, axes = plt.subplots(2, 2, figsize=(11.0, 7.6))

    # (a) tanh map: y = tanh(x) on x in [-1, 1]
    ax = axes[0, 0]
    y_tanh_line = np.tanh(x_line)
    y_tanh_grid = np.tanh(x_demo)
    ax.plot(x_line, y_tanh_line, color=TEAL, lw=1.6,
            label=r"$y = \tanh(x)$")
    ax.scatter(x_demo, np.full_like(x_demo, -1.18),
               marker="|", color=NAVY, s=90, lw=1.4)
    ax.scatter(np.full_like(y_tanh_grid, 1.10), y_tanh_grid,
               marker="_", color=TEAL, s=90, lw=1.4)
    for j in range(0, Ngrid + 1, 3):
        ax.plot([x_demo[j], x_demo[j], 1.10],
                [-1.18, y_tanh_grid[j], y_tanh_grid[j]],
                color=TEAL, lw=0.4, ls=":", alpha=0.45)
    ax.axhline(0, color="0.85", lw=0.4)
    ax.axvline(0, color="0.85", lw=0.4)
    ax.set_xlim(-1.20, 1.20)
    ax.set_ylim(-1.30, 1.20)
    ax.set_xlabel(r"computational coordinate $x$")
    ax.set_ylabel(r"physical $y$")
    ax.set_title(r"(a) tanh map: $x_j$ (bottom) $\to y_j$ (right)")
    ax.legend(loc="upper left", frameon=False, fontsize=9)

    # (b) algebraic semi-infinite map: y = ell(1+x)/(1-x)
    ax = axes[0, 1]
    y_alg_line = ELL * (1 + x_line) / (1 - x_line)
    visible_line = y_alg_line <= Y_LIM_ALG
    ax.plot(x_line[visible_line], y_alg_line[visible_line],
            color=CORAL, lw=1.6, label=r"$y = \ell(1+x)/(1-x)$")
    x_finite = x_demo[:-1]
    y_alg_grid = ELL * (1 + x_finite) / (1 - x_finite)
    visible_grid = y_alg_grid <= Y_LIM_ALG
    ax.scatter(x_demo, np.full_like(x_demo, -0.7),
               marker="|", color=NAVY, s=90, lw=1.4)
    ax.scatter(np.full_like(y_alg_grid[visible_grid], 1.10),
               y_alg_grid[visible_grid],
               marker="_", color=CORAL, s=90, lw=1.4)
    for j, xj in enumerate(x_finite):
        if y_alg_grid[j] <= Y_LIM_ALG and j % 2 == 0:
            ax.plot([xj, xj, 1.10],
                    [-0.7, y_alg_grid[j], y_alg_grid[j]],
                    color=CORAL, lw=0.4, ls=":", alpha=0.45)
    n_off = int(np.sum(~visible_grid)) + 1   # +1 for x=1 endpoint at infinity
    ax.annotate(rf"$y \to \infty$ ({n_off} ticks beyond view)",
                xy=(1.0, Y_LIM_ALG), xytext=(0.10, Y_LIM_ALG - 1.6),
                fontsize=8, color=CORAL,
                arrowprops=dict(arrowstyle="->", color=CORAL, lw=0.6))
    ax.axhline(0, color="0.85", lw=0.4)
    ax.axvline(0, color="0.85", lw=0.4)
    ax.set_xlim(-1.20, 1.20)
    ax.set_ylim(-1.6, Y_LIM_ALG + 1.0)
    ax.set_xlabel(r"computational coordinate $x$")
    ax.set_ylabel(r"physical $y$")
    ax.set_title(r"(b) algebraic map ($\ell = 2$): clusters near $y = 0$")
    ax.legend(loc="upper left", frameon=False, fontsize=9)

    # (c) First-derivative convergence
    ax = axes[1, 0]
    ax.semilogy(Ns, err_alg_1 + 1e-18, "-o", color=CORAL,
                label="algebraic")
    ax.semilogy(Ns, err_tanh_1 + 1e-18, "-s", color=TEAL,
                label="tanh")
    ax.set_xlabel(r"$N$")
    ax.set_ylabel(r"$\|u'_N - u'\|_\infty$")
    ax.set_title(r"(c) First derivative: $D_y\, u$")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(frameon=False, fontsize=9)

    # (d) Second-derivative convergence -- the N^2-lag panel
    ax = axes[1, 1]
    ax.semilogy(Ns, err_alg_2 + 1e-18, "--o", color=CORAL, mfc="none",
                label="algebraic")
    ax.semilogy(Ns, err_tanh_2 + 1e-18, "--s", color=TEAL, mfc="none",
                label="tanh")
    ax.set_xlabel(r"$N$")
    ax.set_ylabel(r"$\|u''_N - u''\|_\infty$")
    ax.set_title(r"(d) Second derivative: $D_y^2\, u$ (N$^2$-lag)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(frameon=False, fontsize=9)

    fig.tight_layout()
    save_fig(fig, OUTPUT_DIR, "map1d_toolkit")
    plt.close(fig)

    print(f"[19.3] saved figure to {OUTPUT_DIR / 'map1d_toolkit.pdf'}")
    for N, a1, a2, t1, t2 in zip(Ns, err_alg_1, err_alg_2, err_tanh_1, err_tanh_2):
        print(f"  N={N:3d}  alg_D1={a1:.3e}  alg_D2={a2:.3e}  "
              f"tanh_D1={t1:.3e}  tanh_D2={t2:.3e}")


if __name__ == "__main__":
    make_figure()
