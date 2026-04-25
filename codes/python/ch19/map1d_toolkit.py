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

    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.6))

    # grid clustering
    ax = axes[0]
    Ngrid = 24
    Dx_demo, x_demo = cheb_matrix(Ngrid)
    y_alg = algebraic_semi_infinite(2.0).grid(x_demo)
    y_tanh = tanh_map().grid(x_demo)
    ax.plot(x_demo, np.full_like(x_demo, 0.0), "|", color=NAVY, ms=14,
            label="Chebyshev-GL $x_j$")
    # only show finite y for the algebraic map
    y_alg_plot = np.clip(y_alg, 0, 12)
    ax.plot(y_alg_plot, np.full_like(y_alg_plot, -0.55), "|", color=CORAL,
            ms=14, label=r"algebraic $y_j$ ($ell=2$)")
    ax.plot(y_tanh, np.full_like(y_tanh, -1.1), "|", color=TEAL, ms=14,
            label="tanh $y_j$")
    ax.set_xlim(-1.2, 12.2)
    ax.set_ylim(-1.5, 0.6)
    ax.set_xlabel("coordinate value")
    ax.set_title(f"(a) Grid clustering at $N = {Ngrid}$")
    ax.legend(loc="upper right", frameon=False, fontsize=9)
    ax.set_yticks([])

    ax = axes[1]
    ax.semilogy(Ns, err_alg_1 + 1e-18, "-o", color=CORAL,
                label=r"algebraic: $\|u'\|$")
    ax.semilogy(Ns, err_alg_2 + 1e-18, "--o", color=CORAL,
                mfc="none", label=r"algebraic: $\|u''\|$")
    ax.semilogy(Ns, err_tanh_1 + 1e-18, "-s", color=TEAL,
                label=r"$\tanh$: $\|u'\|$")
    ax.semilogy(Ns, err_tanh_2 + 1e-18, "--s", color=TEAL, mfc="none",
                label=r"$\tanh$: $\|u''\|$")
    ax.set_xlabel(r"$N$")
    ax.set_ylabel("max error of mapped derivative")
    ax.set_title("(b) Validation of mapped derivative formulas")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(frameon=False, fontsize=8)

    fig.tight_layout()
    save_fig(fig, OUTPUT_DIR, "map1d_toolkit")
    plt.close(fig)

    print(f"[19.3] saved figure to {OUTPUT_DIR / 'map1d_toolkit.pdf'}")
    for N, a1, a2, t1, t2 in zip(Ns, err_alg_1, err_alg_2, err_tanh_1, err_tanh_2):
        print(f"  N={N:3d}  alg_D1={a1:.3e}  alg_D2={a2:.3e}  "
              f"tanh_D1={t1:.3e}  tanh_D2={t2:.3e}")


if __name__ == "__main__":
    make_figure()
