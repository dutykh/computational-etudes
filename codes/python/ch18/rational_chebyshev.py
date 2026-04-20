#!/usr/bin/env python3
"""
rational_chebyshev.py

Chapter 18: Linear Spectral Eigenproblems --- shared utility.

Rational Chebyshev basis for the infinite interval (-inf, +inf) via the
algebraic map
                    x = L * t / sqrt(1 - t^2),
with t on the Chebyshev-Gauss-Lobatto grid t_j = cos(pi j / N) in [-1, 1].
The endpoints t = +-1 map to x = +-inf, and the interior grid gives
(N - 1) collocation points on the real line.

The chain-rule factor is
                    dt/dx = (1 - t^2)^(3/2) / L,
so that if D_t is the (N+1) x (N+1) Chebyshev differentiation matrix in t,
the differentiation matrix in x on the interior nodes is
                    D_x = diag((1 - t^2)^(3/2) / L) * D_t
restricted to rows and columns 1:N.

For the second derivative we apply the chain rule again:
    d^2/dx^2 = (dt/dx)^2 * d^2/dt^2 + (d^2 t / d x^2) * d/dt,
where d^2 t / d x^2 = d/dx (dt/dx) = (dt/dx) * d/dt (dt/dx)
                    = -(3 t (1 - t^2)^2) / L^2.

The parameter L is the "map parameter": choosing L to match the
characteristic length of the target eigenfunction (e.g. L ~ few for the
harmonic oscillator) controls resolution. Too small an L crowds the
interior grid; too large an L puts most grid points where the function
has already decayed.

This utility is used by Etude 18.4 (infinite-interval tax / harmonic
oscillator) and Etude 18.6 (bound states plus continuum).

Author: Dr. Denys Dutykh (Khalifa University of Science and Technology)
Part of "Computational Etudes: A Spectral Approach"
https://github.com/dutykh/computational-etudes
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

_CH07_DIR = Path(__file__).resolve().parent.parent / "ch07"
if str(_CH07_DIR) not in sys.path:
    sys.path.insert(0, str(_CH07_DIR))
from cheb_matrix import cheb_matrix  # noqa: E402


def rational_chebyshev_grid(N: int, L: float = 4.0):
    """Return (x_interior, t_full, D_t_full) for the TB_n algebraic map.

    Parameters
    ----------
    N : int
        Number of intervals on the underlying Chebyshev-Lobatto t-grid
        (so there are N+1 t-nodes, of which 2 are the endpoints and
        the remaining N-1 are the interior nodes we use for collocation).
    L : float
        Map parameter; larger L stretches the interior grid further
        along the real line.

    Returns
    -------
    x_int : ndarray, shape (N-1,)
        Interior collocation points on the real line, from +x_max to
        -x_max (matching the Chebyshev ordering which descends from
        t = +1 to t = -1).
    t_full : ndarray, shape (N+1,)
        Full Chebyshev-Lobatto t-grid including the endpoints.
    D_t_full : ndarray, shape (N+1, N+1)
        Chebyshev differentiation matrix in t on the full grid.
    """
    D_t, t = cheb_matrix(N)
    # avoid divide-by-zero at endpoints (x = inf); we only use the interior
    interior = slice(1, N)
    t_int = t[interior]
    x_int = L * t_int / np.sqrt(1.0 - t_int ** 2)
    return x_int, t, D_t


def rational_chebyshev_derivative_matrices(N: int, L: float = 4.0):
    """Return (D1_x, D2_x, x_int) with first and second derivative matrices
    on the interior rational Chebyshev grid.

    The derivatives use the chain rule
        D1_x = (dt/dx) * D_t
        D2_x = (dt/dx)^2 * D_t^2 + (d^2 t / d x^2) * D_t

    where
        dt/dx        = (1 - t^2)^(3/2) / L
        d^2 t / d x^2 = -3 t (1 - t^2)^2 / L^2.

    Returns
    -------
    D1_x : ndarray, shape (N-1, N-1)
    D2_x : ndarray, shape (N-1, N-1)
    x_int : ndarray, shape (N-1,)
    """
    x_int, t, D_t = rational_chebyshev_grid(N, L)

    # chain-rule factors on the interior nodes
    t_int = t[1:N]
    one_minus_t2 = 1.0 - t_int ** 2
    dt_dx = (one_minus_t2 ** 1.5) / L
    d2t_dx2 = -3.0 * t_int * (one_minus_t2 ** 2) / (L ** 2)

    # restrict D_t to interior rows/cols; boundary rows of D_t contain
    # values at t = +-1 where x = +-inf; any decaying function vanishes
    # there so dropping those rows is consistent with the boundary
    # condition u(+-inf) = 0.
    D_t_int  = D_t[1:N, 1:N]
    D_t2_int = (D_t @ D_t)[1:N, 1:N]

    D1_x = np.diag(dt_dx) @ D_t_int
    D2_x = np.diag(dt_dx ** 2) @ D_t2_int + np.diag(d2t_dx2) @ D_t_int

    return D1_x, D2_x, x_int


if __name__ == "__main__":
    # sanity test: u(x) = exp(-x^2/2) has u_xx = (x^2 - 1) u
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    for N in (16, 32, 64):
        D1, D2, x = rational_chebyshev_derivative_matrices(N, L=4.0)
        u = np.exp(-x ** 2 / 2.0)
        u_xx_num = D2 @ u
        u_xx_ex  = (x ** 2 - 1.0) * u
        err = np.linalg.norm(u_xx_num - u_xx_ex, np.inf)
        print(f"N = {N:3d}  L = 4   max error in u_xx (Gaussian) = {err:.3e}")
