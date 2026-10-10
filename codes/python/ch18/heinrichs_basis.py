#!/usr/bin/env python3
"""
heinrichs_basis.py

Chapter 18: Linear Spectral Eigenproblems --- shared utility.

Heinrichs boundary-adapted bases for second- and fourth-order operators
on [-1, 1]. For a second-order Dirichlet problem u(+-1) = 0, the basis

    phi_j(x) = (1 - x^2) T_j(x),   j = 0, 1, ..., N-2

automatically satisfies the boundary conditions. For a fourth-order
clamped problem u(+-1) = u'(+-1) = 0,

    phi_j(x) = (1 - x^2)^2 T_j(x),   j = 0, 1, ..., N-4

satisfies all four conditions. The pedagogical payoff is conditioning:
the matrix discretising d^2/dx^2 in the Heinrichs basis has condition
number O(N^2) versus O(N^4) for the naive difference basis
T_{2j} - T_0, and similarly O(N^4) vs. O(N^8) for fourth derivatives.

We construct the Heinrichs operator matrix by collocation: evaluate the
basis and its derivatives at interior Chebyshev points, assemble the
differentiated basis, and solve for eigenvalues by standard
LAPACK routines.

Author: Dr. Denys Dutykh (Khalifa University of Science and Technology)
Part of "Computational Etudes: A Spectral Approach"
https://github.com/dutykh/computational-etudes
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

_CH07 = Path(__file__).resolve().parent.parent / "ch07"
if str(_CH07) not in sys.path:
    sys.path.insert(0, str(_CH07))
from cheb_matrix import cheb_matrix  # noqa: E402


def heinrichs_dirichlet_matrix(N: int):
    """Return (A, x_int) for the second-order Dirichlet eigenproblem
            -u'' = lambda u,   u(+-1) = 0
    in the Heinrichs basis phi_j(x) = (1 - x^2) T_j(x). The returned
    matrix A has the property that eigenvalues of A approximate the
    operator spectrum (j pi/2)^2.

    The implementation builds the matrix via Chebyshev collocation on
    the interior grid, using the factorisation u = (1 - x^2) q with q
    an arbitrary polynomial. By the Leibniz rule,

        u'' = (1 - x^2) q'' - 4 x q' - 2 q,

    so on interior Chebyshev nodes
        -u'' = lambda (1 - x^2) q
    becomes the GENERALISED eigenproblem
        (-S D^2 + 4 X D + 2 I) q = lambda S q,
    where S = diag(1 - x_j^2) and X = diag(x_j).

    Returning the equivalent standard problem A = S^{-1} (-S D^2 + 4 X D + 2 I)
    allows a direct comparison of condition number with the
    "naive" restricted second-derivative approach of Etude 18.3.
    """
    D, x = cheb_matrix(N)
    D2 = D @ D
    idx = slice(1, N)
    Di = D[idx, idx]
    D2i = D2[idx, idx]
    xi = x[idx]
    S = np.diag(1.0 - xi ** 2)
    X = np.diag(xi)
    # left matrix L = -S D^2 + 4 X D + 2 I (this is the operator for -u'')
    L = -S @ D2i + 4.0 * X @ Di + 2.0 * np.eye(N - 1)
    # convert to standard form A = S^{-1} L so A q = lambda q
    A = np.linalg.solve(S, L)
    return A, xi


def heinrichs_clamped_matrix(N: int):
    """Return (A, M, x_int) for the fourth-order clamped eigenproblem
            u'''' = lambda u,  u(+-1) = u'(+-1) = 0,
    in the Heinrichs basis phi_j(x) = (1 - x^2)^2 T_j(x). The problem
    is posed as a GENERALISED eigenproblem A q = lambda M q.

    With u = (1 - x^2)^2 q, the Leibniz rule gives
        u^(4) = (1 - x^2)^2 q^(4) - 16 x (1 - x^2) q''' +
                (48 x^2 - 24) q'' + 96 x q' + 24 q.
    Collocating at interior Chebyshev nodes yields
        (S2 D^4 - 16 X S D^3 + (48 X^2 - 24 I) D^2 + 96 X D + 24 I) q
           = lambda S2 q,
    where S = diag(1 - x^2), S2 = diag((1 - x^2)^2), X = diag(x).
    """
    D, x = cheb_matrix(N)
    D2 = D @ D
    D3 = D2 @ D
    D4 = D3 @ D
    idx = slice(1, N)
    Di = D[idx, idx]
    D2i = D2[idx, idx]
    D3i = D3[idx, idx]
    D4i = D4[idx, idx]
    xi = x[idx]
    I_int = np.eye(N - 1)
    S  = np.diag(1.0 - xi ** 2)
    S2 = np.diag((1.0 - xi ** 2) ** 2)
    X  = np.diag(xi)
    X2 = np.diag(xi ** 2)
    A = (S2 @ D4i
         - 16.0 * X @ S @ D3i
         + (48.0 * X2 - 24.0 * I_int) @ D2i
         + 96.0 * X @ Di
         + 24.0 * I_int)
    M = S2
    return A, M, xi


def naive_dirichlet_matrix(N: int):
    """-D^2 restricted to interior nodes. Same as Etude 18.1 / 18.3."""
    D, x = cheb_matrix(N)
    D2 = D @ D
    return -D2[1:N, 1:N], x[1:N]


def naive_clamped_operator(N: int):
    """Return A for u^(4) = lambda u with clamped BCs via boundary
    bordering: D^4 with rows 0, 1, N-1, N replaced by identity/derivative
    conditions. Returns (A, B) pencil."""
    D, x = cheb_matrix(N)
    D4 = (D @ D) @ (D @ D)
    A = D4.copy()
    B = np.eye(N + 1)
    # u(1) = 0 on row 0, u'(1) = 0 on row 1,
    # u'(-1) = 0 on row N-1, u(-1) = 0 on row N
    A[0, :] = 0.0;   A[0, 0] = 1.0;      B[0, :] = 0.0
    A[1, :] = D[0, :];                   B[1, :] = 0.0
    A[N - 1, :] = D[N, :];               B[N - 1, :] = 0.0
    A[N, :] = 0.0;   A[N, N] = 1.0;      B[N, :] = 0.0
    return A, B


if __name__ == "__main__":
    # Sanity test: compare Heinrichs and naive Dirichlet spectra
    print("Condition-number comparison for second-order operator:")
    print(f"{'N':>6} {'kappa(naive)':>16} {'kappa(Heinrichs)':>20}")
    for N in (16, 32, 64, 96):
        A_naive, _ = naive_dirichlet_matrix(N)
        A_hein,  _ = heinrichs_dirichlet_matrix(N)
        k_naive = np.linalg.cond(A_naive)
        k_hein  = np.linalg.cond(A_hein)
        print(f"{N:>6d} {k_naive:>16.3e} {k_hein:>20.3e}")
