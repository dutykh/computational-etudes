#!/usr/bin/env python3
"""
behavioral_bc_laguerre.py

Chapter 20: Spectral Methods on Unbounded Intervals.
Computational Etude 20.7: Boundedness is sometimes enough -- and
sometimes not.

Boyd (1987b) poses the Laguerre-type eigenvalue problem on [0, infty):

    y u'' + (y + 1) u' + lambda u = 0,      u(0) & u(infty)  bounded.

The exact spectrum is  lambda_n = n,  n = 0, 1, 2, ...  and the
eigenfunctions are  u_n(y) = exp(-y) L_n^1(y), with
L_n^1  the first-order associated Laguerre polynomial.

STRATEGY A (naive).  Discretise u on the TL_n grid, impose NO
boundary constraint.  One accepts the unnecessary assumption that
the second linearly-independent solution (which blows up like a
polynomial power of y, not exponentially) will be numerically
filtered out.  Boyd reports that this works, but painfully: for
N = 40 collocation points, the lowest eigenvalue agrees with
lambda_0 = 0 to only a couple of decimal places.

STRATEGY B (behavioural recast).  Change of unknown w = exp(y/2) u(y).
The ODE becomes

    y w'' + w' + [-1/2 - y/4 + lambda] w = 0,

whose blow-up solutions now blow up EXPONENTIALLY, and the
rational-Chebyshev matrix sees a far better-conditioned problem.

The etude follows Boyd's Fig. 17.9 spirit: plot, as a function of N,
the number of eigenvalues that match the exact spectrum to within a
relative tolerance 0.05.  Strategy A has to climb; Strategy B starts
climbing immediately.

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

from unbounded_common import (CORAL, NAVY, TEAL, output_dir_for, save_fig,
                              setup_matplotlib, tl_map_forward, tl_map_fprime,
                              tl_map_fdoubleprime)

OUTPUT_DIR = output_dir_for(SCRIPT_DIR)


def tln_differentiation_matrices(N, L):
    """Return the mapped first and second derivative matrices acting on
    samples of a function at the TL_n grid  y_j = L (1 + x_j) / (1 - x_j),
    x_j Chebyshev-Gauss-Lobatto.  We suppress the right endpoint (y = infty)
    by returning N x N interior matrices.
    """
    Dx, x = cheb_matrix(N)
    with np.errstate(divide="ignore", invalid="ignore"):
        fp = tl_map_fprime(x, L)
        fpp = tl_map_fdoubleprime(x, L)
        Dy = np.diag(1.0 / fp) @ Dx
        Dy2 = np.diag(1.0 / fp ** 2) @ (Dx @ Dx) - np.diag(fpp / fp ** 3) @ Dx
    y = tl_map_forward(x, L)
    # interior: drop y = infty (index 0) but keep y = 0 (index N)
    # reorder to have increasing y: x descending -> y ascending
    # we keep indices 1..N (drop x = 1 = y = infty)
    return y[1:], Dy[1:, 1:], Dy2[1:, 1:]


def solve_strategy_A(N, L):
    """Directly discretise  y u'' + (y + 1) u' + lambda u = 0  on the TL_n
    interior grid; no boundary condition imposed (u(0) and u(infty) are
    assumed to be treated as behavioural).
    """
    y, Dy, Dy2 = tln_differentiation_matrices(N, L)
    Y = np.diag(y)
    A = Y @ Dy2 + (Y + np.eye(len(y))) @ Dy
    eigs = np.linalg.eigvals(-A)  # rearrange: lambda u = - [y D2 + (y+1) D1] u
    return np.sort(eigs.real)


def solve_strategy_B(N, L):
    """Change of unknown  w = exp(y/2) u(y), yielding the rescaled problem
        y w'' + w' + [-1/2 - y/4 + lambda] w = 0.
    We remove the quartic blow-up of the unphysical solutions, and the
    eigenvalues are the same.
    """
    y, Dy, Dy2 = tln_differentiation_matrices(N, L)
    Y = np.diag(y)
    M = np.diag(0.5 + 0.25 * y)
    A = Y @ Dy2 + Dy - M
    eigs = np.linalg.eigvals(-A)
    return np.sort(eigs.real)


def count_good(eigs, tol=0.05):
    """Count eigenvalues that lie within relative tolerance tol of some
    integer n = 0, 1, 2, ..., M.  Each matrix eigenvalue is matched to
    its nearest integer target; integers may be claimed by at most one
    eigenvalue (`greedy closest` matching).
    """
    targets = np.arange(50)
    taken = set()
    good = 0
    for lam in eigs:
        if lam < -0.5:   # well below the positive axis
            continue
        dists = np.abs(lam - targets)
        order = np.argsort(dists)
        for n in order:
            if n in taken:
                continue
            rel_err = abs(lam - n) / max(1.0, n)
            if rel_err < tol:
                taken.add(n)
                good += 1
            break
    return good


def make_figure():
    setup_matplotlib()

    L = 32.0
    Ns = np.array([10, 20, 30, 40, 60, 80, 120])
    good_A = [count_good(solve_strategy_A(N, L)) for N in Ns]
    good_B = [count_good(solve_strategy_B(N, L)) for N in Ns]

    fig, axes = plt.subplots(1, 2, figsize=(10.0, 3.8))

    ax = axes[0]
    ax.plot(Ns, good_A, "-o", color=CORAL, lw=1.2, label="(A) naive")
    ax.plot(Ns, good_B, "-s", color=TEAL, lw=1.2, mfc="none",
            label=r"(B) behavioural, $w = e^{y/2} u$")
    ax.set_xlabel(r"$N$")
    ax.set_ylabel("number of good eigenvalues")
    ax.set_title("(a) How many matrix eigenvalues track $\\lambda_n = n$?")
    ax.grid(True, alpha=0.3)
    ax.legend(frameon=False, fontsize=9)

    ax = axes[1]
    eigs_A = solve_strategy_A(40, L)[:20]
    eigs_B = solve_strategy_B(40, L)[:20]
    ax.plot(np.real(eigs_A), np.imag(eigs_A), "o", color=CORAL,
            label="Strategy A, $N = 40$")
    ax.plot(np.real(eigs_B), np.imag(eigs_B), "s", color=TEAL,
            mfc="none", label="Strategy B, $N = 40$")
    ax.axhline(0, color="gray", lw=0.5)
    ax.set_xlim(-1, 20)
    ax.set_ylim(-3, 3)
    ax.set_xlabel(r"$\mathrm{Re}(\lambda)$")
    ax.set_ylabel(r"$\mathrm{Im}(\lambda)$")
    ax.set_title(r"(b) Spectrum: A is contaminated; B is clean")
    ax.grid(True, alpha=0.3)
    ax.legend(frameon=False, fontsize=9)

    fig.tight_layout()
    save_fig(fig, OUTPUT_DIR, "behavioral_bc_laguerre")
    plt.close(fig)

    print(f"[20.7] saved figure to {OUTPUT_DIR / 'behavioral_bc_laguerre.pdf'}")
    for N, gA, gB in zip(Ns, good_A, good_B):
        print(f"  N={N:3d}  good(A)={gA}  good(B)={gB}")


if __name__ == "__main__":
    make_figure()
