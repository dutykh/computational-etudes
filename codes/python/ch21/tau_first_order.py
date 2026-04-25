#!/usr/bin/env python3
"""
tau_first_order.py

Chapter 21: Special Tricks for the Spectral Researcher
Computational Etude 21.11: tau method on a first-order ODE.

Lanczos's philosophy: instead of approximating the solution of the
exact problem, find the *exact* solution of a *nearby* problem.  The
nearby problem differs from the original only by a single, small
perturbation chosen in a spectrally good basis (Chebyshev T_N).

Original problem (Boyd Eq 21.13):
        u'(x) + u(x) = 0,   u(-1) = 1,
with exact solution  u(x) = exp(-(x + 1)).

A polynomial cannot satisfy this equation identically, but the
*modified* equation
        v'(x) + v(x) = tau T_N(x),   v(-1) = 1,
admits a unique polynomial solution of degree N for each N.  Compared
to standard Galerkin/collocation, the tau method is essentially
'boundary bordering with a Chebyshev test function on the residual'.

We compare two things:
    (1) tau scheme: solve the modified equation; compare v_N against u
    (2) standard pseudospectral: collocate u' + u = 0 at the interior
        Chebyshev points and impose u(-1) = 1.
The two are within a few decimals of each other.  More importantly,
the perturbation amplitude tau decays geometrically with N -- the
sense in which Lanczos's 'small modification' is small.

Author: Dr. Denys Dutykh (Khalifa University of Science and Technology)
Part of "Computational Etudes: A Spectral Approach"
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
sys.path.insert(0, str(SCRIPT_DIR.parent / "ch07"))
from tricks_common import setup_matplotlib, output_dir_for, save_fig, dump_json
from tricks_common import NAVY, CORAL, TEAL
from cheb_matrix import cheb_matrix


def tau_solve(N):
    """Solve  v'(x) + v(x) = tau T_N(x), v(-1) = 1, on the (N+1)-pt
    Chebyshev grid.  Returns the grid x, the solution vector v, and tau.

    The system has N+2 unknowns (v_0, ..., v_N, tau) and is closed by:
       row 1: v(-1) = 1
       rows 2..N+1: v'(x_i) + v(x_i) - tau T_N(x_i) = 0  at interior nodes
    """
    D, x = cheb_matrix(N)
    # Build augmented system [A | -T_N(x)] [v; tau] = b
    L = D + np.eye(N + 1)               # L v = right-hand side (= 0 in the
                                        # original problem)
    TN = np.cos(N * np.arccos(x))       # T_N at the CGL points
    A = np.zeros((N + 2, N + 2))
    b = np.zeros(N + 2)
    # Equation block: residual = 0 at all N+1 grid points
    A[0:N+1, 0:N+1] = L
    A[0:N+1, N+1] = -TN
    # Replace one row with the boundary condition v(-1) = 1
    # The CGL grid is in descending order: x[0] = +1, x[N] = -1.
    A[N+1, :] = 0.0
    A[N+1, N] = 1.0  # v_N corresponds to x = -1
    b[N+1] = 1.0
    sol = np.linalg.solve(A, b)
    return x, sol[0:N+1], sol[N+1]


def standard_pseudospectral(N):
    """Solve u' + u = 0 by Chebyshev pseudospectral collocation at
    interior nodes, with u(-1) = 1 imposed exactly."""
    D, x = cheb_matrix(N)
    L = D + np.eye(N + 1)
    # Replace the row corresponding to x = -1 (the last row) with BC.
    L[N, :] = 0.0
    L[N, N] = 1.0
    rhs = np.zeros(N + 1)
    rhs[N] = 1.0
    return x, np.linalg.solve(L, rhs)


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    p.add_argument("--dump", type=str, default=None)
    args = p.parse_args()

    setup_matplotlib()
    out = output_dir_for(SCRIPT_DIR)

    Ns = np.array([4, 6, 8, 12, 16, 24, 32])
    tau_vals  = np.empty_like(Ns, dtype=float)
    err_tau   = np.empty_like(Ns, dtype=float)
    err_pseu  = np.empty_like(Ns, dtype=float)
    for k, N in enumerate(Ns):
        x, v, tau = tau_solve(N)
        tau_vals[k] = abs(tau)
        u_exact = np.exp(-(x + 1.0))
        err_tau[k] = np.max(np.abs(v - u_exact))
        x2, u_ps = standard_pseudospectral(N)
        err_pseu[k] = np.max(np.abs(u_ps - np.exp(-(x2 + 1.0))))

    # ---- Visual at N = 16 ------------------------------------------
    N_show = 16
    x_show, v_show, tau_show = tau_solve(N_show)
    x_dense = np.linspace(-1, 1, 401)
    u_exact_dense = np.exp(-(x_dense + 1.0))

    # ---- Figure -----------------------------------------------------
    fig = plt.figure(figsize=(11.5, 4.4))
    gs = fig.add_gridspec(1, 2, wspace=0.28)

    # Panel A: tau and error vs N
    ax = fig.add_subplot(gs[0, 0])
    ax.semilogy(Ns, tau_vals, "o-", color=NAVY, markerfacecolor="white",
                markersize=6, linewidth=1.0,
                label=r"$|\tau|$ (size of the Lanczos perturbation)")
    ax.semilogy(Ns, err_tau, "s--", color=TEAL, markerfacecolor="white",
                markersize=5, linewidth=0.9,
                label=r"max $|v_N - u_\mathrm{exact}|$")
    ax.semilogy(Ns, err_pseu, "^:", color=CORAL, markerfacecolor="white",
                markersize=5, linewidth=0.9,
                label=r"max $|u_N^{\mathrm{p\text{-}s}} - u_\mathrm{exact}|$")
    ax.set_xlabel(r"$N$")
    ax.set_ylabel(r"$|\tau|$ and pointwise error")
    ax.set_title("geometric decay of $\\tau$, "
                 "geometric convergence of both methods", fontsize=10)
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(True, alpha=0.25, linewidth=0.4)
    ax.set_ylim(1e-16, 5)

    # Panel B: solution and error at N=16
    ax = fig.add_subplot(gs[0, 1])
    ax.plot(x_dense, u_exact_dense, color=NAVY, linewidth=1.4,
            label=r"$u_\mathrm{exact}(x) = e^{-(x+1)}$")
    ax.plot(x_show, v_show, "o", color=TEAL, markerfacecolor="white",
            markersize=6, label=fr"$v_{{{N_show}}}(x_j)$ tau solution")
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$u(x)$")
    ax.set_title(rf"$N = {N_show}$:  $\tau \approx {tau_show:+.3e}$, "
                 rf"max err $\approx${err_tau[np.where(Ns == N_show)[0][0]]:.1e}",
                 fontsize=10)
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(True, alpha=0.25, linewidth=0.4)

    save_fig(fig, out, "tau_first_order")
    plt.close(fig)

    print(f"[Etude 21.11]  tau method on u' + u = 0,  u(-1) = 1")
    for k, N in enumerate(Ns):
        print(f"  N = {N:2d}:  |tau| = {tau_vals[k]:.3e}, "
              f"err_tau = {err_tau[k]:.3e}, err_p-s = {err_pseu[k]:.3e}")
    print(f"  figure: {out / 'tau_first_order.pdf'}")

    dump_json(args.dump, {
        "Ns": Ns.tolist(),
        "tau": tau_vals.tolist(),
        "err_tau": err_tau.tolist(),
        "err_pseudospectral": err_pseu.tolist(),
    })


if __name__ == "__main__":
    main()
