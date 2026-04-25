#!/usr/bin/env python3
"""
ioakimidis_root.py

Chapter 21: Special Tricks for the Spectral Researcher
Computational Etude 21.7: Ioakimidis' non-iterative root via Chebyshev quadrature.

Suppose f(x) has exactly one simple root rho on [a, b].  Boyd notes
(after Ioakimidis, in unpublished work cited in Boyd 2000) that one
can locate rho non-iteratively by writing the Chebyshev-quadrature
identity

   I(N) = sum_{j=0}^{2N} '' (-1)^j x_j / f(x_j)
   J(N) = sum_{j=0}^{2N} '' (-1)^j     / f(x_j)

where x_j is the (2N+1)-point Chebyshev-Lobatto grid on [a, b] and
the double prime denotes the Clenshaw-Curtis half-weighting at the
endpoints.  Then

   rho_N = I(N) / J(N)

converges to rho geometrically as N grows, with rate determined by
the Bernstein ellipse of f restricted to the relevant interval.

The remarkable property is that NO iteration is needed: the formula
is a direct quotient of two finite sums.  As an alternative, simple
bisection on the same number of evaluations decreases the error only
linearly.

We test on Boyd's example f(x) = sin(x - pi/4) / sqrt(1 + 10 x^2)
on [-1, 1], whose unique root is rho = pi/4 = 0.7853981633974...

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
from tricks_common import setup_matplotlib, output_dir_for, save_fig, dump_json
from tricks_common import NAVY, CORAL, TEAL


def f_test(x):
    return np.sin(x - np.pi / 4.0) / np.sqrt(1.0 + 10.0 * x ** 2)


RHO_EXACT = np.pi / 4.0


def ioakimidis_root(N, a, b):
    """Boyd's Eq (19.36)-(19.37): non-iterative root of f on [a, b]."""
    j = np.arange(2 * N + 1)
    half = 0.5 * ((a + b) - (a - b) * np.cos(j * np.pi / (2.0 * N)))  # Boyd's x_j
    fj = f_test(half)
    sign = (-1.0) ** j
    weight = np.ones_like(j, dtype=float)
    weight[0] = 0.5
    weight[-1] = 0.5
    num = np.sum(weight * sign * half / fj)
    den = np.sum(weight * sign      / fj)
    return num / den


def bisection_root(n_evals, a, b, tol=0.0):
    """Plain bisection using exactly n_evals function calls."""
    fa = f_test(a)
    fb = f_test(b)
    assert fa * fb < 0
    history = []
    for _ in range(n_evals - 2):
        c = 0.5 * (a + b)
        fc = f_test(c)
        history.append(c)
        if fc * fa < 0:
            b = c; fb = fc
        else:
            a = c; fa = fc
    return 0.5 * (a + b)


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    p.add_argument("--dump", type=str, default=None)
    args = p.parse_args()

    setup_matplotlib()
    out = output_dir_for(SCRIPT_DIR)

    a, b = -1.0, 1.0
    N_axis = np.arange(2, 31)
    err_io = np.empty_like(N_axis, dtype=float)
    err_bi = np.empty_like(N_axis, dtype=float)
    for k, N in enumerate(N_axis):
        rho_N = ioakimidis_root(N, a, b)
        err_io[k] = abs(rho_N - RHO_EXACT)
        rho_bi = bisection_root(2 * N + 1, a, b)
        err_bi[k] = abs(rho_bi - RHO_EXACT)

    # ---- Figure: error vs evaluations -------------------------------
    fig = plt.figure(figsize=(7.5, 4.4))
    ax = fig.add_subplot(111)
    n_evals = 2 * N_axis + 1
    ax.semilogy(n_evals, err_io + 1e-18, "o-", color=NAVY, markerfacecolor="white",
                markersize=6, linewidth=1.0, label="Ioakimidis (non-iterative)")
    ax.semilogy(n_evals, err_bi + 1e-18, "s--", color=CORAL, markerfacecolor="white",
                markersize=5, linewidth=0.9, label="bisection (same # of evals)")
    ax.axhline(1e-15, color="gray", linewidth=0.4, alpha=0.5)
    ax.set_xlabel(r"number of $f$ evaluations  ($= 2N+1$)")
    ax.set_ylabel(r"$|\rho_N - \pi/4|$")
    ax.set_title(r"Ioakimidis: geometric convergence, no Newton, "
                 r"on $f(x) = \sin(x-\pi/4)/\sqrt{1+10x^2}$",
                 fontsize=10)
    ax.legend(loc="lower left", fontsize=10)
    ax.grid(True, alpha=0.25, linewidth=0.4)
    ax.set_ylim(1e-16, 1)

    save_fig(fig, out, "ioakimidis_root")
    plt.close(fig)

    print(f"[Etude 21.7]  Ioakimidis non-iterative root")
    print(f"  exact rho = pi/4 = {RHO_EXACT:.16f}")
    for kk in (3, 6, 10, 20, 30):
        rho = ioakimidis_root(kk, a, b)
        print(f"  N = {kk:2d} ({2*kk+1:2d} evals): rho = {rho:.16f},  err = {abs(rho - RHO_EXACT):.3e}")
    print(f"  figure: {out / 'ioakimidis_root.pdf'}")

    dump_json(args.dump, {
        "rho_exact": float(RHO_EXACT),
        "N_axis": N_axis.tolist(),
        "err_ioakimidis": err_io.tolist(),
        "err_bisection": err_bi.tolist(),
    })


if __name__ == "__main__":
    main()
