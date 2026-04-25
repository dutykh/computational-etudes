#!/usr/bin/env python3
"""
singularity_subtraction.py

Chapter 21: Special Tricks for the Spectral Researcher
Computational Etude 21.3: singularity subtraction.

When the solution of a boundary-value problem has known asymptotic
singular behaviour at a corner or endpoint (such as the r^(2/3)
streamfunction singularity at a re-entrant corner of a 2-D Stokes
flow, Boyd Sec 19.3), the standard repair is to write
        u  =  u_singular  +  u_smooth,
expand u_singular ANALYTICALLY in known special functions, and let
the spectral solver approximate only u_smooth.  This decouples the
algebraic decay of the coefficients (driven by u_singular) from the
geometric decay of the actual unknown coefficients (of u_smooth).

To demonstrate the mechanism cleanly without solving a 2-D corner
problem, we use a one-dimensional surrogate.  Consider

        f(x)  =  e^x  +  c (1 + x)^{2/3},     x in [-1, 1],

with c = 0.1.  The first term is entire; the second has a 'corner'-
type singularity at x = -1 (continuous but not smooth -- its first
derivative is unbounded).  Expanded directly in Chebyshev polynomials,
f has coefficients that decay only algebraically (driven by the
singular term), and a max pointwise error that decays correspondingly
slowly.

The trick.  Write
        f(x) - c (1 + x)^{2/3}  =  e^x,
expand the RIGHT-HAND SIDE in Chebyshev polynomials (geometric decay,
machine epsilon at modest N), then re-add the analytic singular part
when reconstructing.  The same N now buys many more decimal places.

This is a one-line example of singularity subtraction's structural
moral: handle the bad part by hand, hand the good part to the
spectral method.

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


CSING = 0.1


def f_full(x):
    return np.exp(x) + CSING * (1.0 + x) ** (2.0 / 3.0)


def f_singular(x):
    return CSING * (1.0 + x) ** (2.0 / 3.0)


def f_smooth(x):
    return np.exp(x)


def cgl(N):
    return np.cos(np.pi * np.arange(N + 1) / N)


def cheb_coeffs(v):
    N = len(v) - 1
    extended = np.concatenate([v, v[N - 1:0:-1]])
    A = np.real(np.fft.fft(extended)) / N
    A[0] *= 0.5
    A[N] *= 0.5
    return A[:N + 1]


def reconstruct(a, t_eval):
    T0 = np.ones_like(t_eval)
    if len(a) == 1:
        return a[0] * T0
    T1 = t_eval.copy()
    val = a[0] * T0 + a[1] * T1
    for n in range(2, len(a)):
        Tk = 2.0 * t_eval * T1 - T0
        val += a[n] * Tk
        T0, T1 = T1, Tk
    return val


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    p.add_argument("--dump", type=str, default=None)
    args = p.parse_args()

    setup_matplotlib()
    out = output_dir_for(SCRIPT_DIR)

    # --- Coefficient panel at fixed N ----------------------------------
    N_show = 80
    x_show = cgl(N_show)
    a_naive = cheb_coeffs(f_full(x_show))
    a_trick = cheb_coeffs(f_smooth(x_show))   # only the smooth part is fit

    # --- Convergence sweep --------------------------------------------
    Ns = np.array([8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256])
    err_naive = np.empty_like(Ns, dtype=float)
    err_trick = np.empty_like(Ns, dtype=float)
    x_eval = np.linspace(-1.0, 1.0, 5001)
    for k, N in enumerate(Ns):
        xg = cgl(N)
        a_n = cheb_coeffs(f_full(xg))
        a_t = cheb_coeffs(f_smooth(xg))
        approx_n = reconstruct(a_n, x_eval)
        approx_t = reconstruct(a_t, x_eval) + f_singular(x_eval)
        err_naive[k] = np.max(np.abs(approx_n - f_full(x_eval)))
        err_trick[k] = np.max(np.abs(approx_t - f_full(x_eval)))

    # ---- Figure ------------------------------------------------------
    fig = plt.figure(figsize=(11.5, 4.4))
    gs = fig.add_gridspec(1, 2, wspace=0.30)

    # Panel A: coefficient decay at N_show
    ax = fig.add_subplot(gs[0, 0])
    n_axis = np.arange(N_show + 1)
    ax.semilogy(n_axis, np.abs(a_naive) + 1e-20, "o-", color=CORAL,
                markerfacecolor="white", markersize=4, linewidth=0.8,
                label=r"naive: $|a_n|$ for $f(x)$")
    ax.semilogy(n_axis, np.abs(a_trick) + 1e-20, "^-", color=TEAL,
                markerfacecolor="white", markersize=4, linewidth=0.8,
                label=r"trick: $|a_n|$ for $f(x) - c (1+x)^{2/3}$")
    ax.set_xlabel(r"Chebyshev degree $n$")
    ax.set_ylabel(r"$|a_n|$")
    ax.set_title(rf"coefficient decay at $N = {N_show}$",
                 fontsize=10)
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(True, alpha=0.25, linewidth=0.4)
    ax.set_ylim(1e-18, 5)

    # Panel B: convergence in N
    ax = fig.add_subplot(gs[0, 1])
    ax.loglog(Ns, err_naive, "o-", color=CORAL, markerfacecolor="white",
              markersize=6, linewidth=1.0,
              label=r"naive: max $|f - f_N|$")
    ax.loglog(Ns, err_trick, "^-", color=TEAL, markerfacecolor="white",
              markersize=6, linewidth=1.0,
              label=r"trick: max $|f - (f - c\,(1+x)^{2/3})_N - c\,(1+x)^{2/3}|$")
    ax.axhline(1e-15, color="gray", linewidth=0.4, alpha=0.5)
    ax.set_xlabel(r"$N$")
    ax.set_ylabel("max pointwise error")
    ax.set_title("subtraction reaches machine $\\epsilon$ at modest $N$",
                 fontsize=10)
    ax.legend(loc="lower left", fontsize=9)
    ax.grid(True, alpha=0.25, linewidth=0.4)
    ax.set_ylim(1e-17, 1)

    save_fig(fig, out, "singularity_subtraction")
    plt.close(fig)

    print(f"[Etude 21.3]  singularity subtraction (1-D surrogate)")
    print(f"  test function: e^x + {CSING} (1+x)^(2/3) on [-1, 1]")
    print(f"  N      naive max err     trick max err")
    for N, en, et in zip(Ns, err_naive, err_trick):
        print(f"  {N:4d}   {en:.3e}      {et:.3e}")
    print(f"  figure: {out / 'singularity_subtraction.pdf'}")

    dump_json(args.dump, {
        "N_show": N_show, "csing": CSING,
        "a_naive": list(map(float, a_naive)),
        "a_trick": list(map(float, a_trick)),
        "Ns": Ns.tolist(),
        "err_naive": err_naive.tolist(),
        "err_trick": err_trick.tolist(),
    })


if __name__ == "__main__":
    main()
