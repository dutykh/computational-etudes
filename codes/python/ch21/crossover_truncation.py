#!/usr/bin/env python3
"""
crossover_truncation.py

Chapter 21: Special Tricks for the Spectral Researcher
Computational Etude 21.4: cross-over truncation.

A pure asymptotic statement is not always a numerical statement.
Boyd's cartoon (his Eq 19.13) makes the point in a single line:
        a_n ~ 10 exp(-n / 3) + 1e-6 / n^5

Asymptotically the algebraic 1/n^5 wins, but it does not win until
n is past the cross-over n_x ~ 120, beyond which most computations
will not even reach.  The correct statement for the practitioner is
not 'a_n decays algebraically' but 'a_n decays exponentially up to
n ~ n_x, then algebraically'.  The two slopes are visible to the eye.

We then back the cartoon with a real example:  consider a function
        f(x) = exp(-(x - 0.3)^2 / 0.01) + 1e-7 (x + 1)^{1/3}
        on [-1, 1].
The first term is entire and contributes geometric coefficient
decay; the second has a weak (algebraic) endpoint singularity at
x = -1 contributing a 1/n^p tail.  At small N the geometric head
dominates; at large N the algebraic tail eventually emerges.

This second panel illustrates that 'cross-over truncation' is a
phenomenon of the real world, not just of textbook cartoons.

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
from tricks_common import NAVY, CORAL, TEAL, PURPLE, ORANGE


def cgl(N):
    return np.cos(np.pi * np.arange(N + 1) / N)


def cheb_coeffs(v):
    N = len(v) - 1
    extended = np.concatenate([v, v[N - 1:0:-1]])
    A = np.real(np.fft.fft(extended)) / N
    A[0] *= 0.5
    A[N] *= 0.5
    return A[:N + 1]


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    p.add_argument("--dump", type=str, default=None)
    args = p.parse_args()

    setup_matplotlib()
    out = output_dir_for(SCRIPT_DIR)

    # ---- Panel A: Boyd's cartoon -------------------------------------
    n_axis = np.arange(1, 401)
    head_geom = 10.0 * np.exp(-n_axis / 3.0)
    tail_alg  = 1.0e-6 / n_axis ** 5
    a_total = head_geom + tail_alg
    # Crossover: smallest n where the algebraic tail exceeds the
    # exponential head.  Compare in log-space to avoid floating-point
    # underflow when both values are tiny.
    log_diff = np.log10(tail_alg + 1e-300) - np.log10(head_geom + 1e-300)
    n_cross = float(n_axis[np.argmin(np.abs(log_diff))])

    # ---- Panel B: real Chebyshev coefficients of a 'two-slope' function -
    # f(x) = exp(-((x - 0.3) / 0.1)^2) + 1e-7 (x + 1)^{1/3}
    def f(x):
        return np.exp(-((x - 0.3) / 0.1) ** 2) + 1e-7 * (x + 1.0) ** (1.0 / 3.0)

    N = 400
    xg = cgl(N)
    a = cheb_coeffs(f(xg))

    # ---- Figure ------------------------------------------------------
    fig = plt.figure(figsize=(11.5, 4.4))
    gs = fig.add_gridspec(1, 2, wspace=0.30)

    # Panel A
    ax = fig.add_subplot(gs[0, 0])
    ax.semilogy(n_axis, head_geom, color=NAVY, linewidth=1.2,
                label=r"$10\,e^{-n/3}$ (geometric head)")
    ax.semilogy(n_axis, tail_alg, color=CORAL, linewidth=1.2,
                label=r"$10^{-6}/n^{5}$ (algebraic tail)")
    ax.semilogy(n_axis, a_total, color=PURPLE, linewidth=1.6,
                label=r"$a_n = $ head $+$ tail")
    ax.axvline(n_cross, color="gray", linestyle="--", linewidth=0.8, alpha=0.6,
               label=rf"$n_\times \approx {int(n_cross)}$")
    ax.set_xlabel(r"index $n$")
    ax.set_ylabel(r"$|a_n|$")
    ax.set_title("Panel A.  Boyd's cartoon: a slope change at $n_\\times$",
                 fontsize=10)
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(True, alpha=0.25, linewidth=0.4)
    ax.set_ylim(1e-22, 50)

    # Panel B
    ax = fig.add_subplot(gs[0, 1])
    ax.semilogy(np.arange(N + 1), np.abs(a) + 1e-300,
                color=NAVY, linewidth=1.0, marker=".",
                markersize=2, markerfacecolor=NAVY,
                label=r"computed $|a_n|$ for $f(x)$")
    # Reference slopes
    n_ref = np.arange(1, N + 1)
    ax.semilogy(n_ref, 0.5 * np.exp(-n_ref / 5.5), color=TEAL, linewidth=0.8,
                linestyle="--", label=r"$\sim e^{-n/5.5}$ (smooth bulk)")
    ax.semilogy(n_ref, 1.5e-6 / n_ref ** (4.0 / 3.0), color=CORAL,
                linewidth=0.8, linestyle="--",
                label=r"$\sim n^{-4/3}$ (endpoint $x = -1$)")
    ax.set_xlabel(r"Chebyshev degree $n$")
    ax.set_ylabel(r"$|a_n|$")
    ax.set_title(r"Panel B.  $f(x) = e^{-((x-0.3)/0.1)^2} + 10^{-7}(x+1)^{1/3}$",
                 fontsize=10)
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(True, alpha=0.25, linewidth=0.4)
    ax.set_ylim(1e-18, 5)

    save_fig(fig, out, "crossover_truncation")
    plt.close(fig)

    print(f"[Etude 21.4]  cross-over truncation")
    print(f"  cartoon n_cross approx {int(n_cross)}")
    print(f"  real example: |a_50|  = {np.abs(a[50]):.3e}")
    print(f"  real example: |a_200| = {np.abs(a[200]):.3e}")
    print(f"  real example: |a_400| = {np.abs(a[400]):.3e}")
    print(f"  figure: {out / 'crossover_truncation.pdf'}")

    dump_json(args.dump, {
        "n_cross_cartoon": float(n_cross),
        "n_axis_cartoon": n_axis.tolist(),
        "a_cartoon": a_total.tolist(),
        "N_real": N,
        "a_real": np.abs(a).tolist(),
    })


if __name__ == "__main__":
    main()
