#!/usr/bin/env python3
"""
L_diagnostic.py

Chapter 20: Spectral Methods on Unbounded Intervals.
Computational Etude 20.9: Read the coefficients before the error.

Boyd's Rule-of-Thumb 16 turns the map parameter  L  from a guess into
a diagnostic:

    ``Plot the coefficients a_j versus degree on a log-linear plot.
    If the graph abruptly flattens at some N, then this implies that
    L is TOO SMALL for the given N, and one should increase L until
    the flattening is postponed to j = N.''

We approximate  f(y) = sech(y)  by the TB_n basis on (-infty, +infty)
for a sweep of L values and plot the |a_j| versus j.  Too-small L
produces an abrupt flattening at low j (the domain-truncation-like
tail becomes visible); as L grows, the flattening is pushed out.  For
too-large L the small-j slope becomes gentle (poles in the Chebyshev
x-plane are drifting toward x = +/- 1).  The valley of good L is
broad but not infinite.

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

from unbounded_common import (CORAL, NAVY, TEAL, PURPLE, ORANGE, output_dir_for,
                              save_fig, setup_matplotlib, dct1_coeffs,
                              tb_map_forward)

OUTPUT_DIR = output_dir_for(SCRIPT_DIR)


def target(y):
    return 1.0 / np.cosh(y)


def tbn_coeffs(N, L):
    _, x = cheb_matrix(N)
    y = tb_map_forward(x, L)
    # x = +/- 1 map to y = +/- infty where f = 0
    fv = np.where(np.abs(x) < 1.0 - 1e-12, target(y), 0.0)
    return dct1_coeffs(fv)


def make_figure():
    setup_matplotlib()
    N = 64
    L_list = [0.5, 1.0, 2.0, 4.0, 8.0, 16.0]
    colours = [CORAL, ORANGE, TEAL, PURPLE, NAVY, "#8B4513"]

    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.0))

    ax = axes[0]
    for L, c in zip(L_list, colours):
        a = np.abs(tbn_coeffs(N, L))
        ax.semilogy(np.arange(N + 1), a + 1e-18, "-o", ms=3,
                    color=c, label=fr"$L = {L}$")
    ax.set_xlabel(r"degree $n$")
    ax.set_ylabel(r"$|a_n|$ of $TB_n$ expansion")
    ax.set_title(r"(a) Coefficient decay of $\mathrm{sech}(y)$, $N = 64$")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(frameon=False, fontsize=9, loc="lower left")
    ax.set_ylim(1e-17, 10)

    ax = axes[1]
    # error versus L for two resolutions
    Ls = np.array([0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 16.0, 24.0])
    for N_ref, col, marker in [(24, CORAL, "o"), (48, TEAL, "s"), (96, NAVY, "^")]:
        errs = []
        for L in Ls:
            a = tbn_coeffs(N_ref, L)
            # rough error estimate: sum of |a_n| for n > N_ref / 2, as proxy for tail
            errs.append(np.sum(np.abs(a[N_ref // 2:])))
        ax.semilogy(Ls, errs, "-" + marker, color=col, mfc="none" if marker != "o" else col,
                    label=fr"$N = {N_ref}$")
    ax.set_xscale("log")
    ax.set_xlabel(r"$L$")
    ax.set_ylabel(r"$\sum_{n > N/2} |a_n|$ (tail size)")
    ax.set_title("(b) Valley of good $L$ broadens with $N$")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(frameon=False, fontsize=9)

    fig.tight_layout()
    save_fig(fig, OUTPUT_DIR, "L_diagnostic")
    plt.close(fig)

    print(f"[20.9] saved figure to {OUTPUT_DIR / 'L_diagnostic.pdf'}")
    for L in L_list:
        a = np.abs(tbn_coeffs(N, L))
        below = np.argmax(a[:N + 1] < 1e-12) if (a < 1e-12).any() else N + 1
        print(f"  L={L:4.1f}  first n with |a_n| < 1e-12: {below}")


if __name__ == "__main__":
    make_figure()
