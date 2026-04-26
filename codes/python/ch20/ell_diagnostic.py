#!/usr/bin/env python3
"""
ell_diagnostic.py

Chapter 20: Spectral Methods on Unbounded Intervals.
Computational Etude 20.9: Read the coefficients before the error.

Boyd's Rule-of-Thumb 16 turns the map parameter  ell  from a guess into
a diagnostic:

    ``Plot the coefficients a_j versus degree on a log-linear plot.
    If the graph abruptly flattens at some N, then this implies that
    ell is TOO SMALL for the given N, and one should increase ell until
    the flattening is postponed to j = N.''

We approximate  f(y) = sech(y)  by the TB_n basis on (-infty, +infty)
for a sweep of ell values and plot the |a_j| versus j.  Too-small ell
produces an abrupt flattening at low j (the domain-truncation-like
tail becomes visible); as ell grows, the flattening is pushed out.  For
too-large ell the small-j slope becomes gentle (poles in the Chebyshev
x-plane are drifting toward x = +/- 1).  The valley of good ell is
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


def tbn_coeffs(N, ell):
    _, x = cheb_matrix(N)
    y = tb_map_forward(x, ell)
    # x = +/- 1 map to y = +/- infty where f = 0
    fv = np.where(np.abs(x) < 1.0 - 1e-12, target(y), 0.0)
    return dct1_coeffs(fv)


def envelope(a, win=3):
    """Rolling-max envelope to suppress the parity zig-zag in |a_n|."""
    out = np.copy(a)
    n = len(a)
    for i in range(n):
        lo = max(0, i - win)
        hi = min(n, i + win + 1)
        out[i] = np.max(a[lo:hi])
    return out


def make_figure():
    setup_matplotlib()
    N = 64
    L_full = [0.5, 1.0, 2.0, 4.0, 8.0, 16.0]
    colours_full = [CORAL, ORANGE, TEAL, PURPLE, NAVY, "#8B4513"]

    # 2x2 layout:
    #   (a) clean three-regime view (small, good, large ell)
    #   (b) full sweep, envelope-only
    #   (c) tail-size vs ell at three N (kept from old panel b)
    #   (d) location of the broad-valley centre vs N
    fig, axes = plt.subplots(2, 2, figsize=(11.0, 8.0))

    # ---- (a) Clean three-regime view (envelopes) -----------------------
    ax = axes[0, 0]
    three = [(0.5, CORAL, "small (early flatten)"),
             (2.0, NAVY, "good (clean descent)"),
             (16.0, TEAL, "large (gentle small-$n$ slope)")]
    for ell, c, label in three:
        a = np.abs(tbn_coeffs(N, ell))
        env = envelope(a, win=3)
        ns = np.arange(N + 1)
        ax.semilogy(ns, env + 1e-18, color=c, lw=1.2,
                    label=fr"$\ell = {ell:g}$ — {label}")
        ax.semilogy(ns[::4], env[::4] + 1e-18, "o", color=c, ms=4,
                    mfc="white")
    ax.set_xlabel(r"degree $n$")
    ax.set_ylabel(r"envelope of $|a_n|$")
    ax.set_title(r"(a) three regimes of $\ell$ (envelope view), $N = 64$")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(frameon=False, fontsize=8, loc="lower left")
    ax.set_ylim(1e-17, 10)

    # ---- (b) Full six-ell sweep, envelope only -------------------------
    ax = axes[0, 1]
    for ell, c in zip(L_full, colours_full):
        a = np.abs(tbn_coeffs(N, ell))
        env = envelope(a, win=3)
        ax.semilogy(np.arange(N + 1), env + 1e-18, color=c, lw=0.9,
                    label=fr"$\ell = {ell:g}$")
    ax.set_xlabel(r"degree $n$")
    ax.set_ylabel(r"envelope of $|a_n|$")
    ax.set_title(r"(b) full $\ell$ sweep at $N = 64$, envelope only")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(frameon=False, fontsize=8, ncol=2, loc="lower left")
    ax.set_ylim(1e-17, 10)

    # ---- (c) Tail size vs ell at three resolutions ---------------------
    ax = axes[1, 0]
    Ls = np.array([0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0,
                   12.0, 16.0, 24.0])
    tails_by_N = {}
    for N_ref, col, marker in [(24, CORAL, "o"), (48, TEAL, "s"),
                                (96, NAVY, "^")]:
        errs = []
        for ell in Ls:
            a = tbn_coeffs(N_ref, ell)
            errs.append(np.sum(np.abs(a[N_ref // 2:])))
        errs = np.array(errs)
        tails_by_N[N_ref] = errs
        ax.semilogy(Ls, errs, "-" + marker, color=col,
                    mfc="none" if marker != "o" else col,
                    label=fr"$N = {N_ref}$")
    ax.set_xscale("log")
    ax.set_xlabel(r"$\ell$")
    ax.set_ylabel(r"$\sum_{n > N/2} |a_n|$ (tail size)")
    ax.set_title("(c) Valley of good $\\ell$ broadens with $N$")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(frameon=False, fontsize=9)

    # ---- (d) NEW: location of the broad-valley centre vs N -----------
    # For each resolution, take the geometric mean of the ell values
    # whose tail size is within a factor of 3 of the minimum -- this
    # gives a stable estimate of the valley centre, less sensitive to
    # the parity zig-zag than argmin alone.
    ax = axes[1, 1]
    Ns_extra = [24, 32, 48, 64, 96, 128]
    centres, widths = [], []
    for N_ref in Ns_extra:
        errs_N = np.array([np.sum(np.abs(tbn_coeffs(N_ref, ell)[N_ref // 2:]))
                           for ell in Ls])
        emin = errs_N.min()
        good = errs_N <= 3.0 * emin
        centre = np.exp(np.mean(np.log(Ls[good])))
        width = Ls[good].max() / Ls[good].min()
        centres.append(centre)
        widths.append(width)
    ax.loglog(Ns_extra, centres, "-o", color=NAVY, lw=1.1,
              mfc="white", ms=6, label=r"valley centre $\bar{\ell}^*(N)$")
    ax.loglog(Ns_extra, widths, "-s", color=TEAL, lw=1.0,
              mfc="white", ms=5,
              label=r"valley width $\ell_{\max}/\ell_{\min}$ (factor)")
    ax.set_xlabel(r"$N$")
    ax.set_ylabel("value (log-log)")
    ax.set_title(r"(d) broad-valley centre and width vs $N$")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(frameon=False, fontsize=9)

    fig.tight_layout()
    save_fig(fig, OUTPUT_DIR, "ell_diagnostic")
    plt.close(fig)

    print(f"[20.9] saved figure to {OUTPUT_DIR / 'ell_diagnostic.pdf'}")
    for ell in L_full:
        a = np.abs(tbn_coeffs(N, ell))
        below = np.argmax(a[:N + 1] < 1e-12) if (a < 1e-12).any() else N + 1
        print(f"  ell={ell:4.1f}  first n with |a_n| < 1e-12: {below}")


if __name__ == "__main__":
    make_figure()
