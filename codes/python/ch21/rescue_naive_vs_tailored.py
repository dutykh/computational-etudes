#!/usr/bin/env python3
"""
rescue_naive_vs_tailored.py

Chapter 21: Special Tricks for the Spectral Researcher
Computational Etude 21.1: a rescue story in miniature.

A small parable for the chapter.  We approximate the same exponentially
decaying pulse f(y) = sech(y) on a *truncated* finite interval [-L, L]
with the *same* Chebyshev resolution N = 32, but choose two very
different L:

  NAIVE:   L = 20.  An honest, conservative choice -- 'make the window
           large enough to contain everything'.  The pulse sech(y)
           occupies a width of order 6 around y = 0, so 88% of the
           Chebyshev grid lies outside the support of the function.
           That is wasted resolution, and it costs us decimals.

  TAILORED: L = 8.  Just large enough that sech(L) ~ 6.7e-4 sits below
           the truncation noise we expect.  Now the Chebyshev points
           live in the support of f, and the polynomial sees a problem
           it can resolve.

This is the chapter's first lesson: the truncation length L is part of
the discretisation, not a cosmetic detail.  At N = 48 the naive choice
gives a maximum pointwise error of order 10^{-1}; the tailored choice
gives 10^{-4}.  Three decimals, recovered by editing one line of code,
is the chapter's opening sermon.

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


def f(y):
    return 1.0 / np.cosh(y)


def cgl(N):
    """Chebyshev-Gauss-Lobatto nodes on [-1, 1] in descending order."""
    return np.cos(np.pi * np.arange(N + 1) / N)


def cheb_coeffs(values_on_cgl):
    """DCT-I coefficients of samples on the (N+1)-point CGL grid."""
    N = len(values_on_cgl) - 1
    extended = np.concatenate([values_on_cgl, values_on_cgl[N - 1:0:-1]])
    A = np.real(np.fft.fft(extended)) / N
    A[0] *= 0.5
    A[N] *= 0.5
    return A[:N + 1]


def reconstruct(a, t_eval):
    """Clenshaw evaluation of a Chebyshev series at arbitrary points t in [-1,1]."""
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


def truncated_chebyshev_error(N, L, y_ref):
    """Build the Chebyshev approximation of f on [-L, L] and return
    the absolute error on the finite reference window y_ref, plus the
    Chebyshev coefficient magnitudes on the truncated interval."""
    t_grid = cgl(N)
    y_grid = L * t_grid               # CGL points pulled back to [-L, L]
    a = cheb_coeffs(f(y_grid))
    # Reconstruct at y_ref by mapping y_ref -> t = y_ref / L.
    t_ref = y_ref / L
    inside = np.abs(t_ref) <= 1.0
    approx = np.full_like(y_ref, np.nan)
    approx[inside] = reconstruct(a, t_ref[inside])
    err = np.abs(approx - f(y_ref))
    err[~inside] = np.nan
    return a, approx, err


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    p.add_argument("--dump", type=str, default=None)
    args = p.parse_args()

    setup_matplotlib()
    out = output_dir_for(SCRIPT_DIR)

    N = 48
    L_naive = 30.0
    L_tail  = 8.0

    y_ref = np.linspace(-35.0, 35.0, 6001)
    a_naive, _, err_naive = truncated_chebyshev_error(N, L_naive, y_ref)
    a_tail,  _, err_tail  = truncated_chebyshev_error(N, L_tail,  y_ref)

    # Worst-case error inside *each window* (over the entire support of f).
    inside_naive = np.abs(y_ref) <= L_naive
    inside_tail  = np.abs(y_ref) <= L_tail
    e_naive = float(np.nanmax(err_naive[inside_naive]))
    e_tail  = float(np.nanmax(err_tail[inside_tail]))

    fig = plt.figure(figsize=(11.5, 4.4))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.05, 1.0], wspace=0.28)

    # --- Left panel: function with both grids overlaid -------------------
    ax = fig.add_subplot(gs[0, 0])
    ax.plot(y_ref, f(y_ref), color=NAVY, linewidth=1.4,
            label=r"$f(y)=\mathrm{sech}\,y$")
    # Naive grid: 33 markers spread over [-20, 20]
    y_naive_grid = L_naive * cgl(N)
    ax.scatter(y_naive_grid, f(y_naive_grid), s=22, marker="o",
               edgecolors=CORAL, facecolors="white", linewidths=1.0,
               label=fr"naive grid: Chebyshev on $[-{int(L_naive)},{int(L_naive)}]$")
    # Tailored grid:
    y_tail_grid = L_tail * cgl(N)
    ax.scatter(y_tail_grid, f(y_tail_grid), s=28, marker="^",
               edgecolors=TEAL, facecolors="white", linewidths=1.0,
               label=fr"tailored grid: Chebyshev on $[-{int(L_tail)},{int(L_tail)}]$")
    ax.axvline(0.0, color="gray", linewidth=0.4, alpha=0.4)
    ax.set_xlabel(r"$y$")
    ax.set_ylabel(r"$f(y)$")
    ax.set_title(rf"same $f$, same $N={N}$, two truncation lengths $L$",
                 fontsize=10)
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(True, alpha=0.25, linewidth=0.4)
    ax.set_xlim(-22, 22)
    ax.set_ylim(-0.06, 1.12)

    # --- Right panel: pointwise error (log scale) ------------------------
    ax = fig.add_subplot(gs[0, 1])
    ax.semilogy(y_ref[inside_naive], err_naive[inside_naive] + 1e-18,
                color=CORAL, linewidth=1.0,
                label=rf"naive ($L={int(L_naive)}$): max err $\approx${e_naive:.1e}")
    ax.semilogy(y_ref[inside_tail], err_tail[inside_tail] + 1e-18,
                color=TEAL, linewidth=1.0,
                label=rf"tailored ($L={int(L_tail)}$): max err $\approx${e_tail:.1e}")
    ax.axhline(1e-15, color="gray", linewidth=0.4, alpha=0.5)
    ax.set_xlabel(r"$y$")
    ax.set_ylabel(r"$|f(y)-f_N(y)|$")
    ax.set_title("pointwise error: a single tuned $L$ buys decimals", fontsize=10)
    ax.legend(loc="lower center", fontsize=9)
    ax.grid(True, alpha=0.25, linewidth=0.4)
    ax.set_xlim(-22, 22)
    ax.set_ylim(1e-16, 1e0)

    save_fig(fig, out, "rescue_naive_vs_tailored")
    plt.close(fig)

    print(f"[Etude 21.1]  rescue story")
    print(f"  N = {N}, f(y) = sech(y), |f(L)| = sech(L) ~ truncation noise floor")
    print(f"  naive    L = {L_naive:5.2f}: max err = {e_naive:.3e}")
    print(f"  tailored L = {L_tail:5.2f}: max err = {e_tail:.3e}")
    print(f"  ratio (tailored / naive) = {e_tail / e_naive:.2e}")
    print(f"  figure: {out / 'rescue_naive_vs_tailored.pdf'}")

    dump_json(args.dump, {
        "N": N,
        "L_naive": L_naive, "L_tailored": L_tail,
        "err_naive": e_naive, "err_tailored": e_tail,
    })


if __name__ == "__main__":
    main()
