#!/usr/bin/env python3
"""
eigen_benchmark_oscillator.py

Chapter 18: Linear Spectral Eigenproblems
Computational Etude 18.4: the infinite-interval tax.

Solves the quantum harmonic oscillator
        u_xx + (lambda - x^2) u = 0,   |u| -> 0 as |x| -> inf,
using the rational_chebyshev utility (algebraic map x = L t / sqrt(1-t^2))
at two resolutions (N = 16 and N = 32) and for three map parameters.

Exact eigenvalues are lambda_j = 2 j + 1 for j = 0, 1, 2, ... with
Hermite-function eigenfunctions H_j(x) exp(-x^2/2).

Reproduces Boyd (2000) Figs 7.4 and 7.6. The pedagogical contrast with
Etude 18.3 is deliberate: for the bounded benchmark, N = 16 gives 6
trusted eigenvalues; for the unbounded benchmark, the same N buys fewer
good modes --- the "infinite-interval tax".

We also scan L in {2, 4, 8} to show that map-parameter choice is a
first-class tuning knob: too small an L crowds the grid near x = 0,
too large an L wastes nodes where the wavefunction has already decayed.

Author: Dr. Denys Dutykh (Khalifa University of Science and Technology)
Part of "Computational Etudes: A Spectral Approach"
https://github.com/dutykh/computational-etudes
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
from rational_chebyshev import rational_chebyshev_derivative_matrices  # noqa: E402

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["CMU Serif", "DejaVu Serif"],
    "mathtext.fontset": "cm",
    "font.size": 11,
    "axes.linewidth": 0.8,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "figure.dpi": 150,
    "savefig.dpi": 300,
})

NAVY, SKY, CORAL = "#142D6E", "#7896D2", "#E74C3C"
TEAL, PURPLE, ORANGE = "#16A085", "#8E44AD", "#E67E22"

OUTPUT_DIR = SCRIPT_DIR.parent.parent.parent / "textbook" / "figures" / "ch18" / "python"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def solve_oscillator(N: int, L: float):
    """Return sorted positive-real spectrum of H = -d^2/dx^2 + x^2."""
    _, D2_x, x = rational_chebyshev_derivative_matrices(N, L)
    H = -D2_x + np.diag(x ** 2)
    lam = np.linalg.eigvals(H)
    lam = lam.real[np.abs(lam.imag) < 1e-6]
    lam = lam[lam > 0]
    return np.sort(lam)


def count_good(err: np.ndarray, tol: float = 1e-2) -> int:
    failures = np.where(err >= tol)[0]
    return failures[0] if failures.size > 0 else err.size


def make_figure(Ns=(16, 32), L_scan=(2.0, 4.0, 8.0), L_best: float = 4.0):
    # --- Panel A: L = L_best, two resolutions, error vs mode ---
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.5))

    markers = ["o", "s"]
    colors = [NAVY, CORAL]
    ax = axes[0]
    counts_by_N = {}
    for N, m, c in zip(Ns, markers, colors):
        lam = solve_oscillator(N, L_best)
        j = np.arange(lam.size)
        lam_exact = 2 * j + 1
        err = np.abs(lam - lam_exact)
        err_plot = np.maximum(err, 1e-17)
        ax.semilogy(j + 1, err_plot, marker=m, linestyle="-",
                    color=c, markerfacecolor="white",
                    markersize=5, linewidth=0.6, markeredgewidth=1.0,
                    label=f"$N = {N}$, $L = {L_best}$")
        counts_by_N[int(N)] = int(count_good(err, 1e-2))
    ax.axhline(1e-2, color="k", linestyle="--", lw=0.6, alpha=0.5)
    ax.set_xlabel("mode number $j$")
    ax.set_ylabel(r"$|\lambda^{\mathrm{num}}_j - \lambda^{\mathrm{exact}}_j|$")
    ax.set_title("infinite-interval tax: oscillator spectrum", fontsize=10)
    ax.set_ylim(1e-17, 1e4)
    ax.grid(True, which="both", alpha=0.25, linewidth=0.4)
    ax.legend(loc="lower right", fontsize=9)

    # --- Panel B: L scan at fixed N = 32 ---
    ax = axes[1]
    colors_B = [SKY, NAVY, PURPLE]
    markers_B = ["^", "o", "v"]
    counts_by_L = {}
    N_scan = 32
    for L, m, c in zip(L_scan, markers_B, colors_B):
        lam = solve_oscillator(N_scan, L)
        j = np.arange(lam.size)
        lam_exact = 2 * j + 1
        err = np.abs(lam - lam_exact)
        err_plot = np.maximum(err, 1e-17)
        ax.semilogy(j + 1, err_plot, marker=m, linestyle="-",
                    color=c, markerfacecolor="white",
                    markersize=5, linewidth=0.6, markeredgewidth=1.0,
                    label=f"$L = {L}$")
        counts_by_L[float(L)] = int(count_good(err, 1e-2))
    ax.axhline(1e-2, color="k", linestyle="--", lw=0.6, alpha=0.5)
    ax.set_xlabel("mode number $j$")
    ax.set_ylabel(r"$|\lambda^{\mathrm{num}}_j - \lambda^{\mathrm{exact}}_j|$")
    ax.set_title(rf"$L$-scan at $N = {N_scan}$: map-parameter sensitivity", fontsize=10)
    ax.set_ylim(1e-17, 1e4)
    ax.grid(True, which="both", alpha=0.25, linewidth=0.4)
    ax.legend(loc="lower right", fontsize=9)

    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "eigen_benchmark_oscillator.pdf", bbox_inches="tight", pad_inches=0.05)
    fig.savefig(OUTPUT_DIR / "eigen_benchmark_oscillator.png", bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)

    spectra = {
        f"N{N}_L{L_best}": solve_oscillator(N, L_best).tolist() for N in Ns
    }
    spectra.update({
        f"N{N_scan}_L{L}": solve_oscillator(N_scan, L).tolist() for L in L_scan
    })
    return {
        "Ns":              list(Ns),
        "L_best":          L_best,
        "L_scan":          list(L_scan),
        "counts_by_N":     counts_by_N,
        "counts_by_L":     counts_by_L,
        "spectra":         spectra,
    }


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--dump", type=str, default=None)
    args = p.parse_args()

    r = make_figure()
    print("[Etude 18.4]  harmonic oscillator on (-inf, +inf)")
    print(f"  L = {r['L_best']} scan over N:")
    for N in r["Ns"]:
        print(f"    N = {N:2d}:  {r['counts_by_N'][N]} good eigenvalues (|err| < 1e-2)")
    print(f"  N = 32 scan over L:")
    for L in r["L_scan"]:
        print(f"    L = {L}:  {r['counts_by_L'][L]} good eigenvalues")
    print(f"  figure: {OUTPUT_DIR / 'eigen_benchmark_oscillator.pdf'}")
    if args.dump:
        Path(args.dump).write_text(json.dumps(r, indent=2))
        print(f"  dumped: {args.dump}")


if __name__ == "__main__":
    main()
