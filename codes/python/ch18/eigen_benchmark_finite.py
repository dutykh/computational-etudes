#!/usr/bin/env python3
"""
eigen_benchmark_finite.py

Chapter 18: Linear Spectral Eigenproblems
Computational Etude 18.3: Finite-interval benchmark and the N/2 rule.

Repeats the Dirichlet Laplacian calculation at N = 16, 32, 64 and plots
the absolute error |lambda_j^num - lambda_j^exact| on a semilog scale
against mode number j. Reproduces Boyd (2000) Fig 7.1 and Fig 7.3.

At each N, the lowest N/2 eigenvalues are accurate to within a few
percent (Boyd's Rule-of-Thumb 8); beyond that threshold the error
climbs exponentially until it saturates near the eigenvalue magnitude
itself (random, meaningless numbers).

Interpretation of the figure:
  - The flat roundoff floor at 1e-13 for low modes is the combined effect
    of floating-point accumulation in eig() and the sigma-scaled
    condition number of the interior Laplacian.
  - The steep rise starts near j ~ N/2 (vertical dividing lines in the
    plot) and levels off at O(lambda_j) for the unresolved modes.
  - Doubling N shifts the "cliff" rightward by roughly N, buying an
    equal number of additional trusted eigenvalues --- EXACTLY the
    behaviour Boyd's heuristic predicts.

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
sys.path.insert(0, str(SCRIPT_DIR.parent / "ch07"))
from cheb_matrix import cheb_matrix  # noqa: E402

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


def solve_dirichlet(N: int) -> np.ndarray:
    """Return sorted interior-block spectrum of u_xx + lam u = 0, u(+-1) = 0."""
    D, _ = cheb_matrix(N)
    D2 = D @ D
    A = -D2[1:N, 1:N]
    lam = np.linalg.eigvals(A).real
    return np.sort(lam)


def count_good(err: np.ndarray, tol: float = 1e-2) -> int:
    """Return the LARGEST j such that every error up to j is below tol."""
    below = err < tol
    if not below.any():
        return 0
    # first failure index
    failures = np.where(~below)[0]
    return failures[0] if failures.size > 0 else err.size


def make_figure(Ns=(16, 32, 64)):
    colours = [NAVY, CORAL, TEAL]
    markers = ["o", "s", "^"]
    counts = {}

    fig, ax = plt.subplots(1, 1, figsize=(7.6, 4.6))

    for N, c, m in zip(Ns, colours, markers):
        lam = solve_dirichlet(N)
        j = np.arange(1, lam.size + 1)
        lam_exact = (j * np.pi / 2) ** 2
        err = np.abs(lam - lam_exact)
        err_plot = np.maximum(err, 1e-17)
        ax.semilogy(j, err_plot, marker=m, linestyle="-", color=c,
                    markerfacecolor="white", markeredgewidth=1.0,
                    markersize=5, linewidth=0.6,
                    label=f"$N = {N}$")
        # vertical cue at N/2
        ax.axvline(N / 2, color=c, linestyle=":", linewidth=1.4, alpha=0.85)
        counts[N] = count_good(err, tol=1e-2)

    ax.axhline(1e-2, color="k", linestyle="--", linewidth=0.6, alpha=0.5)
    ax.text(1.5, 2e-2, "0.01 tolerance", fontsize=8, color="k", alpha=0.7)
    ax.set_xlabel("mode number $j$")
    ax.set_ylabel(r"$|\lambda^{\mathrm{num}}_j - \lambda^{\mathrm{exact}}_j|$")
    ax.set_title("Dirichlet Laplacian: errors at three resolutions",
                 fontsize=11)
    ax.set_xlim(0, max(Ns) + 2)
    ax.set_ylim(1e-17, 1e6)
    ax.grid(True, which="both", alpha=0.25, linewidth=0.4)
    ax.legend(loc="lower right", fontsize=9, framealpha=0.9)

    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "eigen_benchmark_finite.pdf", bbox_inches="tight", pad_inches=0.05)
    fig.savefig(OUTPUT_DIR / "eigen_benchmark_finite.png", bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)

    spectra = {int(N): solve_dirichlet(N).tolist() for N in Ns}
    return {
        "Ns":            list(Ns),
        "spectra":       spectra,
        "good_counts":   {int(N): int(k) for N, k in counts.items()},
    }


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--dump", type=str, default=None)
    args = p.parse_args()

    r = make_figure()
    print(f"[Etude 18.3]")
    for N in r["Ns"]:
        print(f"  N = {N:2d}:  {r['good_counts'][N]} eigenvalues with |err| < 1e-2   (heuristic N/2 = {N//2})")
    print(f"  figure: {OUTPUT_DIR / 'eigen_benchmark_finite.pdf'}")
    if args.dump:
        Path(args.dump).write_text(json.dumps(r, indent=2))
        print(f"  dumped: {args.dump}")


if __name__ == "__main__":
    main()
