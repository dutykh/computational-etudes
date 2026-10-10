#!/usr/bin/env python3
"""
eigen_drift_diagnostic.py

Chapter 18: Linear Spectral Eigenproblems
Computational Etude 18.5: build a spectrum lie detector.

Applies the drift-with-N diagnostic (spectrum_verify utility) to two
problems studied earlier in the chapter:

  (A) Dirichlet Laplacian u_xx + lambda u = 0 on (-1, 1), at N1 = 32
      and N2 = 48. Exact spectrum, well-ordered modes, ordinal matching
      sufficient.

  (B) Harmonic oscillator u_xx + (lambda - x^2) u = 0 on the real line,
      via rational_chebyshev with L = 4, at N1 = 32 and N2 = 48. Exact
      spectrum (Hermite), decay at infinity introduces the "infinite-
      interval tax".

For each problem we plot the reciprocal scaled drift 1/delta vs. mode
number j on a semilog axis (following Boyd Fig 7.7): trusted modes
appear at the TOP, suspect modes at the bottom. The cliff between the
two is the boundary of the trustworthy spectrum.

Reproduces the pedagogical lesson of Boyd Chapter 7.5: both the
finite-interval and unbounded problems show a sharp drop in reciprocal
drift ratio near the N/2 threshold.

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
sys.path.insert(0, str(SCRIPT_DIR.parent / "ch07"))
from cheb_matrix import cheb_matrix  # noqa: E402
from rational_chebyshev import rational_chebyshev_derivative_matrices  # noqa: E402
from spectrum_verify import verify_spectrum  # noqa: E402

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

NAVY, SKY, CORAL, TEAL = "#142D6E", "#7896D2", "#E74C3C", "#16A085"

OUTPUT_DIR = SCRIPT_DIR.parent.parent.parent / "textbook" / "figures" / "ch18" / "python"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def solve_laplacian(N: int) -> np.ndarray:
    D, _ = cheb_matrix(N)
    A = -(D @ D)[1:N, 1:N]
    return np.sort(np.linalg.eigvals(A).real)


def solve_oscillator(N: int, ell: float = 4.0) -> np.ndarray:
    _, D2, x = rational_chebyshev_derivative_matrices(N, ell)
    H = -D2 + np.diag(x ** 2)
    lam = np.linalg.eigvals(H)
    lam = lam.real[np.abs(lam.imag) < 1e-6]
    lam = lam[lam > 0]
    return np.sort(lam)


def plot_reciprocal_drift(ax, report, title, tol):
    j = np.arange(1, report.lam1.size + 1)
    # reciprocal drift: trusted modes at top
    inv_ord = 1.0 / np.maximum(report.delta_ordinal, 1e-18)
    inv_nst = 1.0 / np.maximum(report.delta_nearest, 1e-18)
    ax.semilogy(j, inv_ord, "o", color=NAVY, markersize=4,
                markerfacecolor="white", markeredgewidth=1.0,
                label=r"$1/\delta_{j,\mathrm{ordinal}}$")
    ax.semilogy(j, inv_nst, "x", color=CORAL, markersize=5,
                markeredgewidth=1.0,
                label=r"$1/\delta_{j,\mathrm{nearest}}$")
    ax.axhline(1.0 / tol, color="k", linestyle="--", linewidth=0.6,
               alpha=0.5)
    ax.text(1, 2.0 / tol, fr"$1/$tol $= 1/{tol:.0e}$",
            fontsize=8, color="k", alpha=0.7)
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("mode number $j$")
    ax.set_ylabel(r"$1 / \delta_j$  (trusted = high)")
    ax.grid(True, which="both", alpha=0.25, linewidth=0.4)
    ax.legend(loc="upper right", fontsize=9)


def make_figure(N1: int = 32, N2: int = 48, tol: float = 1e-3):
    # Problem A: Dirichlet Laplacian
    lam1_A = solve_laplacian(N1)
    lam2_A = solve_laplacian(N2)
    rep_A = verify_spectrum(lam1_A, lam2_A, tol=tol)

    # Problem B: Oscillator on the real line
    lam1_B = solve_oscillator(N1)
    lam2_B = solve_oscillator(N2)
    rep_B = verify_spectrum(lam1_B, lam2_B, tol=tol)

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.5))
    plot_reciprocal_drift(axes[0], rep_A,
        rf"Dirichlet Laplacian,  $N_1 = {N1}$, $N_2 = {N2}$", tol)
    plot_reciprocal_drift(axes[1], rep_B,
        rf"Harmonic oscillator (real line),  $N_1 = {N1}$, $N_2 = {N2}$", tol)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "eigen_drift_diagnostic.pdf", bbox_inches="tight", pad_inches=0.05)
    fig.savefig(OUTPUT_DIR / "eigen_drift_diagnostic.png", bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)

    return {
        "N1": N1, "N2": N2, "tol": tol,
        "laplacian": {
            "n_trusted":       int(rep_A.n_trusted),
            "lam1":            rep_A.lam1.tolist(),
            "lam2":            rep_A.lam2.tolist(),
            "delta_ordinal":   rep_A.delta_ordinal.tolist(),
            "delta_nearest":   rep_A.delta_nearest.tolist(),
            "trusted":         rep_A.trusted.astype(int).tolist(),
        },
        "oscillator": {
            "n_trusted":       int(rep_B.n_trusted),
            "lam1":            rep_B.lam1.tolist(),
            "lam2":            rep_B.lam2.tolist(),
            "delta_ordinal":   rep_B.delta_ordinal.tolist(),
            "delta_nearest":   rep_B.delta_nearest.tolist(),
            "trusted":         rep_B.trusted.astype(int).tolist(),
        },
    }


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--N1",  type=int,   default=32)
    p.add_argument("--N2",  type=int,   default=48)
    p.add_argument("--tol", type=float, default=1e-3)
    p.add_argument("--dump", type=str,  default=None)
    args = p.parse_args()

    r = make_figure(N1=args.N1, N2=args.N2, tol=args.tol)
    print(f"[Etude 18.5]  drift diagnostic at N1 = {args.N1}, N2 = {args.N2}, tol = {args.tol}")
    print(f"  Dirichlet Laplacian : trusted = {r['laplacian']['n_trusted']} of {len(r['laplacian']['lam1'])}")
    print(f"  Harmonic oscillator  : trusted = {r['oscillator']['n_trusted']} of {len(r['oscillator']['lam1'])}")
    print(f"  figure: {OUTPUT_DIR / 'eigen_drift_diagnostic.pdf'}")
    if args.dump:
        Path(args.dump).write_text(json.dumps(r, indent=2))
        print(f"  dumped: {args.dump}")


if __name__ == "__main__":
    main()
