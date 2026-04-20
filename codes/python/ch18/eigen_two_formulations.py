#!/usr/bin/env python3
"""
eigen_two_formulations.py

Chapter 18: Linear Spectral Eigenproblems
Computational Etude 18.2: Assemble the pencil.

Same eigenproblem as Etude 18.1 --- the Dirichlet Laplacian
u_xx + lambda u = 0, u(+-1) = 0 --- solved two ways:

  (A) Boundary bordering: full (N+1)x(N+1) matrix pencil A v = lambda B v,
      where A is D^2 with rows 1 and N+1 replaced by identity rows, and
      B is the identity with zeros in those two rows. The two singular
      rows of B produce two infinite eigenvalues, which are filtered out.

  (B) Basis recombination: the Dirichlet condition u(+-1) = 0 is built
      into the basis by restricting to interior grid points. The result
      is the regular matrix eigenproblem -D^2_int u = lambda u.

Both formulations recover the same finite spectrum. We display the
spectra side by side and the sparsity pattern of A in each approach.

This is the pedagogical companion to Etude 18.1: two different
roads to the same matrix eigenvalue problem. The lesson is that the
matrix-pencil formulation is NOT an algebraic nuisance to be
multiplied away; it is part of the mathematical structure of the
operator problem.

Author: Dr. Denys Dutykh
        Mathematics Department
        Khalifa University of Science and Technology
        Abu Dhabi, UAE

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
import scipy.linalg

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


def exact_eigenvalue(j: int) -> float:
    return (j * np.pi / 2.0) ** 2


def boundary_bordering(N: int):
    """Full (N+1)x(N+1) pencil A v = lambda B v with rows 1 and N+1
    replaced by identity rows. Returns (A, B, x)."""
    D, x = cheb_matrix(N)
    D2 = D @ D
    A = -D2.copy()                      # -u_xx on all rows
    B = np.eye(N + 1)
    # impose u(1) = 0 at row 0 and u(-1) = 0 at row N
    A[0, :] = 0.0; A[0, 0] = 1.0; B[0, :] = 0.0
    A[N, :] = 0.0; A[N, N] = 1.0; B[N, :] = 0.0
    return A, B, x


def basis_recombination(N: int):
    """Interior-block regular problem -D^2_int u = lambda u.
    Returns (A_int, x)."""
    D, x = cheb_matrix(N)
    D2 = D @ D
    A_int = -D2[1:N, 1:N]
    return A_int, x


def solve_pencil(A, B):
    """Generalised eigensolve, filter infinite eigenvalues, sort ascending."""
    lam = scipy.linalg.eigvals(A, B)
    lam = lam[np.isfinite(lam)]
    lam = np.sort(lam.real[np.abs(lam.imag) < 1e-6])
    return lam


def solve_regular(A):
    lam = np.linalg.eigvals(A)
    lam = np.sort(lam.real[np.abs(lam.imag) < 1e-6])
    return lam


def make_figure(N: int = 16):
    A_pen, B_pen, x = boundary_bordering(N)
    A_reg, _         = basis_recombination(N)

    lam_pen = solve_pencil(A_pen, B_pen)         # N-1 finite eigenvalues
    lam_reg = solve_regular(A_reg)               # N-1 eigenvalues

    # Align shapes (both should return N-1 finite modes)
    m = min(len(lam_pen), len(lam_reg))
    lam_pen = lam_pen[:m]
    lam_reg = lam_reg[:m]

    j = np.arange(1, m + 1)
    lam_exact = (j * np.pi / 2.0) ** 2

    err_pen = np.abs(lam_pen - lam_exact)
    err_reg = np.abs(lam_reg - lam_exact)
    diff    = np.abs(lam_pen - lam_reg)         # cross-formulation agreement

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2))

    ax = axes[0]
    ax.semilogy(j, np.maximum(err_pen, 1e-17), "o", color=NAVY,
                markersize=6, markerfacecolor="white", markeredgewidth=1.1,
                label="boundary bordering (pencil)")
    ax.semilogy(j, np.maximum(err_reg, 1e-17), "s", color=CORAL,
                markersize=5, markerfacecolor="white", markeredgewidth=1.1,
                label="basis recombination (regular)")
    ax.semilogy(j, np.maximum(diff, 1e-17), "x", color=TEAL,
                markersize=6, markeredgewidth=1.1,
                label=r"|pencil $-$ regular|")
    ax.axhline(1e-2, color="k", ls="--", lw=0.6, alpha=0.5)
    ax.set_xlabel("mode number $j$")
    ax.set_ylabel(r"$|\lambda^{\mathrm{num}}_j - \lambda^{\mathrm{exact}}_j|$")
    ax.set_title(rf"both formulations, $N={N}$", fontsize=10)
    ax.set_ylim(1e-17, 1e5)
    ax.grid(True, alpha=0.25, linewidth=0.4, which="both")
    ax.legend(loc="lower right", fontsize=8, framealpha=0.9)

    # sparsity-style heatmaps of |A|
    ax = axes[1]
    M = np.log10(np.maximum(np.abs(A_pen), 1e-20))
    im = ax.imshow(M, cmap="Blues", aspect="equal")
    ax.set_title(r"$\log_{10}|A|$ for boundary-bordered pencil", fontsize=10)
    ax.set_xlabel("column index"); ax.set_ylabel("row index")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "eigen_two_formulations.pdf", bbox_inches="tight", pad_inches=0.05)
    fig.savefig(OUTPUT_DIR / "eigen_two_formulations.png", bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)

    return {
        "N": N,
        "lam_pencil":         lam_pen.tolist(),
        "lam_regular":        lam_reg.tolist(),
        "lam_exact":          lam_exact.tolist(),
        "max_abs_agreement":  float(diff.max()),
        "max_err_pencil":     float(err_pen.max()),
        "max_err_regular":    float(err_reg.max()),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--N", type=int, default=16)
    parser.add_argument("--dump", type=str, default=None)
    args = parser.parse_args()

    r = make_figure(N=args.N)
    print(f"[Etude 18.2]  N = {args.N}")
    print(f"  max |lam_pencil - lam_regular| = {r['max_abs_agreement']:.3e}")
    print(f"  max |lam_pencil  - exact|      = {r['max_err_pencil']:.3e}")
    print(f"  max |lam_regular - exact|      = {r['max_err_regular']:.3e}")
    print(f"  figure: {OUTPUT_DIR / 'eigen_two_formulations.pdf'}")
    if args.dump:
        Path(args.dump).write_text(json.dumps(r, indent=2))
        print(f"  dumped: {args.dump}")


if __name__ == "__main__":
    main()
