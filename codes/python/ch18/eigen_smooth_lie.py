#!/usr/bin/env python3
"""
eigen_smooth_lie.py

Chapter 18: Linear Spectral Eigenproblems
Computational Etude 18.1: The smooth lie.

Solves the benchmark eigenvalue problem

    u_xx + lambda u = 0,   u(+-1) = 0

using Chebyshev pseudospectral collocation at N = 16. The exact eigenvalues
are lambda_j = (j pi / 2)^2 with eigenfunctions cos(j pi x / 2) for odd j
and sin(j pi x / 2) for even j.

The etude reproduces the central pedagogical lesson of Boyd (2000) Fig 7.2:
the N=16 collocation recovers low modes to machine precision but returns
high modes that LOOK smooth and oscillatory, yet are quantitatively wrong.
We plot side by side:

    (a) mode j = 3          --- exact and numerical agree visually
    (b) near-top mode       --- numerical is smooth but unrelated to exact

The figure is the book's opening shock: smoothness is not enough.

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

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR.parent / "ch07"))
from cheb_matrix import cheb_matrix  # noqa: E402

# --- publication-quality matplotlib configuration ---
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
    """Exact eigenvalue lambda_j = (j pi / 2)^2 of u_xx + lam u = 0, u(+-1)=0."""
    return (j * np.pi / 2.0) ** 2


def exact_eigenfunction(x: np.ndarray, j: int) -> np.ndarray:
    """cos(j pi x / 2) for j odd, sin(j pi x / 2) for j even."""
    if j % 2 == 1:
        return np.cos(j * np.pi * x / 2.0)
    return np.sin(j * np.pi * x / 2.0)


def solve(N: int):
    """Assemble and solve the interior Laplacian eigenproblem -D^2_int u = lambda u.

    The interior block of -D^2 has size (N-1) x (N-1) and its eigenvalues
    approximate (j pi / 2)^2 for j = 1, 2, ... up to resolution.
    """
    D, x = cheb_matrix(N)
    D2 = D @ D
    A = -D2[1:N, 1:N]
    eigvals, eigvecs = np.linalg.eig(A)
    # sort by real part (they are essentially real for this self-adjoint problem)
    order = np.argsort(eigvals.real)
    eigvals = eigvals[order].real
    eigvecs = eigvecs[:, order].real
    # sign-fix: make the first interior component positive (stabilises plots/JSON)
    for k in range(eigvecs.shape[1]):
        if eigvecs[0, k] < 0:
            eigvecs[:, k] *= -1
    return eigvals, eigvecs, x


def pad_with_boundary_zeros(v_int: np.ndarray) -> np.ndarray:
    """Wrap an interior eigenvector with u(+-1) = 0 boundary values."""
    return np.concatenate(([0.0], v_int, [0.0]))


def normalise_max(u: np.ndarray) -> np.ndarray:
    """Normalise so max|u| = 1 and the global max is positive."""
    m = np.max(np.abs(u))
    if m == 0:
        return u
    u = u / m
    # orient so the first nonzero extremum is positive
    idx = np.argmax(np.abs(u))
    if u[idx] < 0:
        u = -u
    return u


def make_figure(N: int = 16, low_mode: int = 3, high_mode: int = 15):
    eigvals, eigvecs, x = solve(N)
    # The numerical modes are the j-th sorted eigenvectors (j = 1 <-> index 0).
    # Boyd's Fig 7.2 shows mode 3 (top) and mode 15 (bottom) at N=16.
    # With N=16 interior = 15, the index for "mode 15" is 14 (zero-based).
    num_low_idx  = low_mode - 1
    num_high_idx = high_mode - 1

    # pad with boundary zeros for a cleaner curve
    u_low_num  = pad_with_boundary_zeros(eigvecs[:, num_low_idx])
    u_high_num = pad_with_boundary_zeros(eigvecs[:, num_high_idx])
    u_low_num  = normalise_max(u_low_num)
    u_high_num = normalise_max(u_high_num)

    xfine = np.linspace(-1, 1, 801)
    u_low_exact  = normalise_max(exact_eigenfunction(xfine, low_mode))
    u_high_exact = normalise_max(exact_eigenfunction(xfine, high_mode))

    # eigenvalue comparison
    lam_num_low    = eigvals[num_low_idx]
    lam_num_high   = eigvals[num_high_idx]
    lam_exact_low  = exact_eigenvalue(low_mode)
    lam_exact_high = exact_eigenvalue(high_mode)

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2))

    ax = axes[0]
    ax.plot(xfine, u_low_exact, color=NAVY, lw=1.2, label=f"exact mode {low_mode}")
    ax.plot(x, u_low_num, "o", color=CORAL, markersize=4.5,
            markerfacecolor="white", markeredgewidth=1.0,
            label=f"numerical (N={N})")
    ax.set_title(
        rf"mode $j={low_mode}$: $\lambda_{{\mathrm{{num}}}}={lam_num_low:.4f}$,"
        rf" $\lambda_{{\mathrm{{exact}}}}={lam_exact_low:.4f}$",
        fontsize=10,
    )
    ax.set_xlabel("$x$")
    ax.set_ylabel("$u(x)$")
    ax.set_xlim(-1, 1)
    ax.set_ylim(-1.15, 1.15)
    ax.grid(True, alpha=0.25, linewidth=0.4)
    ax.legend(loc="lower center", fontsize=9, framealpha=0.9)

    ax = axes[1]
    ax.plot(xfine, u_high_exact, color=NAVY, lw=1.2, label=f"exact mode {high_mode}")
    ax.plot(x, u_high_num, "o-", color=CORAL, markersize=4.5, lw=0.8,
            markerfacecolor="white", markeredgewidth=1.0,
            label=f"numerical (N={N})")
    ax.set_title(
        rf"mode $j={high_mode}$: $\lambda_{{\mathrm{{num}}}}={lam_num_high:.4f}$,"
        rf" $\lambda_{{\mathrm{{exact}}}}={lam_exact_high:.4f}$",
        fontsize=10,
    )
    ax.set_xlabel("$x$")
    ax.set_ylabel("$u(x)$")
    ax.set_xlim(-1, 1)
    ax.set_ylim(-1.15, 1.15)
    ax.grid(True, alpha=0.25, linewidth=0.4)
    ax.legend(loc="lower center", fontsize=9, framealpha=0.9)

    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "eigen_smooth_lie.pdf", bbox_inches="tight", pad_inches=0.05)
    fig.savefig(OUTPUT_DIR / "eigen_smooth_lie.png", bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)

    return {
        "N": N,
        "low_mode": low_mode,
        "high_mode": high_mode,
        "eigvals_sorted": eigvals.tolist(),
        "lam_num_low": float(lam_num_low),
        "lam_num_high": float(lam_num_high),
        "lam_exact_low": float(lam_exact_low),
        "lam_exact_high": float(lam_exact_high),
        "abs_err_low": float(abs(lam_num_low - lam_exact_low)),
        "abs_err_high": float(abs(lam_num_high - lam_exact_high)),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--N", type=int, default=16)
    parser.add_argument("--dump", type=str, default=None,
                        help="Path to JSON file capturing numerical results.")
    args = parser.parse_args()

    results = make_figure(N=args.N)

    print(f"[Etude 18.1]  N = {args.N}")
    print(f"  mode {results['low_mode']:2d}:  lambda_num = {results['lam_num_low']:.6e}"
          f"   exact = {results['lam_exact_low']:.6e}"
          f"   |err| = {results['abs_err_low']:.2e}")
    print(f"  mode {results['high_mode']:2d}:  lambda_num = {results['lam_num_high']:.6e}"
          f"   exact = {results['lam_exact_high']:.6e}"
          f"   |err| = {results['abs_err_high']:.2e}")
    print(f"  figure: {OUTPUT_DIR / 'eigen_smooth_lie.pdf'}")

    if args.dump:
        Path(args.dump).write_text(json.dumps(results, indent=2))
        print(f"  dumped:  {args.dump}")


if __name__ == "__main__":
    main()
