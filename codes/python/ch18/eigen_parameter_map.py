#!/usr/bin/env python3
"""
eigen_parameter_map.py

Chapter 18: Linear Spectral Eigenproblems
Computational Etude 18.10: map the neutral curve safely.

Two-parameter Sturm-Liouville model
        -u'' + (alpha + beta x^2) u = lambda u,   u(+-1) = 0,
over the parameter rectangle (alpha, beta) in [-2, 2] x [0, 8].
The goal is to compute the leading eigenvalue as a function of
(alpha, beta) while demonstrating the chapter's recommended workflow:

  (1) QR solve on a coarse grid of parameters (GLOBAL solver, guaranteed
      to return all eigenvalues);
  (2) inverse iteration fill-in on a finer grid (LOCAL solver, cheaper
      but seeded by the coarse result).

We also demonstrate the pathology of naive local tracking: starting a
local iteration with a sub-optimal initial shift can miss the true
lowest eigenvalue if an avoided crossing intervenes. In this model
the first and second eigenvalues do not actually cross (they are
both even/odd modes of the same symmetric potential), so we fabricate
a crossing by pairing a NAIVE ordinal-match local solver that always
expects the current lowest to stay lowest against a SAFE anchored
solver that re-runs QR at every coarse node.

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
import scipy.linalg

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR.parent / "ch07"))
from cheb_matrix import cheb_matrix  # noqa: E402

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["CMU Serif", "DejaVu Serif"],
    "mathtext.fontset": "cm",
    "font.size": 11, "axes.linewidth": 0.8,
    "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "figure.dpi": 150, "savefig.dpi": 300,
})
NAVY, SKY, CORAL, TEAL, PURPLE = "#142D6E", "#7896D2", "#E74C3C", "#16A085", "#8E44AD"
OUTPUT_DIR = SCRIPT_DIR.parent.parent.parent / "textbook" / "figures" / "ch18" / "python"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def build_operator(N: int, alpha: float, beta: float) -> np.ndarray:
    D, x = cheb_matrix(N)
    D2 = D @ D
    A = -D2[1:N, 1:N] + np.diag(alpha + beta * x[1:N] ** 2)
    return A


def qr_lowest_two(N: int, alpha: float, beta: float):
    A = build_operator(N, alpha, beta)
    lam = np.sort(np.linalg.eigvals(A).real)
    return lam[0], lam[1]


def inverse_iteration_shift(N: int, alpha: float, beta: float,
                            shift: float, max_iter: int = 20):
    A = build_operator(N, alpha, beta)
    n = A.shape[0]
    rng = np.random.default_rng(0)
    v = rng.standard_normal(n); v /= np.linalg.norm(v)
    lu, piv = scipy.linalg.lu_factor(A - shift * np.eye(n))
    lam = 0.0
    for _ in range(max_iter):
        w = scipy.linalg.lu_solve((lu, piv), v)
        v = w / np.linalg.norm(w)
        lam = v @ (A @ v)
    return lam


def make_figure(N: int = 40,
                alpha_coarse=np.linspace(-2.0, 2.0, 9),
                beta_coarse=np.linspace(0.0, 8.0, 9),
                beta_fine_k: int = 40):
    """Compute the leading eigenvalue surface on a coarse (alpha, beta) grid
    via QR, then a finer beta-only slice at fixed alpha via inverse
    iteration anchored on coarse QR."""
    Aa, Bb = np.meshgrid(alpha_coarse, beta_coarse, indexing="ij")
    lambda1 = np.zeros_like(Aa)
    lambda2 = np.zeros_like(Aa)
    for i, a in enumerate(alpha_coarse):
        for j, b in enumerate(beta_coarse):
            l1, l2 = qr_lowest_two(N, a, b)
            lambda1[i, j] = l1
            lambda2[i, j] = l2

    # Finer beta scan at fixed alpha = 0
    alpha_fixed = 0.0
    beta_fine = np.linspace(0.0, 8.0, beta_fine_k)
    # Anchor: every 5th beta is a full QR solve; others are inverse iteration
    lam1_fine_safe = np.zeros_like(beta_fine)
    lam1_fine_naive = np.zeros_like(beta_fine)
    # "safe": re-anchor with QR at each beta
    for k, b in enumerate(beta_fine):
        lam1_fine_safe[k], _ = qr_lowest_two(N, alpha_fixed, b)

    # "naive": continue from previous shift, expecting the SECOND mode tracks
    # (simulating a user who tracks the wrong branch); we use
    # shift from PREVIOUS beta, seeded from beta=0 SECOND eigenvalue
    l1_0, l2_0 = qr_lowest_two(N, alpha_fixed, beta_fine[0])
    shift = l2_0                # DELIBERATELY start on the wrong branch
    lam1_fine_naive[0] = l2_0
    for k in range(1, len(beta_fine)):
        lam = inverse_iteration_shift(N, alpha_fixed, beta_fine[k], shift)
        lam1_fine_naive[k] = lam
        shift = lam             # continue from previous converged value

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.5))

    ax = axes[0]
    c = ax.contourf(Aa, Bb, lambda1, levels=18, cmap="viridis")
    cs = ax.contour(Aa, Bb, lambda1, levels=8, colors="white", linewidths=0.5)
    ax.clabel(cs, inline=True, fontsize=7)
    fig.colorbar(c, ax=ax, label=r"$\lambda_1$")
    ax.scatter(Aa.ravel(), Bb.ravel(), s=14, c="white",
               edgecolors="black", linewidths=0.6)
    ax.set_xlabel(r"$\alpha$")
    ax.set_ylabel(r"$\beta$")
    ax.set_title(r"coarse QR grid ($9 \times 9$): lowest eigenvalue",
                 fontsize=10)

    ax = axes[1]
    ax.plot(beta_fine, lam1_fine_safe, "o-", color=NAVY, markersize=4,
            markerfacecolor="white", markeredgewidth=1.0, linewidth=0.8,
            label="SAFE: QR re-anchor each step")
    ax.plot(beta_fine, lam1_fine_naive, "x-", color=CORAL, markersize=5,
            linewidth=0.8,
            label="NAIVE: inverse iter seeded from wrong branch")
    ax.set_xlabel(r"$\beta$  (with $\alpha = 0$)")
    ax.set_ylabel("eigenvalue")
    ax.set_title("eigenvalue tracking: safe vs. naive local iteration",
                 fontsize=10)
    ax.grid(True, alpha=0.25, linewidth=0.4)
    ax.legend(loc="upper left", fontsize=9)

    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "eigen_parameter_map.pdf", bbox_inches="tight", pad_inches=0.05)
    fig.savefig(OUTPUT_DIR / "eigen_parameter_map.png", bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)

    return {
        "N": N,
        "alpha_grid": alpha_coarse.tolist(),
        "beta_grid":  beta_coarse.tolist(),
        "lambda1_coarse": lambda1.tolist(),
        "beta_fine":       beta_fine.tolist(),
        "lam1_fine_safe":  lam1_fine_safe.tolist(),
        "lam1_fine_naive": lam1_fine_naive.tolist(),
    }


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--dump", type=str, default=None)
    args = p.parse_args()

    r = make_figure()
    print(f"[Etude 18.10]  parameter map for -u'' + (alpha + beta x^2) u = lam u")
    print(f"  coarse grid: {len(r['alpha_grid'])} x {len(r['beta_grid'])}  = {len(r['alpha_grid'])*len(r['beta_grid'])} QR solves")
    print(f"  safe tracker: starts at {r['lam1_fine_safe'][0]:.4f}, ends at {r['lam1_fine_safe'][-1]:.4f}")
    print(f"  naive tracker: starts at {r['lam1_fine_naive'][0]:.4f}, ends at {r['lam1_fine_naive'][-1]:.4f}")
    print(f"  figure: {OUTPUT_DIR / 'eigen_parameter_map.pdf'}")
    if args.dump:
        Path(args.dump).write_text(json.dumps(r, indent=2))
        print(f"  dumped: {args.dump}")


if __name__ == "__main__":
    main()
