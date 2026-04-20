#!/usr/bin/env python3
"""
eigen_power_inverse.py

Chapter 18: Linear Spectral Eigenproblems
Computational Etude 18.9: one mode at a time.

Applies the power method and inverse iteration (with shift) to the
interior-block Dirichlet Laplacian -D^2_int at N = 32. The target
spectrum is lambda_j = (j pi / 2)^2 for j = 1, ..., 31.

Three panels:
  (A) Power method convergence to the LARGEST eigenvalue. The
      iteration converges GEOMETRICALLY at rate |lambda_{N-1}/lambda_N|.
  (B) Inverse iteration with three different shifts recovers three
      different interior eigenvalues. Each shift locks onto the
      eigenvalue closest to it.
  (C) Cautionary case: a deliberately-chosen shift exactly between two
      eigenvalues makes inverse iteration indecisive (slow convergence).

Pedagogical takeaway: the power method is cheap per step but finds
only the extreme eigenvalue; inverse iteration is more expensive per
step (requires a linear solve with (A - mu I)) but selects the mode
closest to the shift. Both can MISS eigenvalues if used carelessly.

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


def build_matrix(N: int = 32):
    D, _ = cheb_matrix(N)
    return -(D @ D)[1:N, 1:N]


def power_method(A, max_iter: int = 80, seed: int = 0):
    rng = np.random.default_rng(seed)
    n = A.shape[0]
    v = rng.standard_normal(n); v /= np.linalg.norm(v)
    history = []
    for k in range(max_iter):
        w = A @ v
        lam = v @ w
        history.append(lam)
        v = w / np.linalg.norm(w)
    return np.asarray(history), v


def inverse_iteration(A, shift: float, max_iter: int = 40, seed: int = 1):
    rng = np.random.default_rng(seed)
    n = A.shape[0]
    v = rng.standard_normal(n); v /= np.linalg.norm(v)
    history = []
    # LU factorisation once, re-used each iteration
    lu, piv = scipy.linalg.lu_factor(A - shift * np.eye(n))
    for k in range(max_iter):
        w = scipy.linalg.lu_solve((lu, piv), v)
        # Rayleigh quotient for A (not the shifted system)
        lam = v @ (A @ v)
        history.append(lam)
        v = w / np.linalg.norm(w)
    # final Rayleigh quotient
    final = v @ (A @ v)
    history.append(final)
    return np.asarray(history), v


def make_figure(N: int = 32):
    A = build_matrix(N)
    # reference spectrum
    lam_ref = np.sort(np.linalg.eigvals(A).real)

    # Panel A: power method to largest
    hist_pow, _ = power_method(A, max_iter=80)
    lam_max = lam_ref[-1]

    # Panel B: inverse iteration with three different shifts
    shifts = [5.0, 90.0, 250.0]
    histories = []
    targets = []
    for mu in shifts:
        hist, _ = inverse_iteration(A, shift=mu, max_iter=30)
        # identify which eigenvalue this converged to (nearest to final)
        final = hist[-1]
        target = lam_ref[np.argmin(np.abs(lam_ref - final))]
        histories.append(hist)
        targets.append(target)

    # Panel C: shift exactly between two eigenvalues
    # pick two adjacent eigenvalues and shift to their midpoint
    lam_a, lam_b = lam_ref[3], lam_ref[4]
    mu_bad = 0.5 * (lam_a + lam_b)
    hist_bad, _ = inverse_iteration(A, shift=mu_bad, max_iter=30)

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2))

    ax = axes[0]
    ax.semilogy(np.abs(hist_pow - lam_max), "o-", color=NAVY,
                markersize=4, markerfacecolor="white", markeredgewidth=0.9,
                linewidth=0.8)
    ax.axhline(1e-14, color="k", linestyle=":", linewidth=0.5, alpha=0.4)
    ax.set_xlabel("iteration $k$")
    ax.set_ylabel(r"$|\mu^{(k)} - \lambda_{\max}|$")
    ax.set_title(f"power method  →  $\\lambda_{{\\max}} = {lam_max:.3f}$", fontsize=10)
    ax.grid(True, which="both", alpha=0.25, linewidth=0.4)
    ax.set_ylim(1e-18, 1e4)

    ax = axes[1]
    colors_B = [NAVY, CORAL, TEAL]
    for mu, hist, target, c in zip(shifts, histories, targets, colors_B):
        err = np.abs(hist - target)
        ax.semilogy(err, "o-", color=c, markersize=4,
                    markerfacecolor="white", markeredgewidth=0.9,
                    linewidth=0.8,
                    label=fr"shift $\mu = {mu}$  →  $\lambda = {target:.3f}$")
    ax.axhline(1e-14, color="k", linestyle=":", linewidth=0.5, alpha=0.4)
    ax.set_xlabel("iteration $k$")
    ax.set_ylabel(r"$|\mu^{(k)} - \lambda_{\rm target}|$")
    ax.set_title("inverse iteration:  three shifts, three modes", fontsize=10)
    ax.grid(True, which="both", alpha=0.25, linewidth=0.4)
    ax.set_ylim(1e-18, 1e4)
    ax.legend(loc="upper right", fontsize=8)

    ax = axes[2]
    # distance to both candidates
    err_a = np.abs(hist_bad - lam_a)
    err_b = np.abs(hist_bad - lam_b)
    ax.semilogy(err_a, "o-", color=NAVY, markersize=4,
                markerfacecolor="white", markeredgewidth=0.9, linewidth=0.8,
                label=fr"dist to $\lambda = {lam_a:.3f}$")
    ax.semilogy(err_b, "s-", color=CORAL, markersize=4,
                markerfacecolor="white", markeredgewidth=0.9, linewidth=0.8,
                label=fr"dist to $\lambda = {lam_b:.3f}$")
    ax.set_xlabel("iteration $k$")
    ax.set_ylabel(r"$|\mu^{(k)} - \lambda|$")
    ax.set_title(f"cautionary: shift $\\mu = {mu_bad:.3f}$ between two modes", fontsize=10)
    ax.grid(True, which="both", alpha=0.25, linewidth=0.4)
    ax.legend(loc="upper right", fontsize=8)

    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "eigen_power_inverse.pdf", bbox_inches="tight", pad_inches=0.05)
    fig.savefig(OUTPUT_DIR / "eigen_power_inverse.png", bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)

    return {
        "N": N,
        "lam_max":   float(lam_max),
        "final_power_method": float(hist_pow[-1]),
        "shifts":    shifts,
        "targets":   [float(t) for t in targets],
        "bad_shift": float(mu_bad),
        "bad_pair":  [float(lam_a), float(lam_b)],
        "final_bad": float(hist_bad[-1]),
    }


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--dump", type=str, default=None)
    args = p.parse_args()

    r = make_figure()
    print(f"[Etude 18.9]  power and inverse iteration, N = {r['N']}")
    print(f"  power method final -> {r['final_power_method']:.6f} (exact lam_max = {r['lam_max']:.6f})")
    for mu, t in zip(r["shifts"], r["targets"]):
        print(f"  inverse iter, shift {mu}  ->  lambda = {t:.6f}")
    print(f"  cautionary: shift={r['bad_shift']:.3f} between {r['bad_pair'][0]:.3f} and {r['bad_pair'][1]:.3f}")
    print(f"  figure: {OUTPUT_DIR / 'eigen_power_inverse.pdf'}")
    if args.dump:
        Path(args.dump).write_text(json.dumps(r, indent=2))
        print(f"  dumped: {args.dump}")


if __name__ == "__main__":
    main()
