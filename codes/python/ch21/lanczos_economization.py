#!/usr/bin/env python3
"""
lanczos_economization.py

Chapter 21: Special Tricks for the Spectral Researcher
Computational Etude 21.6: Lanczos economization for an expensive determinant.

Suppose D(lambda) = det(A(lambda))  is the determinant of a parameter-
dependent matrix.  In real applications -- nonlinear eigenvalue problems
arising from boundary integral methods, transfer-matrix reformulations
of high-order BVPs, or fluid-dynamic stability calculations -- A is
typically several hundred to a few thousand on a side, and each LU
factorisation is O(M^3) in cost.  The toy problem below uses only
M = 20 to keep the demonstration fast, but the *number* of LU calls
saved is the figure of merit, not the wall-clock time.

Strategy A (naive): scan a fine grid of K_dense lambda values for sign
changes of D, then refine each bracket by 30 bisection steps.  Cost:
K_dense + 30 * (number of roots) LU calls.

Strategy B (Lanczos economization): sample D at K = 17 Chebyshev
extrema points on [a, b], build the degree-(K-1) Chebyshev interpolant,
and find its roots via the colleague-matrix eigenvalue problem (free).
Cost: K = 17 LU calls.  No bracketing, no Newton, no iteration.

Setup.  We pick a deliberately expensive A(lambda):

        A(lambda) = T_{40}  +  diag(lambda - mu_j),

where T_{40} is the 40-by-40 Trefethen second-derivative differentiation
matrix on Chebyshev points and mu_j = 1, 4, 9, ..., 400 are 'shifts'
arranged so that D(lambda) = det A vanishes at each mu_j.  The roots
of D are therefore the unperturbed values mu_j -- on the chosen
interval [0.5, 30] there are five of them (at 1, 4, 9, 16, 25).

Compare:
  (A) plain root-bracketing on a fine sampling: O(K) LUs to locate
      sign changes, then a few extra to refine each root by bisection.
  (B) Lanczos economization: K = 17 LUs, then near-instantaneous
      polynomial root-finding.

Both should recover the same five roots; (B) does it with fewer LUs
and cleaner arithmetic.

Author: Dr. Denys Dutykh (Khalifa University of Science and Technology)
Part of "Computational Etudes: A Spectral Approach"
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
from tricks_common import setup_matplotlib, output_dir_for, save_fig, dump_json
from tricks_common import NAVY, CORAL, TEAL


M = 20
SHIFTS = np.arange(1, M + 1).astype(float) ** 2     # 1, 4, 9, ..., M^2


def _make_T2():
    """Build a small (M x M) symmetric perturbation matrix once."""
    rng = np.random.default_rng(0)
    Q, _ = np.linalg.qr(rng.standard_normal((M, M)))
    diag = 0.05 * rng.standard_normal(M)
    return Q @ np.diag(diag) @ Q.T


_T2 = _make_T2()


def expensive_D(lam):
    """Determinant of A(lam) = T2 + diag(lam - mu_j) (truncated to M shifts).
    The 'expensive' work is the LU factorisation inside slogdet for an
    M-by-M dense matrix.  At M = 20 the wall time is microseconds, but
    the etude's *cost ratio* (LUs naive vs LUs Lanczos) is what matters,
    not absolute wall time."""
    A = _T2 + np.diag(lam - SHIFTS[:M])
    sign, logabs = np.linalg.slogdet(A)
    return sign * np.exp(logabs)


def roots_chebyshev_surrogate(K, a, b):
    """Sample D at K = K Chebyshev-extrema points, build a Chebyshev
    interpolant, and find its roots via the colleague matrix.  Returns
    (roots in [a, b], list of (lam_j, D(lam_j)) sample pairs)."""
    j = np.arange(K)
    t = np.cos(np.pi * j / (K - 1))                  # Chebyshev-Lobatto on [-1,1]
    lam = 0.5 * (a + b) + 0.5 * (b - a) * t          # mapped to [a, b]
    D_samples = np.array([expensive_D(L) for L in lam])

    # DCT-I to obtain Chebyshev coefficients on the (K)-pt extrema grid.
    extended = np.concatenate([D_samples, D_samples[K - 2:0:-1]])
    A = np.real(np.fft.fft(extended)) / (K - 1)
    A[0] *= 0.5
    A[K - 1] *= 0.5
    coeffs = A[:K]

    # Build the colleague (companion-like) matrix for Chebyshev coefficients.
    n = K - 1
    C = np.zeros((n, n))
    for i in range(n - 1):
        C[i, i + 1] = 0.5
        C[i + 1, i] = 0.5
    C[0, 1] = 1.0       # special row in colleague matrix
    last = -coeffs[:n] / (2.0 * coeffs[n])
    last[n - 2] += 0.5
    C[n - 1, :] = last

    eigs = np.linalg.eigvals(C)
    real_in_window = np.array([float(np.real(z))
                               for z in eigs
                               if abs(z.imag) < 1e-6 and -1.0 <= z.real <= 1.0])
    real_in_window.sort()
    roots = 0.5 * (a + b) + 0.5 * (b - a) * real_in_window
    return roots, lam, D_samples


def roots_naive_brackets(K_dense, a, b):
    """Sample D at K_dense uniformly-spaced points, find sign changes,
    then refine each by 30 bisection steps."""
    lam = np.linspace(a, b, K_dense)
    D = np.array([expensive_D(L) for L in lam])
    roots = []
    for k in range(len(lam) - 1):
        if D[k] * D[k + 1] < 0:
            l, r = lam[k], lam[k + 1]
            fl = D[k]
            for _ in range(30):
                mid = 0.5 * (l + r)
                fm = expensive_D(mid)
                if fl * fm <= 0:
                    r = mid
                else:
                    l, r, fl = mid, r, fm
            roots.append(0.5 * (l + r))
    return np.array(roots), lam, D


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    p.add_argument("--dump", type=str, default=None)
    args = p.parse_args()

    setup_matplotlib()
    out = output_dir_for(SCRIPT_DIR)

    a, b = 0.5, 30.0

    # --- Strategy A: naive bracketing ---------------------------------
    K_dense = 60
    t0 = time.perf_counter()
    roots_A, lam_A, D_A = roots_naive_brackets(K_dense, a, b)
    t_A = time.perf_counter() - t0
    n_LUs_A = K_dense + 30 * len(roots_A)

    # --- Strategy B: Lanczos economization ----------------------------
    K = 17
    t0 = time.perf_counter()
    roots_B, lam_B, D_B = roots_chebyshev_surrogate(K, a, b)
    t_B = time.perf_counter() - t0
    n_LUs_B = K

    # The 'exact' perturbed roots: the unperturbed mu_j shifted by T2, found
    # by 60 Newton steps from each unperturbed mu_j as initial guess.
    seeds = SHIFTS[(SHIFTS >= a) & (SHIFTS <= b)]
    exact_roots = []
    for seed in seeds:
        x = float(seed)
        for _ in range(80):
            d  = expensive_D(x)
            dp = (expensive_D(x + 1e-7) - expensive_D(x - 1e-7)) / 2e-7
            if abs(dp) < 1e-30: break
            x -= d / dp
        exact_roots.append(x)
    exact_roots = np.array(sorted(exact_roots))

    err_A = np.array([min(np.abs(r - exact_roots)) for r in sorted(roots_A)])
    err_B = np.array([min(np.abs(r - exact_roots)) for r in sorted(roots_B)])

    # ---- Figure -----------------------------------------------------
    fig = plt.figure(figsize=(11.5, 4.4))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.05, 1.0], wspace=0.30)

    # Panel A: D(lambda) and detected roots
    ax = fig.add_subplot(gs[0, 0])
    lam_dense = np.linspace(a, b, 600)
    D_dense = np.array([expensive_D(L) for L in lam_dense])
    ax.plot(lam_dense, D_dense, color=NAVY, linewidth=1.0,
            label=r"$D(\lambda) = \det A(\lambda)$")
    ax.axhline(0.0, color="gray", linewidth=0.4, alpha=0.5)
    ax.scatter(roots_A, np.zeros_like(roots_A), s=80, marker="x",
               color=CORAL, linewidth=1.4,
               label=f"naive: {len(roots_A)} roots ({n_LUs_A} LUs)")
    ax.scatter(roots_B, np.zeros_like(roots_B) + 0.05 * np.max(np.abs(D_dense)),
               s=80, marker="o", facecolors="none", edgecolors=TEAL,
               linewidth=1.4,
               label=f"Lanczos: {len(roots_B)} roots ({n_LUs_B} LUs)")
    ax.scatter(lam_B, D_B, s=20, color=TEAL,
               label=f"Lanczos sample points ($K = {K}$)")
    ax.set_xlabel(r"$\lambda$")
    ax.set_ylabel(r"$D(\lambda)$")
    ax.set_title(r"the determinant and the two strategies' detected roots",
                 fontsize=10)
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(True, alpha=0.25, linewidth=0.4)
    ax.set_xlim(a - 1, b + 1)

    # Panel B: cost vs accuracy
    ax = fig.add_subplot(gs[0, 1])
    width = 0.35
    methods = ["naive\nbracketing\n+ 30 bisects", f"Lanczos\nsurrogate\n($K={K}$)"]
    counts = [n_LUs_A, n_LUs_B]
    accs   = [float(np.max(err_A)), float(np.max(err_B))]
    ax.bar([0, 1], counts, width=width, color=[CORAL, TEAL],
           edgecolor="black", linewidth=0.6)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(methods, fontsize=9)
    ax.set_ylabel(r"\# expensive LU calls", fontsize=10)
    ax.set_title(r"cost: 17 vs %d LU factorisations" % n_LUs_A,
                 fontsize=10)
    for i, c in enumerate(counts):
        ax.text(i, c + 5, f"{c} LUs", ha="center", fontsize=9)
        ax.text(i, c * 0.5, f"max root err\n{accs[i]:.1e}",
                ha="center", fontsize=8.5, color="white",
                fontweight="bold")
    ax.set_ylim(0, max(counts) * 1.15)

    save_fig(fig, out, "lanczos_economization")
    plt.close(fig)

    print(f"[Etude 21.6]  Lanczos economization for an expensive determinant")
    print(f"  exact roots in [{a}, {b}]: {exact_roots}")
    print(f"  naive  : {len(roots_A)} roots, {n_LUs_A} LUs, max err = {np.max(err_A):.3e}, "
          f"wall-clock = {t_A * 1000:.1f} ms")
    print(f"  Lanczos: {len(roots_B)} roots, {n_LUs_B} LUs, max err = {np.max(err_B):.3e}, "
          f"wall-clock = {t_B * 1000:.1f} ms")
    print(f"  figure: {out / 'lanczos_economization.pdf'}")

    dump_json(args.dump, {
        "exact_roots": exact_roots.tolist(),
        "roots_naive": roots_A.tolist(),
        "roots_lanczos": roots_B.tolist(),
        "n_LUs_naive": int(n_LUs_A),
        "n_LUs_lanczos": int(n_LUs_B),
        "max_err_naive": float(np.max(err_A)),
        "max_err_lanczos": float(np.max(err_B)),
    })


if __name__ == "__main__":
    main()
