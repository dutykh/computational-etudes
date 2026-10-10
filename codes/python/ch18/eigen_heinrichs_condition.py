#!/usr/bin/env python3
"""
eigen_heinrichs_condition.py

Chapter 18: Linear Spectral Eigenproblems
Computational Etude 18.8: condition-number surgery.

For the fourth-order clamped eigenproblem
        u'''' = lambda u,   u(+-1) = u'(+-1) = 0,
we compare two discretisations:

  (naive) boundary-bordered D^4 pencil (as in Ch 14 Etude 14.2)
  (Heinrichs) u = (1 - x^2)^2 q substitution; solve the generalised
              problem A q = lambda S2 q directly

and record:
  (a) condition number of the matrices vs N,
  (b) drift-with-N of selected low eigenvalues.

Exact eigenvalues satisfy cos(2 beta) cosh(2 beta) = 1 with
lambda = beta^4. Lowest four: lambda_1 approx 500.564,
lambda_2 approx 3803.537, lambda_3 approx 14617.6, lambda_4 approx 39943.8.

The empirical finding matches Boyd (2000) Fig 7.9: the Heinrichs-basis
matrix has condition number O(N^4), whereas the naive D^4 boundary-
bordered pencil has condition number O(N^8).

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
sys.path.insert(0, str(SCRIPT_DIR))
from heinrichs_basis import (
    heinrichs_clamped_matrix, naive_clamped_operator
)  # noqa: E402

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["CMU Serif", "DejaVu Serif"],
    "mathtext.fontset": "cm",
    "font.size": 11, "axes.linewidth": 0.8,
    "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "figure.dpi": 150, "savefig.dpi": 300,
})
NAVY, CORAL, TEAL = "#142D6E", "#E74C3C", "#16A085"
OUTPUT_DIR = SCRIPT_DIR.parent.parent.parent / "textbook" / "figures" / "ch18" / "python"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def exact_clamped_beta(n_modes: int = 6) -> np.ndarray:
    """Return the first n_modes positive roots beta of cos(2b)cosh(2b) = 1."""
    def f(b): return np.cos(2 * b) * np.cosh(2 * b) - 1.0
    # asymptotic: beta_j -> (2j+1) pi / 4
    # bracket search around these asymptotes
    betas = []
    for j in range(1, n_modes + 1):
        b_asym = (2 * j + 1) * np.pi / 4.0
        lo, hi = b_asym - 0.5, b_asym + 0.5
        flo, fhi = f(lo), f(hi)
        if flo * fhi < 0:
            # bisection
            for _ in range(200):
                mid = 0.5 * (lo + hi)
                fmid = f(mid)
                if flo * fmid < 0:
                    hi = mid; fhi = fmid
                else:
                    lo = mid; flo = fmid
                if hi - lo < 1e-14:
                    break
            betas.append(0.5 * (lo + hi))
    return np.asarray(betas)


def cond_and_spectrum_naive(N: int):
    """Condition of the bordered naive pencil and the smallest four finite
    real eigenvalues."""
    A, B = naive_clamped_operator(N)
    lam = scipy.linalg.eigvals(A, B)
    lam = lam[np.isfinite(lam)]
    lam = np.sort(lam.real[np.abs(lam.imag) < 1e-6])
    lam = lam[lam > 1.0]   # drop near-zero eigenvalues of B-singular rows
    # condition number via generalised eigvals doesn't have a standard
    # single number; use the condition of A (ignoring infinite eigenvalues)
    kappa = np.linalg.cond(A)
    return kappa, lam[:4]


def cond_and_spectrum_heinrichs(N: int):
    A, M, _ = heinrichs_clamped_matrix(N)
    lam = scipy.linalg.eigvals(A, M)
    lam = lam[np.isfinite(lam)]
    lam = np.sort(lam.real[np.abs(lam.imag) < 1e-6])
    lam = lam[lam > 1.0]
    # standardised form for condition:
    A_std = np.linalg.solve(M, A)
    kappa = np.linalg.cond(A_std)
    return kappa, lam[:4]


def make_figure(Ns=(12, 16, 24, 32, 48, 64, 96)):
    betas = exact_clamped_beta(6)
    lam_exact = betas ** 4

    kappa_naive = []; kappa_hein = []
    drift_naive = []; drift_hein = []
    mode_track = 0   # track the lowest eigenvalue

    # reference at high N to define "drift"
    _, lam_ref = cond_and_spectrum_heinrichs(max(Ns))
    lam1_ref = lam_ref[mode_track]

    for N in Ns:
        kn, lam_n = cond_and_spectrum_naive(N)
        kh, lam_h = cond_and_spectrum_heinrichs(N)
        kappa_naive.append(kn); kappa_hein.append(kh)
        drift_naive.append(abs(lam_n[mode_track] - lam_exact[mode_track]) if lam_n.size > mode_track else np.nan)
        drift_hein.append( abs(lam_h[mode_track] - lam_exact[mode_track]) if lam_h.size > mode_track else np.nan)

    Ns = np.asarray(Ns, dtype=float)
    kappa_naive = np.asarray(kappa_naive)
    kappa_hein  = np.asarray(kappa_hein)
    drift_naive = np.asarray(drift_naive)
    drift_hein  = np.asarray(drift_hein)

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.5))

    ax = axes[0]
    ax.loglog(Ns, kappa_naive, "o-", color=NAVY, markersize=6,
              markerfacecolor="white", markeredgewidth=1.1, linewidth=0.8,
              label="naive ($D^4$ bordered)")
    ax.loglog(Ns, kappa_hein, "s-", color=CORAL, markersize=5,
              markerfacecolor="white", markeredgewidth=1.1, linewidth=0.8,
              label="Heinrichs $(1-x^2)^2 T_j$")
    # reference slopes
    Nref = Ns
    ax.loglog(Nref, kappa_naive[0] * (Nref / Nref[0]) ** 8,
              "--", color=NAVY, linewidth=0.6, alpha=0.7,
              label=r"$N^8$")
    ax.loglog(Nref, kappa_hein[0] * (Nref / Nref[0]) ** 4,
              "--", color=CORAL, linewidth=0.6, alpha=0.7,
              label=r"$N^4$")
    ax.set_xlabel("$N$")
    ax.set_ylabel(r"$\kappa$ (condition number)")
    ax.set_title("conditioning of the fourth-derivative matrix", fontsize=10)
    ax.grid(True, which="both", alpha=0.25, linewidth=0.4)
    ax.legend(loc="upper left", fontsize=9)

    ax = axes[1]
    ax.loglog(Ns, np.maximum(drift_naive, 1e-17), "o-", color=NAVY,
              markersize=6, markerfacecolor="white", markeredgewidth=1.1,
              linewidth=0.8, label=fr"naive: $|\lambda_1 - {lam_exact[0]:.4f}|$")
    ax.loglog(Ns, np.maximum(drift_hein, 1e-17), "s-", color=CORAL,
              markersize=5, markerfacecolor="white", markeredgewidth=1.1,
              linewidth=0.8, label=fr"Heinrichs: $|\lambda_1 - {lam_exact[0]:.4f}|$")
    ax.set_xlabel("$N$")
    ax.set_ylabel("error in first eigenvalue")
    ax.set_title("first-eigenvalue error vs. $N$", fontsize=10)
    ax.grid(True, which="both", alpha=0.25, linewidth=0.4)
    ax.legend(loc="upper right", fontsize=9)

    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "eigen_heinrichs_condition.pdf", bbox_inches="tight", pad_inches=0.05)
    fig.savefig(OUTPUT_DIR / "eigen_heinrichs_condition.png", bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)

    return {
        "Ns":           Ns.tolist(),
        "kappa_naive":  kappa_naive.tolist(),
        "kappa_heinrichs": kappa_hein.tolist(),
        "drift_naive":  drift_naive.tolist(),
        "drift_heinrichs": drift_hein.tolist(),
        "lam_exact_first4": lam_exact[:4].tolist(),
    }


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--dump", type=str, default=None)
    args = p.parse_args()

    r = make_figure()
    print(f"[Etude 18.8]  Heinrichs condition-number surgery")
    print(f"  exact lambda_1..4 = {[round(v, 3) for v in r['lam_exact_first4']]}")
    print(f"  N, kappa(naive), kappa(Heinrichs):")
    for N, kn, kh in zip(r["Ns"], r["kappa_naive"], r["kappa_heinrichs"]):
        print(f"    N={int(N):3d}   kappa_naive={kn:.2e}   kappa_Heinrichs={kh:.2e}")
    print(f"  figure: {OUTPUT_DIR / 'eigen_heinrichs_condition.pdf'}")
    if args.dump:
        Path(args.dump).write_text(json.dumps(r, indent=2))
        print(f"  dumped: {args.dump}")


if __name__ == "__main__":
    main()
