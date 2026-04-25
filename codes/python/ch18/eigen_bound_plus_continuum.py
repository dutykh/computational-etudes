#!/usr/bin/env python3
"""
eigen_bound_plus_continuum.py

Chapter 18: Linear Spectral Eigenproblems
Computational Etude 18.6: bound states versus continuum.

Solves the Schrödinger equation
        -u''(x) + V(x) u(x) = E u(x),   |u| -> 0 as |x| -> inf
with the Pöschl-Teller potential
        V(x) = -nu(nu + 1) sech^2(x).
For integer nu, there are exactly floor(nu) + 1 bound states at
        E_j = -(nu - j)^2,   j = 0, 1, ..., floor(nu),
plus a continuous spectrum at E >= 0.

The etude uses the rational_chebyshev utility on the real line at two
resolutions (N = 60 and N = 96), then applies the spectrum_verify
drift diagnostic. The discrete bound states should be stable under
refinement; the continuum modes should wander.

Reproduces the pedagogical lesson of Boyd (2000) Fig 7.7: for an
operator with a finite number of discrete modes plus continuum, the
drift test sharply isolates the trusted subset.

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
from spectrum_verify import verify_spectrum  # noqa: E402

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["CMU Serif", "DejaVu Serif"],
    "mathtext.fontset": "cm",
    "font.size": 11, "axes.linewidth": 0.8,
    "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "figure.dpi": 150, "savefig.dpi": 300,
})
NAVY, SKY, CORAL, TEAL = "#142D6E", "#7896D2", "#E74C3C", "#16A085"
OUTPUT_DIR = SCRIPT_DIR.parent.parent.parent / "textbook" / "figures" / "ch18" / "python"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def solve_poschl_teller(N: int, nu: float, ell: float = 6.0):
    """Return the sorted spectrum of H = -d^2/dx^2 - nu(nu+1) sech^2(x)."""
    _, D2, x = rational_chebyshev_derivative_matrices(N, ell)
    V = -nu * (nu + 1.0) / np.cosh(x) ** 2
    H = -D2 + np.diag(V)
    lam = np.linalg.eigvals(H)
    lam = lam.real[np.abs(lam.imag) < 1e-6]
    return np.sort(lam)


def make_figure(N1: int = 60, N2: int = 96, nu: float = 4.0, tol: float = 1e-3):
    lam1 = solve_poschl_teller(N1, nu)
    lam2 = solve_poschl_teller(N2, nu)

    # expected bound states for integer nu = 4: E = -16, -9, -4, -1 plus E_0 = 0 (threshold)
    expected_bound = -np.arange(int(nu), -1, -1) ** 2  # [-16, -9, -4, -1, 0]
    # Drop the zero-energy threshold (marginal bound state); keep strictly negative
    expected_bound = expected_bound[expected_bound < -0.5]

    report = verify_spectrum(lam1, lam2, tol=tol)

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.5))

    # Panel A: spectra at two resolutions (scatter of low modes)
    ax = axes[0]
    n_show = 25
    j = np.arange(1, n_show + 1)
    ax.plot(j, lam1[:n_show], "o", color=NAVY, markersize=6,
            markerfacecolor="white", markeredgewidth=1.1,
            label=f"$N_1 = {N1}$")
    ax.plot(j, lam2[:n_show], "x", color=CORAL, markersize=7,
            markeredgewidth=1.1,
            label=f"$N_2 = {N2}$")
    for E in expected_bound:
        ax.axhline(E, color=TEAL, linestyle="--", linewidth=0.6, alpha=0.6)
    ax.axhline(0.0, color="k", linestyle="-", linewidth=0.6, alpha=0.4)
    ax.set_xlabel("mode number $j$")
    ax.set_ylabel("eigenvalue $E$")
    ax.set_title(
        rf"Pöschl-Teller $V = -{int(nu)}\cdot{int(nu)+1}\,\mathrm{{sech}}^2 x$"
        rf" (exact: {expected_bound.size} bound states)",
        fontsize=10)
    ax.grid(True, alpha=0.25, linewidth=0.4)
    ax.set_ylim(min(-20, lam1[:n_show].min() - 1), max(30, lam1[:n_show].max() + 2))
    ax.legend(loc="upper left", fontsize=9)

    # Panel B: reciprocal drift ratio (trusted at top)
    ax = axes[1]
    j = np.arange(1, report.lam1.size + 1)
    inv_ord = 1.0 / np.maximum(report.delta_ordinal, 1e-18)
    inv_nst = 1.0 / np.maximum(report.delta_nearest, 1e-18)
    ax.semilogy(j, inv_ord, "o", color=NAVY, markersize=4,
                markerfacecolor="white", markeredgewidth=1.0,
                label=r"$1/\delta_\mathrm{ordinal}$")
    ax.semilogy(j, inv_nst, "x", color=CORAL, markersize=5,
                markeredgewidth=1.0,
                label=r"$1/\delta_\mathrm{nearest}$")
    ax.axhline(1.0 / tol, color="k", linestyle="--", linewidth=0.6, alpha=0.5)
    ax.set_xlim(0, 30)
    ax.set_xlabel("mode number $j$")
    ax.set_ylabel(r"$1/\delta_j$  (trusted = high)")
    ax.set_title(f"drift diagnostic:  trusted = {report.n_trusted}", fontsize=10)
    ax.grid(True, which="both", alpha=0.25, linewidth=0.4)
    ax.legend(loc="upper right", fontsize=9)

    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "eigen_bound_plus_continuum.pdf", bbox_inches="tight", pad_inches=0.05)
    fig.savefig(OUTPUT_DIR / "eigen_bound_plus_continuum.png", bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)

    return {
        "N1": N1, "N2": N2, "nu": nu, "tol": tol,
        "expected_bound":  expected_bound.tolist(),
        "n_trusted":       int(report.n_trusted),
        "lam1_low":        lam1[:25].tolist(),
        "lam2_low":        lam2[:25].tolist(),
        "trusted_first_k": report.trusted[:10].astype(int).tolist(),
    }


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--N1",  type=int,   default=60)
    p.add_argument("--N2",  type=int,   default=96)
    p.add_argument("--nu",  type=float, default=4.0)
    p.add_argument("--tol", type=float, default=1e-3)
    p.add_argument("--dump", type=str,  default=None)
    args = p.parse_args()

    r = make_figure(N1=args.N1, N2=args.N2, nu=args.nu, tol=args.tol)
    print(f"[Etude 18.6]  Pöschl-Teller with nu = {args.nu}")
    print(f"  expected bound states: {r['expected_bound']}")
    print(f"  trusted (drift < {args.tol}): {r['n_trusted']} of {len(r['lam1_low'])*'?'}")
    print(f"  lowest {len(r['lam1_low'][:6])} at N1: {[round(v,4) for v in r['lam1_low'][:6]]}")
    print(f"  lowest {len(r['lam2_low'][:6])} at N2: {[round(v,4) for v in r['lam2_low'][:6]]}")
    print(f"  figure: {OUTPUT_DIR / 'eigen_bound_plus_continuum.pdf'}")
    if args.dump:
        Path(args.dump).write_text(json.dumps(r, indent=2))
        print(f"  dumped: {args.dump}")


if __name__ == "__main__":
    main()
