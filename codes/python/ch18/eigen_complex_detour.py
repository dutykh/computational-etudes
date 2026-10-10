#!/usr/bin/env python3
"""
eigen_complex_detour.py

Chapter 18: Linear Spectral Eigenproblems
Computational Etude 18.11: detour around the singularity.

Singular Sturm-Liouville problem (Boyd 2000, Eq. 7.56)
        u''(y) + [1/y - lambda] u(y) = 0,   u(-1) = u(+1) = 0,
with the singularity at y = 0 on the interior of the integration
interval. The solution and its derivatives are singular at y = 0,
and eigenvalues are complex-valued.

The parabolic complex-plane mapping
        y(x) = x + i Delta (x^2 - 1),   Delta > 0,
moves the integration path OFF the real axis into the lower half
y-plane, provided the branch cut of log(y) is in the UPPER half-
plane (as required by physical causality for the associated
fluid-dynamics problem). The transformed problem is solved by
standard Chebyshev collocation in x; the computed eigenvalues are
complex but well-defined and well-resolved.

We scan Delta in {0.25, 0.50, 0.75} at N = 20 and 40, and reproduce
Boyd Table 7.2 (first six eigenvalues at N = 40, Delta = 0.5).

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
    "font.size": 11, "axes.linewidth": 0.8,
    "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "figure.dpi": 150, "savefig.dpi": 300,
})
NAVY, SKY, CORAL, TEAL, PURPLE = "#142D6E", "#7896D2", "#E74C3C", "#16A085", "#8E44AD"
OUTPUT_DIR = SCRIPT_DIR.parent.parent.parent / "textbook" / "figures" / "ch18" / "python"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def solve_detoured(N: int, Delta: float):
    """Return the sorted-by-real-part eigenvalues of the detoured
    singular SL problem."""
    D, x_real = cheb_matrix(N)
    # map y(x) = x + i Delta (x^2 - 1); at x = +-1, y = +-1 (real); interior
    # contour dips below the real axis
    y = x_real + 1j * Delta * (x_real ** 2 - 1.0)
    dy_dx = 1.0 + 2j * Delta * x_real
    d2y_dx2 = 2.0 * 1j * Delta * np.ones_like(x_real)
    # chain rule: d/dy = (dy/dx)^(-1) d/dx
    #             d^2/dy^2 = (dy/dx)^(-2) d^2/dx^2 - (d^2 y / d x^2) (dy/dx)^(-3) d/dx
    # Operator: u_yy + (1/y) u = lambda u
    # On interior grid (imposing u(+-1) = 0 via restriction)
    idx = slice(1, N)
    D_int  = D[idx, idx]
    D2 = D @ D
    D2_int = D2[idx, idx]
    dydx_i  = dy_dx[idx]
    d2ydx2_i = d2y_dx2[idx]
    y_i      = y[idx]
    # operator L = d^2/dy^2 + 1/y
    L = (np.diag(dydx_i ** (-2)) @ D2_int
         - np.diag(d2ydx2_i * dydx_i ** (-3)) @ D_int
         + np.diag(1.0 / y_i))
    # eigenvalue problem: L u = lambda u
    lam = np.linalg.eigvals(L)
    # sort by magnitude (smallest |lambda| = physically meaningful bound states)
    idx_sort = np.argsort(np.abs(lam))
    return lam[idx_sort], y, y_i


def make_figure(Delta_scan=(0.25, 0.50, 0.75), N_demo: int = 40, N_low: int = 20):
    # reference at N = 40, Delta = 0.5
    lam_ref, y_ref, _ = solve_detoured(N_demo, 0.5)
    first6_ref = lam_ref[:6]

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.5))

    # Panel A: complex-plane contour for Delta = 0.5
    ax = axes[0]
    for Delta, c, ls in zip(Delta_scan, [NAVY, CORAL, TEAL], ["-", "--", ":"]):
        x_dense = np.linspace(-1, 1, 401)
        y_dense = x_dense + 1j * Delta * (x_dense ** 2 - 1.0)
        ax.plot(y_dense.real, y_dense.imag, ls, color=c, linewidth=1.2,
                label=fr"$\Delta = {Delta}$")
    # mark endpoints and singularity
    ax.axhline(0, color="k", linewidth=0.4, alpha=0.5)
    ax.axvline(0, color="k", linewidth=0.4, alpha=0.5)
    ax.plot([0], [0], "x", color=PURPLE, markersize=10, markeredgewidth=2,
            label="singularity $y = 0$")
    ax.plot([-1, 1], [0, 0], "o", color="k", markerfacecolor="white",
            markersize=8, markeredgewidth=1.5)
    ax.set_xlabel("Re $y$")
    ax.set_ylabel("Im $y$")
    ax.set_title("parabolic detour $y(x) = x + i\\Delta(x^2 - 1)$", fontsize=10)
    ax.legend(loc="lower right", fontsize=9)
    ax.grid(True, alpha=0.25, linewidth=0.4)
    ax.set_xlim(-1.2, 1.2)
    ax.set_ylim(-1.0, 0.3)

    # Panel B: complex eigenvalues at Delta scan
    ax = axes[1]
    markers = ["o", "s", "^"]
    for Delta, c, m in zip(Delta_scan, [NAVY, CORAL, TEAL], markers):
        lam, _, _ = solve_detoured(N_demo, Delta)
        first_few = lam[:7]
        ax.scatter(first_few.real, first_few.imag, s=45,
                   edgecolors=c, facecolors="white", linewidths=1.2,
                   marker=m, label=fr"$\Delta = {Delta}$, $N = {N_demo}$")
    ax.axhline(0, color="k", linewidth=0.4, alpha=0.5)
    ax.axvline(0, color="k", linewidth=0.4, alpha=0.5)
    ax.set_xlabel("Re $\\lambda$")
    ax.set_ylabel("Im $\\lambda$")
    ax.set_title("first eigenvalues in the complex plane", fontsize=10)
    ax.legend(loc="lower left", fontsize=9)
    ax.grid(True, alpha=0.25, linewidth=0.4)

    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "eigen_complex_detour.pdf", bbox_inches="tight", pad_inches=0.05)
    fig.savefig(OUTPUT_DIR / "eigen_complex_detour.png", bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)

    # sensitivity table: first 6 modes at three Delta values
    table = {}
    for Delta in Delta_scan:
        lam, _, _ = solve_detoured(N_demo, Delta)
        table[float(Delta)] = [(float(z.real), float(z.imag)) for z in lam[:6]]
    # low-N reference
    lam_low, _, _ = solve_detoured(N_low, 0.5)
    low_first6 = [(float(z.real), float(z.imag)) for z in lam_low[:6]]

    return {
        "N_demo": N_demo, "N_low": N_low,
        "Delta_scan": list(Delta_scan),
        "first6_ref_Delta0p5": [(float(z.real), float(z.imag)) for z in first6_ref],
        "low_N_first6":        low_first6,
        "sensitivity":         table,
    }


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--dump", type=str, default=None)
    args = p.parse_args()

    r = make_figure()
    print(f"[Etude 18.11]  complex-plane detour")
    print(f"  first 6 eigenvalues at N = {r['N_demo']}, Delta = 0.5:")
    for zr, zi in r["first6_ref_Delta0p5"]:
        print(f"    {zr:+.6f} {('+' if zi >= 0 else '-')} i {abs(zi):.6f}")
    print(f"  low-N (N={r['N_low']}, Delta=0.5) check, first 6:")
    for zr, zi in r["low_N_first6"]:
        print(f"    {zr:+.6f} {('+' if zi >= 0 else '-')} i {abs(zi):.6f}")
    print(f"  figure: {OUTPUT_DIR / 'eigen_complex_detour.pdf'}")
    if args.dump:
        rjson = dict(r)
        rjson["sensitivity"] = {str(k): v for k, v in r["sensitivity"].items()}
        Path(args.dump).write_text(json.dumps(rjson, indent=2))
        print(f"  dumped: {args.dump}")


if __name__ == "__main__":
    main()
