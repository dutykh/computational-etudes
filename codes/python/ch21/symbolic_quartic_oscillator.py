#!/usr/bin/env python3
"""
symbolic_quartic_oscillator.py

Chapter 21: Special Tricks for the Spectral Researcher
Computational Etude 21.10: the quartic oscillator on an infinite interval,
                           solved symbolically.

The eigenvalue problem
        u_yy + (E - y^4) u = 0,    |u| -> 0  as |y| -> infty,
is the classic 'quartic oscillator' of quantum mechanics.  We attack it
in three stages:

1. Map the line to (-1, 1) via the rational-Chebyshev coordinate
        y = ell x / sqrt(1 - x^2).
   In x-space the equation becomes (Boyd Eq 20.13)
        (1-x^2)^4 {(1-x^2) u_xx - 3 x u_x} / ell^2
        + {(1-x^2)^2 E - ell^4 x^4} u  =  0.

2. Build in the boundary condition u(+/-1) = 0 by writing
        u(x) = (1 - x^2) (a_1 + a_2 x^2 + a_3 x^4 + a_4 x^6 + a_5 x^8).
   This is 'basis recombination' (Boyd Sec 20.2): a five-term ansatz
   that automatically satisfies the BC.

3. Demand orthogonality of the residual against the test functions
        {1, x^2, x^4, x^6, x^8}
   on x in [-1, 1].  This produces a 5x5 LINEAR system whose entries
   depend linearly on E.  The vanishing of its determinant is a
   degree-5 polynomial in E -- the 'secular determinant' D(E).

We compare the resulting D(E) to Boyd Eq 20.16 (after normalising the
constant term to 1) and the resulting eigenvalues to the high-precision
reference of Bender and Orszag (Adv Math Methods, 1978, p 523):

        E_0 = 1.060362,  E_2 = 7.45570,  E_4 = 16.2618,
        E_6 = 26.528,    E_8 = 37.92.

The lower three eigenvalues are accurate to 1-2 decimals; the upper two
are spectacularly inaccurate.  This is Boyd's enduring symbolic lesson:
trust the LOWER spectrum first.

Author: Dr. Denys Dutykh (Khalifa University of Science and Technology)
Part of "Computational Etudes: A Spectral Approach"
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import sympy as sp

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
from tricks_common import setup_matplotlib, output_dir_for, save_fig, dump_json
from tricks_common import NAVY, CORAL, TEAL


ELL = sp.Integer(2)


def build_secular_determinant():
    x, E = sp.symbols("x E", real=True)
    ell = ELL
    a = sp.symbols("a1:6")  # (a_1, a_2, a_3, a_4, a_5)

    u = (1 - x ** 2) * sum(a[k] * x ** (2 * k) for k in range(5))
    u_x = sp.diff(u, x)
    u_xx = sp.diff(u, x, 2)

    # Boyd Eq 20.13 (multiplied through by (1-x^2)^2 to clear the y^4 pole).
    R = (1 - x ** 2) ** 4 * ((1 - x ** 2) * u_xx - 3 * x * u_x) / ell ** 2 \
        + ((1 - x ** 2) ** 2 * E - ell ** 4 * x ** 4) * u

    # Galerkin moments with test functions x^{2j}, j = 0, 1, 2, 3, 4.
    eqs = [sp.integrate(x ** (2 * j) * R, (x, -1, 1)) for j in range(5)]

    # Build the 5x5 matrix in the unknowns (a_1, ..., a_5).
    M = sp.Matrix([[sp.simplify(eqs[i].coeff(a[j])) for j in range(5)] for i in range(5)])
    D_expr = sp.Poly(sp.simplify(M.det()), E)
    return D_expr, M, x, E


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    p.add_argument("--dump", type=str, default=None)
    args = p.parse_args()

    setup_matplotlib()
    out = output_dir_for(SCRIPT_DIR)

    D_poly, M, x, E = build_secular_determinant()
    coeffs_E = D_poly.all_coeffs()                      # high-degree first
    print("Symbolic secular determinant (high degree -> low):")
    for k, c in enumerate(coeffs_E):
        print(f"   E^{len(coeffs_E)-1-k}: {sp.nsimplify(c)} = {float(c):.6e}")
    norm_const = float(coeffs_E[-1])
    print("\nBoyd Eq 20.16 form  D(E) / D(0) = 1 + c_1 E + c_2 E^2 + ...:")
    for k, c in enumerate(coeffs_E[::-1]):
        print(f"   E^{k}: {float(c) / norm_const:.6e}")

    # Find numerical roots of the polynomial in E.
    coeffs_num = [float(c) for c in coeffs_E]
    roots = np.roots(coeffs_num)
    roots_real = sorted([float(np.real(r)) for r in roots if abs(np.imag(r)) < 1e-6])

    bender_orszag = [1.060362, 7.45570, 16.2618, 26.528, 37.92]   # n = 0, 2, 4, 6, 8
    print("\nEigenvalues from secular determinant vs Bender-Orszag (n = 0, 2, 4, 6, 8):")
    for r, ref in zip(roots_real, bender_orszag):
        rel = (r - ref) / ref
        print(f"   E_num = {r:11.4f}    E_exact = {ref:9.4f}    rel err = {rel:+.2e}")

    # ---- Figure ------------------------------------------------------
    fig = plt.figure(figsize=(11.5, 4.4))
    gs = fig.add_gridspec(1, 2, wspace=0.30)

    # Panel A: log10 |D(E)| on the relevant window [0, 50] -- zeros
    # appear as deep V-shaped dips, sharply localised, easy to read off.
    # The five Bender-Orszag references are placed as a "rug" of teal
    # circles on the x-axis; coral X markers locate the numerical roots
    # the symbolic determinant actually has in this window.  The lone
    # rogue root (the wholly-worthless E_8 ~ 720) is annotated with
    # an arrow rather than allowed to stretch the axis.
    ax = fig.add_subplot(gs[0, 0])
    E_window = (0.0, 50.0)
    Es = np.linspace(E_window[0], E_window[1], 2000)
    Ds = np.polyval(coeffs_num, Es)
    log_absD = np.log10(np.abs(Ds) + 1e-30)
    ax.plot(Es, log_absD, color=NAVY, linewidth=1.0, label=r"$\log_{10}|D(E)|$")

    # Roots and references inside the window
    in_window = lambda v: E_window[0] <= v <= E_window[1]
    roots_in = [r for r in roots_real if in_window(r)]
    refs_in  = [r for r in bender_orszag if in_window(r)]
    roots_out = [r for r in roots_real if not in_window(r)]

    # Rug strip just below the curve range
    rug_y = log_absD.min() - 0.6
    ax.scatter(refs_in, [rug_y] * len(refs_in), s=70, marker="o",
               facecolors="none", edgecolors=TEAL, linewidth=1.2,
               label=r"Bender-Orszag $E_n$ (reference)")
    ax.scatter(roots_in, [rug_y] * len(roots_in), s=80, marker="x",
               color=CORAL, linewidth=1.5,
               label="numerical roots of $D$")
    # vertical guide lines at references
    for ref in refs_in:
        ax.axvline(ref, color=TEAL, linewidth=0.6, alpha=0.35, linestyle=":")

    # Annotate the out-of-window rogue roots (the failure modes).
    # There are two: E_6 ~ 64 (way past its reference 26.5) and E_8 ~ 716
    # (way past its reference 37.9). Both are mirages; both deserve mention.
    if roots_out:
        labels = [r"$E_6$", r"$E_8$"]
        text_parts = [fr"{labels[i]} $\approx$ {r:.0f}"
                      for i, r in enumerate(roots_out[:2])]
        annot_text = ", ".join(text_parts) + " (numerical mirages, see right panel)"
        ax.annotate(annot_text,
                    xy=(E_window[1] - 0.5, rug_y),
                    xytext=(E_window[1] - 28, rug_y - 1.6),
                    color=CORAL, fontsize=8,
                    arrowprops=dict(arrowstyle="->", color=CORAL,
                                    lw=0.8, shrinkA=2, shrinkB=2))

    ax.set_xlabel(r"eigenvalue $E$")
    ax.set_ylabel(r"$\log_{10}|D(E)|$")
    ax.set_title(r"secular determinant zeros on the lower window $E\in[0,50]$",
                 fontsize=10)
    ax.set_xlim(*E_window)
    ax.set_ylim(rug_y - 2.2, log_absD.max() + 0.5)
    ax.legend(loc="upper left", fontsize=9)
    ax.grid(True, alpha=0.25, linewidth=0.4)

    # Panel B: relative error vs mode index
    ax = fig.add_subplot(gs[0, 1])
    ns = [0, 2, 4, 6, 8]
    rels = np.abs([(r - ref) / ref for r, ref in zip(roots_real, bender_orszag)])
    ax.semilogy(ns, rels, "o-", color=NAVY, markerfacecolor="white",
                markersize=8, linewidth=1.0,
                label=r"$|E_n^{\mathrm{sym}} - E_n| / E_n$")
    ax.axhline(0.01, color=TEAL, linewidth=0.4, alpha=0.5,
               label="1% threshold (lower spectrum target)")
    ax.set_xlabel(r"physical mode index $n$ (even modes only)")
    ax.set_ylabel("relative error")
    ax.set_title("trust the lower spectrum: bottom 3 modes good, top 2 unusable",
                 fontsize=10)
    ax.legend(loc="upper left", fontsize=9)
    ax.grid(True, alpha=0.25, linewidth=0.4)
    ax.set_ylim(1e-4, 1e2)

    save_fig(fig, out, "symbolic_quartic_oscillator")
    plt.close(fig)

    print(f"  figure: {out / 'symbolic_quartic_oscillator.pdf'}")

    dump_json(args.dump, {
        "ell": float(ELL),
        "secular_coeffs": [float(c) for c in coeffs_E],
        "roots": roots_real,
        "bender_orszag": bender_orszag,
        "rel_err": rels.tolist(),
    })


if __name__ == "__main__":
    main()
