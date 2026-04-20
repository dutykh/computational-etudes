#!/usr/bin/env python3
"""
rK1_composed_map.py

Chapter 20: Spectral Methods on Unbounded Intervals.
Computational Etude 20.10: One global expansion instead of two local
formulas.

The modified Bessel function  r K_1(r),  r in [0, infty),  has a
logarithmic singularity at r = 0 and decays like exp(-r) at infinity.
Library software typically patches together TWO formulas: a power
series in r (with log terms) for small r, and an asymptotic series
in 1/r for large r.  We want a SINGLE global spectral expansion.

Boyd's trick (Sec 17.12, Example 2) is a composed map:

    r = arcsinh(exp y),        y in (-infty, +infty).

For large y, r ~ exp(y)/2, so the map is roughly linear in log-space
(good for the exponential decay at infinity).  For small y (i.e.,
r close to 0), r ~ exp(y), so the map resolves the logarithmic
singularity by using exponentially fine spacing near r = 0.

We then apply the TB_n rational-Chebyshev basis in the y-coordinate.
The result is geometric (subgeometrically so) convergence for  f(r) =
r K_1(r).

Reference: Abramowitz-Stegun series for K_1(r).  For small r,
    K_1(r) = 1/r + (r/2) log(r/2) + O(r);
    r K_1(r) = 1 + (r^2/2) log(r/2) + O(r^2).
For large r, K_1(r) ~ sqrt(pi / (2 r)) exp(-r) (1 + 3/(8r) + ...).

Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.special import k1

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR.parent / "ch07"))
from cheb_matrix import cheb_matrix  # noqa: E402

from unbounded_common import (CORAL, NAVY, TEAL, output_dir_for, save_fig,
                              setup_matplotlib, dct1_coeffs, cheb_eval,
                              tb_map_forward, tb_map_inverse)

OUTPUT_DIR = output_dir_for(SCRIPT_DIR)


def r_of_y(y):
    """Composed map  r = arcsinh(exp y)."""
    return np.arcsinh(np.exp(y))


def y_of_r(r):
    """Inverse:  y = log(sinh(r))."""
    # safe for large r: sinh(r) ~ exp(r)/2
    return np.log(np.sinh(r))


def target_rK1(r):
    """f(r) = r K_1(r).  Use scipy.special.k1 for the reference."""
    r = np.asarray(r)
    with np.errstate(invalid="ignore"):
        return r * k1(r)


def tbn_approx(N, L):
    """Chebyshev coefficients of  f(r(y))  on the TB_n grid in y.

    At x = +/- 1, y = +/- infty.
    At x = 1  (y = +infty): r = +infty,  r K_1(r) -> 0;  put 0.
    At x = -1 (y = -infty): r = 0^+, r K_1(r) -> 1;  put 1.
    """
    _, x = cheb_matrix(N)
    # interior grid: x in (-1, 1), y finite
    interior = np.abs(x) < 1.0 - 1e-12
    y = np.full_like(x, 0.0)
    y[interior] = tb_map_forward(x[interior], L)
    r = r_of_y(y)
    fv = np.where(interior, target_rK1(r), 0.0)
    # fix left endpoint: x = -1 corresponds to y = -infty -> r = 0 -> f = 1
    fv[-1] = 1.0                  # cheb_matrix returns x descending, so x[-1] = -1
    fv[0] = 0.0                   # x[0] = 1 -> y = +infty -> f = 0
    return dct1_coeffs(fv), interior


def evaluate(coeffs, r, L):
    y = y_of_r(r)
    x = tb_map_inverse(y, L)
    return cheb_eval(coeffs, x, len(coeffs) - 1)


def make_figure():
    setup_matplotlib()
    r_fine = np.logspace(-3, 1.6, 400)
    truth = target_rK1(r_fine)

    Ns = np.array([8, 12, 16, 20, 24, 32, 48, 64])
    L = 4.0
    errs = []
    for N in Ns:
        a, _ = tbn_approx(N, L)
        approx = evaluate(a, r_fine, L)
        errs.append(np.max(np.abs(approx - truth)))
    errs = np.array(errs)

    # Boyd's (approximate) coefficients for r K_1(r) expansion at small-N pre-analysis
    # from his Eq 17.68; just for the plot caption
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 3.8))

    ax = axes[0]
    ax.semilogx(r_fine, truth, color=NAVY, lw=1.2, label=r"exact $r K_1(r)$")
    a, _ = tbn_approx(16, L)
    ax.semilogx(r_fine, evaluate(a, r_fine, L), "--", color=CORAL, lw=1.0,
                label="$N = 16$ global $TB_n$ approximation")
    ax.set_xlabel(r"$r$")
    ax.set_ylabel(r"$f(r)$")
    ax.set_title(r"(a) Global expansion of $r K_1(r)$")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(frameon=False, fontsize=9)

    ax = axes[1]
    ax.semilogy(Ns, errs + 1e-18, "-o", color=TEAL, lw=1.1)
    ax.set_xlabel(r"$N$")
    ax.set_ylabel(r"$\max_r |f - f_N|$")
    ax.set_title(fr"(b) Subgeometric descent, $L = {L}$")
    ax.grid(True, which="both", alpha=0.3)

    fig.tight_layout()
    save_fig(fig, OUTPUT_DIR, "rK1_composed_map")
    plt.close(fig)

    print(f"[20.10] saved figure to {OUTPUT_DIR / 'rK1_composed_map.pdf'}")
    for N, e in zip(Ns, errs):
        print(f"  N={N:3d}  err={e:.3e}")


if __name__ == "__main__":
    make_figure()
