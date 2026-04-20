#!/usr/bin/env python3
"""
unbounded_common.py

Shared utilities for Chapter 20 (Spectral Methods on Unbounded Intervals).

Provides:
  * matplotlib setup and consistent colour palette;
  * output directory resolution;
  * the algebraic map y = L (x / sqrt(1 - x^2))   and its inverse
       x = y / sqrt(L^2 + y^2),
    which is Boyd's rational-Chebyshev mapping (TB_n) on (-infty, +infty);
  * the semi-infinite cotangent map  y = L cot^2(t/2)  and its inverse,
    which is Boyd's rational-Chebyshev mapping (TL_n) on [0, +infty);
  * utilities for Chebyshev-Gauss-Lobatto grids and DCT-I coefficients.

Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
Part of "Computational Etudes: A Spectral Approach"
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

NAVY, SKY, CORAL = "#142D6E", "#7896D2", "#E74C3C"
TEAL, PURPLE, ORANGE = "#16A085", "#8E44AD", "#E67E22"
GOLD, OLIVE = "#D4A017", "#6B8E23"

RC = {
    "font.family": "serif",
    "font.serif": ["CMU Serif", "DejaVu Serif"],
    "mathtext.fontset": "cm",
    "font.size": 11,
    "axes.linewidth": 0.8,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "figure.dpi": 150,
    "savefig.dpi": 300,
}


def setup_matplotlib():
    plt.rcParams.update(RC)


def output_dir_for(script_path: Path) -> Path:
    p = script_path.resolve()
    if p.is_file():
        p = p.parent
    while p != p.parent and not (p / "textbook").is_dir():
        p = p.parent
    out = p / "textbook" / "figures" / "ch20" / "python"
    out.mkdir(parents=True, exist_ok=True)
    return out


def save_fig(fig, out_dir: Path, stem: str) -> None:
    for ext in ("pdf", "png"):
        fig.savefig(out_dir / f"{stem}.{ext}", bbox_inches="tight")


# -----------------------------------------------------------------------------
# Chebyshev utilities (local copy to keep ch20 scripts self-contained).
# -----------------------------------------------------------------------------
def cheb_grid(N):
    """Chebyshev-Gauss-Lobatto nodes in descending order (Trefethen convention)."""
    return np.cos(np.pi * np.arange(N + 1) / N)


def dct1_coeffs(v):
    """DCT-I coefficients of samples v on Chebyshev-Gauss-Lobatto grid.

    v[k] = sum_n a_n T_n(x_k) exact at the N+1 collocation points.
    """
    N = len(v) - 1
    V = np.concatenate([v, v[N - 1:0:-1]])
    A = np.real(np.fft.fft(V)) / N
    A[0] *= 0.5
    A[N] *= 0.5
    return A[:N + 1]


def cheb_eval(a, xfine, N):
    """Clenshaw evaluation of sum_n a_n T_n(x) at xfine."""
    T0 = np.ones_like(xfine)
    T1 = xfine.copy()
    val = a[0] * T0 + (a[1] if N >= 1 else 0.0) * T1
    for n in range(2, N + 1):
        Tk = 2.0 * xfine * T1 - T0
        val += a[n] * Tk
        T0, T1 = T1, Tk
    return val


# -----------------------------------------------------------------------------
# The TB_n rational-Chebyshev map on (-infty, +infty).  Boyd's Eq 17.38a.
# -----------------------------------------------------------------------------
def tb_map_forward(x, L):
    """x in [-1, 1] -> y in (-infty, +infty)."""
    return L * x / np.sqrt(1.0 - x ** 2)


def tb_map_inverse(y, L):
    """y in (-infty, +infty) -> x in [-1, 1]."""
    return y / np.sqrt(L ** 2 + y ** 2)


def tb_map_fprime(x, L):
    """dy/dx for the TB_n map, evaluated at x-grid points."""
    return L / (1.0 - x ** 2) ** 1.5


def tb_map_fdoubleprime(x, L):
    """d2y/dx2 for the TB_n map."""
    return 3.0 * L * x / (1.0 - x ** 2) ** 2.5


# -----------------------------------------------------------------------------
# The TL_n rational-Chebyshev map on [0, +infty).  Boyd's Eq 17.61a.
# -----------------------------------------------------------------------------
def tl_map_forward(x, L):
    """x in [-1, 1] -> y in [0, +infty).  y = L(1+x)/(1-x)."""
    return L * (1.0 + x) / (1.0 - x)


def tl_map_inverse(y, L):
    """y in [0, +infty) -> x in [-1, 1].  x = (y-L)/(y+L)."""
    return (y - L) / (y + L)


def tl_map_fprime(x, L):
    return 2.0 * L / (1.0 - x) ** 2


def tl_map_fdoubleprime(x, L):
    return 4.0 * L / (1.0 - x) ** 3
