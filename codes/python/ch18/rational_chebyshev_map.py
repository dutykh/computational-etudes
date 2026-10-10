#!/usr/bin/env python3
"""
rational_chebyshev_map.py

Chapter 18: Linear Spectral Eigenproblems --- illustrative figure for
the rational Chebyshev TB_n basis.

Three panels visualise the algebraic map and the resulting basis:

  (a) The map x = ell * t / sqrt(1 - t^2) plotted as x vs t on [-1, 1] for
      ell in {1, 2, 4, 8}, showing how the choice of ell stretches the
      interior toward the real line and how t -> +-1 maps to x -> +- inf.

  (b) The corresponding collocation grids on the real line at fixed N = 24
      for the same four ell values: each grid is drawn as a row of tick
      marks so the reader can see the clustering shift directly. With sech
      profiles superimposed, the panel shows the trade-off: small ell
      crowds the grid near x = 0; large ell wastes nodes in the decayed
      tails.

  (c) The first five rational Chebyshev basis functions TB_n(x), n = 0..4,
      at ell = 2 over x in [-10, 10]. Each TB_n asymptotes to a constant
      at +- infinity, which is precisely why the basis is suitable for
      functions that decay algebraically or asymptote.

Author: Dr. Denys Dutykh (Khalifa University of Science and Technology)
Part of "Computational Etudes: A Spectral Approach"
https://github.com/dutykh/computational-etudes
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
from rational_chebyshev import rational_chebyshev_grid  # noqa: E402

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

NAVY    = "#142D6E"
SKY     = "#7896D2"
CORAL   = "#E74C3C"
TEAL    = "#16A085"
ORANGE  = "#E67E22"
PURPLE  = "#8E44AD"
GOLD    = "#D4A017"

OUTPUT_DIR = SCRIPT_DIR.parent.parent.parent / "textbook" / "figures" / "ch18" / "python"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def TB_n(n: int, x: np.ndarray, ell: float) -> np.ndarray:
    """Rational Chebyshev basis: TB_n(x) = T_n( x / sqrt(ell^2 + x^2) )."""
    t = x / np.sqrt(ell ** 2 + x ** 2)
    if n == 0:
        return np.ones_like(t)
    if n == 1:
        return t
    Tk_minus2 = np.ones_like(t)
    Tk_minus1 = t
    Tk = None
    for _ in range(2, n + 1):
        Tk = 2 * t * Tk_minus1 - Tk_minus2
        Tk_minus2, Tk_minus1 = Tk_minus1, Tk
    return Tk


def make_figure(ell_values=(1.0, 2.0, 4.0, 8.0), N_grid: int = 24,
                ell_basis: float = 2.0):
    fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.0))

    # ----- Panel (a): the map x(t) for several ell values -----
    ax = axes[0]
    t_dense = np.linspace(-0.99, 0.99, 401)
    palette_a = [NAVY, CORAL, TEAL, ORANGE]
    for ell, color in zip(ell_values, palette_a):
        x_dense = ell * t_dense / np.sqrt(1.0 - t_dense ** 2)
        ax.plot(t_dense, x_dense, color=color, linewidth=1.4,
                label=fr"$\ell = {ell:g}$")
    ax.axhline(0.0, color="k", linewidth=0.4, alpha=0.5)
    ax.axvline(-1.0, color="gray", linestyle=":", linewidth=0.6)
    ax.axvline(+1.0, color="gray", linestyle=":", linewidth=0.6)
    ax.set_xlabel(r"computational coordinate $t$")
    ax.set_ylabel(r"physical coordinate $x$")
    ax.set_xlim(-1.05, 1.05)
    ax.set_ylim(-15, 15)
    ax.set_title(r"(a) the map $x = \ell\, t / \sqrt{1 - t^2}$",
                 fontsize=10)
    ax.grid(True, alpha=0.25, linewidth=0.4)
    ax.legend(loc="upper left", fontsize=9, frameon=False)

    # ----- Panel (b): grid clustering at N = N_grid for same ell values --
    ax = axes[1]
    # background reference: sech^2 (a generic decaying target)
    x_ref = np.linspace(-15, 15, 601)
    profile = 0.4 * np.exp(-(x_ref / 4.0) ** 2)
    ax.plot(x_ref, profile + 0.45, color="gray", linewidth=0.8, alpha=0.4)
    ax.text(8.5, 0.88, r"target $e^{-(x/4)^2}$", fontsize=8,
            color="gray", alpha=0.7)
    # one row of ticks per ell
    row_y = np.linspace(0.32, 0.04, len(ell_values))
    for (ell, color, y_row) in zip(ell_values, palette_a, row_y):
        x_int, _, _ = rational_chebyshev_grid(N_grid, ell)
        # only show |x| <= 15 to keep the panel readable
        x_show = x_int[np.abs(x_int) <= 15.0]
        ax.scatter(x_show, np.full_like(x_show, y_row),
                   marker="|", s=55, color=color,
                   linewidths=1.1, label=fr"$\ell = {ell:g}$")
    ax.set_xlim(-15, 15)
    ax.set_ylim(0.0, 1.0)
    ax.set_xlabel(r"physical coordinate $x$")
    ax.set_yticks([])
    ax.set_title(rf"(b) collocation grids at $N = {N_grid}$ for several $\ell$",
                 fontsize=10)
    ax.legend(loc="lower right", fontsize=9, frameon=False, ncol=2)
    ax.grid(True, axis="x", alpha=0.25, linewidth=0.4)

    # ----- Panel (c): first five basis functions at ell = ell_basis ------
    ax = axes[2]
    x_basis = np.linspace(-10, 10, 801)
    palette_c = [NAVY, CORAL, TEAL, PURPLE, GOLD]
    for n, color in zip(range(5), palette_c):
        ax.plot(x_basis, TB_n(n, x_basis, ell_basis),
                color=color, linewidth=1.2,
                label=fr"$\mathrm{{TB}}_{n}$")
    ax.axhline(0.0, color="k", linewidth=0.4, alpha=0.5)
    ax.axhline(+1.0, color="gray", linestyle=":", linewidth=0.5, alpha=0.6)
    ax.axhline(-1.0, color="gray", linestyle=":", linewidth=0.5, alpha=0.6)
    ax.set_xlabel(r"physical coordinate $x$")
    ax.set_ylabel(r"$\mathrm{TB}_n(x)$")
    ax.set_xlim(-10, 10)
    ax.set_ylim(-1.2, 1.25)
    ax.set_title(rf"(c) basis functions $\mathrm{{TB}}_0,\dots,\mathrm{{TB}}_4$ at $\ell = {ell_basis:g}$",
                 fontsize=10)
    ax.grid(True, alpha=0.25, linewidth=0.4)
    ax.legend(loc="lower right", fontsize=9, frameon=False, ncol=5,
              columnspacing=0.8, handlelength=1.4)

    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "rational_chebyshev_map.pdf",
                bbox_inches="tight", pad_inches=0.05)
    fig.savefig(OUTPUT_DIR / "rational_chebyshev_map.png",
                bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    return {
        "ell_values":  list(ell_values),
        "N_grid":      N_grid,
        "ell_basis":   ell_basis,
    }


if __name__ == "__main__":
    r = make_figure()
    print("[Etude 18.4 / illustrative figure]  rational Chebyshev map")
    print(f"  ell values shown: {r['ell_values']}")
    print(f"  N = {r['N_grid']} for grid panel; ell = {r['ell_basis']} for basis panel")
    print(f"  figure: {OUTPUT_DIR / 'rational_chebyshev_map.pdf'}")
