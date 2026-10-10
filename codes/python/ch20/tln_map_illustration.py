#!/usr/bin/env python3
"""
tln_map_illustration.py

Chapter 20: Spectral Methods on Unbounded Intervals --- illustrative figure
for the rational Chebyshev TL_n basis on the semi-infinite interval [0, infty).

Three panels visualise the algebraic semi-infinite map and the resulting
basis:

  (a) The map y = ell * (1 + x) / (1 - x) plotted as y vs x on [-1, 1) for
      ell in {1, 2, 4, 8}, showing how the singularity at x = 1 carries
      x in (-1, 1) onto the half-line y in [0, +infty), and how ell shifts
      where the grid concentrates.

  (b) The corresponding collocation grids on [0, +infty) at fixed N = 24
      for the same four ell values, as horizontal tick rows over a generic
      sech-decay target. Smaller ell crowds the grid near y = 0; larger
      ell pushes nodes outward.

  (c) The first five basis functions TL_n(y), n = 0..4, at ell = 2 over
      y in [0, 20]. Each TL_n asymptotes to a constant at y = +infty,
      which is exactly the right qualitative behaviour for problems with
      algebraic decay or asymptote-to-a-constant.

The map parameter is named `ell` (Greek script ell) rather than `L`
because chapter 20 also reserves `L` for the truncation half-width
[-L, L] of the domain-truncation methods of section 20.1.

Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
Part of "Computational Etudes: A Spectral Approach"
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
from unbounded_common import (NAVY, CORAL, TEAL, ORANGE, PURPLE, GOLD,  # noqa: E402
                              cheb_grid, output_dir_for, save_fig,
                              setup_matplotlib, tl_map_forward)

OUTPUT_DIR = output_dir_for(SCRIPT_DIR)


def TL_n(n: int, y: np.ndarray, ell: float) -> np.ndarray:
    """Rational Chebyshev TL_n basis: TL_n(y) = T_n( (y - ell)/(y + ell) )."""
    t = (y - ell) / (y + ell)
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


def tln_coeffs(N: int, ell: float):
    """Approximate TL_n expansion coefficients of sech(y) on [0, infty).

    Sample sech(y) at the TL_n collocation grid (interior CGL nodes
    pulled forward through y = ell (1+x)/(1-x)), then take a DCT-I to
    recover the coefficients of T_n(x) = TL_n(y(x)).
    """
    # We use N+1 grid points x_j = cos(pi j / N), j = 0..N, with the
    # endpoints x = +/-1 corresponding to y = +infty (where sech = 0)
    # and y = 0 respectively.
    j = np.arange(N + 1)
    x = np.cos(np.pi * j / N)
    fv = np.zeros_like(x)
    interior = np.abs(x) < 1.0 - 1e-12
    y = ell * (1.0 + x[interior]) / (1.0 - x[interior])
    fv[interior] = 1.0 / np.cosh(y)
    # Endpoints: x = -1 -> y = 0, sech(0) = 1; x = +1 -> y = inf, sech = 0
    fv[x < -1.0 + 1e-12] = 1.0
    extended = np.concatenate([fv, fv[N - 1:0:-1]])
    A = np.real(np.fft.fft(extended)) / N
    A[0] *= 0.5
    A[N] *= 0.5
    return A[:N + 1]


def make_figure(ell_values=(1.0, 2.0, 4.0, 8.0), N_grid: int = 24,
                ell_basis: float = 2.0):
    setup_matplotlib()
    fig, axes = plt.subplots(2, 2, figsize=(11.0, 8.0))
    palette_a = [NAVY, CORAL, TEAL, ORANGE]

    # ----- (a) The map y(x) for several ell -------------------------
    ax = axes[0, 0]
    x_dense = np.linspace(-1.0, 0.985, 401)
    for ell, color in zip(ell_values, palette_a):
        y_dense = ell * (1.0 + x_dense) / (1.0 - x_dense)
        ax.plot(x_dense, y_dense, color=color, linewidth=1.4,
                label=fr"$\ell = {ell:g}$")
    ax.axhline(0.0, color="k", linewidth=0.4, alpha=0.5)
    ax.axvline(+1.0, color="gray", linestyle=":", linewidth=0.6)
    ax.set_xlabel(r"computational coordinate $x$")
    ax.set_ylabel(r"physical coordinate $y$")
    ax.set_xlim(-1.05, 1.05)
    ax.set_ylim(0, 60)
    ax.set_title(r"(a) the map $y = \ell\, (1 + x) / (1 - x)$",
                 fontsize=10)
    ax.grid(True, alpha=0.25, linewidth=0.4)
    ax.legend(loc="upper left", fontsize=9, frameon=False)

    # ----- (b) Collocation grids at N = N_grid -----------------------
    ax = axes[0, 1]
    y_ref = np.linspace(0, 30, 601)
    profile = 0.4 / np.cosh(y_ref / 4.0)
    ax.plot(y_ref, profile + 0.45, color="gray", linewidth=0.8, alpha=0.4)
    ax.text(18, 0.88, r"target $\mathrm{sech}(y/4)$", fontsize=8,
            color="gray", alpha=0.7)
    row_y = np.linspace(0.32, 0.04, len(ell_values))
    x_int = cheb_grid(N_grid)[1:-1]      # interior collocation in t
    for ell, color, y_row in zip(ell_values, palette_a, row_y):
        y_int = tl_map_forward(x_int, ell)
        y_show = y_int[y_int <= 30.0]
        ax.scatter(y_show, np.full_like(y_show, y_row),
                   marker="|", s=55, color=color,
                   linewidths=1.1, label=fr"$\ell = {ell:g}$")
    ax.set_xlim(0, 30)
    ax.set_ylim(0.0, 1.0)
    ax.set_xlabel(r"physical coordinate $y$")
    ax.set_yticks([])
    ax.set_title(rf"(b) collocation grids at $N = {N_grid}$ for several $\ell$",
                 fontsize=10)
    ax.legend(loc="upper right", fontsize=9, frameon=False, ncol=2)
    ax.grid(True, axis="x", alpha=0.25, linewidth=0.4)

    # ----- (c) First five basis functions at ell = ell_basis -------
    ax = axes[1, 0]
    y_basis = np.linspace(0.001, 20, 801)
    palette_c = [NAVY, CORAL, TEAL, PURPLE, GOLD]
    for n, color in zip(range(5), palette_c):
        ax.plot(y_basis, TL_n(n, y_basis, ell_basis),
                color=color, linewidth=1.2,
                label=fr"$\mathrm{{TL}}_{n}$")
    ax.axhline(0.0, color="k", linewidth=0.4, alpha=0.5)
    ax.axhline(+1.0, color="gray", linestyle=":", linewidth=0.5, alpha=0.6)
    ax.axhline(-1.0, color="gray", linestyle=":", linewidth=0.5, alpha=0.6)
    ax.set_xlabel(r"physical coordinate $y$")
    ax.set_ylabel(r"$\mathrm{TL}_n(y)$")
    ax.set_xlim(0, 20)
    ax.set_ylim(-1.2, 1.25)
    ax.set_title(rf"(c) basis functions $\mathrm{{TL}}_0,\dots,\mathrm{{TL}}_4$ at $\ell = {ell_basis:g}$",
                 fontsize=10)
    ax.grid(True, alpha=0.25, linewidth=0.4)
    ax.legend(loc="lower right", fontsize=9, frameon=False, ncol=5,
              columnspacing=0.8, handlelength=1.4)

    # ----- (d) NEW: coefficient diagnostic for sech(y) ---------------
    ax = axes[1, 1]
    N_diag = 48
    diag_palette = [CORAL, NAVY, TEAL]
    for ell, color in zip([0.5, 2.0, 8.0], diag_palette):
        a = np.abs(tln_coeffs(N_diag, ell))
        ax.semilogy(np.arange(N_diag + 1), a + 1e-18, "-o", ms=3,
                    color=color, lw=1.0,
                    label=fr"$\ell = {ell:g}$")
    ax.set_xlabel(r"degree $n$")
    ax.set_ylabel(r"$|a_n|$ of $\mathrm{TL}_n$ expansion")
    ax.set_title(r"(d) coefficient diagnostic: $\mathrm{TL}_n$ expansion of "
                 r"$\mathrm{sech}(y)$", fontsize=10)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="upper right", fontsize=9, frameon=False)
    ax.set_ylim(1e-17, 10)

    fig.tight_layout()
    save_fig(fig, OUTPUT_DIR, "tln_map_illustration")
    plt.close(fig)
    return {"ell_values": list(ell_values), "N_grid": N_grid,
            "ell_basis": ell_basis}


if __name__ == "__main__":
    r = make_figure()
    print("[Etude 20.6 / illustrative figure]  rational Chebyshev TL_n map")
    print(f"  ell values shown: {r['ell_values']}")
    print(f"  N = {r['N_grid']} for grid panel; ell = {r['ell_basis']} for basis panel")
    print(f"  figure: {OUTPUT_DIR / 'tln_map_illustration.pdf'}")
