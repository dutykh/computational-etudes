#!/usr/bin/env python3
"""
weideman_cloot_map.py

Chapter 20: Spectral Methods on Unbounded Intervals --- illustrative figure
for the Weideman--Cloot (1990) sinh map.

The Weideman--Cloot map
    y = sinh(ell * t),   t in [-pi, pi]
turns Fourier pseudospectral methods on the periodic interval [-pi, pi]
into a hybrid map+truncation method on the unbounded line. Because
y_max = sinh(ell * pi) grows EXPONENTIALLY with ell, even a modest
logarithmic increase in ell drives the domain-truncation error
geometrically to zero. The trade-off is sensitivity: too-large ell
crowds the grid into the steeply growing tails and high-frequency
content is lost.

Three panels:

  (a) The map y = sinh(ell t) plotted as y vs t on [-pi, pi] for ell in
      {0.5, 1.0, 1.5, 2.0}, with the y_max = sinh(ell pi) reach annotated.

  (b) The corresponding uniform-in-t grids transferred to the y-line at
      fixed N = 64 for the same four ell values, plotted as horizontal
      tick rows. The visible compression near t = 0 reflects the y'(0) =
      ell, which is small; the wild stretch at t -> +- pi is the
      exponential reach.

  (c) The target sech(y) shown both directly on [-y_max, y_max] (where
      Fourier truncation must capture the full real-line decay) and
      transported back to t-space (where the same function becomes
      smooth and bounded with a sharp central peak that Fourier
      naturally resolves).

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
from unbounded_common import (NAVY, CORAL, TEAL, ORANGE,  # noqa: E402
                              output_dir_for, save_fig, setup_matplotlib)

OUTPUT_DIR = output_dir_for(SCRIPT_DIR)


def make_figure(ell_values=(0.5, 1.0, 1.5, 2.0), N_grid: int = 64,
                ell_demo: float = 1.0):
    setup_matplotlib()
    fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.0))
    palette_a = [NAVY, CORAL, TEAL, ORANGE]

    # ----- Panel (a): the map y(t) for several ell values -----
    ax = axes[0]
    t_dense = np.linspace(-np.pi, np.pi, 401)
    for ell, color in zip(ell_values, palette_a):
        y_dense = np.sinh(ell * t_dense)
        y_max = np.sinh(ell * np.pi)
        ax.plot(t_dense, y_dense, color=color, linewidth=1.4,
                label=fr"$\ell = {ell:g}$, $y_{{\max}} \approx {y_max:.1f}$")
    ax.axhline(0.0, color="k", linewidth=0.4, alpha=0.5)
    ax.axvline(-np.pi, color="gray", linestyle=":", linewidth=0.6)
    ax.axvline(+np.pi, color="gray", linestyle=":", linewidth=0.6)
    ax.set_xlabel(r"computational coordinate $t \in [-\pi, \pi]$")
    ax.set_ylabel(r"physical coordinate $y$")
    ax.set_xlim(-np.pi - 0.1, np.pi + 0.1)
    ax.set_ylim(-30, 30)
    ax.set_title(r"(a) the map $y = \sinh(\ell\, t)$",
                 fontsize=10)
    ax.grid(True, alpha=0.25, linewidth=0.4)
    ax.legend(loc="upper left", fontsize=8, frameon=False)

    # ----- Panel (b): grid clustering at N = N_grid for same ell values --
    ax = axes[1]
    # Backdrop: sech(y) at ell_demo, plotted in y space
    y_axis = np.linspace(-30, 30, 601)
    profile = 0.4 / np.cosh(y_axis / 1.0)
    ax.plot(y_axis, profile + 0.45, color="gray", linewidth=0.8, alpha=0.4)
    ax.text(15, 0.88, r"target $\mathrm{sech}(y)$", fontsize=8,
            color="gray", alpha=0.7)
    row_y = np.linspace(0.32, 0.04, len(ell_values))
    # uniform Fourier grid in t
    j = np.arange(N_grid)
    t_nodes = -np.pi + 2 * np.pi * j / N_grid
    for ell, color, y_row in zip(ell_values, palette_a, row_y):
        y_nodes = np.sinh(ell * t_nodes)
        y_show = y_nodes[np.abs(y_nodes) <= 30.0]
        ax.scatter(y_show, np.full_like(y_show, y_row),
                   marker="|", s=55, color=color,
                   linewidths=1.1, label=fr"$\ell = {ell:g}$")
    ax.set_xlim(-30, 30)
    ax.set_ylim(0.0, 1.0)
    ax.set_xlabel(r"physical coordinate $y$")
    ax.set_yticks([])
    ax.set_title(rf"(b) sinh-mapped grids at $N = {N_grid}$ for several $\ell$",
                 fontsize=10)
    ax.legend(loc="lower right", fontsize=9, frameon=False, ncol=2)
    ax.grid(True, axis="x", alpha=0.25, linewidth=0.4)

    # ----- Panel (c): sech(y) in y-space vs in t-space at ell = ell_demo --
    ax = axes[2]
    # Direct plot of sech(y) in y-space
    y_max = np.sinh(ell_demo * np.pi)
    y_phys = np.linspace(-y_max, y_max, 801)
    sech_in_y = 1.0 / np.cosh(y_phys)
    ax.plot(y_phys, sech_in_y, color=NAVY, linewidth=1.3,
            label=r"$\mathrm{sech}(y)$ in $y$-space")
    # Same function in t-space: f(t) = sech(sinh(ell t)) -- bounded, smooth
    t_axis = np.linspace(-np.pi, np.pi, 801)
    sech_in_t = 1.0 / np.cosh(np.sinh(ell_demo * t_axis))
    # Plot in t-space using the right-hand x-axis (we'll scale so both fit
    # nicely). The reader sees: same function, two coordinate views.
    # Map t into the same y-range visually for comparison: scale t into
    # [-y_max, y_max] linearly (it's just a re-parametrisation of the x
    # axis from t to scaled-t for plotting).
    t_scaled = t_axis * (y_max / np.pi)
    ax.plot(t_scaled, sech_in_t, color=CORAL, linewidth=1.3,
            linestyle="--",
            label=r"same function in $t$-space (scaled)")
    ax.axhline(0.0, color="k", linewidth=0.4, alpha=0.5)
    ax.set_xlabel(r"physical coordinate $y$ (or scaled $t$)")
    ax.set_ylabel(r"function value")
    ax.set_xlim(-y_max, y_max)
    ax.set_ylim(-0.1, 1.1)
    ax.set_title(rf"(c) $\mathrm{{sech}}(y)$ in two coordinate frames at $\ell = {ell_demo:g}$",
                 fontsize=10)
    ax.grid(True, alpha=0.25, linewidth=0.4)
    ax.legend(loc="upper right", fontsize=9, frameon=False)

    fig.tight_layout()
    save_fig(fig, OUTPUT_DIR, "weideman_cloot_map")
    plt.close(fig)
    return {"ell_values": list(ell_values),
            "N_grid": N_grid, "ell_demo": ell_demo,
            "y_max_per_ell": [float(np.sinh(e * np.pi)) for e in ell_values]}


if __name__ == "__main__":
    r = make_figure()
    print("[Section 20.12 / illustrative figure]  Weideman--Cloot sinh map")
    print(f"  ell values shown: {r['ell_values']}")
    print(f"  y_max = sinh(ell * pi) per ell: "
          f"{[round(y, 2) for y in r['y_max_per_ell']]}")
    print(f"  figure: {OUTPUT_DIR / 'weideman_cloot_map.pdf'}")
