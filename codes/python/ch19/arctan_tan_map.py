#!/usr/bin/env python3
"""
arctan_tan_map.py

Chapter 19: Coordinate Transformations
Illustration of the period-2 pi arctan/tan map (eq. 19.2),

    y = 2 * arctan( ell * tan(x / 2) ),    x, y in [-pi, pi].

Two panels: (a) the map curves y(x) for several values of the width
parameter ell, with the identity (ell = 1) drawn as a dashed grey line
for reference; (b) the images of a uniform N = 24 x-grid under the
map, plotted as tick marks on five stacked horizontal bands -- one per
ell -- making visible the textbook claim that ell < 1 clusters the
grid near y = 0, while ell > 1 clusters it near y = +/- pi.

Author: Dr. Denys Dutykh
        Mathematics Department
        Khalifa University of Science and Technology
        Abu Dhabi, UAE

Part of "Computational Etudes: A Spectral Approach"
https://github.com/dutykh/computational-etudes
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from map_common import (CORAL, NAVY, ORANGE, PURPLE, TEAL,
                        output_dir_for, save_fig, setup_matplotlib)

SCRIPT_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = output_dir_for(SCRIPT_DIR)

EPS_END = 1.0e-12  # safety margin so tan(x/2) does not blow up at x = +/- pi


def arctan_tan_map(x, ell):
    xc = np.clip(x, -np.pi + EPS_END, np.pi - EPS_END)
    return 2.0 * np.arctan(ell * np.tan(xc / 2.0))


# Five ell values + their colours.  ell = 1 (identity) is drawn dashed in grey.
ELLS = [0.1, 0.3, 1.0, 3.0, 10.0]
COLOURS = {0.1: CORAL, 0.3: TEAL, 1.0: "0.45", 3.0: ORANGE, 10.0: PURPLE}
LINESTYLES = {0.1: "-", 0.3: "-", 1.0: "--", 3.0: "-", 10.0: "-"}


def make_figure():
    setup_matplotlib()

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.4))

    # --- Panel (a): the map curves y(x) ------------------------------------
    ax = axes[0]
    x_line = np.linspace(-np.pi + EPS_END, np.pi - EPS_END, 1001)
    for ell in ELLS:
        ax.plot(x_line, arctan_tan_map(x_line, ell),
                color=COLOURS[ell], lw=1.6, ls=LINESTYLES[ell],
                label=fr"$\ell = {ell:g}$")
    # Light reference grid lines at +/- pi/2 and 0
    ax.axhline(0.0, color="0.7", lw=0.5)
    ax.axvline(0.0, color="0.7", lw=0.5)
    ticks = [-np.pi, -np.pi / 2, 0, np.pi / 2, np.pi]
    tick_labels = [r"$-\pi$", r"$-\pi/2$", r"$0$", r"$\pi/2$", r"$\pi$"]
    ax.set_xticks(ticks)
    ax.set_xticklabels(tick_labels)
    ax.set_yticks(ticks)
    ax.set_yticklabels(tick_labels)
    ax.set_xlim(-np.pi, np.pi)
    ax.set_ylim(-np.pi, np.pi)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, alpha=0.25)
    ax.set_xlabel(r"computational coordinate $x$")
    ax.set_ylabel(r"physical coordinate $y$")
    ax.set_title(r"(a) the map  $y = 2\,\arctan(\ell\,\tan(x/2))$")
    ax.legend(loc="upper left", frameon=False, fontsize=10)

    # --- Panel (b): clustering of a uniform N = 24 x-grid ------------------
    ax = axes[1]
    N = 24
    x_uniform = -np.pi + 2.0 * np.pi * (np.arange(N) + 0.5) / N
    # midpoint convention to keep both endpoints away from the singularity
    bands_y = np.arange(len(ELLS), 0, -1)  # 5, 4, 3, 2, 1 (top down)
    for band_y, ell in zip(bands_y, ELLS):
        y_mapped = arctan_tan_map(x_uniform, ell)
        col = COLOURS[ell]
        # base line
        ax.hlines(band_y, -np.pi, np.pi, color="0.85", lw=0.6, zorder=1)
        # tick marks at the mapped node positions
        ax.vlines(y_mapped, band_y - 0.18, band_y + 0.18,
                  color=col, lw=1.5, zorder=2)
        ax.text(-np.pi - 0.18, band_y, fr"$\ell = {ell:g}$",
                ha="right", va="center", fontsize=10, color=NAVY)
    ax.set_xlim(-np.pi, np.pi)
    ax.set_ylim(0.4, len(ELLS) + 0.6)
    ax.set_xticks(ticks)
    ax.set_xticklabels(tick_labels)
    ax.set_yticks([])
    ax.set_xlabel(r"physical coordinate $y$")
    ax.set_title(fr"(b) image of a uniform $N = {N}$ $x$-grid")

    fig.tight_layout()
    save_fig(fig, OUTPUT_DIR, "arctan_tan_map")
    plt.close(fig)

    print(f"[19.0] saved figure to {OUTPUT_DIR / 'arctan_tan_map.pdf'}")
    # Quick numerical sanity check: tightest spacing for each ell.
    for ell in ELLS:
        y_mapped = np.sort(arctan_tan_map(x_uniform, ell))
        gaps = np.diff(y_mapped)
        print(f"  ell={ell:5.2f}  min gap={gaps.min():.4f}  "
              f"max gap={gaps.max():.4f}  ratio={gaps.max()/gaps.min():.2f}")


if __name__ == "__main__":
    make_figure()
