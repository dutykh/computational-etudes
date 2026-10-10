#!/usr/bin/env python3
"""
map_common.py

Shared plotting & style utilities for Chapter 19 (Coordinate Transformations).

Keeps the colour palette, matplotlib rc configuration, and the output
directory logic consistent across the eight etudes of Chapter 19.

Author: Dr. Denys Dutykh
        Mathematics Department
        Khalifa University of Science and Technology
        Abu Dhabi, UAE

Part of the book "Computational Etudes: A Spectral Approach"
https://github.com/dutykh/computational-etudes
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt

NAVY, SKY, CORAL = "#142D6E", "#7896D2", "#E74C3C"
TEAL, PURPLE, ORANGE = "#16A085", "#8E44AD", "#E67E22"
GOLD = "#D4A017"

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
    """Return the canonical figure output directory for a ch19 script.

    Accepts either the script's Path (with a .py suffix) or the directory
    containing it.  Walks up until the repository root (containing a
    ``textbook/`` sub-directory) is found.
    """
    p = script_path.resolve()
    if p.is_file():
        p = p.parent
    while p != p.parent and not (p / "textbook").is_dir():
        p = p.parent
    out = p / "textbook" / "figures" / "ch19" / "python"
    out.mkdir(parents=True, exist_ok=True)
    return out


def save_fig(fig, out_dir: Path, stem: str) -> None:
    for ext in ("pdf", "png"):
        fig.savefig(out_dir / f"{stem}.{ext}", bbox_inches="tight")
