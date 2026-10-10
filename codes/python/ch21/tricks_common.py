#!/usr/bin/env python3
"""
tricks_common.py

Shared utilities for Chapter 21 (Special Tricks for the Spectral Researcher).

Provides:
  * matplotlib setup and the consistent NAVY/CORAL/TEAL/PURPLE palette;
  * output directory resolution to textbook/figures/ch21/python/;
  * tiny helpers (save_fig, dump_json) used by every etude in this chapter.

Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
Part of "Computational Etudes: A Spectral Approach"
"""
from __future__ import annotations

import json
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
    out = p / "textbook" / "figures" / "ch21" / "python"
    out.mkdir(parents=True, exist_ok=True)
    return out


def save_fig(fig, out_dir: Path, stem: str) -> None:
    for ext in ("pdf", "png"):
        fig.savefig(out_dir / f"{stem}.{ext}", bbox_inches="tight")


def dump_json(path: str | None, payload: dict) -> None:
    if path is None:
        return
    with open(path, "w") as f:
        json.dump(payload, f, indent=2, default=_json_default)


def _json_default(obj):
    if isinstance(obj, (np.floating, np.integer)):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, complex):
        return {"re": obj.real, "im": obj.imag}
    raise TypeError(f"Cannot JSON-serialise {type(obj).__name__}")
