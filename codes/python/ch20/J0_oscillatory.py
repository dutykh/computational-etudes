#!/usr/bin/env python3
"""
J0_oscillatory.py

Chapter 20: Spectral Methods on Unbounded Intervals.
Computational Etude 20.11: Teach the basis the oscillation first.

The Bessel function  J_0(y),  y in [0, +infty),  decays only as
1 / sqrt(y) and continues to oscillate at a non-decaying carrier
frequency of unity.  A naive TL_n-type rational-Chebyshev expansion
would have to resolve thousands of wavelengths to achieve even
modest accuracy, because truncating the basis at any finite N leaves
an O(1) oscillation error in the tail.

Boyd's (1987b) cure is ASYMPTOTIC AUGMENTATION: write

    sqrt(1 + y) J_0(y) ~ a(y) cos(y - pi/4 - phi(y))

and expand  a(y)  and  phi(y)  -- both smoothly asymptoting to
constants -- as TL_n series.  The basis now "knows" the oscillation
and only needs to spell out the slowly-varying amplitude and phase.

In this etude we do a slightly simplified version of Boyd's (1987b)
construction: we represent

    sqrt(1 + y) J_0(y) ~  sum_{n=0}^{N} a_n TL_n(y) cos(y - pi/4)
                         + sum_{n=0}^{N} b_n TL_n(y) sin(y - pi/4),

fit the coefficients {a_n, b_n} by least squares on a fine sample,
and compare with the naive  sqrt(1 + y) J_0(y) ~ sum_n c_n TL_n(y)
fit of the same degree.  The figure reproduces Boyd's Fig 17.11 in
qualitative spirit.

Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.special import j0

from unbounded_common import (CORAL, NAVY, TEAL, output_dir_for, save_fig,
                              setup_matplotlib, tl_map_inverse)

SCRIPT_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = output_dir_for(SCRIPT_DIR)


def target(y):
    return np.sqrt(1.0 + y) * j0(y)


def TL_n(y, n, L):
    """TL_n(y; L) = T_n((y - L)/(y + L))."""
    x = tl_map_inverse(y, L)
    if n == 0:
        return np.ones_like(y)
    if n == 1:
        return x
    # Clenshaw
    T0 = np.ones_like(y)
    T1 = x.copy()
    for k in range(2, n + 1):
        T2 = 2.0 * x * T1 - T0
        T0, T1 = T1, T2
    return T1


def build_TL_block(ys, N, L):
    """Design matrix of TL_0, ..., TL_N at sample points ys."""
    M = np.zeros((len(ys), N + 1))
    x = tl_map_inverse(ys, L)
    M[:, 0] = 1.0
    if N >= 1:
        M[:, 1] = x
    T0 = M[:, 0]; T1 = M[:, 1] if N >= 1 else None
    for k in range(2, N + 1):
        M[:, k] = 2.0 * x * M[:, k - 1] - M[:, k - 2]
    return M


def naive_fit(y_sample, N, L):
    """Fit g(y) = sum_n c_n TL_n(y)."""
    f_sample = target(y_sample)
    M = build_TL_block(y_sample, N, L)
    c, *_ = np.linalg.lstsq(M, f_sample, rcond=None)
    return c


def augmented_fit(y_sample, N, L):
    """Fit g(y) = cos(y - pi/4) sum a_n TL_n + sin(y - pi/4) sum b_n TL_n."""
    f_sample = target(y_sample)
    TL_block = build_TL_block(y_sample, N, L)
    C = np.cos(y_sample - np.pi / 4.0)
    S = np.sin(y_sample - np.pi / 4.0)
    design = np.hstack([TL_block * C[:, None], TL_block * S[:, None]])
    ab, *_ = np.linalg.lstsq(design, f_sample, rcond=None)
    a = ab[:N + 1]; b = ab[N + 1:]
    return a, b


def naive_eval(c, y, L):
    return build_TL_block(y, len(c) - 1, L) @ c


def aug_eval(a, b, y, L):
    TL_block = build_TL_block(y, len(a) - 1, L)
    return (TL_block @ a) * np.cos(y - np.pi / 4.0) + (TL_block @ b) * np.sin(y - np.pi / 4.0)


def make_figure():
    setup_matplotlib()
    L = 4.0
    y_fine = np.linspace(0.01, 50.0, 8001)
    truth = target(y_fine)

    # dense sample for least-squares fit
    y_sample = np.linspace(0.01, 80.0, 2001)

    Ns = np.array([4, 6, 8, 10, 15, 20, 30, 40])
    err_naive, err_aug = [], []
    for N in Ns:
        c = naive_fit(y_sample, N, L)
        err_naive.append(np.max(np.abs(naive_eval(c, y_fine, L) - truth)))
        a, b = augmented_fit(y_sample, N, L)
        err_aug.append(np.max(np.abs(aug_eval(a, b, y_fine, L) - truth)))
    err_naive = np.array(err_naive); err_aug = np.array(err_aug)

    fig, axes = plt.subplots(1, 2, figsize=(10.0, 3.8))

    ax = axes[0]
    ax.plot(y_fine, truth, color=NAVY, lw=1.0, label=r"$\sqrt{1+y}\,J_0(y)$")
    a, b = augmented_fit(y_sample, 15, L)
    # amplitude and phase functions
    amp = build_TL_block(y_fine, 15, L) @ a
    phi = build_TL_block(y_fine, 15, L) @ b
    ax.plot(y_fine, np.abs(amp), "--", color=CORAL, lw=1.0, label=r"$|a(y)|$")
    ax.plot(y_fine, np.abs(phi), ":", color=TEAL, lw=1.0, label=r"$|\phi(y)|$")
    ax.set_xlim(0, 20)
    ax.set_xlabel(r"$y$")
    ax.set_title(r"(a) $\sqrt{1+y}\,J_0$ with amp-phase decomposition")
    ax.legend(frameon=False, fontsize=9)
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    ax.semilogy(Ns, err_naive + 1e-18, "-o", color=CORAL, lw=1.1,
                label=r"naive $\sum c_n TL_n$")
    ax.semilogy(Ns, err_aug + 1e-18, "-s", color=TEAL, lw=1.1, mfc="none",
                label=r"augmented $\cos/\sin$ carrier")
    ax.set_xlabel(r"$N$")
    ax.set_ylabel(r"$\|g - g_N\|_\infty$ on $[0, 50]$")
    ax.set_title("(b) Asymptotic augmentation wins")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(frameon=False, fontsize=9)

    fig.tight_layout()
    save_fig(fig, OUTPUT_DIR, "J0_oscillatory")
    plt.close(fig)

    print(f"[20.11] saved figure to {OUTPUT_DIR / 'J0_oscillatory.pdf'}")
    for N, n, a in zip(Ns, err_naive, err_aug):
        print(f"  N={N:3d}  naive={n:.3e}  aug={a:.3e}")


if __name__ == "__main__":
    make_figure()
