#!/usr/bin/env python3
"""
pseudospectra_demo.py
Chapter 17: Stability of Spectral Methods

Computational Etude 17.5: Pseudospectra of Normal vs Nonnormal Matrices.

Two-panel figure comparing pseudospectra:
  (a) Normal diagonal matrix with eigenvalues -1, -2, ..., -16.
      Concentric circles around each eigenvalue.
  (b) Nonnormal upper triangular matrix with same diagonal but A[i,i+1] = 10.
      Contours bulge into the right half-plane, showing transient growth
      potential.

Pseudospectra computed via resolvent norm:
  ||(zI - A)^{-1}|| = 1 / sigma_min(zI - A)

Generates Figure 17.5: pseudospectra_comparison.pdf

Author: Dr. Denys Dutykh
        Mathematics Department
        Khalifa University of Science and Technology
        Abu Dhabi, UAE

Part of the book "Computational Etudes: A Spectral Approach"
https://github.com/dutykh/computational-etudes
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from pathlib import Path

rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Computer Modern Roman', 'CMU Serif', 'DejaVu Serif']
rcParams['mathtext.fontset'] = 'cm'
rcParams['axes.labelsize'] = 11
rcParams['axes.titlesize'] = 12
rcParams['xtick.labelsize'] = 10
rcParams['ytick.labelsize'] = 10
rcParams['legend.fontsize'] = 9
rcParams['figure.dpi'] = 150
rcParams['savefig.dpi'] = 300
rcParams['text.usetex'] = False
rcParams['axes.linewidth'] = 0.8
rcParams['xtick.major.width'] = 0.8
rcParams['ytick.major.width'] = 0.8

NAVY = '#142D6E'
SKY = '#7896D2'
CORAL = '#E74C3C'
TEAL = '#16A085'

SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / '..' / '..' / '..' / 'textbook' / 'figures' / 'ch17' / 'python'


def compute_pseudospectra(A, x_range, y_range, n_grid=200):
    """Compute log10(resolvent norm) on a grid.

    resolvent_norm(z) = ||(zI - A)^{-1}|| = 1 / sigma_min(zI - A)
    """
    x = np.linspace(x_range[0], x_range[1], n_grid)
    y = np.linspace(y_range[0], y_range[1], n_grid)
    X, Y = np.meshgrid(x, y)

    n = A.shape[0]
    log_resolvent = np.zeros_like(X)

    for i in range(n_grid):
        for j in range(n_grid):
            z = X[i, j] + 1j * Y[i, j]
            M = z * np.eye(n) - A
            smin = np.linalg.svd(M, compute_uv=False)[-1]
            if smin > 0:
                log_resolvent[i, j] = np.log10(1.0 / smin)
            else:
                log_resolvent[i, j] = 15.0  # effectively infinity

    return X, Y, log_resolvent


def main():
    N = 16

    # Normal matrix: diagonal with eigenvalues -1, -2, ..., -N
    A_normal = np.diag(-np.arange(1.0, N + 1))

    # Nonnormal matrix: same diagonal, upper triangular with A[i,i+1] = 10
    A_nonnormal = A_normal.copy()
    for i in range(N - 1):
        A_nonnormal[i, i + 1] = 10.0

    eigenvalues = -np.arange(1.0, N + 1)

    print("Computing pseudospectra of normal matrix...")
    x_range_n = (-18, 4)
    y_range_n = (-8, 8)
    X_n, Y_n, log_res_n = compute_pseudospectra(
        A_normal, x_range_n, y_range_n, n_grid=250)

    print("Computing pseudospectra of nonnormal matrix...")
    x_range_nn = (-20, 6)
    y_range_nn = (-10, 10)
    X_nn, Y_nn, log_res_nn = compute_pseudospectra(
        A_nonnormal, x_range_nn, y_range_nn, n_grid=250)

    # Contour levels: epsilon^{-1} = 10^1, 10^2, 10^3
    # i.e., log10(resolvent_norm) = 1, 2, 3
    levels = [1, 2, 3]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    # --- Panel (a): Normal matrix ---
    cs1 = ax1.contour(X_n, Y_n, log_res_n, levels=levels,
                       colors=[SKY, NAVY, '#0A1A40'],
                       linewidths=[1.0, 1.2, 1.5])
    ax1.clabel(cs1, fmt={1: r'$\varepsilon^{-1}=10$',
                          2: r'$10^2$',
                          3: r'$10^3$'},
               fontsize=8, inline=True)
    ax1.plot(eigenvalues, np.zeros(N), 'o', color=CORAL, markersize=5,
             markerfacecolor=CORAL, markeredgecolor='white',
             markeredgewidth=0.5, zorder=5, label='eigenvalues')
    ax1.axhline(0, color='gray', linewidth=0.5, zorder=0)
    ax1.axvline(0, color='gray', linewidth=0.5, zorder=0)
    ax1.set_xlabel(r'Re$(z)$')
    ax1.set_ylabel(r'Im$(z)$')
    ax1.set_title(r'(a) Normal matrix ($A = \mathrm{diag}(-1, \ldots, -16)$)')
    ax1.legend(loc='upper left', fontsize=8)
    ax1.grid(True, alpha=0.15, linewidth=0.5)
    ax1.set_aspect('equal')

    # --- Panel (b): Nonnormal matrix ---
    cs2 = ax2.contour(X_nn, Y_nn, log_res_nn, levels=levels,
                       colors=[SKY, NAVY, '#0A1A40'],
                       linewidths=[1.0, 1.2, 1.5])
    ax2.clabel(cs2, fmt={1: r'$\varepsilon^{-1}=10$',
                          2: r'$10^2$',
                          3: r'$10^3$'},
               fontsize=8, inline=True)
    ax2.plot(eigenvalues, np.zeros(N), 'o', color=CORAL, markersize=5,
             markerfacecolor=CORAL, markeredgecolor='white',
             markeredgewidth=0.5, zorder=5, label='eigenvalues')
    ax2.axhline(0, color='gray', linewidth=0.5, zorder=0)
    ax2.axvline(0, color='gray', linewidth=0.5, zorder=0)
    ax2.set_xlabel(r'Re$(z)$')
    ax2.set_ylabel(r'Im$(z)$')
    ax2.set_title(r'(b) Nonnormal ($A_{i,i+1} = 10$, same diagonal)')
    ax2.legend(loc='upper left', fontsize=8)
    ax2.grid(True, alpha=0.15, linewidth=0.5)
    ax2.set_aspect('equal')

    # Highlight the bulge into RHP
    ax2.annotate('pseudospectral\nbulge into RHP',
                 xy=(2, 3), fontsize=9, color=CORAL, ha='center',
                 style='italic')

    plt.tight_layout()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_DIR / 'pseudospectra_comparison.pdf',
                bbox_inches='tight')
    fig.savefig(OUTPUT_DIR / 'pseudospectra_comparison.png',
                bbox_inches='tight')
    print(f"Figure saved to {OUTPUT_DIR / 'pseudospectra_comparison.pdf'}")

    plt.show()


if __name__ == '__main__':
    main()
