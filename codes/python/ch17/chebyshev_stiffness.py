#!/usr/bin/env python3
"""
chebyshev_stiffness.py
Chapter 17: Stability of Spectral Methods

Computational Etude 17.4: Chebyshev Stiffness and Eigenvalue Structure.

Three figures:
  (1) cheb_eigenvalues.pdf -- Eigenvalues of Chebyshev D^2 (Dirichlet BCs)
      for N=32 plotted on the real axis.  Physical eigenvalues (close to
      -n^2*pi^2/4) marked with circles; outliers in a different color.
  (2) cheb_eigenvectors.pdf -- Two panels showing a physical eigenvector
      (smooth sinusoidal, mode N//4) and the outlier eigenvector (localised
      near x = +/- 1).
  (3) cheb_cfl_scaling.pdf -- Log-log plot of spectral radius rho(D2) vs N,
      showing slope = 4 and reference line 0.048*N^4.

Generates Figures 17.4a,b,c

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
AMBER = '#E67E22'

SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / '..' / '..' / '..' / 'textbook' / 'figures' / 'ch17' / 'python'


def cheb_diff_matrix(N):
    """Chebyshev differentiation matrix of size (N+1) x (N+1)."""
    if N == 0:
        return np.array([[0.0]]), np.array([1.0])
    x = np.cos(np.pi * np.arange(N + 1) / N)
    c = np.ones(N + 1)
    c[0] = 2.0
    c[N] = 2.0
    c *= (-1.0) ** np.arange(N + 1)
    X = np.outer(x, np.ones(N + 1))
    dX = X - X.T + np.eye(N + 1)
    D = np.outer(c, 1.0 / c) / dX
    D -= np.diag(D.sum(axis=1))
    return D, x


def figure_eigenvalues():
    """Plot eigenvalues of Chebyshev D^2 for N=32.

    Two-panel figure: top shows the full range (outliers visible),
    bottom is a close-up of the physical eigenvalues with the
    predicted values -n^2*pi^2/4 overlaid.
    """
    N = 32
    D, x = cheb_diff_matrix(N)
    D2 = D @ D
    D2_int = D2[1:N, 1:N]

    eigs = np.sort(np.real(np.linalg.eigvals(D2_int)))

    # Predicted physical eigenvalues: -(n*pi/2)^2 for n = 1, 2, ...
    n_phys = np.arange(1, N)
    phys_approx = -(n_phys * np.pi / 2.0) ** 2

    # Separate physical from outliers using the predicted range
    # Physical eigenvalues cluster near -(n*pi/2)^2 for n up to ~2N/pi
    threshold = -(N * np.pi / 2.0) ** 2 * 0.5
    phys_mask = eigs > threshold
    out_mask = ~phys_mask

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 5), height_ratios=[1, 1.2])

    # --- Top panel: full range ---
    ax1.plot(eigs[phys_mask], np.zeros(np.sum(phys_mask)), 'o',
             color=NAVY, markersize=6, label='physical modes', zorder=5)
    if np.any(out_mask):
        ax1.plot(eigs[out_mask], np.zeros(np.sum(out_mask)), 's',
                 color=CORAL, markersize=7, label='outlier modes', zorder=5)

    ax1.axhline(0, color='gray', linewidth=0.5, zorder=0)
    ax1.set_xlabel(r'$\lambda$')
    ax1.set_title(r'Eigenvalues of Chebyshev $D^2$ with Dirichlet BCs ($N = 32$)')
    ax1.set_yticks([])
    ax1.set_ylim(-0.3, 0.3)
    ax1.legend(loc='upper left', fontsize=9)
    ax1.grid(True, axis='x', alpha=0.3, linewidth=0.5)

    # Annotate the worst outlier
    ax1.annotate(fr'$\lambda_{{\max}} \approx {eigs[0]:.0f}$',
                 xy=(eigs[0], 0), xytext=(eigs[0] * 0.7, 0.2),
                 fontsize=9, color=CORAL,
                 arrowprops=dict(arrowstyle='->', color=CORAL, lw=1.2))

    # --- Bottom panel: close-up of physical eigenvalues ---
    # Show only eigenvalues in the physical range
    close_range = max(phys_approx.min() * 1.2, eigs[phys_mask].min() * 1.1)
    phys_eigs = eigs[eigs > close_range]

    ax2.plot(phys_eigs, np.zeros(len(phys_eigs)), 'o',
             color=NAVY, markersize=7, label='computed', zorder=5)

    # Overlay predicted values
    n_show = min(len(phys_eigs), 15)
    ax2.plot(phys_approx[:n_show], np.full(n_show, 0.015), 'v',
             color=TEAL, markersize=6, alpha=0.8,
             label=r'predicted $-n^2\pi^2/4$', zorder=4)

    ax2.axhline(0, color='gray', linewidth=0.5, zorder=0)
    ax2.set_xlabel(r'$\lambda$')
    ax2.set_title(r'Close-up: physical eigenvalues vs.\ predicted $-n^2\pi^2/4$')
    ax2.set_ylim(-0.06, 0.06)
    ax2.set_yticks([])
    ax2.set_xlim(close_range, 10)
    ax2.legend(loc='upper left', fontsize=9)
    ax2.grid(True, axis='x', alpha=0.3, linewidth=0.5)

    plt.tight_layout()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_DIR / 'cheb_eigenvalues.pdf', bbox_inches='tight')
    fig.savefig(OUTPUT_DIR / 'cheb_eigenvalues.png', bbox_inches='tight')
    print(f"Figure saved to {OUTPUT_DIR / 'cheb_eigenvalues.pdf'}")
    return eigs


def figure_eigenvectors():
    """Plot physical and outlier eigenvectors."""
    N = 32
    D, x = cheb_diff_matrix(N)
    D2 = D @ D
    D2_int = D2[1:N, 1:N]
    x_int = x[1:N]

    eigs, vecs = np.linalg.eig(D2_int)
    idx = np.argsort(np.real(eigs))
    eigs_sorted = np.real(eigs[idx])
    vecs_sorted = vecs[:, idx]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

    # Physical eigenvector: mode N//4 = 8
    mode_idx = N // 4 - 1  # 0-indexed, mode 8
    v_phys = vecs_sorted[:, mode_idx]
    # Normalise and fix sign
    v_phys = v_phys / np.max(np.abs(v_phys))
    if v_phys[len(v_phys) // 4] < 0:
        v_phys = -v_phys

    # Full vector including BCs
    v_full = np.zeros(N + 1)
    v_full[1:N] = v_phys

    ax1.plot(x, v_full, '-', color=NAVY, lw=1.8)
    ax1.plot(x_int, v_phys, 'o', color=NAVY, markersize=4,
             markerfacecolor='white', markeredgewidth=1.0)
    ax1.set_xlabel(r'$x$')
    ax1.set_ylabel(r'$v(x)$')
    ax1.set_title(fr'(a) Physical eigenvector (mode {N // 4},'
                  fr' $\lambda = {eigs_sorted[mode_idx]:.1f}$)')
    ax1.grid(True, alpha=0.3, linewidth=0.5)
    ax1.set_xlim(-1.05, 1.05)

    # Outlier eigenvector: largest eigenvalue (most negative)
    v_out = vecs_sorted[:, 0]  # most negative eigenvalue
    v_out = v_out / np.max(np.abs(v_out))
    # Fix sign so we see localisation clearly
    v_full_out = np.zeros(N + 1)
    v_full_out[1:N] = np.real(v_out)

    ax2.plot(x, v_full_out, '-', color=CORAL, lw=1.8)
    ax2.plot(x_int, np.real(v_out), 'o', color=CORAL, markersize=4,
             markerfacecolor='white', markeredgewidth=1.0)
    ax2.set_xlabel(r'$x$')
    ax2.set_ylabel(r'$v(x)$')
    ax2.set_title(fr'(b) Outlier eigenvector ($\lambda = {eigs_sorted[0]:.0f}$)')
    ax2.grid(True, alpha=0.3, linewidth=0.5)
    ax2.set_xlim(-1.05, 1.05)

    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / 'cheb_eigenvectors.pdf', bbox_inches='tight')
    fig.savefig(OUTPUT_DIR / 'cheb_eigenvectors.png', bbox_inches='tight')
    print(f"Figure saved to {OUTPUT_DIR / 'cheb_eigenvectors.pdf'}")


def figure_cfl_scaling():
    """Log-log plot of spectral radius of D^2 vs N."""
    N_values = [8, 12, 16, 20, 24, 32, 48, 64]
    rho_values = []

    for N in N_values:
        D, x = cheb_diff_matrix(N)
        D2 = D @ D
        D2_int = D2[1:N, 1:N]
        eigs = np.linalg.eigvals(D2_int)
        rho = np.max(np.abs(eigs))
        rho_values.append(rho)
        print(f"N = {N:3d}, rho(D2) = {rho:.2e}, "
              f"rho/N^4 = {rho / N ** 4:.4f}")

    rho_values = np.array(rho_values)
    N_arr = np.array(N_values, dtype=float)

    fig, ax = plt.subplots(figsize=(7, 5))

    ax.loglog(N_arr, rho_values, 'o-', color=NAVY, markersize=7, lw=1.5,
              label=r'$\rho(\tilde{D}_2)$')

    # Reference line: 0.048 * N^4
    N_ref = np.linspace(6, 80, 100)
    ax.loglog(N_ref, 0.048 * N_ref ** 4, '--', color=CORAL, lw=1.2,
              label=r'$0.048\, N^4$')

    ax.set_xlabel(r'$N$')
    ax.set_ylabel(r'$\rho(\tilde{D}_2)$')
    ax.set_title(r'Spectral radius of Chebyshev $D^2$ (Dirichlet BCs)')
    ax.legend(loc='upper left', fontsize=10)
    ax.grid(True, which='both', alpha=0.3, linewidth=0.5)

    # Slope annotation
    ax.text(30, 1e5, r'slope $= 4$', fontsize=11, color=CORAL,
            rotation=40, ha='center')

    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / 'cheb_cfl_scaling.pdf', bbox_inches='tight')
    fig.savefig(OUTPUT_DIR / 'cheb_cfl_scaling.png', bbox_inches='tight')
    print(f"Figure saved to {OUTPUT_DIR / 'cheb_cfl_scaling.pdf'}")


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    figure_eigenvalues()
    figure_eigenvectors()
    figure_cfl_scaling()
    plt.show()


if __name__ == '__main__':
    main()
