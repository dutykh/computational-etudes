#!/usr/bin/env python3
"""
aliasing_demo.py

Visualization of the aliasing phenomenon (Theorem 2: Poisson summation).

When a function is sampled at N equispaced points, frequencies that differ
by multiples of N become indistinguishable. High-frequency content "folds"
onto low frequencies, corrupting the discrete Fourier transform.

The aliasing formula:
    f_tilde_k = sum_{j=-infty}^{infty} f_hat_{k + jN}

This demo shows:
    1. A function with high-frequency content sampled coarsely
    2. The "folding" diagram showing how frequencies overlap
    3. Comparison of true spectrum vs aliased spectrum

For smooth functions, aliasing is harmless because high-frequency
coefficients are negligible. For non-smooth functions, aliasing
introduces significant errors.

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
from matplotlib.patches import FancyArrowPatch
from pathlib import Path

# -----------------------------------------------------------------------------
# Publication-quality matplotlib configuration
# -----------------------------------------------------------------------------
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Computer Modern Roman', 'CMU Serif', 'DejaVu Serif']
rcParams['mathtext.fontset'] = 'cm'
rcParams['axes.labelsize'] = 11
rcParams['axes.titlesize'] = 12
rcParams['xtick.labelsize'] = 10
rcParams['ytick.labelsize'] = 10
rcParams['legend.fontsize'] = 10
rcParams['figure.dpi'] = 150
rcParams['savefig.dpi'] = 300
rcParams['text.usetex'] = False
rcParams['axes.linewidth'] = 0.8

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
NAVY = '#142D6E'
SKY = '#7896D2'
CORAL = '#E74C3C'
TEAL = '#16A085'
PURPLE = '#9B59B6'

SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / '..' / '..' / '..' / 'textbook' / 'figures' / 'ch06' / 'python'
OUTPUT_FILE = OUTPUT_DIR / 'aliasing_visualization.pdf'


def main():
    """Create a three-panel figure demonstrating aliasing."""

    fig, axes = plt.subplots(1, 3, figsize=(12, 4))

    # -------------------------------------------------------------------------
    # Panel 1: Physical space - coarse sampling misses oscillations
    # -------------------------------------------------------------------------
    ax1 = axes[0]

    # A function with both low and high frequency content
    def f(x):
        return np.sin(x) + 0.5 * np.sin(10 * x)

    x_fine = np.linspace(0, 2 * np.pi, 500)
    N_coarse = 12  # Coarse sampling

    x_coarse = np.linspace(0, 2 * np.pi, N_coarse, endpoint=False)

    # Plot fine function
    ax1.plot(x_fine, f(x_fine), '-', color=NAVY, linewidth=1.5,
             label='True function')

    # Plot coarse samples
    ax1.plot(x_coarse, f(x_coarse), 'o', color=CORAL, markersize=8,
             label=f'Sampled ($N={N_coarse}$)')

    # Interpolant through coarse samples (trigonometric)
    # The aliased version sees sin(10x) as sin(10x - 12x) = sin(-2x)
    # So the perceived function is sin(x) + 0.5*sin(-2x) = sin(x) - 0.5*sin(2x)
    def f_aliased(x):
        return np.sin(x) - 0.5 * np.sin(2 * x)

    ax1.plot(x_fine, f_aliased(x_fine), '--', color=TEAL, linewidth=1.5,
             label='Aliased interpolant')

    ax1.set_xlabel(r'$x$')
    ax1.set_ylabel(r'$f(x)$')
    ax1.set_title('Coarse Sampling Misses High Frequencies')
    ax1.legend(loc='upper right', fontsize=9)
    ax1.set_xlim(0, 2 * np.pi)
    ax1.set_ylim(-1.8, 1.8)

    # -------------------------------------------------------------------------
    # Panel 2: Frequency folding schematic
    # -------------------------------------------------------------------------
    ax2 = axes[1]

    N = 12  # Number of sample points
    nyquist = N // 2  # Nyquist frequency

    # Draw the frequency axis
    k_range = np.arange(-20, 21)
    ax2.axhline(y=0, color='gray', linewidth=0.5)

    # Mark the Nyquist boundaries
    ax2.axvline(x=-nyquist, color='gray', linestyle='--', linewidth=1)
    ax2.axvline(x=nyquist, color='gray', linestyle='--', linewidth=1)

    # Shade the resolved region
    ax2.axvspan(-nyquist, nyquist, alpha=0.2, color=TEAL)
    ax2.text(0, 0.8, 'Resolved\nfrequencies', ha='center', va='center',
             fontsize=10, color=TEAL)

    # Show folding arrows for a specific frequency (k=10)
    k_original = 10

    # k=10 folds to k=10-12=-2
    k_aliased = k_original - N
    ax2.annotate('', xy=(k_aliased, 0.3), xytext=(k_original, 0.3),
                 arrowprops=dict(arrowstyle='->', color=CORAL, lw=2))
    ax2.plot(k_original, 0.3, 'o', color=CORAL, markersize=8)
    ax2.text(k_original, 0.45, f'$k={k_original}$', ha='center', fontsize=10,
             color=CORAL)
    ax2.text((k_original + k_aliased) / 2, 0.5, 'folds to', ha='center',
             fontsize=9, color=CORAL)

    # k=-10 folds to k=-10+12=2
    k_original2 = -10
    k_aliased2 = k_original2 + N
    ax2.annotate('', xy=(k_aliased2, -0.3), xytext=(k_original2, -0.3),
                 arrowprops=dict(arrowstyle='->', color=PURPLE, lw=2))
    ax2.plot(k_original2, -0.3, 'o', color=PURPLE, markersize=8)
    ax2.text(k_original2, -0.45, f'$k={k_original2}$', ha='center', fontsize=10,
             color=PURPLE)

    # Labels for Nyquist boundaries
    ax2.text(-nyquist, -0.7, f'$-N/2$\n$={-nyquist}$', ha='center',
             fontsize=10, color='gray')
    ax2.text(nyquist, -0.7, f'$N/2$\n$={nyquist}$', ha='center',
             fontsize=10, color='gray')

    ax2.set_xlabel(r'Wavenumber $k$')
    ax2.set_title(f'Frequency Folding ($N={N}$)')
    ax2.set_xlim(-18, 18)
    ax2.set_ylim(-1, 1)
    ax2.set_yticks([])

    # -------------------------------------------------------------------------
    # Panel 3: True vs aliased spectrum
    # -------------------------------------------------------------------------
    ax3 = axes[2]

    # A function with significant high-frequency content
    def g(x):
        return 1.0 / (1.5 + np.cos(x))

    # Compute "true" spectrum with many points
    N_fine = 128
    x_fine = np.linspace(0, 2 * np.pi, N_fine, endpoint=False)
    g_fine = g(x_fine)
    g_hat_fine = np.fft.fft(g_fine) / N_fine

    # Compute aliased spectrum with few points
    N_coarse = 16
    x_coarse = np.linspace(0, 2 * np.pi, N_coarse, endpoint=False)
    g_coarse = g(x_coarse)
    g_hat_coarse = np.fft.fft(g_coarse) / N_coarse

    # Plot true spectrum (positive frequencies only)
    k_fine = np.arange(N_fine // 2 + 1)
    ax3.semilogy(k_fine, np.abs(g_hat_fine[:N_fine // 2 + 1]), 'o-',
                 color=NAVY, markersize=3, linewidth=1,
                 label='True spectrum ($N=128$)')

    # Plot aliased spectrum
    k_coarse = np.arange(N_coarse // 2 + 1)
    ax3.semilogy(k_coarse, np.abs(g_hat_coarse[:N_coarse // 2 + 1]), 's-',
                 color=CORAL, markersize=6, linewidth=1.5,
                 label=f'Aliased spectrum ($N={N_coarse}$)')

    # Mark Nyquist frequency
    ax3.axvline(x=N_coarse // 2, color='gray', linestyle='--', linewidth=1)
    ax3.text(N_coarse // 2 + 1, 1e-1, f'Nyquist\n($N/2={N_coarse // 2}$)',
             fontsize=9, color='gray')

    ax3.set_xlabel(r'Wavenumber $k$')
    ax3.set_ylabel(r'$|\hat{g}_k|$')
    ax3.set_title('True vs. Aliased Spectrum')
    ax3.legend(loc='upper right', fontsize=9)
    ax3.set_xlim(-1, 50)
    ax3.set_ylim(1e-10, 1e1)

    # Clean styling for all panels
    for ax in axes:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_color('#444444')
        ax.spines['bottom'].set_color('#444444')
        ax.tick_params(colors='#444444')

    plt.tight_layout()

    # Save figure
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_FILE, bbox_inches='tight', pad_inches=0.05)
    fig.savefig(OUTPUT_FILE.with_suffix('.png'), bbox_inches='tight',
                pad_inches=0.05, dpi=300)

    print(f'Figure saved to: {OUTPUT_FILE.resolve()}')

    plt.close(fig)


if __name__ == '__main__':
    main()
