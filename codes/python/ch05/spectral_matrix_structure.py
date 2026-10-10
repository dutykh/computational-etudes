#!/usr/bin/env python3
"""
spectral_matrix_structure.py

Visualizes the structure of the periodic spectral differentiation matrix,
showing its Toeplitz/circulant structure and the characteristic cotangent
decay of the off-diagonal entries.

Key Properties Illustrated:
    - Skew-symmetry: D^T = -D
    - Toeplitz structure: D[i,j] depends only on (i-j)
    - Circulant (due to periodicity)
    - Dense: all off-diagonal entries nonzero
    - Diagonal entries are zero

Author: Dr. Denys Dutykh
        Mathematics Department
        Khalifa University of Science and Technology
        Abu Dhabi, UAE

Part of the book "Computational Études: A Spectral Approach"
https://github.com/dutykh/computational-etudes
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.colors import TwoSlopeNorm
from pathlib import Path

# Import our modules
import sys
sys.path.insert(0, str(Path(__file__).parent))
from spectral_matrix_periodic import spectral_diff_periodic

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

SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / '..' / '..' / '..' / 'textbook' / 'figures' / 'ch05' / 'python'
OUTPUT_FILE = OUTPUT_DIR / 'spectral_matrix_structure.pdf'


# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
def main():
    N = 16  # Number of grid points

    D, x = spectral_diff_periodic(N)

    # Create figure with two panels
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))

    # ==========================================================================
    # Left panel: Heatmap of the matrix
    # ==========================================================================
    # Use diverging colormap centered at 0
    vmax = np.max(np.abs(D))
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    im = ax1.imshow(D, cmap='RdBu_r', norm=norm, interpolation='nearest',
                     aspect='equal')

    # Add colorbar
    cbar = plt.colorbar(im, ax=ax1, shrink=0.8)
    cbar.set_label('Matrix entry value', fontsize=10)

    ax1.set_xlabel('Column index $j$', fontsize=11)
    ax1.set_ylabel('Row index $i$', fontsize=11)
    ax1.set_title(f'Spectral Matrix $D$ ($N = {N}$)', fontsize=12)

    # Add annotation about skew-symmetry
    ax1.text(0.02, 0.98, r'$D^T = -D$', transform=ax1.transAxes,
             fontsize=11, fontweight='bold', color=NAVY,
             verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    # Add ticks
    ax1.set_xticks([0, N//2, N-1])
    ax1.set_xticklabels(['0', f'{N//2}', f'{N-1}'])
    ax1.set_yticks([0, N//2, N-1])
    ax1.set_yticklabels(['0', f'{N//2}', f'{N-1}'])

    # ==========================================================================
    # Right panel: First row profile (showing the cotangent decay)
    # ==========================================================================
    row = D[0, :]  # First row
    j_indices = np.arange(N)

    # Shift indices to show the pattern centered around 0
    j_shifted = np.where(j_indices <= N//2, j_indices, j_indices - N)

    # Sort for nicer plotting
    sort_idx = np.argsort(j_shifted)
    j_plot = j_shifted[sort_idx]
    row_plot = row[sort_idx]

    ax2.stem(j_plot, row_plot, linefmt=NAVY, markerfmt='o', basefmt='gray')

    # Add theoretical cotangent curve
    j_theory = np.linspace(-N/2 + 0.01, N/2 - 0.01, 500)
    j_theory = j_theory[np.abs(j_theory) > 0.1]  # Avoid singularity at 0
    D_theory = 0.5 * ((-1) ** np.round(j_theory)) / np.tan(j_theory * np.pi / N)
    ax2.plot(j_theory, D_theory, '--', color=CORAL, alpha=0.7, linewidth=1.5,
             label=r'$\frac{1}{2}(-1)^k \cot\left(\frac{k\pi}{N}\right)$')

    ax2.axhline(y=0, color='gray', linewidth=0.5, linestyle='-')
    ax2.axvline(x=0, color='gray', linewidth=0.5, linestyle='-')

    ax2.set_xlabel(r'Column offset $k = j - i$', fontsize=11)
    ax2.set_ylabel(r'Matrix entry $D_{0,k}$', fontsize=11)
    ax2.set_title('First Row Profile (Toeplitz Structure)', fontsize=12)
    ax2.set_xlim(-N/2 - 0.5, N/2 + 0.5)

    ax2.legend(loc='upper right', fontsize=9)

    # Clean styling
    for ax in [ax1, ax2]:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_color('#444444')
        ax.spines['bottom'].set_color('#444444')

    ax2.yaxis.grid(True, linestyle='-', alpha=0.2, color='#888888')
    ax2.xaxis.grid(True, linestyle='-', alpha=0.2, color='#888888')
    ax2.set_axisbelow(True)

    plt.tight_layout()

    # Save figure
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_FILE, bbox_inches='tight', pad_inches=0.05)
    fig.savefig(OUTPUT_FILE.with_suffix('.png'), bbox_inches='tight',
                pad_inches=0.05, dpi=300)

    print(f'Figure saved to: {OUTPUT_FILE.resolve()}')

    # Print matrix properties
    print(f"\nSpectral Matrix Properties (N = {N}):")
    print("-" * 50)
    print(f"  Skew-symmetry error ||D + D^T||_max: {np.max(np.abs(D + D.T)):.2e}")
    print(f"  Diagonal sum: {np.sum(np.diag(D)):.2e} (should be 0)")
    print(f"  Max entry: {np.max(D):.4f}")
    print(f"  Min entry: {np.min(D):.4f}")
    print("-" * 50)

    plt.close(fig)


if __name__ == '__main__':
    main()
