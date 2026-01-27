#!/usr/bin/env python3
"""
fd_matrix_bandwidth.py

Visualizes the sparsity patterns of finite difference and spectral
differentiation matrices, showing how the bandwidth expands as the
order of accuracy increases.

Key Message:
    - 2nd order FD: Tridiagonal (bandwidth 3)
    - 4th order FD: Pentadiagonal (bandwidth 5)
    - Spectral: Dense (bandwidth N)

    This illustrates that spectral methods are the limit case of finite
    differences as the stencil width extends to the full domain.

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
from pathlib import Path

# Import our modules
import sys
sys.path.insert(0, str(Path(__file__).parent))
from fdweights import fdweights
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
# Book color scheme
NAVY = '#142D6E'
SKY = '#7896D2'
CORAL = '#E74C3C'
TEAL = '#16A085'

# Output path
SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / '..' / '..' / '..' / 'textbook' / 'figures' / 'ch05' / 'python'
OUTPUT_FILE = OUTPUT_DIR / 'fd_matrix_bandwidth.pdf'


# -----------------------------------------------------------------------------
# Finite difference matrix construction
# -----------------------------------------------------------------------------
def fd_diff_periodic(N, order):
    """
    Construct a periodic finite difference differentiation matrix.

    Parameters
    ----------
    N : int
        Number of grid points.
    order : int
        Order of accuracy (2, 4, 6, ...).

    Returns
    -------
    D : ndarray of shape (N, N)
        The finite difference differentiation matrix.
    """
    h = 2 * np.pi / N

    # Determine stencil size
    stencil_half = order // 2
    stencil_size = order + 1

    # Create stencil nodes centered at 0
    stencil_nodes = h * np.arange(-stencil_half, stencil_half + 1)

    # Compute FD weights for first derivative at center
    weights = fdweights(0, stencil_nodes, 1)

    # Build circulant matrix
    D = np.zeros((N, N))
    for i in range(N):
        for k, w in enumerate(weights):
            j = (i + k - stencil_half) % N  # Periodic wrap-around
            D[i, j] += w

    return D


# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
def main():
    N = 20  # Number of grid points

    # Create the three matrices
    D_fd2 = fd_diff_periodic(N, order=2)   # 2nd order: tridiagonal
    D_fd4 = fd_diff_periodic(N, order=4)   # 4th order: pentadiagonal
    D_spec, _ = spectral_diff_periodic(N)  # Spectral: dense

    # Create figure with three panels
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.5))

    matrices = [D_fd2, D_fd4, D_spec]
    titles = ['2nd Order FD', '4th Order FD', 'Spectral']
    bandwidths = ['Bandwidth = 3', 'Bandwidth = 5', f'Dense ($N \\times N$)']

    for ax, D, title, bw in zip(axes, matrices, titles, bandwidths):
        # Create binary sparsity pattern
        # Mark entries as nonzero if |D[i,j]| > threshold
        threshold = 1e-14
        sparsity = np.abs(D) > threshold

        # Plot using imshow
        ax.imshow(sparsity, cmap='Blues', interpolation='nearest',
                  aspect='equal', vmin=0, vmax=1)

        ax.set_title(title, fontsize=12, fontweight='bold')
        # Combine column label and bandwidth info
        ax.set_xlabel(f'Column index $j$\n{bw}', fontsize=10)

        # Only first panel gets y-label
        if ax == axes[0]:
            ax.set_ylabel('Row index $i$', fontsize=10)

        # Reduce tick density
        ax.set_xticks([0, N//2, N-1])
        ax.set_xticklabels(['0', f'{N//2}', f'{N-1}'])
        ax.set_yticks([0, N//2, N-1])
        ax.set_yticklabels(['0', f'{N//2}', f'{N-1}'])

    # Add overall title
    fig.suptitle('Finite Difference Matrix Sparsity: From Local to Global',
                 fontsize=13, y=1.02)

    plt.tight_layout()

    # Save figure
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_FILE, bbox_inches='tight', pad_inches=0.1)
    fig.savefig(OUTPUT_FILE.with_suffix('.png'), bbox_inches='tight',
                pad_inches=0.1, dpi=300)

    print(f'Figure saved to: {OUTPUT_FILE.resolve()}')

    # Print matrix statistics
    print("\nMatrix Statistics (N = 20):")
    print("-" * 50)
    for D, title in zip(matrices, titles):
        nnz = np.sum(np.abs(D) > 1e-14)
        density = nnz / D.size * 100
        print(f"{title:15s}: {nnz:4d} nonzeros ({density:5.1f}% dense)")
    print("-" * 50)

    plt.close(fig)


if __name__ == '__main__':
    main()
