#!/usr/bin/env python3
"""
periodic_cardinal_functions.py

Visualizes periodic cardinal functions (discrete Dirichlet kernel) for equispaced
nodes on a periodic domain.

The periodic cardinal function is:
    phi_j(x) = sin(N(x - x_j)/2) / (N sin((x - x_j)/2))

These functions satisfy phi_j(x_k) = delta_{jk} (Kronecker delta), making them
the "cardinal functions" for trigonometric interpolation on periodic domains.
The interpolant is:
    p(x) = sum_{j=0}^{N-1} u_j * phi_j(x)

Unlike Lagrange basis polynomials on non-periodic domains, the periodic cardinal
functions are bounded and exhibit damped oscillations away from their peak,
reflecting the sin(x)/x-like structure of the discrete Dirichlet kernel.

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
rcParams['text.usetex'] = False  # Set to True if LaTeX is available
rcParams['axes.linewidth'] = 0.8
rcParams['xtick.major.width'] = 0.8
rcParams['ytick.major.width'] = 0.8

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
N = 16  # Number of grid points (should be even)
N_FINE = 1000  # Number of points for smooth plotting

# Book color scheme
NAVY = '#142D6E'
SKY = '#7896D2'
CORAL = '#E74C3C'
TEAL = '#16A085'
PURPLE = '#8E44AD'
ORANGE = '#E67E22'

# Output path
SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / '..' / '..' / '..' / 'textbook' / 'figures' / 'ch05' / 'python'
OUTPUT_FILE = OUTPUT_DIR / 'periodic_cardinal_functions.pdf'


# -----------------------------------------------------------------------------
# Periodic cardinal function computation
# -----------------------------------------------------------------------------
def periodic_cardinal(x, x_j, N):
    """
    Compute the periodic cardinal function phi_j(x) centered at x_j.

    phi_j(x) = sin(N(x - x_j)/2) / (N sin((x - x_j)/2))

    This is the discrete Dirichlet kernel, satisfying phi_j(x_k) = delta_{jk}
    for equispaced nodes x_k = 2*pi*k/N on [0, 2*pi).

    Parameters
    ----------
    x : array_like
        Points at which to evaluate the cardinal function
    x_j : float
        Center point (node location)
    N : int
        Number of grid points

    Returns
    -------
    phi : ndarray
        Values of phi_j at x
    """
    x = np.atleast_1d(x)
    theta = (x - x_j) / 2.0

    # Handle the singularity at theta = 0 (x = x_j)
    # Use L'Hopital's rule: lim_{theta->0} sin(N*theta)/(N*sin(theta)) = 1
    phi = np.zeros_like(x, dtype=float)

    # Find points where theta is effectively zero (at the node x_j)
    small = np.abs(np.sin(theta)) < 1e-14

    # For non-singular points
    phi[~small] = np.sin(N * theta[~small]) / (N * np.sin(theta[~small]))

    # For singular points (at x_j and periodic copies)
    phi[small] = 1.0

    return phi


# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
def main():
    # Grid points on [0, 2*pi)
    h = 2 * np.pi / N
    x_nodes = h * np.arange(N)

    # Fine grid for plotting
    x_fine = np.linspace(0, 2 * np.pi, N_FINE)

    # Create figure
    fig, ax = plt.subplots(figsize=(9, 4.5))

    # Indices of cardinal functions to plot (spread across the domain)
    indices_to_plot = [0, 4, 8, 12]
    colors = [CORAL, TEAL, SKY, PURPLE]
    linestyles = ['-', '-', '-', '-']

    # Plot cardinal functions
    for j, color in zip(indices_to_plot, colors):
        phi_j = periodic_cardinal(x_fine, x_nodes[j], N)
        ax.plot(x_fine, phi_j, color=color, linewidth=1.5,
                label=rf'$\phi_{{{j}}}(x)$')

    # Mark all nodes on x-axis
    ax.plot(x_nodes, np.zeros_like(x_nodes), 'o', color=NAVY, markersize=5,
            zorder=5, label='Nodes')

    # Highlight the peak values at each plotted cardinal function's node
    for j, color in zip(indices_to_plot, colors):
        ax.plot(x_nodes[j], 1.0, 'o', color=color, markersize=8,
                markeredgecolor='white', markeredgewidth=1.5, zorder=6)

    # Reference lines
    ax.axhline(y=0, color='#888888', linewidth=0.5, linestyle='-')
    ax.axhline(y=1, color='#888888', linewidth=0.5, linestyle='--', alpha=0.5)

    # Axis labels and title
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$\phi_j(x)$')
    ax.set_title(rf'Periodic Cardinal Functions ($N = {N}$)', fontsize=11)

    # Set x-axis ticks at multiples of pi/2
    ax.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi])
    ax.set_xticklabels([r'$0$', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'])

    ax.set_xlim(0, 2 * np.pi)
    ax.set_ylim(-0.25, 1.15)

    # Legend
    ax.legend(loc='upper right', fontsize=9, ncol=2)

    # Add annotation about cardinal property
    ax.annotate(r'$\phi_j(x_k) = \delta_{jk}$',
                xy=(x_nodes[4], 1.0), xytext=(x_nodes[4] + 0.8, 0.85),
                fontsize=10, color='gray',
                arrowprops=dict(arrowstyle='->', color='gray', lw=0.8))

    # Clean styling
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color('#444444')
    ax.spines['bottom'].set_color('#444444')
    ax.tick_params(colors='#444444')
    ax.yaxis.grid(True, linestyle='-', alpha=0.2, color='#888888')
    ax.xaxis.grid(True, linestyle='-', alpha=0.2, color='#888888')
    ax.set_axisbelow(True)

    plt.tight_layout()

    # Save figure
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_FILE, bbox_inches='tight', pad_inches=0.05)
    fig.savefig(OUTPUT_FILE.with_suffix('.png'), bbox_inches='tight',
                pad_inches=0.05, dpi=300)

    print(f'Figure saved to: {OUTPUT_FILE.resolve()}')

    # Verify cardinal property
    print(f"\nVerification of cardinal property for N = {N}:")
    print("-" * 50)
    print(f"{'j':>4} {'phi_j(x_j)':>12} {'max |phi_j(x_k)|, k!=j':>25}")
    print("-" * 50)
    for j in range(N):
        phi_at_nodes = periodic_cardinal(x_nodes, x_nodes[j], N)
        phi_at_own = phi_at_nodes[j]
        phi_at_others = np.delete(phi_at_nodes, j)
        max_at_others = np.max(np.abs(phi_at_others))
        print(f"{j:>4d} {phi_at_own:>12.6f} {max_at_others:>25.2e}")
    print("-" * 50)

    # Print maximum oscillation amplitude
    print(f"\nMaximum |phi_j(x)| away from x_j:")
    for j in indices_to_plot:
        phi_j = periodic_cardinal(x_fine, x_nodes[j], N)
        # Find points away from the node
        away_from_node = np.abs(x_fine - x_nodes[j]) > 0.1
        if np.any(away_from_node):
            max_away = np.max(np.abs(phi_j[away_from_node]))
            print(f"  phi_{j}: max |phi_{j}(x)| = {max_away:.4f}")

    plt.close(fig)


if __name__ == '__main__':
    main()
