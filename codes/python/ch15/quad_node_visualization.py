#!/usr/bin/env python3
"""
quad_node_visualization.py
Chapter 15: Quadrature in Spectral Methods

Computational Etude 15.1: Quadrature Node Distributions.

Visualises three classical sets of quadrature nodes on [-1, 1] for n = 32
points, displayed on parallel horizontal lines:

    1. Newton--Cotes (equispaced) nodes: x_j = -1 + 2j/n, j = 0, ..., n.
    2. Gauss--Legendre nodes: zeros of P_{n+1}(x).
    3. Clenshaw--Curtis (Chebyshev) nodes: x_j = cos(j pi / n), j = 0, ..., n.

The figure highlights the fundamental difference in node clustering near
the endpoints: equispaced nodes remain uniformly distributed, whereas both
Gauss--Legendre and Clenshaw--Curtis nodes cluster near x = +-1. This
clustering is essential for stable high-order polynomial interpolation and
for the convergence of spectral quadrature rules.

Reference: L. N. Trefethen, "Is Gauss Quadrature Better than
           Clenshaw--Curtis?", SIAM Review 50(1):67--87, 2008, Figure 1.

Generates Figure 15.1: node_visualization.pdf

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
rcParams['text.usetex'] = False
rcParams['axes.linewidth'] = 0.8
rcParams['xtick.major.width'] = 0.8
rcParams['ytick.major.width'] = 0.8

# -----------------------------------------------------------------------------
# Book color scheme
# -----------------------------------------------------------------------------
NAVY = '#142D6E'
SKY = '#7896D2'
CORAL = '#E74C3C'
TEAL = '#16A085'
PURPLE = '#8E44AD'
ORANGE = '#E67E22'

# Output paths
SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / '..' / '..' / '..' / 'textbook' / 'figures' / 'ch15' / 'python'


def newton_cotes_nodes(n):
    """Return n+1 equispaced nodes on [-1, 1]."""
    return np.linspace(-1, 1, n + 1)


def gauss_legendre_nodes(n):
    """Return n+1 Gauss--Legendre nodes on [-1, 1]."""
    nodes, _ = np.polynomial.legendre.leggauss(n + 1)
    return nodes


def clenshaw_curtis_nodes(n):
    """Return n+1 Clenshaw--Curtis (Chebyshev extremal) nodes on [-1, 1]."""
    j = np.arange(n + 1)
    return np.cos(j * np.pi / n)


def main():
    """Generate the node visualisation figure."""
    n = 32  # Number of intervals (n+1 points)

    # Compute the three node sets
    nc_nodes = newton_cotes_nodes(n)
    gl_nodes = gauss_legendre_nodes(n)
    cc_nodes = clenshaw_curtis_nodes(n)

    # Vertical positions for the three rows
    y_nc = 3.0
    y_gl = 2.0
    y_cc = 1.0

    # Create figure
    fig, ax = plt.subplots(figsize=(8, 3))

    # Plot horizontal base lines
    for y in [y_nc, y_gl, y_cc]:
        ax.plot([-1, 1], [y, y], color='gray', linewidth=0.5, zorder=1)

    # Plot nodes as filled circles
    ms = 5  # Marker size

    ax.plot(nc_nodes, np.full_like(nc_nodes, y_nc), 'o',
            color=CORAL, markersize=ms, markeredgecolor=CORAL,
            markeredgewidth=0.5, zorder=2)

    ax.plot(gl_nodes, np.full_like(gl_nodes, y_gl), 'o',
            color=NAVY, markersize=ms, markeredgecolor=NAVY,
            markeredgewidth=0.5, zorder=2)

    ax.plot(cc_nodes, np.full_like(cc_nodes, y_cc), 'o',
            color=TEAL, markersize=ms, markeredgecolor=TEAL,
            markeredgewidth=0.5, zorder=2)

    # Labels on the left
    label_x = -1.25
    ax.text(label_x, y_nc, r'Newton$\mathrm{-}$Cotes',
            fontsize=11, va='center', ha='right', color=CORAL, fontweight='bold')
    ax.text(label_x, y_gl, r'Gauss$\mathrm{-}$Legendre',
            fontsize=11, va='center', ha='right', color=NAVY, fontweight='bold')
    ax.text(label_x, y_cc, r'Clenshaw$\mathrm{-}$Curtis',
            fontsize=11, va='center', ha='right', color=TEAL, fontweight='bold')

    # Mark the interval endpoints
    for y in [y_nc, y_gl, y_cc]:
        ax.plot([-1, 1], [y, y], '|', color='black', markersize=8,
                markeredgewidth=1.0, zorder=3)

    # Axis formatting
    ax.set_xlim(-1.8, 1.15)
    ax.set_ylim(0.4, 3.6)
    ax.set_xlabel(r'$x$')
    ax.set_xticks([-1, -0.5, 0, 0.5, 1])

    # Remove y-axis and top/right spines
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    # Title
    ax.set_title(r'Quadrature nodes on $[-1, 1]$ with $n = 32$', fontsize=12)

    plt.tight_layout()

    # Save output
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_DIR / 'node_visualization.pdf', bbox_inches='tight')
    fig.savefig(OUTPUT_DIR / 'node_visualization.png', bbox_inches='tight')
    print(f"Figure saved to {OUTPUT_DIR / 'node_visualization.pdf'}")

    plt.show()


if __name__ == '__main__':
    main()
