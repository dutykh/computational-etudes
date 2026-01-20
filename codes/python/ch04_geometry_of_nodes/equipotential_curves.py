#!/usr/bin/env python3
"""
equipotential_curves.py

Visualizes the potential theory explanation for polynomial interpolation
convergence. Based on Trefethen's Spectral Methods in MATLAB, Program 10.

The potential function associated with a set of interpolation nodes is:
    φ(z) = -∫_{-1}^{1} μ(x) ln|z - x| dx

where μ(x) is the node density. The equipotential curves φ(z) = const
determine the region of convergence for polynomial interpolation.

For equispaced nodes:
    μ(x) = 1/2 (uniform density)
    Equipotentials are ellipse-like curves symmetric about [-1,1]

For Chebyshev density:
    μ(x) = 1 / (π√(1-x²))
    The potential simplifies to: φ(z) = log|z + √(z² - 1)| - log 2
    Equipotentials are Bernstein ellipses with foci at ±1

The Runge function has poles at z = ±0.2i. The largest equipotential that
doesn't enclose these poles determines the convergence region.

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
N_GRID = 500  # Grid resolution

# Book color scheme
NAVY = '#142D6E'
SKY = '#7896D2'
CORAL = '#E74C3C'
TEAL = '#16A085'

# Output path
SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / '..' / '..' / '..' / 'textbook' / 'figures' / 'ch04' / 'python'
OUTPUT_FILE = OUTPUT_DIR / 'equipotential_curves.pdf'


# -----------------------------------------------------------------------------
# Potential functions
# -----------------------------------------------------------------------------
def potential_uniform(z, N_nodes=64):
    """
    Potential for uniform (equispaced) node distribution.

    For uniform density μ(x) = 1/2 on [-1,1], the potential is computed
    numerically as the average of log distances to equispaced nodes.

    This approximates: φ(z) = -∫_{-1}^{1} (1/2) ln|z - x| dx
    """
    x_nodes = np.linspace(-1, 1, N_nodes)
    phi = np.zeros_like(z, dtype=float)

    for x_k in x_nodes:
        phi += np.log(np.abs(z - x_k) + 1e-15)

    return phi / N_nodes


def potential_chebyshev(z):
    """
    Potential for Chebyshev node distribution.

    For Chebyshev density μ(x) = 1/(π√(1-x²)), the potential has
    the closed form:
        φ(z) = log|z + √(z² - 1)| - log 2

    The equipotential curves are Bernstein ellipses with foci at ±1.
    An ellipse with parameter ρ > 1 has semi-axes:
        a = (ρ + ρ⁻¹)/2,  b = (ρ - ρ⁻¹)/2
    """
    # Compute z + √(z² - 1) carefully for numerical stability
    sqrt_term = np.sqrt(z**2 - 1 + 0j)
    # Choose the branch that gives |result| > 1 for z outside [-1,1]
    result = z + sqrt_term
    # For points where |z + sqrt| < 1, use the other branch
    mask = np.abs(result) < 1
    result[mask] = z[mask] - sqrt_term[mask]

    return np.log(np.abs(result)) - np.log(2)


# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
def main():
    # Create complex grid
    x = np.linspace(-1.5, 1.5, N_GRID)
    y = np.linspace(-1.0, 1.0, N_GRID)
    X, Y = np.meshgrid(x, y)
    Z = X + 1j * Y

    # Compute potentials
    phi_uniform = potential_uniform(Z)
    phi_chebyshev = potential_chebyshev(Z)

    # Contour levels
    levels_uniform = np.linspace(-0.6, 0.6, 13)
    levels_chebyshev = np.log([1.1, 1.2, 1.4, 1.6, 2.0, 2.5, 3.0]) - np.log(2)

    # Runge function singularities
    pole_y = 0.2  # Poles at ±0.2i

    # Create two-panel figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

    # Left panel: Uniform (equispaced) potential
    cs1 = ax1.contour(X, Y, phi_uniform, levels=levels_uniform,
                       colors=[SKY], linewidths=0.8)
    ax1.clabel(cs1, inline=True, fontsize=7, fmt='%.2f')

    # Draw the interval [-1,1]
    ax1.plot([-1, 1], [0, 0], color=NAVY, linewidth=2)
    ax1.plot([-1, 1], [0, 0], 'o', color=NAVY, markersize=6)

    # Mark singularities
    ax1.plot([0, 0], [pole_y, -pole_y], 'x', color=CORAL, markersize=10,
             markeredgewidth=2, label=f'Poles at $\\pm{pole_y}i$')

    ax1.set_xlabel(r'$\mathrm{Re}(z)$')
    ax1.set_ylabel(r'$\mathrm{Im}(z)$')
    ax1.set_title('Equispaced Nodes (Uniform Density)', fontsize=11)
    ax1.set_xlim(-1.5, 1.5)
    ax1.set_ylim(-1, 1)
    ax1.set_aspect('equal')
    ax1.legend(loc='upper right', fontsize=9)

    # Right panel: Chebyshev potential
    cs2 = ax2.contour(X, Y, phi_chebyshev.real, levels=levels_chebyshev,
                       colors=[TEAL], linewidths=0.8)
    # Label with ρ values
    rho_labels = {lev: f'ρ={np.exp(lev + np.log(2)):.1f}' for lev in levels_chebyshev}
    ax2.clabel(cs2, inline=True, fontsize=7,
               fmt=lambda x: rho_labels.get(x, f'{x:.2f}'))

    # Draw the interval [-1,1]
    ax2.plot([-1, 1], [0, 0], color=NAVY, linewidth=2)
    ax2.plot([-1, 1], [0, 0], 'o', color=NAVY, markersize=6)

    # Mark singularities
    ax2.plot([0, 0], [pole_y, -pole_y], 'x', color=CORAL, markersize=10,
             markeredgewidth=2, label=f'Poles at $\\pm{pole_y}i$')

    # Draw the critical Bernstein ellipse (passes through poles)
    # For z = 0.2i: ρ = |z + √(z²-1)| = |0.2i + √(-0.04-1)| = |0.2i + √(-1.04)|
    rho_critical = np.abs(0.2j + np.sqrt(-1.04 + 0j))
    theta = np.linspace(0, 2*np.pi, 200)
    a_crit = (rho_critical + 1/rho_critical) / 2
    b_crit = (rho_critical - 1/rho_critical) / 2
    ax2.plot(a_crit * np.cos(theta), b_crit * np.sin(theta), '--',
             color=CORAL, linewidth=1.5, label=f'Critical ellipse (ρ≈{rho_critical:.2f})')

    ax2.set_xlabel(r'$\mathrm{Re}(z)$')
    ax2.set_ylabel(r'$\mathrm{Im}(z)$')
    ax2.set_title('Chebyshev Nodes (Bernstein Ellipses)', fontsize=11)
    ax2.set_xlim(-1.5, 1.5)
    ax2.set_ylim(-1, 1)
    ax2.set_aspect('equal')
    ax2.legend(loc='upper right', fontsize=9)

    # Clean styling for both panels
    for ax in [ax1, ax2]:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_color('#444444')
        ax.spines['bottom'].set_color('#444444')
        ax.tick_params(colors='#444444')
        ax.axhline(y=0, color='#CCCCCC', linewidth=0.5, zorder=0)
        ax.axvline(x=0, color='#CCCCCC', linewidth=0.5, zorder=0)

    plt.tight_layout()

    # Save figure
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_FILE, bbox_inches='tight', pad_inches=0.05)
    fig.savefig(OUTPUT_FILE.with_suffix('.png'), bbox_inches='tight',
                pad_inches=0.05, dpi=300)

    print(f'Figure saved to: {OUTPUT_FILE.resolve()}')

    # Print analysis
    print("\nPotential Theory Analysis:")
    print("-" * 60)
    print(f"Runge function poles at z = ±{pole_y}i")
    print(f"\nChebyshev (Bernstein ellipse) analysis:")
    print(f"  Critical ρ = |{pole_y}i + √({pole_y}i² - 1)| = {rho_critical:.4f}")
    print(f"  Convergence rate: O(ρ^{{-N}}) = O({1/rho_critical:.4f}^N)")
    print(f"\nEquispaced analysis:")
    print(f"  The poles lie inside the smallest equipotential curve")
    print(f"  that encloses [-1,1], so interpolation diverges.")
    print("-" * 60)

    plt.close(fig)


if __name__ == '__main__':
    main()
