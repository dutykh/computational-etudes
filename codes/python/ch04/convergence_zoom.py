#!/usr/bin/env python3
"""
convergence_zoom.py

Zoomed view of Chebyshev interpolation convergence for the Runge function,
without the divergent equispaced curve.

This allows better visualization of the geometric convergence rate O(ρ^{-N})
where ρ ≈ 1.22 for the Runge function.

The Runge function f(x) = 1/(1 + 25x²) has poles at ±0.2i, giving:
    ρ = |0.2i + √((0.2i)² - 1)| ≈ 1.22

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
N_FINE = 2000  # Number of points for error evaluation

# Book color scheme
NAVY = '#142D6E'
SKY = '#7896D2'
CORAL = '#E74C3C'
TEAL = '#16A085'

# Output path
SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / '..' / '..' / '..' / 'textbook' / 'figures' / 'ch04' / 'python'
OUTPUT_FILE = OUTPUT_DIR / 'convergence_zoom.pdf'


# -----------------------------------------------------------------------------
# Functions
# -----------------------------------------------------------------------------
def runge(x):
    """The Runge function: f(x) = 1 / (1 + 25x²)"""
    return 1.0 / (1.0 + 25.0 * x**2)


def chebyshev_nodes(N):
    """Chebyshev-Gauss-Lobatto nodes on [-1, 1]."""
    j = np.arange(N + 1)
    return np.cos(j * np.pi / N)


def lagrange_interpolate(x_nodes, f_nodes, x_eval):
    """Evaluate the Lagrange interpolating polynomial at points x_eval."""
    n = len(x_nodes)
    x_eval = np.atleast_1d(x_eval)
    p_eval = np.zeros_like(x_eval)

    for k in range(n):
        L_k = np.ones_like(x_eval)
        for j in range(n):
            if j != k:
                L_k *= (x_eval - x_nodes[j]) / (x_nodes[k] - x_nodes[j])
        p_eval += f_nodes[k] * L_k

    return p_eval


def max_interpolation_error(node_func, N, n_eval=N_FINE):
    """Compute max |f - p_N| for given node distribution."""
    x_eval = np.linspace(-1, 1, n_eval)
    x_nodes = node_func(N)
    f_nodes = runge(x_nodes)
    p_interp = lagrange_interpolate(x_nodes, f_nodes, x_eval)
    return np.max(np.abs(runge(x_eval) - p_interp))


# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
def main():
    # Range of polynomial degrees (extended for better view)
    N_values = np.arange(2, 71, 2)

    # Compute errors for Chebyshev only
    errors_cheb = []

    for N in N_values:
        errors_cheb.append(max_interpolation_error(chebyshev_nodes, N))

    errors_cheb = np.array(errors_cheb)

    # Theoretical convergence parameter
    # The Runge function has poles at ±0.2i
    # ρ = |0.2i + √((0.2i)² - 1)|
    rho = np.abs(0.2j + np.sqrt(-1.04 + 0j))

    # Create figure
    fig, ax = plt.subplots(figsize=(7, 5))

    # Plot Chebyshev errors
    ax.semilogy(N_values, errors_cheb, 's-', color=TEAL, linewidth=1.5,
                markersize=5, label='Chebyshev interpolation')

    # Add theoretical convergence rate
    N_theory = np.linspace(10, 70, 100)
    # Fit the constant using data from N >= 20
    idx_fit = N_values >= 20
    C_fit = np.mean(errors_cheb[idx_fit] * rho**N_values[idx_fit])
    ax.semilogy(N_theory, C_fit * rho**(-N_theory), '--', color=NAVY,
                alpha=0.7, linewidth=2,
                label=rf'Theory: $C \cdot \rho^{{-N}}$, $\rho \approx {rho:.2f}$')

    # Add machine precision line
    ax.axhline(y=1e-15, color='#888888', linewidth=0.8, linestyle=':',
               label='Machine precision')

    # Styling
    ax.set_xlabel(r'Polynomial degree $N$')
    ax.set_ylabel(r'Maximum interpolation error $\|f - p_N\|_\infty$')
    ax.set_title('Chebyshev Interpolation Convergence (Runge Function)', fontsize=11)
    ax.set_xlim(0, 72)
    ax.set_ylim(1e-16, 1e1)

    ax.legend(loc='upper right', frameon=True, fancybox=False,
              edgecolor='none', facecolor='white', framealpha=0.9)

    # Add annotation for convergence rate
    ax.annotate(rf'$O(\rho^{{-N}})$ with $\rho \approx {rho:.2f}$',
                xy=(50, errors_cheb[-11]), xytext=(35, 1e-6),
                fontsize=11, color=TEAL,
                arrowprops=dict(arrowstyle='->', color=TEAL, lw=1.0))

    # Add formula box
    textstr = r'$f(x) = \frac{1}{1 + 25x^2}$' + '\n' + r'Poles at $z = \pm 0.2i$'
    props = dict(boxstyle='round', facecolor='white', alpha=0.9,
                 edgecolor=NAVY, linewidth=0.5)
    ax.text(0.05, 0.05, textstr, transform=ax.transAxes, fontsize=11,
            verticalalignment='bottom', bbox=props)

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

    # Print summary table
    print("\nChebyshev Convergence for Runge Function:")
    print("-" * 50)
    print(f"{'N':>4} {'Error':>15} {'Theoretical':>15}")
    print("-" * 50)
    for i, N in enumerate(N_values[::5]):  # Print every 5th value
        idx = i * 5
        err = errors_cheb[idx]
        theory = C_fit * rho**(-N)
        print(f"{N:>4d} {err:>15.2e} {theory:>15.2e}")
    print("-" * 50)

    # Analysis
    print("\nAnalysis:")
    print(f"  Runge function poles at z = ±0.2i")
    print(f"  Theoretical ρ = {rho:.6f}")

    # Estimate ρ from the data
    log_errors = np.log(errors_cheb[10:])  # Skip initial transient
    N_fit = N_values[10:]
    slope, intercept = np.polyfit(N_fit, log_errors, 1)
    rho_estimated = np.exp(-slope)
    print(f"  Estimated from data: ρ ≈ {rho_estimated:.6f}")
    print(f"  Error reduction per degree: factor of {rho_estimated:.3f}")

    plt.close(fig)


if __name__ == '__main__':
    main()
