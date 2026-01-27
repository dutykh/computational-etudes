#!/usr/bin/env python3
"""
cheb_convergence.py

Demonstrates spectral convergence rates for four functions of increasing
smoothness. Shows the fundamental relationship between function regularity
and convergence rate of spectral differentiation.

Test functions (different from Trefethen's standard examples):
1. |x|^(5/2) - Third derivative of bounded variation (algebraic: ~N^{-2.5})
2. exp(-1/(1-x^2)) - Smooth bump function, C^∞ but not analytic (superalgebraic)
3. tanh(5x) - Analytic in [-1,1], poles at ±iπ/10 (exponential)
4. x^8 - Polynomial of degree 8 (exact for N ≥ 8)

This script generates Figure 6.5 and Table 6.2 for Chapter 6.

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

# Import Chebyshev matrix function
import sys
sys.path.insert(0, str(Path(__file__).parent))
from cheb_matrix import cheb_matrix

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
OUTPUT_DIR = SCRIPT_DIR / '..' / '..' / '..' / 'textbook' / 'figures' / 'ch06' / 'python'


# -----------------------------------------------------------------------------
# Test functions and their exact derivatives
# -----------------------------------------------------------------------------

# Function 1: |x|^(5/2) - C^2, 3rd derivative in BV
def f1(x):
    return np.abs(x) ** 2.5

def f1_prime(x):
    # d/dx |x|^(5/2) = (5/2)|x|^(3/2) * sign(x)
    return 2.5 * np.abs(x) ** 1.5 * np.sign(x)


# Function 2: exp(-1/(1-x^2)) - Smooth bump, C^∞ but not analytic at boundaries
def f2(x):
    result = np.zeros_like(x)
    mask = np.abs(x) < 1
    result[mask] = np.exp(-1.0 / (1.0 - x[mask]**2))
    return result

def f2_prime(x):
    result = np.zeros_like(x)
    mask = np.abs(x) < 1
    inner = 1.0 - x[mask]**2
    result[mask] = np.exp(-1.0 / inner) * (-2.0 * x[mask]) / inner**2
    return result


# Function 3: tanh(5x) - Analytic, poles at x = ±iπ/10
def f3(x):
    return np.tanh(5.0 * x)

def f3_prime(x):
    return 5.0 / np.cosh(5.0 * x)**2


# Function 4: x^8 - Polynomial degree 8
def f4(x):
    return x ** 8

def f4_prime(x):
    return 8.0 * x ** 7


# List of test functions
TEST_FUNCTIONS = [
    {'name': r'$|x|^{5/2}$', 'f': f1, 'fp': f1_prime,
     'smoothness': r'$C^2$, 3rd deriv in BV', 'rate': r'$O(N^{-2.5})$',
     'color': CORAL},
    {'name': r'$e^{-1/(1-x^2)}$', 'f': f2, 'fp': f2_prime,
     'smoothness': r'$C^\infty$, not analytic', 'rate': r'faster than any $N^{-k}$',
     'color': TEAL},
    {'name': r'$\tanh(5x)$', 'f': f3, 'fp': f3_prime,
     'smoothness': 'Analytic', 'rate': r'$O(\rho^{-N})$',
     'color': PURPLE},
    {'name': r'$x^8$', 'f': f4, 'fp': f4_prime,
     'smoothness': 'Polynomial', 'rate': r'Exact for $N \geq 8$',
     'color': ORANGE},
]


def main():
    """Create convergence waterfall figure."""

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Range of N values
    N_values = np.array([4, 6, 8, 10, 12, 16, 20, 24, 32, 40, 48, 64, 80, 96, 128])

    # Create 2x2 figure
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))

    # Store errors for table
    error_data = {tf['name']: [] for tf in TEST_FUNCTIONS}

    for idx, (ax, tf) in enumerate(zip(axes.flat, TEST_FUNCTIONS)):
        errors = []

        for N in N_values:
            D, x = cheb_matrix(N)
            v = tf['f'](x)
            w = D @ v
            w_exact = tf['fp'](x)

            # Compute max error (avoiding boundary issues for bump function)
            if idx == 1:  # Bump function
                interior = (x > -0.99) & (x < 0.99)
                if np.any(interior):
                    error = np.max(np.abs(w[interior] - w_exact[interior]))
                else:
                    error = np.max(np.abs(w - w_exact))
            else:
                error = np.max(np.abs(w - w_exact))

            errors.append(error if error > 1e-16 else 1e-16)
            error_data[tf['name']].append((N, error))

        errors = np.array(errors)

        # Plot convergence
        ax.semilogy(N_values, errors, 'o-', color=tf['color'], linewidth=1.5,
                    markersize=5, markeredgecolor='white', markeredgewidth=0.5)

        # Add reference lines for theoretical rates
        if idx == 0:  # |x|^(5/2) - algebraic O(N^-2.5)
            N_ref = N_values[N_values >= 16]
            err_ref = errors[N_values >= 16]
            if len(err_ref) > 0:
                C = err_ref[0] * N_ref[0]**2.5
                ax.semilogy(N_ref, C * N_ref**(-2.5), '--', color='gray',
                            linewidth=1, alpha=0.7, label=r'$O(N^{-2.5})$')
        elif idx == 2:  # tanh(5x) - exponential
            # rho = distance to nearest singularity = pi/10 ≈ 0.314
            # So convergence rate ~ (1 + pi/10)^(-N) is too slow
            # Actually for tanh, poles at ±iπ/10 give exponential convergence
            pass
        elif idx == 3:  # x^8 - exact for N >= 8
            ax.axhline(1e-14, color='gray', linewidth=1, linestyle='--', alpha=0.7)

        ax.set_xlabel(r'$N$')
        ax.set_ylabel('Max Error')
        ax.set_title(f'{tf["name"]}: {tf["smoothness"]}', fontsize=11)
        ax.set_xlim(0, 135)
        ax.set_ylim(1e-16, 1e2)
        ax.grid(True, which='both', alpha=0.3)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # Add convergence rate annotation
        ax.text(0.95, 0.95, tf['rate'], transform=ax.transAxes,
                fontsize=9, verticalalignment='top', horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    # Main title
    fig.suptitle('Spectral Convergence: Four Functions of Increasing Smoothness',
                 fontsize=13, y=0.98)

    plt.tight_layout(rect=[0, 0, 1, 0.96])

    # Save figure
    output_file = OUTPUT_DIR / 'convergence_waterfall.pdf'
    fig.savefig(output_file, bbox_inches='tight', pad_inches=0.05)
    fig.savefig(output_file.with_suffix('.png'), bbox_inches='tight',
                pad_inches=0.05, dpi=300)

    print(f'Figure saved to: {output_file.resolve()}')
    plt.close(fig)

    # Generate convergence table
    print("\n" + "=" * 80)
    print("Convergence Table: Chebyshev Spectral Differentiation")
    print("=" * 80)

    # Table header
    print(f"{'N':>6}", end='')
    for tf in TEST_FUNCTIONS:
        name_short = tf['name'].replace('$', '').replace('\\', '').replace('{', '').replace('}', '')
        print(f"  {name_short:>14}", end='')
    print()
    print("-" * 80)

    # Table rows
    for i, N in enumerate(N_values):
        print(f"{N:>6d}", end='')
        for tf in TEST_FUNCTIONS:
            _, err = error_data[tf['name']][i]
            if err < 1e-14:
                print(f"  {'< 1e-14':>14}", end='')
            else:
                print(f"  {err:>14.2e}", end='')
        print()

    print("=" * 80)

    # Summary
    print("\nSummary:")
    print("-" * 80)
    print(f"{'Function':<25} {'Smoothness':<25} {'Expected Rate':<20}")
    print("-" * 80)
    for tf in TEST_FUNCTIONS:
        name_clean = tf['name'].replace('$', '')
        smooth_clean = tf['smoothness'].replace('$', '').replace('\\', '')
        rate_clean = tf['rate'].replace('$', '').replace('\\', '')
        print(f"{name_clean:<25} {smooth_clean:<25} {rate_clean:<20}")
    print("-" * 80)


if __name__ == '__main__':
    main()
