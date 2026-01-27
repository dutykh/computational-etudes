#!/usr/bin/env python3
"""
lebesgue_random_nodes.py

Computational experiment: Lebesgue constant for random nodes on [-1, 1].

This script investigates the growth of the Lebesgue constant when interpolation
nodes are chosen randomly (uniformly distributed) on [-1, 1]. Through Monte Carlo
simulation, we discover the statistical behavior and derive an empirical
asymptotic formula.

The experiment:
    1. For each N, generate N+1 random points uniformly on [-1, 1]
    2. Compute the Lebesgue constant for these random nodes
    3. Repeat M times to obtain statistics (mean, std, min, max)
    4. Fit an empirical asymptotic formula to the averaged data

Key finding: Random nodes exhibit exponential growth similar to equispaced nodes,
demonstrating that the special structure of Chebyshev points (endpoint clustering)
is essential for stable interpolation.

Author: Dr. Denys Dutykh
        Mathematics Department
        Khalifa University of Science and Technology
        Abu Dhabi, UAE

Part of the book "Computational Etudes: A Spectral Approach"
https://github.com/dutykh/computational-etudes

Created: January 2026
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from pathlib import Path
from scipy.optimize import curve_fit

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
N_FINE = 2000           # Number of points for Lebesgue function evaluation
M_TRIALS = 200          # Number of Monte Carlo trials per N
N_MAX = 30              # Maximum polynomial degree
SEED = 42               # Random seed for reproducibility

# Book color scheme
NAVY = '#142D6E'
SKY = '#7896D2'
CORAL = '#E74C3C'
TEAL = '#16A085'
GOLD = '#F39C12'        # New color for random nodes

# Output path
SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / '..' / '..' / '..' / 'textbook' / 'figures' / 'ch04' / 'python'
OUTPUT_FILE = OUTPUT_DIR / 'lebesgue_random_nodes.pdf'


# -----------------------------------------------------------------------------
# Node generation
# -----------------------------------------------------------------------------
def random_nodes(N, rng=None):
    """
    Generate N+1 random nodes, sorted, uniformly distributed on [-1, 1].

    Parameters
    ----------
    N : int
        Polynomial degree (generates N+1 nodes)
    rng : numpy.random.Generator, optional
        Random number generator for reproducibility

    Returns
    -------
    ndarray
        Sorted array of N+1 random nodes in [-1, 1]
    """
    if rng is None:
        rng = np.random.default_rng()
    return np.sort(rng.uniform(-1, 1, N + 1))


def equispaced_nodes(N):
    """Equispaced nodes on [-1, 1]."""
    return np.linspace(-1, 1, N + 1)


def chebyshev_nodes(N):
    """Chebyshev-Gauss-Lobatto nodes on [-1, 1]."""
    j = np.arange(N + 1)
    return np.cos(j * np.pi / N)


# -----------------------------------------------------------------------------
# Lebesgue function computation
# -----------------------------------------------------------------------------
def lagrange_basis(x_nodes, k, x_eval):
    """
    Compute the k-th Lagrange basis polynomial at points x_eval.

    Parameters
    ----------
    x_nodes : ndarray
        Interpolation nodes
    k : int
        Index of the basis polynomial
    x_eval : ndarray
        Evaluation points

    Returns
    -------
    ndarray
        Values of L_k(x) at x_eval
    """
    n = len(x_nodes)
    x_eval = np.atleast_1d(x_eval)
    L_k = np.ones_like(x_eval, dtype=float)

    for j in range(n):
        if j != k:
            L_k *= (x_eval - x_nodes[j]) / (x_nodes[k] - x_nodes[j])

    return L_k


def lebesgue_function(x_nodes, x_eval):
    """
    Compute the Lebesgue function: Lambda_N(x) = sum_{k=0}^{N} |L_k(x)|

    Parameters
    ----------
    x_nodes : ndarray
        Interpolation nodes
    x_eval : ndarray
        Evaluation points

    Returns
    -------
    ndarray
        Values of the Lebesgue function at x_eval
    """
    N = len(x_nodes) - 1
    Lambda = np.zeros_like(x_eval, dtype=float)

    for k in range(N + 1):
        Lambda += np.abs(lagrange_basis(x_nodes, k, x_eval))

    return Lambda


def lebesgue_constant(x_nodes, n_eval=2000):
    """
    Compute the Lebesgue constant: Lambda_N = max_{x in [-1,1]} Lambda_N(x)

    Parameters
    ----------
    x_nodes : ndarray
        Interpolation nodes
    n_eval : int
        Number of evaluation points for finding the maximum

    Returns
    -------
    float
        The Lebesgue constant
    """
    x_eval = np.linspace(-1, 1, n_eval)
    Lambda = lebesgue_function(x_nodes, x_eval)
    return np.max(Lambda)


# -----------------------------------------------------------------------------
# Monte Carlo experiment
# -----------------------------------------------------------------------------
def monte_carlo_lebesgue(N, M=100, n_eval=2000, rng=None):
    """
    Compute M samples of Lebesgue constant for random nodes.

    Parameters
    ----------
    N : int
        Polynomial degree
    M : int
        Number of Monte Carlo trials
    n_eval : int
        Number of evaluation points for Lebesgue constant
    rng : numpy.random.Generator, optional
        Random number generator

    Returns
    -------
    ndarray
        Array of M Lebesgue constant samples
    """
    if rng is None:
        rng = np.random.default_rng()

    samples = np.zeros(M)
    for m in range(M):
        x_rand = random_nodes(N, rng)
        samples[m] = lebesgue_constant(x_rand, n_eval)

    return samples


# -----------------------------------------------------------------------------
# Empirical fitting
# -----------------------------------------------------------------------------
def exponential_model(N, a, b, c):
    """Model: Lambda ~ a * exp(b * N) / N^c"""
    return a * np.exp(b * N) / np.power(N, c)


def fit_exponential(N_values, Lambda_values):
    """
    Fit an exponential model to the Lebesgue constant data.

    Model: Lambda_N ~ a * exp(b * N) / N^c

    Returns fitted parameters (a, b, c) and the covariance matrix.
    """
    # Initial guess
    p0 = [0.5, 0.5, 0.5]

    try:
        popt, pcov = curve_fit(exponential_model, N_values, Lambda_values,
                               p0=p0, maxfev=10000)
        return popt, pcov
    except RuntimeError:
        # Fallback: simple exponential fit in log space
        log_Lambda = np.log(Lambda_values)
        coeffs = np.polyfit(N_values, log_Lambda, 1)
        # Lambda ~ exp(coeffs[1]) * exp(coeffs[0] * N)
        return [np.exp(coeffs[1]), coeffs[0], 0.0], None


# -----------------------------------------------------------------------------
# Main computation and visualization
# -----------------------------------------------------------------------------
def main():
    print("=" * 70)
    print("Computational Experiment: Lebesgue Constant for Random Nodes")
    print("=" * 70)

    # Set random seed for reproducibility
    rng = np.random.default_rng(SEED)

    # Arrays to store results
    N_values = np.arange(2, N_MAX + 1)
    n_N = len(N_values)

    # Statistics arrays
    Lambda_mean = np.zeros(n_N)
    Lambda_std = np.zeros(n_N)
    Lambda_min = np.zeros(n_N)
    Lambda_max = np.zeros(n_N)
    Lambda_median = np.zeros(n_N)
    Lambda_p25 = np.zeros(n_N)
    Lambda_p75 = np.zeros(n_N)

    # Also compute equispaced for comparison
    Lambda_equi = np.zeros(n_N)
    Lambda_cheb = np.zeros(n_N)

    print(f"\nRunning Monte Carlo simulation with M = {M_TRIALS} trials per N...")
    print(f"N ranges from 2 to {N_MAX}")
    print()

    for i, N in enumerate(N_values):
        # Monte Carlo for random nodes
        samples = monte_carlo_lebesgue(N, M=M_TRIALS, n_eval=N_FINE, rng=rng)

        Lambda_mean[i] = np.mean(samples)
        Lambda_std[i] = np.std(samples)
        Lambda_min[i] = np.min(samples)
        Lambda_max[i] = np.max(samples)
        Lambda_median[i] = np.median(samples)
        Lambda_p25[i] = np.percentile(samples, 25)
        Lambda_p75[i] = np.percentile(samples, 75)

        # Deterministic distributions
        Lambda_equi[i] = lebesgue_constant(equispaced_nodes(N), N_FINE)
        Lambda_cheb[i] = lebesgue_constant(chebyshev_nodes(N), N_FINE)

        if (i + 1) % 5 == 0 or N == N_MAX:
            print(f"  N = {N:2d}: Lambda_mean = {Lambda_mean[i]:12.4f}, "
                  f"std = {Lambda_std[i]:10.4f}, "
                  f"range = [{Lambda_min[i]:.2f}, {Lambda_max[i]:.2f}]")

    # =========================================================================
    # Fit empirical formula
    # =========================================================================
    print("\n" + "-" * 70)
    print("Fitting empirical asymptotic formula...")

    # Use data for N >= 5 to avoid small-N effects
    mask = N_values >= 5
    popt, pcov = fit_exponential(N_values[mask], Lambda_mean[mask])
    a_fit, b_fit, c_fit = popt

    print(f"\nEmpirical formula: Lambda_N ~ {a_fit:.4f} * exp({b_fit:.4f} * N) / N^{c_fit:.4f}")

    # Simplified formula (if c is small, we can ignore it)
    if abs(c_fit) < 0.3:
        # Refit with c = 0
        log_Lambda = np.log(Lambda_mean[mask])
        coeffs = np.polyfit(N_values[mask], log_Lambda, 1)
        b_simple = coeffs[0]
        a_simple = np.exp(coeffs[1])
        print(f"Simplified formula: Lambda_N ~ {a_simple:.4f} * exp({b_simple:.4f} * N)")
        print(f"                  or Lambda_N ~ {a_simple:.4f} * {np.exp(b_simple):.4f}^N")

    # Compare growth rate with equispaced
    log_equi = np.log(Lambda_equi[mask])
    coeffs_equi = np.polyfit(N_values[mask], log_equi, 1)
    b_equi = coeffs_equi[0]

    print(f"\nFor comparison:")
    print(f"  Random nodes growth rate:    exp({b_fit:.4f} * N) ~ {np.exp(b_fit):.4f}^N")
    print(f"  Equispaced nodes growth rate: exp({b_equi:.4f} * N) ~ {np.exp(b_equi):.4f}^N")

    # =========================================================================
    # Create figure
    # =========================================================================
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

    # -------------------------------------------------------------------------
    # Left panel: Growth of Lebesgue constant with confidence bands
    # -------------------------------------------------------------------------
    # Shaded region for random nodes (min to max)
    ax1.fill_between(N_values, Lambda_min, Lambda_max,
                     color=GOLD, alpha=0.2, label='Random (min-max range)')

    # Shaded region for interquartile range
    ax1.fill_between(N_values, Lambda_p25, Lambda_p75,
                     color=GOLD, alpha=0.3, label='Random (25th-75th percentile)')

    # Mean line for random nodes
    ax1.semilogy(N_values, Lambda_mean, 'o-', color=GOLD, linewidth=1.5,
                 markersize=4, label='Random (mean)', zorder=5)

    # Equispaced for comparison
    ax1.semilogy(N_values, Lambda_equi, 's--', color=CORAL, linewidth=1.2,
                 markersize=3, alpha=0.7, label='Equispaced')

    # Chebyshev for reference
    ax1.semilogy(N_values, Lambda_cheb, '^--', color=SKY, linewidth=1.2,
                 markersize=3, alpha=0.7, label='Chebyshev')

    # Fitted curve
    N_fit = np.linspace(5, N_MAX, 100)
    Lambda_fit = exponential_model(N_fit, *popt)
    ax1.semilogy(N_fit, Lambda_fit, '-', color=NAVY, linewidth=1.5,
                 alpha=0.8, label=f'Fit: $\\Lambda_N \\sim {a_fit:.2f} e^{{{b_fit:.3f}N}}$')

    ax1.set_xlabel(r'$N$ (polynomial degree)')
    ax1.set_ylabel(r'Lebesgue constant $\Lambda_N$')
    ax1.set_title('Growth of Lebesgue Constant', fontsize=11)
    ax1.set_xlim(0, N_MAX + 2)
    ax1.legend(loc='upper left', fontsize=8)

    # Add text annotations for growth rates
    ax1.text(N_MAX - 5, Lambda_mean[-1] * 2, r'$O(e^{bN})$',
             fontsize=10, color=GOLD)
    ax1.text(N_MAX - 5, Lambda_cheb[-1] * 1.5, r'$O(\ln N)$',
             fontsize=10, color=SKY)

    # -------------------------------------------------------------------------
    # Right panel: Distribution at fixed N (log-scale histogram)
    # -------------------------------------------------------------------------
    N_hist = 15  # Fixed N for histogram
    M_hist = 500
    samples_hist = monte_carlo_lebesgue(N_hist, M=M_hist, n_eval=N_FINE,
                                        rng=np.random.default_rng(SEED + 1))

    # Compute statistics
    mean_val = np.mean(samples_hist)
    median_val = np.median(samples_hist)
    equi_val = lebesgue_constant(equispaced_nodes(N_hist), N_FINE)
    cheb_val = lebesgue_constant(chebyshev_nodes(N_hist), N_FINE)

    # Use log-scale histogram to show the spread across orders of magnitude
    log_samples = np.log10(samples_hist)

    # Create histogram with log-spaced bins
    bins = np.linspace(np.floor(log_samples.min()), np.ceil(log_samples.max()), 35)
    ax2.hist(log_samples, bins=bins, density=False, color=GOLD, alpha=0.7,
             edgecolor='white', linewidth=0.5)

    # Add vertical lines for statistics (on log scale)
    ax2.axvline(np.log10(median_val), color=TEAL, linewidth=2, linestyle='--',
                label=f'Median = $10^{{{np.log10(median_val):.1f}}}$')
    ax2.axvline(np.log10(mean_val), color=NAVY, linewidth=2, linestyle='-',
                label=f'Mean = $10^{{{np.log10(mean_val):.1f}}}$')
    ax2.axvline(np.log10(equi_val), color=CORAL, linewidth=2, linestyle=':',
                label=f'Equispaced = $10^{{{np.log10(equi_val):.1f}}}$')
    ax2.axvline(np.log10(cheb_val), color=SKY, linewidth=2, linestyle='-.',
                label=f'Chebyshev = {cheb_val:.1f}')

    ax2.set_xlabel(r'$\log_{10}(\Lambda_{15})$')
    ax2.set_ylabel('Count')
    ax2.set_title(f'Distribution of $\\Lambda_N$ for $N = {N_hist}$ ({M_hist} trials)',
                  fontsize=11)
    ax2.legend(loc='upper right', fontsize=8)

    # Add annotation about the spread
    spread = np.log10(samples_hist.max()) - np.log10(samples_hist.min())
    ax2.text(0.05, 0.95, f'Spread: {spread:.1f} orders\nof magnitude',
             transform=ax2.transAxes, fontsize=9, ha='left', va='top',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    # -------------------------------------------------------------------------
    # Clean styling for both panels
    # -------------------------------------------------------------------------
    for ax in [ax1, ax2]:
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

    print(f'\nFigure saved to: {OUTPUT_FILE.resolve()}')

    # =========================================================================
    # Print summary table
    # =========================================================================
    print("\n" + "=" * 70)
    print("Summary Statistics")
    print("=" * 70)
    print(f"{'N':>4} {'Mean':>12} {'Std':>12} {'Min':>12} {'Max':>12} {'Equispaced':>12}")
    print("-" * 70)
    for i in range(0, n_N, 3):  # Print every 3rd value
        N = N_values[i]
        print(f"{N:>4d} {Lambda_mean[i]:>12.4f} {Lambda_std[i]:>12.4f} "
              f"{Lambda_min[i]:>12.4f} {Lambda_max[i]:>12.4f} {Lambda_equi[i]:>12.4f}")
    print("-" * 70)

    print("\n" + "=" * 70)
    print("EMPIRICAL ASYMPTOTIC FORMULA")
    print("=" * 70)
    print(f"\nFor random nodes uniformly distributed on [-1, 1]:")
    print(f"\n    Lambda_N ~ {a_fit:.4f} * exp({b_fit:.4f} * N) / N^{c_fit:.4f}")
    print(f"\nor approximately:")
    print(f"\n    Lambda_N ~ {a_fit:.4f} * {np.exp(b_fit):.4f}^N")
    print(f"\nThis confirms exponential growth similar to equispaced nodes,")
    print(f"demonstrating that random node placement provides no advantage")
    print(f"over structured distributions without the clustering property")
    print(f"of Chebyshev points near the interval endpoints.")
    print("=" * 70)

    plt.close(fig)

    return {
        'N_values': N_values,
        'Lambda_mean': Lambda_mean,
        'Lambda_std': Lambda_std,
        'fit_params': popt,
        'Lambda_equi': Lambda_equi,
        'Lambda_cheb': Lambda_cheb
    }


if __name__ == '__main__':
    results = main()
