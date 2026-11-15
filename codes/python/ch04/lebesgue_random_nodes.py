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
    3. Repeat M times to obtain statistics (median, quartiles, min, max)
    4. Fit an empirical exponential law to the median

Why the median: two independently drawn nodes can fall arbitrarily close to
each other, and a single such draw makes the Lebesgue constant astronomically
large. The M values therefore have a very heavy right tail, and their sample
mean is dominated by the worst draws: it jumps by orders of magnitude from one
N to the next. The median and the quartiles describe the typical draw and do
not feel the tail, so the empirical law Lambda_N ~ a * exp(b * N) is fitted to
the median, by linear least squares on its logarithm.

Key finding: over the degrees computed, the median Lebesgue constant of random
nodes grows exponentially and faster than that of equispaced nodes, so that the
special structure of Chebyshev points (endpoint clustering) is essential for
stable interpolation.

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
def exponential_model(N, a, b):
    """Model: Lambda ~ a * exp(b * N)"""
    return a * np.exp(b * N)


def fit_exponential(N_values, Lambda_values):
    """
    Fit the exponential law Lambda_N ~ a * exp(b * N) to positive data.

    The fit is the least-squares straight line through the points
    (N, ln Lambda_N): every decade of Lambda_N weighs the same, no initial
    guess is needed, and the fitted curve is positive by construction.

    Returns the fitted parameters (a, b) and the coefficient of
    determination R^2 of the straight line in the (N, ln Lambda_N) plane.
    """
    log_Lambda = np.log(Lambda_values)
    b, log_a = np.polyfit(N_values, log_Lambda, 1)
    residual = log_Lambda - (log_a + b * N_values)
    r_squared = 1.0 - np.sum(residual**2) / np.sum((log_Lambda - log_Lambda.mean())**2)
    return np.exp(log_a), b, r_squared


def format_power_of_ten(value):
    """Label text '$10^{p}$' for a positive value, with p to one decimal."""
    return f'$10^{{{np.log10(value):.1f}}}$'


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
        Lambda_min[i] = np.min(samples)
        Lambda_max[i] = np.max(samples)
        Lambda_median[i] = np.median(samples)
        Lambda_p25[i] = np.percentile(samples, 25)
        Lambda_p75[i] = np.percentile(samples, 75)

        # Deterministic distributions
        Lambda_equi[i] = lebesgue_constant(equispaced_nodes(N), N_FINE)
        Lambda_cheb[i] = lebesgue_constant(chebyshev_nodes(N), N_FINE)

        if (i + 1) % 5 == 0 or N == N_MAX:
            print(f"  N = {N:2d}: Lambda_median = {Lambda_median[i]:10.4e}, "
                  f"IQR = [{Lambda_p25[i]:.2e}, {Lambda_p75[i]:.2e}], "
                  f"range = [{Lambda_min[i]:.2e}, {Lambda_max[i]:.2e}]")

    # =========================================================================
    # Fit the empirical exponential law to the medians
    # =========================================================================
    print("\n" + "-" * 70)
    print("Fitting the empirical law Lambda_N ~ a * exp(b * N) to the medians...")

    # Use data for N >= 5 to avoid small-N effects
    mask = N_values >= 5
    a_fit, b_fit, r2_fit = fit_exponential(N_values[mask], Lambda_median[mask])

    # The same straight-line fit for equispaced nodes, for comparison
    a_equi, b_equi, r2_equi = fit_exponential(N_values[mask], Lambda_equi[mask])

    print(f"\nRandom nodes (median): Lambda_N ~ {a_fit:.4f} * exp({b_fit:.4f} * N)"
          f"  (R^2 = {r2_fit:.4f} in the (N, ln Lambda_N) plane)")
    print(f"Equispaced nodes:      Lambda_N ~ {a_equi:.4f} * exp({b_equi:.4f} * N)"
          f"  (R^2 = {r2_equi:.4f} in the (N, ln Lambda_N) plane)")

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

    # Median line for random nodes
    ax1.semilogy(N_values, Lambda_median, 'o-', color=GOLD, linewidth=1.5,
                 markersize=4, label='Random (median)', zorder=5)

    # Equispaced for comparison
    ax1.semilogy(N_values, Lambda_equi, 's--', color=CORAL, linewidth=1.2,
                 markersize=3, alpha=0.7, label='Equispaced')

    # Chebyshev for reference
    ax1.semilogy(N_values, Lambda_cheb, '^--', color=SKY, linewidth=1.2,
                 markersize=3, alpha=0.7, label='Chebyshev')

    # Exponential law fitted to the median
    N_fit = np.linspace(5, N_MAX, 100)
    Lambda_fit = exponential_model(N_fit, a_fit, b_fit)
    fit_label = f'$\\Lambda_N \\approx {a_fit:.3g}\\,e^{{{b_fit:.3f}N}}$ (fit to median)'
    ax1.semilogy(N_fit, Lambda_fit, '-', color=NAVY, linewidth=1.5,
                 alpha=0.8, label=fit_label)

    ax1.set_xlabel(r'$N$ (polynomial degree)')
    ax1.set_ylabel(r'Lebesgue constant $\Lambda_N$')
    ax1.set_title('Growth of Lebesgue Constant', fontsize=11)
    ax1.set_xlim(0, N_MAX + 2)
    ax1.legend(loc='upper left', fontsize=8)

    # Growth law of the Chebyshev constant
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
    counts, _, _ = ax2.hist(log_samples, bins=bins, density=False, color=GOLD,
                            alpha=0.7, edgecolor='white', linewidth=0.5)

    # The bars fill the lower 62 % of the panel and the vertical lines stop at
    # 72 %, which leaves the top of the panel to the legend
    ax2.set_ylim(0, counts.max() / 0.62)
    ax2.set_xlim(0, bins[-1] + 0.5)
    line_top = 0.72

    # Add vertical lines for statistics (on log scale)
    ax2.axvline(np.log10(median_val), ymax=line_top, color=TEAL, linewidth=2,
                linestyle='--', label=f'Median = {format_power_of_ten(median_val)}')
    ax2.axvline(np.log10(mean_val), ymax=line_top, color=NAVY, linewidth=2,
                linestyle='-', label=f'Mean = {format_power_of_ten(mean_val)}')
    ax2.axvline(np.log10(equi_val), ymax=line_top, color=CORAL, linewidth=2,
                linestyle=':', label=f'Equispaced = {format_power_of_ten(equi_val)}')
    ax2.axvline(np.log10(cheb_val), ymax=line_top, color=SKY, linewidth=2,
                linestyle='-.', label=f'Chebyshev = {cheb_val:.1f}')

    ax2.set_xlabel(r'$\log_{10}(\Lambda_{15})$')
    ax2.set_ylabel('Count')
    ax2.set_title(f'Distribution of $\\Lambda_N$ for $N = {N_hist}$ ({M_hist} trials)',
                  fontsize=11)

    # Legend, titled with the spread of the samples
    spread = np.log10(samples_hist.max()) - np.log10(samples_hist.min())
    ax2.legend(loc='upper center', ncol=2, fontsize=8, title_fontsize=8,
               title=f'Spread of the samples: {spread:.1f} orders of magnitude')

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
    print("\n" + "=" * 78)
    print("Summary Statistics (medians over the draws; the mean is shown for contrast)")
    print("=" * 78)
    print(f"{'N':>4} {'Median':>12} {'(mean)':>12} {'Min':>12} {'Max':>12} "
          f"{'Equispaced':>12}")
    print("-" * 78)
    for i in range(0, n_N, 3):  # Print every 3rd value
        N = N_values[i]
        print(f"{N:>4d} {Lambda_median[i]:>12.4e} {Lambda_mean[i]:>12.4e} "
              f"{Lambda_min[i]:>12.4e} {Lambda_max[i]:>12.4e} {Lambda_equi[i]:>12.4e}")
    print("-" * 78)

    print("\n" + "=" * 78)
    print("EMPIRICAL GROWTH LAW OF THE MEDIAN")
    print("=" * 78)
    print(f"\nFor random nodes uniformly distributed on [-1, 1],"
          f" N from {N_values[mask][0]} to {N_MAX}:")
    print(f"\n    Lambda_N ~ {a_fit:.4f} * exp({b_fit:.4f} * N)"
          f" = {a_fit:.4f} * {np.exp(b_fit):.4f}^N")
    print(f"\nFor equispaced nodes, the same fit gives"
          f" Lambda_N ~ {a_equi:.4f} * {np.exp(b_equi):.4f}^N.")
    print(f"For Chebyshev nodes, Lambda_{N_MAX} = {Lambda_cheb[-1]:.4f}.")
    print("=" * 78)

    plt.close(fig)

    return {
        'N_values': N_values,
        'Lambda_median': Lambda_median,
        'Lambda_p25': Lambda_p25,
        'Lambda_p75': Lambda_p75,
        'Lambda_min': Lambda_min,
        'Lambda_max': Lambda_max,
        'Lambda_mean': Lambda_mean,
        'fit_params': (a_fit, b_fit),
        'fit_params_equispaced': (a_equi, b_equi),
        'Lambda_equi': Lambda_equi,
        'Lambda_cheb': Lambda_cheb
    }


if __name__ == '__main__':
    results = main()
