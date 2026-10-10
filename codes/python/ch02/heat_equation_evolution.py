#!/usr/bin/env python3
"""
heat_equation_evolution.py

Visualizes the time evolution of the heat equation solution using the
truncated Fourier series derived in Chapter 2. The initial condition is
a triangle wave, which demonstrates how higher frequency modes decay
faster than lower ones.

Heat equation (periodic BCs on [0, 2π]):
    u_t = u_xx

Analytical solution:
    u(x,t) = a_0 + Σ (a_n cos(nx) + b_n sin(nx)) exp(-n² t)

Triangle wave initial condition:
    f(x) = π - |x - π|  on [0, 2π]

Fourier coefficients (computed analytically):
    a_0 = π/2
    a_n = (4/π) / n²  for odd n, 0 for even n
    b_n = 0           (symmetric function)

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
N_MODES = 50           # Number of Fourier modes in truncation
N_POINTS = 1000        # Spatial grid points for smooth curves
TIMES = [0.0, 0.005, 0.02, 0.05, 0.1, 0.3]  # Time snapshots

# Book color scheme
NAVY = '#142D6E'       # Primary navy from book
SKY = '#7896D2'        # Secondary sky blue from book

# Output path (relative to this script's location)
SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / '..' / '..' / '..' / 'textbook' / 'figures' / 'ch02' / 'python'
OUTPUT_FILE = OUTPUT_DIR / 'heat_evolution.pdf'


# -----------------------------------------------------------------------------
# Fourier coefficients for the triangle wave
# -----------------------------------------------------------------------------
def triangle_fourier_coefficients(n_modes):
    """
    Compute Fourier coefficients for the triangle wave f(x) = π - |x - π|.

    The Fourier series is:
        f(x) = π/2 + (4/π) Σ_{k=1}^∞ cos((2k-1)x) / (2k-1)²

    Returns:
        a0: constant term (π/2)
        a_n: array of cosine coefficients (index n corresponds to cos(nx))
        b_n: array of sine coefficients (all zeros for symmetric function)
    """
    a0 = np.pi / 2

    a_n = np.zeros(n_modes + 1)  # a_n[n] = coefficient of cos(nx)
    b_n = np.zeros(n_modes + 1)  # b_n[n] = coefficient of sin(nx)

    for n in range(1, n_modes + 1):
        if n % 2 == 1:  # odd n only
            a_n[n] = 4.0 / (np.pi * n**2)

    return a0, a_n, b_n


# -----------------------------------------------------------------------------
# Evaluate the truncated Fourier series solution
# -----------------------------------------------------------------------------
def heat_solution(x, t, a0, a_n, b_n):
    """
    Evaluate the heat equation solution at positions x and time t.

    u(x,t) = a_0 + Σ_{n=1}^N (a_n cos(nx) + b_n sin(nx)) exp(-n² t)

    Parameters:
        x: spatial positions (array)
        t: time (scalar)
        a0, a_n, b_n: Fourier coefficients

    Returns:
        u: solution values at positions x
    """
    u = np.full_like(x, a0, dtype=float)

    n_modes = len(a_n) - 1
    for n in range(1, n_modes + 1):
        decay = np.exp(-n**2 * t)
        u += (a_n[n] * np.cos(n * x) + b_n[n] * np.sin(n * x)) * decay

    return u


# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
def main():
    # Compute Fourier coefficients
    a0, a_n, b_n = triangle_fourier_coefficients(N_MODES)

    # Spatial grid
    x = np.linspace(0, 2 * np.pi, N_POINTS)

    # Create figure with golden ratio proportions
    fig, ax = plt.subplots(figsize=(6, 3.7))

    # Elegant color palette: navy to sky blue gradient
    n_times = len(TIMES)
    colors = []
    for i in range(n_times):
        # Interpolate from navy to a lighter blue
        t_frac = i / (n_times - 1) if n_times > 1 else 0
        # Navy RGB: (20, 45, 110) -> lighter (120, 150, 210)
        r = int(20 + t_frac * 100) / 255
        g = int(45 + t_frac * 105) / 255
        b = int(110 + t_frac * 100) / 255
        colors.append((r, g, b))

    # Plot solution at each time with varying line styles for clarity
    line_widths = np.linspace(2.0, 1.2, n_times)

    for i, t in enumerate(TIMES):
        u = heat_solution(x, t, a0, a_n, b_n)
        if t == 0:
            label = r'$t = 0$'
        elif t < 0.01:
            label = f'$t = {t:.3f}$'
        elif t < 0.1:
            label = f'$t = {t:.2f}$'
        else:
            label = f'$t = {t:.1f}$'
        ax.plot(x, u, color=colors[i], linewidth=line_widths[i], label=label)

    # Elegant axis formatting
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$u(x,t)$')
    ax.set_xlim(0, 2 * np.pi)
    ax.set_ylim(-0.1, 3.4)

    # Custom x-ticks with π notation
    ax.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi])
    ax.set_xticklabels([r'$0$', r'$\frac{\pi}{2}$', r'$\pi$',
                        r'$\frac{3\pi}{2}$', r'$2\pi$'])

    # Clean legend - centered at top
    legend = ax.legend(loc='upper center', frameon=True, fancybox=False,
                       edgecolor='none', facecolor='white', framealpha=0.9,
                       ncol=3, columnspacing=1.0, bbox_to_anchor=(0.5, 1.0))

    # Minimalist spine styling
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color('#444444')
    ax.spines['bottom'].set_color('#444444')
    ax.tick_params(colors='#444444')

    # Add subtle grid
    ax.yaxis.grid(True, linestyle='-', alpha=0.2, color='#888888')
    ax.set_axisbelow(True)

    plt.tight_layout()

    # Ensure output directory exists and save
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Save as PDF (vector) and PNG (for preview)
    fig.savefig(OUTPUT_FILE, bbox_inches='tight', pad_inches=0.05)
    fig.savefig(OUTPUT_FILE.with_suffix('.png'), bbox_inches='tight',
                pad_inches=0.05, dpi=300)

    print(f'Figure saved to: {OUTPUT_FILE.resolve()}')

    plt.close(fig)


if __name__ == '__main__':
    main()
