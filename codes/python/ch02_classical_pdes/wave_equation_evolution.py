#!/usr/bin/env python3
"""
wave_equation_evolution.py

Visualizes the time evolution of the wave equation solution using the
truncated Fourier sine series derived in Chapter 2. The initial condition is
a plucked string (triangular displacement), which demonstrates standing wave
oscillations.

Wave equation (Dirichlet BCs on [0, L]):
    u_tt = c² u_xx

Analytical solution:
    u(x,t) = Σ (a_n cos(ω_n t) + b_n sin(ω_n t)) sin(nπx/L)

where ω_n = cnπ/L

Plucked string initial condition (triangle):
    f(x) = (2h/L)x           for 0 ≤ x ≤ L/2
    f(x) = 2h(1 - x/L)       for L/2 ≤ x ≤ L
    g(x) = 0                 (initially at rest)

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
L = np.pi              # Domain length [0, L]
C = 1.0                # Wave speed
H = 1.0                # Height of plucked string at center

# Period of oscillation: T = 2L/c
T_PERIOD = 2 * L / C

# Time snapshots: show half a period of oscillation
TIMES = [0.0, T_PERIOD/8, T_PERIOD/4, 3*T_PERIOD/8, T_PERIOD/2]

# Book color scheme
NAVY = '#142D6E'       # Primary navy from book
SKY = '#7896D2'        # Secondary sky blue from book

# Output path (relative to this script's location)
SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / '..' / '..' / '..' / 'book' / 'figures' / 'ch02' / 'python'
OUTPUT_FILE = OUTPUT_DIR / 'wave_evolution.pdf'


# -----------------------------------------------------------------------------
# Fourier coefficients for the plucked string (triangle)
# -----------------------------------------------------------------------------
def plucked_string_coefficients(n_modes, L, h):
    """
    Compute Fourier sine coefficients for the plucked string initial condition.

    f(x) = (2h/L)x           for 0 ≤ x ≤ L/2
    f(x) = 2h(1 - x/L)       for L/2 ≤ x ≤ L

    The Fourier sine series is:
        f(x) = Σ a_n sin(nπx/L)

    where a_n = (8h/(n²π²)) sin(nπ/2)
              = 8h/(n²π²) for n = 1, 5, 9, ...  (n ≡ 1 mod 4)
              = -8h/(n²π²) for n = 3, 7, 11, ... (n ≡ 3 mod 4)
              = 0 for even n

    Returns:
        a_n: array of sine coefficients (index n corresponds to sin(nπx/L))
        b_n: array of time derivative coefficients (all zeros since g(x) = 0)
    """
    a_n = np.zeros(n_modes + 1)  # a_n[n] = coefficient of sin(nπx/L)
    b_n = np.zeros(n_modes + 1)  # b_n[n] = 0 since initially at rest

    for n in range(1, n_modes + 1):
        if n % 2 == 1:  # odd n only
            # sin(nπ/2) = 1 for n = 1, 5, 9, ... and -1 for n = 3, 7, 11, ...
            sign = 1 if (n % 4 == 1) else -1
            a_n[n] = sign * 8.0 * h / (n**2 * np.pi**2)

    return a_n, b_n


# -----------------------------------------------------------------------------
# Evaluate the truncated Fourier series solution
# -----------------------------------------------------------------------------
def wave_solution(x, t, a_n, b_n, L, c):
    """
    Evaluate the wave equation solution at positions x and time t.

    u(x,t) = Σ_{n=1}^N (a_n cos(ω_n t) + b_n sin(ω_n t)) sin(nπx/L)

    where ω_n = cnπ/L

    Parameters:
        x: spatial positions (array)
        t: time (scalar)
        a_n, b_n: Fourier coefficients
        L: domain length
        c: wave speed

    Returns:
        u: solution values at positions x
    """
    u = np.zeros_like(x, dtype=float)

    n_modes = len(a_n) - 1
    for n in range(1, n_modes + 1):
        omega_n = c * n * np.pi / L
        spatial = np.sin(n * np.pi * x / L)
        temporal = a_n[n] * np.cos(omega_n * t) + b_n[n] * np.sin(omega_n * t)
        u += temporal * spatial

    return u


# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
def main():
    # Compute Fourier coefficients
    a_n, b_n = plucked_string_coefficients(N_MODES, L, H)

    # Spatial grid
    x = np.linspace(0, L, N_POINTS)

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

    # Plot solution at each time with varying line widths
    line_widths = np.linspace(2.0, 1.2, n_times)

    for i, t in enumerate(TIMES):
        u = wave_solution(x, t, a_n, b_n, L, C)
        # Create time label as fraction of period
        if t == 0:
            label = r'$t = 0$'
        else:
            # Express as fraction of T
            frac = t / T_PERIOD
            if abs(frac - 1/8) < 0.001:
                label = r'$t = T/8$'
            elif abs(frac - 1/4) < 0.001:
                label = r'$t = T/4$'
            elif abs(frac - 3/8) < 0.001:
                label = r'$t = 3T/8$'
            elif abs(frac - 1/2) < 0.001:
                label = r'$t = T/2$'
            else:
                label = f'$t = {t:.3f}$'
        ax.plot(x, u, color=colors[i], linewidth=line_widths[i], label=label)

    # Elegant axis formatting
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$u(x,t)$')
    ax.set_xlim(0, L)
    ax.set_ylim(-1.2, 1.2)

    # Custom x-ticks with π notation
    ax.set_xticks([0, L/4, L/2, 3*L/4, L])
    ax.set_xticklabels([r'$0$', r'$\frac{\pi}{4}$', r'$\frac{\pi}{2}$',
                        r'$\frac{3\pi}{4}$', r'$\pi$'])

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

    # Add horizontal line at y=0 for reference
    ax.axhline(y=0, color='#888888', linewidth=0.5, linestyle='-', alpha=0.3)

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
