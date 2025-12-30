#!/usr/bin/env python3
"""
wave_equation_waterfall.py

Waterfall (3D surface) plot showing the complete space-time evolution of the
wave equation solution. This visualization demonstrates the standing wave
oscillations of a plucked string over time.

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
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LinearSegmentedColormap
from pathlib import Path

# -----------------------------------------------------------------------------
# Publication-quality matplotlib configuration
# -----------------------------------------------------------------------------
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Computer Modern Roman', 'CMU Serif', 'DejaVu Serif']
rcParams['mathtext.fontset'] = 'cm'
rcParams['axes.labelsize'] = 11
rcParams['axes.titlesize'] = 12
rcParams['xtick.labelsize'] = 9
rcParams['ytick.labelsize'] = 9
rcParams['legend.fontsize'] = 10
rcParams['figure.dpi'] = 150
rcParams['savefig.dpi'] = 300
rcParams['text.usetex'] = False
rcParams['axes.linewidth'] = 0.8

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
N_MODES = 50           # Number of Fourier modes in truncation
NX = 150               # Spatial grid points
NT = 100               # Time grid points
L = np.pi              # Domain length [0, L]
C = 1.0                # Wave speed
H = 1.0                # Height of plucked string at center

# Period of oscillation: T = 2L/c
T_PERIOD = 2 * L / C
T_MAX = T_PERIOD       # Show one full period

# Book color scheme
NAVY = '#142D6E'
SKY = '#7896D2'

# Output path
SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / '..' / '..' / '..' / 'book' / 'figures' / 'ch02' / 'python'
OUTPUT_FILE = OUTPUT_DIR / 'wave_waterfall.pdf'


# -----------------------------------------------------------------------------
# Fourier coefficients for the plucked string
# -----------------------------------------------------------------------------
def plucked_string_coefficients(n_modes, L, h):
    """
    Compute Fourier sine coefficients for the plucked string initial condition.
    """
    a_n = np.zeros(n_modes + 1)
    b_n = np.zeros(n_modes + 1)

    for n in range(1, n_modes + 1):
        if n % 2 == 1:  # odd n only
            sign = 1 if (n % 4 == 1) else -1
            a_n[n] = sign * 8.0 * h / (n**2 * np.pi**2)

    return a_n, b_n


# -----------------------------------------------------------------------------
# Evaluate the wave equation solution
# -----------------------------------------------------------------------------
def wave_solution(x, t, a_n, b_n, L, c):
    """
    Evaluate the wave equation solution at positions x and time t.
    u(x,t) = Σ_{n=1}^N (a_n cos(ω_n t) + b_n sin(ω_n t)) sin(nπx/L)
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

    # Create 2D grid for space-time
    x = np.linspace(0, L, NX)
    t = np.linspace(0, T_MAX, NT)
    X, T = np.meshgrid(x, t)

    # Evaluate solution on grid
    U = np.zeros_like(X)
    for i, ti in enumerate(t):
        U[i, :] = wave_solution(x, ti, a_n, b_n, L, C)

    # Create figure with 3D axes
    fig = plt.figure(figsize=(7, 5))
    ax = fig.add_subplot(111, projection='3d')

    # Create custom colormap (navy to white to sky)
    n_colors = 256
    navy_rgb = np.array([20, 45, 110]) / 255
    white_rgb = np.array([1, 1, 1])
    sky_rgb = np.array([120, 150, 210]) / 255

    colors_list = []
    for i in range(n_colors // 2):
        t_val = i / (n_colors // 2 - 1)
        colors_list.append(tuple(navy_rgb * (1 - t_val) + white_rgb * t_val))
    for i in range(n_colors // 2):
        t_val = i / (n_colors // 2 - 1)
        colors_list.append(tuple(white_rgb * (1 - t_val) + sky_rgb * t_val))

    custom_cmap = LinearSegmentedColormap.from_list('navy_sky', colors_list)

    # Plot the surface
    surf = ax.plot_surface(X, T, U, cmap=custom_cmap, linewidth=0.1,
                           antialiased=True, alpha=0.9, edgecolor='gray',
                           rcount=50, ccount=50)

    # Add wireframe for better depth perception
    ax.plot_wireframe(X, T, U, color='gray', linewidth=0.2, alpha=0.3,
                      rstride=5, cstride=5)

    # Labels
    ax.set_xlabel(r'$x$', labelpad=8)
    ax.set_ylabel(r'$t$', labelpad=8)
    ax.set_zlabel(r'$u(x,t)$', labelpad=8)

    # Custom x-ticks with π notation
    ax.set_xticks([0, L/4, L/2, 3*L/4, L])
    ax.set_xticklabels([r'$0$', r'$\frac{\pi}{4}$', r'$\frac{\pi}{2}$',
                        r'$\frac{3\pi}{4}$', r'$\pi$'])

    # Custom t-ticks as fractions of period
    ax.set_yticks([0, T_PERIOD/4, T_PERIOD/2, 3*T_PERIOD/4, T_PERIOD])
    ax.set_yticklabels([r'$0$', r'$\frac{T}{4}$', r'$\frac{T}{2}$',
                        r'$\frac{3T}{4}$', r'$T$'])

    # Set viewing angle for best visualization
    ax.view_init(elev=25, azim=-60)

    # Adjust layout
    plt.tight_layout()

    # Ensure output directory exists and save
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    fig.savefig(OUTPUT_FILE, bbox_inches='tight', pad_inches=0.1)
    fig.savefig(OUTPUT_FILE.with_suffix('.png'), bbox_inches='tight',
                pad_inches=0.1, dpi=300)

    print(f'Figure saved to: {OUTPUT_FILE.resolve()}')

    plt.close(fig)


if __name__ == '__main__':
    main()
