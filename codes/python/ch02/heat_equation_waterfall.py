#!/usr/bin/env python3
"""
heat_equation_waterfall.py

Waterfall (3D surface) plot showing the complete space-time evolution of the
heat equation solution. This visualization demonstrates how the initial
triangle wave smooths out over time as higher frequency modes decay faster.

Heat equation (periodic BCs on [0, 2π]):
    u_t = u_xx

Analytical solution:
    u(x,t) = a_0 + Σ (a_n cos(nx) + b_n sin(nx)) exp(-n² t)

Triangle wave initial condition:
    f(x) = π - |x - π|  on [0, 2π]

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
T_MAX = 0.3            # Maximum time

# Book color scheme
NAVY = '#142D6E'
SKY = '#7896D2'

# Output path
SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / '..' / '..' / '..' / 'textbook' / 'figures' / 'ch02' / 'python'
OUTPUT_FILE = OUTPUT_DIR / 'heat_waterfall.pdf'


# -----------------------------------------------------------------------------
# Fourier coefficients for the triangle wave
# -----------------------------------------------------------------------------
def triangle_fourier_coefficients(n_modes):
    """
    Compute Fourier coefficients for the triangle wave f(x) = π - |x - π|.
    """
    a0 = np.pi / 2
    a_n = np.zeros(n_modes + 1)
    b_n = np.zeros(n_modes + 1)

    for n in range(1, n_modes + 1):
        if n % 2 == 1:  # odd n only
            a_n[n] = 4.0 / (np.pi * n**2)

    return a0, a_n, b_n


# -----------------------------------------------------------------------------
# Evaluate the heat equation solution
# -----------------------------------------------------------------------------
def heat_solution(x, t, a0, a_n, b_n):
    """
    Evaluate the heat equation solution at positions x and time t.
    u(x,t) = a_0 + Σ_{n=1}^N (a_n cos(nx) + b_n sin(nx)) exp(-n² t)
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

    # Create 2D grid for space-time
    x = np.linspace(0, 2 * np.pi, NX)
    t = np.linspace(0, T_MAX, NT)
    X, T = np.meshgrid(x, t)

    # Evaluate solution on grid
    U = np.zeros_like(X)
    for i, ti in enumerate(t):
        U[i, :] = heat_solution(x, ti, a0, a_n, b_n)

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
    ax.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi])
    ax.set_xticklabels([r'$0$', r'$\frac{\pi}{2}$', r'$\pi$',
                        r'$\frac{3\pi}{2}$', r'$2\pi$'])

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
