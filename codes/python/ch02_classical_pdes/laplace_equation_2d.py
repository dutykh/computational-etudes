#!/usr/bin/env python3
"""
laplace_equation_2d.py

Visualizes the solution of the Laplace equation in a periodic strip using the
truncated Fourier series derived in Chapter 2. The boundary condition is a
smooth function demonstrating how different Fourier modes decay toward the
interior.

Laplace equation in strip D = [0, 2π] × [0, 1]:
    u_xx + u_yy = 0

Boundary conditions:
    u(x, 0) = f(x) = sin(x) + 0.5*sin(3x)   (prescribed at bottom)
    u(x, 1) = 0                              (zero at top)
    u periodic in x with period 2π

Analytical solution:
    u(x,y) = a_0(1-y) + Σ [a_n cos(nx) + b_n sin(nx)] sinh(n(1-y))/sinh(n)

where a_n, b_n are Fourier coefficients of f(x).

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
NX = 300               # Grid points in x direction
NY = 150               # Grid points in y direction

# Book color scheme
NAVY = '#142D6E'       # Primary navy from book
SKY = '#7896D2'        # Secondary sky blue from book

# Output path (relative to this script's location)
SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / '..' / '..' / '..' / 'book' / 'figures' / 'ch02' / 'python'
OUTPUT_FILE = OUTPUT_DIR / 'laplace_solution.pdf'


# -----------------------------------------------------------------------------
# Boundary data and Fourier coefficients
# -----------------------------------------------------------------------------
def boundary_function(x):
    """
    Boundary data at y = 0: f(x) = sin(x) + 0.5*sin(3x)
    """
    return np.sin(x) + 0.5 * np.sin(3 * x)


def boundary_fourier_coefficients(n_modes):
    """
    Compute Fourier coefficients for f(x) = sin(x) + 0.5*sin(3x).

    For this specific function:
        a_0 = 0 (integral of sin over [0, 2π] is zero)
        a_n = 0 for all n (no cosine terms)
        b_1 = 1, b_3 = 0.5, b_n = 0 otherwise

    Returns:
        a0: constant term (0 for this f)
        a_n: array of cosine coefficients
        b_n: array of sine coefficients
    """
    a0 = 0.0
    a_n = np.zeros(n_modes + 1)
    b_n = np.zeros(n_modes + 1)

    # f(x) = sin(x) + 0.5*sin(3x)
    b_n[1] = 1.0
    b_n[3] = 0.5

    return a0, a_n, b_n


# -----------------------------------------------------------------------------
# Evaluate the Laplace equation solution
# -----------------------------------------------------------------------------
def laplace_solution(x, y, a0, a_n, b_n):
    """
    Evaluate the Laplace equation solution in the strip.

    u(x,y) = a_0(1-y) + Σ_{n=1}^N [a_n cos(nx) + b_n sin(nx)] sinh(n(1-y))/sinh(n)

    Parameters:
        x: x coordinates (2D array or 1D)
        y: y coordinates (2D array or 1D)
        a0, a_n, b_n: Fourier coefficients

    Returns:
        u: solution values (same shape as x, y)
    """
    # Start with constant mode
    u = a0 * (1 - y)

    n_modes = len(a_n) - 1
    for n in range(1, n_modes + 1):
        # Skip if both coefficients are zero
        if abs(a_n[n]) < 1e-15 and abs(b_n[n]) < 1e-15:
            continue

        # y-dependent factor: sinh(n(1-y)) / sinh(n)
        # Handle numerical stability for large n
        sinh_n = np.sinh(n)
        if sinh_n > 1e-10:
            y_factor = np.sinh(n * (1 - y)) / sinh_n
        else:
            y_factor = np.zeros_like(y)

        # Add contribution from this mode
        u += (a_n[n] * np.cos(n * x) + b_n[n] * np.sin(n * x)) * y_factor

    return u


# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
def main():
    # Compute Fourier coefficients
    a0, a_n, b_n = boundary_fourier_coefficients(N_MODES)

    # Create 2D grid
    x = np.linspace(0, 2 * np.pi, NX)
    y = np.linspace(0, 1, NY)
    X, Y = np.meshgrid(x, y)

    # Evaluate solution on grid
    U = laplace_solution(X, Y, a0, a_n, b_n)

    # Create figure
    fig, ax = plt.subplots(figsize=(6, 3.5))

    # Create a custom colormap based on book colors (navy to white to sky)
    # For harmonic functions, use a diverging colormap centered at 0
    colors_neg = [(20/255, 45/255, 110/255), (1, 1, 1)]  # Navy to white
    colors_pos = [(1, 1, 1), (120/255, 150/255, 210/255)]  # White to sky

    # Create custom diverging colormap
    n_colors = 256
    navy_rgb = np.array([20, 45, 110]) / 255
    white_rgb = np.array([1, 1, 1])
    sky_rgb = np.array([120, 150, 210]) / 255

    colors_list = []
    for i in range(n_colors // 2):
        t = i / (n_colors // 2 - 1)
        colors_list.append(tuple(navy_rgb * (1 - t) + white_rgb * t))
    for i in range(n_colors // 2):
        t = i / (n_colors // 2 - 1)
        colors_list.append(tuple(white_rgb * (1 - t) + sky_rgb * t))

    custom_cmap = LinearSegmentedColormap.from_list('navy_sky', colors_list)

    # Determine symmetric color limits
    vmax = np.max(np.abs(U))
    vmin = -vmax

    # Filled contour plot
    levels = np.linspace(vmin, vmax, 25)
    cf = ax.contourf(X, Y, U, levels=levels, cmap=custom_cmap, extend='both')

    # Add contour lines for emphasis
    cs = ax.contour(X, Y, U, levels=np.linspace(vmin, vmax, 13),
                    colors='#333333', linewidths=0.3, alpha=0.5)

    # Colorbar
    cbar = fig.colorbar(cf, ax=ax, shrink=0.9, pad=0.02)
    cbar.set_label(r'$u(x,y)$', fontsize=11)
    cbar.ax.tick_params(labelsize=9)

    # Axis labels
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$y$')

    # Custom x-ticks with π notation
    ax.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi])
    ax.set_xticklabels([r'$0$', r'$\frac{\pi}{2}$', r'$\pi$',
                        r'$\frac{3\pi}{2}$', r'$2\pi$'])

    # y-ticks
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])

    # Set aspect ratio to make the strip look appropriate
    ax.set_aspect('auto')

    # Minimalist spine styling
    ax.spines['top'].set_visible(True)
    ax.spines['right'].set_visible(True)
    for spine in ax.spines.values():
        spine.set_color('#444444')
        spine.set_linewidth(0.5)
    ax.tick_params(colors='#444444')

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
