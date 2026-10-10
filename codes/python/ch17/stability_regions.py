#!/usr/bin/env python3
"""
stability_regions.py
Chapter 17: Stability of Spectral Methods

Computational Etude 17.2: Stability Regions and Spectral Eigenvalue Overlay.

Two figures:
  (1) stability_regions.pdf -- 2x3 grid showing stability regions of six
      classical time-stepping methods: Forward Euler, RK4, Leapfrog (1st
      order), Adams-Bashforth 2, Crank-Nicolson, Backward Euler.
  (2) spectra_overlay.pdf -- 1x2 figure overlaying spectral eigenvalues
      on the stability regions of RK4 (imaginary eigenvalues of Fourier D)
      and Forward Euler (real eigenvalues of Chebyshev D^2).

Generates Figures 17.2a,b: stability_regions.pdf, spectra_overlay.pdf

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

rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Computer Modern Roman', 'CMU Serif', 'DejaVu Serif']
rcParams['mathtext.fontset'] = 'cm'
rcParams['axes.labelsize'] = 11
rcParams['axes.titlesize'] = 12
rcParams['xtick.labelsize'] = 10
rcParams['ytick.labelsize'] = 10
rcParams['legend.fontsize'] = 9
rcParams['figure.dpi'] = 150
rcParams['savefig.dpi'] = 300
rcParams['text.usetex'] = False
rcParams['axes.linewidth'] = 0.8
rcParams['xtick.major.width'] = 0.8
rcParams['ytick.major.width'] = 0.8

NAVY = '#142D6E'
SKY = '#7896D2'
CORAL = '#E74C3C'
TEAL = '#16A085'
AMBER = '#E67E22'
LIGHT_NAVY = '#C5D0E8'

SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / '..' / '..' / '..' / 'textbook' / 'figures' / 'ch17' / 'python'


def stability_function_fe(z):
    """Forward Euler: g(z) = 1 + z."""
    return 1.0 + z


def stability_function_rk4(z):
    """RK4: g(z) = 1 + z + z^2/2 + z^3/6 + z^4/24."""
    return 1.0 + z + z ** 2 / 2.0 + z ** 3 / 6.0 + z ** 4 / 24.0


def stability_function_be(z):
    """Backward Euler: g(z) = 1/(1 - z)."""
    return 1.0 / (1.0 - z)


def stability_function_cn(z):
    """Crank-Nicolson: g(z) = (1 + z/2)/(1 - z/2)."""
    return (1.0 + z / 2.0) / (1.0 - z / 2.0)


def ab2_boundary(n_pts=500):
    """Adams-Bashforth 2: g^2 - g - z*(3g/2 - 1/2) = 0.
    On |g|=1, g = e^{i*theta}, solve for z(theta)."""
    theta = np.linspace(0, 2 * np.pi, n_pts)
    g = np.exp(1j * theta)
    z = (g ** 2 - g) / (1.5 * g - 0.5)
    return z.real, z.imag


def leapfrog_boundary(n_pts=500):
    """Leapfrog (midpoint): g^2 - 1 - 2*z*g = 0.
    On |g|=1, g = e^{i*theta}, z = (g^2 - 1)/(2g)."""
    theta = np.linspace(0, 2 * np.pi, n_pts)
    g = np.exp(1j * theta)
    z = (g ** 2 - 1.0) / (2.0 * g)
    return z.real, z.imag


def plot_stability_region(ax, method, xlim, ylim, title, n_grid=400):
    """Plot stability region for a given method."""
    x = np.linspace(xlim[0], xlim[1], n_grid)
    y = np.linspace(ylim[0], ylim[1], n_grid)
    X, Y = np.meshgrid(x, y)
    Z = X + 1j * Y

    if method == 'fe':
        G = np.abs(stability_function_fe(Z))
        ax.contourf(X, Y, G, levels=[0, 1], colors=[LIGHT_NAVY], alpha=0.5)
        ax.contour(X, Y, G, levels=[1], colors=[NAVY], linewidths=1.5)
    elif method == 'rk4':
        G = np.abs(stability_function_rk4(Z))
        ax.contourf(X, Y, G, levels=[0, 1], colors=[LIGHT_NAVY], alpha=0.5)
        ax.contour(X, Y, G, levels=[1], colors=[NAVY], linewidths=1.5)
    elif method == 'be':
        G = np.abs(stability_function_be(Z))
        # Stable region is OUTSIDE the circle |z-1| = 1, i.e., |g| <= 1
        # For BE, |g(z)| <= 1 means |1/(1-z)| <= 1, i.e., |1-z| >= 1
        # Shade the stable region (complement of disk)
        ax.contourf(X, Y, G, levels=[0, 1], colors=[LIGHT_NAVY], alpha=0.5)
        ax.contour(X, Y, G, levels=[1], colors=[NAVY], linewidths=1.5)
    elif method == 'cn':
        G = np.abs(stability_function_cn(Z))
        ax.contourf(X, Y, G, levels=[0, 1], colors=[LIGHT_NAVY], alpha=0.5)
        ax.contour(X, Y, G, levels=[1], colors=[NAVY], linewidths=1.5)
    elif method == 'ab2':
        bx, by = ab2_boundary(1000)
        ax.fill(bx, by, color=LIGHT_NAVY, alpha=0.5)
        ax.plot(bx, by, color=NAVY, linewidth=1.5)
    elif method == 'leapfrog':
        bx, by = leapfrog_boundary(1000)
        ax.fill(bx, by, color=LIGHT_NAVY, alpha=0.5)
        ax.plot(bx, by, color=NAVY, linewidth=1.5)

    ax.axhline(0, color='gray', linewidth=0.5, zorder=0)
    ax.axvline(0, color='gray', linewidth=0.5, zorder=0)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_xlabel(r'Re$(z)$')
    ax.set_ylabel(r'Im$(z)$')
    ax.set_title(title)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.2, linewidth=0.5)


def cheb_diff_matrix(N):
    """Chebyshev differentiation matrix."""
    if N == 0:
        return np.array([[0.0]]), np.array([1.0])
    x = np.cos(np.pi * np.arange(N + 1) / N)
    c = np.ones(N + 1)
    c[0] = 2.0
    c[N] = 2.0
    c *= (-1.0) ** np.arange(N + 1)
    X = np.outer(x, np.ones(N + 1))
    dX = X - X.T + np.eye(N + 1)
    D = np.outer(c, 1.0 / c) / dX
    D -= np.diag(D.sum(axis=1))
    return D, x


def figure_stability_regions():
    """Generate 2x3 grid of stability regions."""
    fig, axes = plt.subplots(2, 3, figsize=(12, 8))

    plot_stability_region(axes[0, 0], 'fe', (-3, 1.5), (-2, 2),
                          '(a) Forward Euler')
    plot_stability_region(axes[0, 1], 'rk4', (-5, 2), (-4, 4),
                          '(b) RK4')
    plot_stability_region(axes[0, 2], 'leapfrog', (-1.5, 1.5), (-1.5, 1.5),
                          '(c) Leapfrog')
    plot_stability_region(axes[1, 0], 'ab2', (-2, 1.5), (-1.5, 1.5),
                          '(d) Adams-Bashforth 2')
    plot_stability_region(axes[1, 1], 'cn', (-5, 1), (-4, 4),
                          '(e) Crank-Nicolson')
    plot_stability_region(axes[1, 2], 'be', (-1.5, 5), (-4, 4),
                          '(f) Backward Euler')

    plt.tight_layout()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_DIR / 'stability_regions.pdf', bbox_inches='tight')
    fig.savefig(OUTPUT_DIR / 'stability_regions.png', bbox_inches='tight')
    print(f"Figure saved to {OUTPUT_DIR / 'stability_regions.pdf'}")


def figure_spectra_overlay():
    """Overlay spectral eigenvalues on stability regions."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5))

    # --- Left: RK4 + Fourier imaginary eigenvalues ---
    n_grid = 400
    x = np.linspace(-5, 2, n_grid)
    y = np.linspace(-4, 4, n_grid)
    X, Y = np.meshgrid(x, y)
    Z = X + 1j * Y
    G = np.abs(stability_function_rk4(Z))
    ax1.contourf(X, Y, G, levels=[0, 1], colors=[LIGHT_NAVY], alpha=0.4)
    ax1.contour(X, Y, G, levels=[1], colors=[NAVY], linewidths=1.5)

    # Fourier eigenvalues for advection: lambda_k = i*k, k = -N/2+1..N/2
    N_f = 32
    k = np.arange(-N_f // 2 + 1, N_f // 2 + 1)
    lam_f = 1j * k  # eigenvalues of -d/dx in Fourier

    # Stable dt
    dt_s = 2.8 / (N_f / 2)
    eigs_s = dt_s * lam_f
    ax1.plot(eigs_s.real, eigs_s.imag, 'o', color=TEAL, markersize=4,
             label=fr'$\Delta t = {dt_s:.4f}$ (stable)', zorder=5)

    # Unstable dt
    dt_u = 3.5 / (N_f / 2)
    eigs_u = dt_u * lam_f
    ax1.plot(eigs_u.real, eigs_u.imag, 's', color=CORAL, markersize=4,
             label=fr'$\Delta t = {dt_u:.4f}$ (unstable)', zorder=5)

    ax1.axhline(0, color='gray', linewidth=0.5, zorder=0)
    ax1.axvline(0, color='gray', linewidth=0.5, zorder=0)
    ax1.set_xlabel(r'Re$(\Delta t \lambda)$')
    ax1.set_ylabel(r'Im$(\Delta t \lambda)$')
    ax1.set_title(r'(a) RK4 + Fourier advection ($N=32$)')
    ax1.set_aspect('equal')
    ax1.legend(loc='lower left', fontsize=8)
    ax1.grid(True, alpha=0.2, linewidth=0.5)

    # --- Right: Forward Euler + Chebyshev D^2 eigenvalues ---
    # Strategy: show the FE disk at a readable scale, place the physical
    # eigenvalues inside it, and use an arrow pointing left to indicate
    # that outlier eigenvalues are far off-screen at dt*lambda ~ -3000.
    N_c = 16
    D, xc = cheb_diff_matrix(N_c)
    D2 = D @ D
    D2_int = D2[1:N_c, 1:N_c]
    eigs_c = np.sort(np.linalg.eigvals(D2_int).real)  # all real, negative

    # Choose dt at the CFL limit: dt = 2 / rho(D2)
    rho = np.abs(eigs_c[0])
    dt_c = 2.0 / rho  # just at the stability boundary

    scaled_eigs = dt_c * eigs_c  # all in [-2, 0]

    # Plot the FE stability disk
    x2 = np.linspace(-3.5, 2, n_grid)
    y2 = np.linspace(-2.5, 2.5, n_grid)
    X2, Y2 = np.meshgrid(x2, y2)
    Z2 = X2 + 1j * Y2
    G2 = np.abs(stability_function_fe(Z2))
    ax2.contourf(X2, Y2, G2, levels=[0, 1], colors=[LIGHT_NAVY], alpha=0.4)
    ax2.contour(X2, Y2, G2, levels=[1], colors=[NAVY], linewidths=1.5)

    # At the CFL limit, all eigenvalues fit inside [-2, 0].
    # But we want to show what happens with a LARGER dt that puts
    # only the physical modes inside.  Use dt = 10 * dt_cfl.
    dt_show = 10 * dt_c
    scaled_show = dt_show * eigs_c

    # Physical modes: those that still fit in the disk (|scaled| <= 2)
    phys_mask = np.abs(scaled_show) <= 2.2
    out_mask = ~phys_mask

    ax2.plot(scaled_show[phys_mask], np.zeros(phys_mask.sum()), 'o',
             color=TEAL, markersize=6, label='physical modes', zorder=5)

    # Show a few outliers that are just outside the disk
    near_out = scaled_show[out_mask & (np.abs(scaled_show) < 6)]
    far_out = scaled_show[out_mask & (np.abs(scaled_show) >= 6)]
    if len(near_out) > 0:
        ax2.plot(near_out, np.zeros(len(near_out)), 's',
                 color=CORAL, markersize=6, zorder=5, label='outlier modes')

    # Arrow pointing left to indicate far outliers off-screen
    worst_val = dt_show * eigs_c[0]
    ax2.annotate('', xy=(-3.4, 0), xytext=(-2.5, 0),
                 arrowprops=dict(arrowstyle='->', color=CORAL, lw=2.5))
    ax2.text(-3.3, 0.5,
             fr'outliers at $\Delta t\lambda \approx {worst_val:.0f}$',
             fontsize=8, color=CORAL, ha='left',
             bbox=dict(facecolor='white', edgecolor='none', alpha=0.8))

    ax2.axhline(0, color='gray', linewidth=0.5, zorder=0)
    ax2.axvline(0, color='gray', linewidth=0.5, zorder=0)
    ax2.set_xlabel(r'Re$(\Delta t \lambda)$')
    ax2.set_ylabel(r'Im$(\Delta t \lambda)$')
    ax2.set_title(r'(b) Forward Euler + Chebyshev $D^2$ ($N=16$)')
    ax2.set_aspect('equal')
    ax2.legend(loc='upper right', fontsize=8)
    ax2.grid(True, alpha=0.2, linewidth=0.5)

    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / 'spectra_overlay.pdf', bbox_inches='tight')
    fig.savefig(OUTPUT_DIR / 'spectra_overlay.png', bbox_inches='tight')
    print(f"Figure saved to {OUTPUT_DIR / 'spectra_overlay.pdf'}")


def main():
    figure_stability_regions()
    figure_spectra_overlay()
    plt.show()


if __name__ == '__main__':
    main()
