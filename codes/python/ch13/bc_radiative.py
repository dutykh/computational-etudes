#!/usr/bin/env python3
"""
bc_radiative.py
Chapter 13: Boundary Conditions in Spectral Methods

Radiative cooling at a boundary: nonlinear Robin BC with Newton iteration.

    u_t = u_xx,         x in [0, 1],  t > 0
    u(0, t) = 1                        (Dirichlet: fixed temperature)
    u_x(1, t) + sigma * u(1, t)^4 = 0  (Nonlinear Robin: radiative cooling)

Initial condition: u(x, 0) = 1

The domain [0, 1] is mapped to the Chebyshev domain [-1, 1] via:
    x = (xi + 1)/2,   D_x = 2 * D_xi

Chebyshev ordering: xi_0 = 1 maps to physical x = 1 (radiating end),
                     xi_N = -1 maps to physical x = 0 (fixed end).

Time stepping: backward Euler with Newton iteration at each step.
The coupled system (PDE + nonlinear BC) is:
    Interior rows:   u - dt * u_xx = u_old      (backward Euler)
    Row 0 (x=1):     D_x[0,:] u + sigma * u[0]^4 = 0  (radiation BC)
    Row N (x=0):     u[N] = 1                           (Dirichlet BC)

Newton iteration solves F(u) = 0 with Jacobian J = dF/du at each step.

Generates Figure 13.6: Radiative cooling profiles and boundary temperature.

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

# Import Chebyshev matrix function from Chapter 7
import sys
sys.path.insert(0, str(Path(__file__).parent.parent / 'ch07'))
from cheb_matrix import cheb_matrix, cheb_second_derivative_matrix

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
OUTPUT_DIR = SCRIPT_DIR / '..' / '..' / '..' / 'textbook' / 'figures' / 'ch13' / 'python'


def solve_radiative_cooling(N=32, sigma=1.0, tmax=5.0, dt=0.005,
                             newton_tol=1e-10, newton_max=20):
    """
    Solve the radiative cooling problem using backward Euler + Newton.

    Domain mapping: physical [0, 1] -> Chebyshev [-1, 1]
        x_phys = (xi + 1)/2,  so D_phys = 2 * D_xi,  D2_phys = 4 * D2_xi

    Chebyshev ordering:
        xi_0 = 1 -> x_phys = 1 (radiating boundary)
        xi_N = -1 -> x_phys = 0 (fixed boundary, u = 1)

    At each time step, we solve F(u^{n+1}) = 0 where:
        F_i = u_i - u_i^old - dt * (D2_phys u)_i    for i = 1, ..., N-1
        F_0 = (D_phys u)_0 + sigma * u_0^4           (radiation BC at x=1)
        F_N = u_N - 1                                  (Dirichlet at x=0)

    Jacobian:
        J_i = (I - dt * D2_phys)_i   for interior rows
        J_0 = D_phys[0, :] + 4*sigma*u_0^3 * e_0^T
        J_N = e_N^T

    Parameters
    ----------
    N : int
        Number of Chebyshev intervals
    sigma : float
        Radiative cooling parameter (Stefan-Boltzmann coefficient)
    tmax : float
        Final time
    dt : float
        Time step
    newton_tol : float
        Newton convergence tolerance
    newton_max : int
        Maximum Newton iterations per step

    Returns
    -------
    x_phys : ndarray
        Physical grid points in [0, 1]
    t_all : list
        All time values
    U_all : list
        Solution at each time step
    """
    # Build Chebyshev matrices on [-1, 1]
    D_xi, xi = cheb_matrix(N)
    D2_xi = D_xi @ D_xi

    # Map to physical domain [0, 1]
    x_phys = (xi + 1.0) / 2.0
    D_phys = 2.0 * D_xi
    D2_phys = 4.0 * D2_xi

    # Identity matrix
    I_mat = np.eye(N + 1)

    # Initial condition: u = 1 everywhere
    u = np.ones(N + 1)

    # Time stepping
    nsteps = int(np.ceil(tmax / dt))

    t_all = [0.0]
    U_all = [u.copy()]

    t = 0.0
    for step in range(nsteps):
        u_old = u.copy()
        t_new = t + dt

        # Newton iteration
        for it in range(newton_max):
            # Residual F(u)
            F = np.zeros(N + 1)

            # Interior: F_i = u_i - u_old_i - dt * (D2_phys @ u)_i
            F[1:N] = u[1:N] - u_old[1:N] - dt * (D2_phys @ u)[1:N]

            # Radiation BC at x_phys = 1 (row 0): D_phys[0,:]*u + sigma*u[0]^4
            F[0] = D_phys[0, :] @ u + sigma * u[0]**4

            # Dirichlet BC at x_phys = 0 (row N): u[N] - 1
            F[N] = u[N] - 1.0

            # Check convergence
            res = np.max(np.abs(F))
            if res < newton_tol:
                break

            # Jacobian
            J = np.zeros((N + 1, N + 1))

            # Interior rows: J = I - dt * D2_phys
            J[1:N, :] = I_mat[1:N, :] - dt * D2_phys[1:N, :]

            # Radiation BC row: J[0,:] = D_phys[0,:]; J[0,0] += 4*sigma*u[0]^3
            J[0, :] = D_phys[0, :]
            J[0, 0] += 4.0 * sigma * u[0]**3

            # Dirichlet row: J[N,:] = e_N
            J[N, :] = 0.0
            J[N, N] = 1.0

            # Newton update
            delta_u = np.linalg.solve(J, -F)
            u += delta_u

        t = t_new
        t_all.append(t)
        U_all.append(u.copy())

    return x_phys, t_all, U_all


def main():
    """Create radiative cooling figure."""

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    N = 32
    tmax = 5.0
    dt = 0.005

    # =========================================================================
    # Run simulations for multiple sigma values
    # =========================================================================
    sigma_values = [0.1, 1.0, 10.0]
    results = {}

    for sigma in sigma_values:
        print(f"Solving radiative cooling with sigma = {sigma}...")
        x_phys, t_all, U_all = solve_radiative_cooling(
            N=N, sigma=sigma, tmax=tmax, dt=dt)
        results[sigma] = (x_phys, t_all, U_all)

    # =========================================================================
    # Figure: 2 panels
    # =========================================================================
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    # ---- Panel 1: Temperature profiles at selected times for sigma=1 ----
    ax1 = axes[0]
    sigma_plot = 1.0
    x_phys, t_all, U_all = results[sigma_plot]

    # Sort by physical x for plotting
    idx_sort = np.argsort(x_phys)
    x_sorted = x_phys[idx_sort]

    t_targets = [0.0, 0.1, 0.5, 1.0, 5.0]
    colors_t = [NAVY, SKY, TEAL, CORAL, PURPLE]

    for t_target, color in zip(t_targets, colors_t):
        # Find closest time index
        idx_t = np.argmin(np.abs(np.array(t_all) - t_target))
        t_actual = t_all[idx_t]
        u_snapshot = U_all[idx_t][idx_sort]

        ax1.plot(x_sorted, u_snapshot, '-', color=color, linewidth=1.5,
                 label=f'$t = {t_actual:.1f}$')

    ax1.set_xlabel(r'$x$')
    ax1.set_ylabel(r'$u(x, t)$')
    ax1.set_title(rf'Temperature Profiles ($\sigma = {sigma_plot:.0f}$)', fontsize=11)
    ax1.legend(loc='lower left', fontsize=9)
    ax1.set_xlim(-0.02, 1.02)
    ax1.set_ylim(0, 1.05)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.yaxis.grid(True, linestyle='-', alpha=0.2, color='#888888')
    ax1.xaxis.grid(True, linestyle='-', alpha=0.2, color='#888888')

    # ---- Panel 2: Boundary temperature u(1,t) vs time for different sigma ----
    ax2 = axes[1]
    colors_sigma = [TEAL, CORAL, PURPLE]

    for sigma, color in zip(sigma_values, colors_sigma):
        x_phys, t_all, U_all = results[sigma]
        # u(1, t) is at x_phys = 1, which is xi = 1, i.e., index 0
        u_boundary = np.array([U_all[i][0] for i in range(len(t_all))])
        ax2.plot(t_all, u_boundary, '-', color=color, linewidth=1.5,
                 label=rf'$\sigma = {sigma}$')

    ax2.set_xlabel(r'$t$')
    ax2.set_ylabel(r'$u(1, t)$')
    ax2.set_title('Radiating Boundary Temperature', fontsize=11)
    ax2.legend(loc='upper right', fontsize=9)
    ax2.set_xlim(0, tmax)
    ax2.set_ylim(0, 1.05)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.yaxis.grid(True, linestyle='-', alpha=0.2, color='#888888')
    ax2.xaxis.grid(True, linestyle='-', alpha=0.2, color='#888888')

    fig.suptitle(r'Radiative Cooling: $u_t = u_{xx}$, $u(0) = 1$, $u_x(1) + \sigma u(1)^4 = 0$',
                 fontsize=13, y=0.98)

    plt.tight_layout(rect=[0, 0, 1, 0.94])

    # Save figure
    output_file = OUTPUT_DIR / 'radiative_cooling.pdf'
    fig.savefig(output_file, bbox_inches='tight', pad_inches=0.05)
    fig.savefig(output_file.with_suffix('.png'), bbox_inches='tight',
                pad_inches=0.05, dpi=300)

    print(f'Figure saved to: {output_file.resolve()}')
    plt.close(fig)

    # Print summary table
    print("\nRadiative Cooling: Steady-State Boundary Temperature u(1, t_final)")
    print("-" * 50)
    print(f"{'sigma':>10} {'u(1, t_final)':>18}")
    print("-" * 50)
    for sigma in sigma_values:
        x_phys, t_all, U_all = results[sigma]
        u_bnd_final = U_all[-1][0]
        print(f"{sigma:>10.1f} {u_bnd_final:>18.10f}")
    print("-" * 50)


if __name__ == '__main__':
    main()
