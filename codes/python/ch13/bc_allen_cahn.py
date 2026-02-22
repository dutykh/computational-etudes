#!/usr/bin/env python3
"""
bc_allen_cahn.py
Chapter 13: Boundary Conditions in Spectral Methods

Allen-Cahn equation with fixed and time-dependent (driven) boundary conditions.

    u_t = eps * u_xx + u - u^3,   x in (-1, 1),   eps = 0.01

Initial condition: u(x, 0) = 0.53 x + 0.47 sin(-1.5 pi x)

Part A -- Fixed BCs: u(-1) = -1, u(1) = 1
    Lifting: w(x) = x satisfies w(-1) = -1, w(1) = 1
    v = u - x satisfies homogeneous Dirichlet BCs
    v_t = eps * v_xx + (v + x) - (v + x)^3
    IMEX: diffusion implicit, reaction explicit

Part B -- Driven BCs: u(-1) = -1, u(1) = 1 + sin^2(t/5)
    Method II (boundary bordering): overwrite boundary values each step

Generates:
    Figure 13.4: Allen-Cahn with fixed BCs (space-time plot)
    Figure 13.5: Allen-Cahn with driven BCs (space-time plot)

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
from scipy.linalg import lu_factor, lu_solve

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


# -----------------------------------------------------------------------------
# Problem parameters
# -----------------------------------------------------------------------------
EPS = 0.01


def initial_condition(x):
    """Initial condition: u(x, 0) = 0.53 x + 0.47 sin(-1.5 pi x)."""
    return 0.53 * x + 0.47 * np.sin(-1.5 * np.pi * x)


# -----------------------------------------------------------------------------
# Part A: Fixed BCs via lifting
# -----------------------------------------------------------------------------
def solve_allen_cahn_fixed(N=80, tmax=60.0, n_snapshots=300):
    """
    Solve Allen-Cahn with fixed BCs u(-1)=-1, u(1)=1 via lifting w(x) = x.

    The variable v = u - x satisfies v(+-1) = 0 and:
        v_t = eps * v_xx + (v + x) - (v + x)^3

    IMEX time stepping:
        Implicit: eps * D2 v^{n+1}
        Explicit: reaction R(v^n) = (v^n + x) - (v^n + x)^3

    Semi-discrete: (I - eps*dt*D2_int) v^{n+1} = v^n + dt * R(v^n)

    Parameters
    ----------
    N : int
        Number of Chebyshev intervals
    tmax : float
        Final time
    n_snapshots : int
        Number of snapshots to store

    Returns
    -------
    x : ndarray
        Chebyshev grid points
    t_save : ndarray
        Snapshot times
    U_save : ndarray
        Solution snapshots (n_snapshots x (N+1))
    """
    D2, D, x = cheb_second_derivative_matrix(N)

    # Time step: stability-limited
    dt = min(0.5 / N**2, 0.01) * EPS * 100  # Scale for eps
    dt = min(dt, 0.01)

    # Interior indices
    D2_int = D2[1:N, 1:N]
    I_int = np.eye(N - 1)

    # Implicit matrix: (I - eps * dt * D2_int)
    A = I_int - EPS * dt * D2_int
    A_lu = lu_factor(A)

    # Initial condition
    u0 = initial_condition(x)
    v = u0 - x  # v = u - w, where w(x) = x
    v_int = v[1:N].copy()

    # Interior grid points
    x_int = x[1:N]

    # Time stepping
    nsteps = int(np.ceil(tmax / dt))
    save_every = max(1, nsteps // n_snapshots)

    # Store full u = v + x
    u_full = np.zeros(N + 1)
    u_full[0] = 1.0   # u(1) = 1
    u_full[N] = -1.0   # u(-1) = -1
    u_full[1:N] = v_int + x_int

    U_save = [u_full.copy()]
    t_save = [0.0]

    t = 0.0
    for n in range(nsteps):
        # Reaction at interior points: R(v) = (v + x) - (v + x)^3
        u_int = v_int + x_int
        reaction = u_int - u_int**3

        # IMEX step: A v^{n+1} = v^n + dt * reaction
        rhs = v_int + dt * reaction
        v_int = lu_solve(A_lu, rhs)

        t += dt

        if (n + 1) % save_every == 0:
            u_full = np.zeros(N + 1)
            u_full[0] = 1.0
            u_full[N] = -1.0
            u_full[1:N] = v_int + x_int
            U_save.append(u_full.copy())
            t_save.append(t)

    return x, np.array(t_save), np.array(U_save)


# -----------------------------------------------------------------------------
# Part B: Driven BCs via Method II
# -----------------------------------------------------------------------------
def solve_allen_cahn_driven(N=80, tmax=60.0, n_snapshots=300):
    """
    Solve Allen-Cahn with driven BC: u(-1)=-1, u(1) = 1 + sin^2(t/5).

    Method II (boundary bordering): at each time step, overwrite boundary
    values of u after the implicit solve.

    IMEX on the full system with boundary rows overwritten:
        Interior: (I - eps*dt*D2) u^{n+1} = u^n + dt * (u^n - (u^n)^3)
        Then set u[0] = g_R(t), u[N] = g_L(t)

    Parameters
    ----------
    N : int
        Number of Chebyshev intervals
    tmax : float
        Final time
    n_snapshots : int
        Number of snapshots to store

    Returns
    -------
    x : ndarray
        Chebyshev grid points
    t_save : ndarray
        Snapshot times
    U_save : ndarray
        Solution snapshots
    """
    D2, D, x = cheb_second_derivative_matrix(N)

    # Time step
    dt = min(0.5 / N**2, 0.01) * EPS * 100
    dt = min(dt, 0.01)

    # Interior system
    D2_int = D2[1:N, 1:N]
    I_int = np.eye(N - 1)

    # Implicit matrix
    A = I_int - EPS * dt * D2_int
    A_lu = lu_factor(A)

    # Boundary-to-interior coupling columns
    # When BCs are non-zero, we need D2[1:N, 0] * u[0] + D2[1:N, N] * u[N]
    d2_col0 = D2[1:N, 0]   # coupling to x[0] = 1 (right boundary)
    d2_colN = D2[1:N, N]   # coupling to x[N] = -1 (left boundary)

    # Boundary functions
    def g_right(t):
        """u(1, t) = 1 + sin^2(t/5)."""
        return 1.0 + np.sin(t / 5.0)**2

    def g_left(t):
        """u(-1, t) = -1."""
        return -1.0

    # Initial condition
    u = initial_condition(x)
    u[0] = g_right(0.0)
    u[N] = g_left(0.0)

    # Time stepping
    nsteps = int(np.ceil(tmax / dt))
    save_every = max(1, nsteps // n_snapshots)

    U_save = [u.copy()]
    t_save = [0.0]

    t = 0.0
    for n in range(nsteps):
        # Explicit reaction at interior points
        u_int = u[1:N]
        reaction = u_int - u_int**3

        # Time-advanced boundary values
        t_new = t + dt
        u_R = g_right(t_new)
        u_L = g_left(t_new)

        # RHS: u_int + dt * reaction - dt * eps * (boundary coupling at new time)
        # The implicit step for the interior is:
        #   (I - eps*dt*D2_int) u_int^{n+1} = u_int^n + dt*R
        #       + eps*dt * (D2[1:N,0]*u_R + D2[1:N,N]*u_L)
        # because we moved boundary terms to the RHS
        rhs = u_int + dt * reaction + EPS * dt * (d2_col0 * u_R + d2_colN * u_L)

        u_int_new = lu_solve(A_lu, rhs)

        # Update full solution
        u[1:N] = u_int_new
        u[0] = u_R
        u[N] = u_L

        t = t_new

        if (n + 1) % save_every == 0:
            U_save.append(u.copy())
            t_save.append(t)

    return x, np.array(t_save), np.array(U_save)


def main():
    """Create Allen-Cahn space-time figures."""

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    N = 80

    # =========================================================================
    # Figure 1: Fixed BCs
    # =========================================================================
    print("Solving Allen-Cahn with fixed BCs...")
    x, t_save, U_save = solve_allen_cahn_fixed(N=N, tmax=60.0, n_snapshots=400)

    fig1, ax1 = plt.subplots(1, 1, figsize=(7, 5))

    # Sort x for proper pcolormesh display (Chebyshev points are cos-ordered)
    idx_sort = np.argsort(x)
    x_sorted = x[idx_sort]
    U_sorted = U_save[:, idx_sort]

    pcm1 = ax1.pcolormesh(x_sorted, t_save, U_sorted, cmap='RdBu_r',
                           shading='gouraud', vmin=-1.0, vmax=1.0)
    fig1.colorbar(pcm1, ax=ax1, shrink=0.8, label=r'$u(x, t)$')
    ax1.set_xlabel(r'$x$')
    ax1.set_ylabel(r'$t$')
    ax1.set_title(r'Allen-Cahn: Fixed BCs $u(\pm 1) = \pm 1$, $\varepsilon = 0.01$',
                  fontsize=12)

    plt.tight_layout()

    output_file1 = OUTPUT_DIR / 'allen_cahn_fixed.pdf'
    fig1.savefig(output_file1, bbox_inches='tight', pad_inches=0.05)
    fig1.savefig(output_file1.with_suffix('.png'), bbox_inches='tight',
                 pad_inches=0.05, dpi=300)

    print(f'Figure 1 saved to: {output_file1.resolve()}')
    plt.close(fig1)

    # =========================================================================
    # Figure 2: Driven BCs
    # =========================================================================
    print("Solving Allen-Cahn with driven BCs...")
    x2, t_save2, U_save2 = solve_allen_cahn_driven(N=N, tmax=60.0, n_snapshots=400)

    fig2, ax2 = plt.subplots(1, 1, figsize=(7, 5))

    idx_sort2 = np.argsort(x2)
    x_sorted2 = x2[idx_sort2]
    U_sorted2 = U_save2[:, idx_sort2]

    pcm2 = ax2.pcolormesh(x_sorted2, t_save2, U_sorted2, cmap='RdBu_r',
                           shading='gouraud', vmin=-1.0, vmax=2.0)
    fig2.colorbar(pcm2, ax=ax2, shrink=0.8, label=r'$u(x, t)$')
    ax2.set_xlabel(r'$x$')
    ax2.set_ylabel(r'$t$')
    ax2.set_title(r'Allen-Cahn: Driven BC $u(1) = 1 + \sin^2(t/5)$, $\varepsilon = 0.01$',
                  fontsize=12)

    plt.tight_layout()

    output_file2 = OUTPUT_DIR / 'allen_cahn_driven.pdf'
    fig2.savefig(output_file2, bbox_inches='tight', pad_inches=0.05)
    fig2.savefig(output_file2.with_suffix('.png'), bbox_inches='tight',
                 pad_inches=0.05, dpi=300)

    print(f'Figure 2 saved to: {output_file2.resolve()}')
    plt.close(fig2)

    # Print summary
    print(f"\nAllen-Cahn simulation parameters:")
    print(f"  N = {N}, eps = {EPS}")
    print(f"  Fixed BCs: {len(t_save)} snapshots, t in [{t_save[0]:.1f}, {t_save[-1]:.1f}]")
    print(f"  Driven BCs: {len(t_save2)} snapshots, t in [{t_save2[0]:.1f}, {t_save2[-1]:.1f}]")


if __name__ == '__main__':
    main()
