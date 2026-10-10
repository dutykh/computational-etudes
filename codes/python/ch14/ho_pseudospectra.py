#!/usr/bin/env python3
"""
ho_pseudospectra.py
Chapter 14: Higher-Order Boundary Value Problems

Computational Etude 14.8: Pseudospectra of the Orr-Sommerfeld Operator.

Computes the epsilon-pseudospectra of the Orr-Sommerfeld operator for
plane Poiseuille flow at the critical Reynolds number R = 5772,
streamwise wavenumber alpha = 1, using N = 100 Chebyshev points.

The Orr-Sommerfeld equation governs the linear stability of parallel
shear flows.  For the base flow U(y) = 1 - y^2 on [-1, 1], the
generalized eigenvalue problem reads A v = lambda B v, where A and B
encode the fourth-order ODE with no-slip boundary conditions
v(+-1) = v'(+-1) = 0.

The Dirichlet conditions are handled by restriction to interior
Chebyshev points.  The derivative (Neumann) conditions v'(+-1) = 0
are imposed by boundary bordering.  The resulting two algebraic
constraints are then eliminated analytically via a projection that
expresses the boundary unknowns in terms of the interior unknowns,
yielding a non-singular reduced system whose matrix C = B_red^{-1} A_red
is the one whose pseudospectra we compute.

The epsilon-pseudospectrum is defined as

    Lambda_eps(C) = { z in C : sigma_min(C - z I) < eps }

where sigma_min denotes the smallest singular value.  For normal
operators the pseudospectra are simply eps-balls around the
eigenvalues.  The Orr-Sommerfeld operator is highly non-normal, so
its pseudospectra extend far from the eigenvalues, revealing extreme
sensitivity to perturbations.

At R = 5772 (the critical Reynolds number for plane Poiseuille flow),
all eigenvalues lie barely in the stable half-plane Im(lambda) < 0.
The contour plot of log10(sigma_min) shows that perturbations of size
as small as 10^{-8} can push eigenvalues across the stability boundary
into Im(lambda) > 0, explaining the subcritical transition to
turbulence observed in pipe-flow experiments.

Reference: Trefethen, "Spectral Methods in MATLAB", Program 40.
           Trefethen & Embree, "Spectra and Pseudospectra" (2005).
           Reddy, Schmid & Henningson, SIAM J. Appl. Math. 53 (1993).

Author: Dr. Denys Dutykh
        Mathematics Department
        Khalifa University of Science and Technology
        Abu Dhabi, UAE

Last modified: 2026-02-22

Part of the book "Computational Etudes: A Spectral Approach"
https://github.com/dutykh/computational-etudes
"""

import numpy as np
from scipy.linalg import svdvals
import matplotlib.pyplot as plt
from matplotlib import rcParams
from pathlib import Path
import time

# Import Chebyshev matrix function from Chapter 7
import sys
sys.path.insert(0, str(Path(__file__).parent.parent / 'ch07'))
from cheb_matrix import cheb_matrix

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
OUTPUT_DIR = SCRIPT_DIR / '..' / '..' / '..' / 'textbook' / 'figures' / 'ch14' / 'python'


# -----------------------------------------------------------------------------
# Build Orr-Sommerfeld matrix C = B^{-1} A via constraint elimination
# -----------------------------------------------------------------------------
def build_orr_sommerfeld_matrix(N, R, alpha):
    """
    Build the matrix C = B^{-1} A for the Orr-Sommerfeld eigenvalue
    problem for plane Poiseuille flow U(y) = 1 - y^2.

    The Orr-Sommerfeld equation, written as a generalized eigenvalue
    problem A v = lambda B v, is:

        A = (D^2 - alpha^2)^2 / (i alpha R)
            - diag(U) (D^2 - alpha^2) + diag(U'')
        B = -(D^2 - alpha^2)

    restricted to interior Chebyshev points (Dirichlet BCs v(+-1) = 0).

    The Neumann conditions v'(+-1) = 0 are imposed by replacing the
    first and last rows with derivative constraints, making B singular.
    These two algebraic constraints are then eliminated analytically:
    the constraint rows express the boundary unknowns (v at interior
    points nearest +-1) as linear combinations of the remaining
    unknowns.  Substituting into the dynamical rows yields a reduced
    system with non-singular B_red, allowing C = B_red^{-1} A_red.

    Parameters
    ----------
    N : int
        Number of Chebyshev intervals (N+1 grid points, N-1 interior).
    R : float
        Reynolds number.
    alpha : float
        Streamwise wavenumber.

    Returns
    -------
    C : ndarray, shape (N-3, N-3)
        Matrix C = B_red^{-1} A_red for pseudospectral analysis.
    eigenvalues : ndarray
        Eigenvalues of C, sorted by imaginary part (descending).
    """
    D, x = cheb_matrix(N)

    # Second and fourth derivative matrices
    D2 = D @ D
    D4 = D2 @ D2

    # Number of interior points
    n = N - 1

    # Interior grid (remove boundary points x_0 = 1, x_N = -1)
    x_int = x[1:N]

    # Base flow: U(y) = 1 - y^2,  U''(y) = -2
    U = np.diag(1.0 - x_int**2)
    U2 = -2.0 * np.eye(n)

    # Interior blocks of D^2 and D^4
    D2_int = D2[1:N, 1:N]
    D4_int = D4[1:N, 1:N]

    # Shift operator: S = D^2 - alpha^2 I
    I_int = np.eye(n)
    S = D2_int - alpha**2 * I_int

    # S^2 = D^4 - 2 alpha^2 D^2 + alpha^4 I
    S2 = D4_int - 2.0 * alpha**2 * D2_int + alpha**4 * I_int

    # Build A and B (boundary-bordering form)
    A = S2 / (1j * alpha * R) - U @ S + U2
    B = -S.copy()

    # Replace first/last rows with Neumann conditions v'(+-1) = 0
    A[0, :] = D[0, 1:N]
    B[0, :] = 0.0
    A[-1, :] = D[N, 1:N]
    B[-1, :] = 0.0

    # --- Constraint elimination ---
    # Rows 0 and n-1 are algebraic constraints: A[0,:] v = 0, A[-1,:] v = 0
    # (since B[0,:] = B[-1,:] = 0, the eigenvalue lambda does not appear).
    #
    # Partition v = [v_0; v_free; v_{n-1}] where v_free = v[1:-1].
    # The constraints give:
    #   A[0,0]  v_0 + A[0,1:-1]  v_free + A[0,-1]  v_{n-1} = 0
    #   A[-1,0] v_0 + A[-1,1:-1] v_free + A[-1,-1] v_{n-1} = 0
    #
    # Solve for [v_0; v_{n-1}] = P @ v_free where P = -CC^{-1} F.

    CC = np.array([[A[0, 0], A[0, -1]],
                   [A[-1, 0], A[-1, -1]]])
    F = np.vstack([A[0, 1:-1], A[-1, 1:-1]])
    P = -np.linalg.solve(CC, F)      # shape (2, n-2)

    # Build the lifting matrix T: v = T @ v_free, shape (n, n-2)
    T = np.zeros((n, n - 2), dtype=complex)
    T[0, :] = P[0, :]                # v_0 row
    T[1:-1, :] = np.eye(n - 2)       # v_free rows
    T[-1, :] = P[1, :]               # v_{n-1} row

    # Reduced system: use the dynamical rows (1 to n-2) of A, B
    A_red = A[1:-1, :] @ T           # shape (n-2, n-2)
    B_red = B[1:-1, :] @ T           # shape (n-2, n-2)

    # Form C = B_red^{-1} A_red (B_red is now non-singular)
    C = np.linalg.solve(B_red, A_red)

    # Eigenvalues
    eigenvalues = np.linalg.eigvals(C)
    idx = np.argsort(-eigenvalues.imag)
    eigenvalues = eigenvalues[idx]

    return C, eigenvalues


# -----------------------------------------------------------------------------
# Pseudospectra computation
# -----------------------------------------------------------------------------
def compute_pseudospectra(C, x_re, x_im):
    """
    Compute log10(sigma_min(C - z I)) on a grid of complex numbers z.

    For each point z = x + iy on the grid, computes the smallest
    singular value of (C - zI).  The level set where this equals eps
    is the boundary of the eps-pseudospectrum.

    Parameters
    ----------
    C : ndarray, shape (M, M)
        The matrix whose pseudospectra are computed.
    x_re : ndarray, shape (nx,)
        Real parts of the grid points.
    x_im : ndarray, shape (ny,)
        Imaginary parts of the grid points.

    Returns
    -------
    sigma_grid : ndarray, shape (ny, nx)
        log10 of the minimum singular value at each grid point.
    """
    nx = len(x_re)
    ny = len(x_im)
    M = C.shape[0]
    I = np.eye(M, dtype=C.dtype)

    sigma_grid = np.zeros((ny, nx))

    total = nx * ny
    count = 0

    for j in range(ny):
        for i in range(nx):
            z = x_re[i] + 1j * x_im[j]
            # Compute smallest singular value of (C - z I)
            sv = svdvals(C - z * I)
            sigma_grid[j, i] = np.log10(max(sv[-1], 1e-20))
            count += 1

        # Progress report every 10 rows
        if (j + 1) % 10 == 0:
            pct = 100.0 * count / total
            print(f"    ... {pct:.0f}% complete ({count}/{total} grid points)")

    return sigma_grid


# -----------------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------------
def plot_pseudospectra(x_re, x_im, sigma_grid, eigenvalues, R, alpha, N):
    """
    Plot the pseudospectra as filled contours with eigenvalues overlaid.

    Parameters
    ----------
    x_re, x_im : ndarray
        Grid coordinates.
    sigma_grid : ndarray
        log10(sigma_min) values on the grid.
    eigenvalues : ndarray
        Eigenvalues of C (complex).
    R : float
        Reynolds number.
    alpha : float
        Streamwise wavenumber.
    N : int
        Number of Chebyshev intervals.

    Returns
    -------
    fig : matplotlib Figure
    """
    fig, ax = plt.subplots(figsize=(8, 6.5))

    # Contour levels: log10(eps) = -1, -2, ..., -8
    levels = np.arange(-8, 0)

    # Filled contours
    cf = ax.contourf(x_re, x_im, sigma_grid, levels=levels, cmap='YlOrRd_r',
                     extend='both')

    # Contour lines (black) for clarity
    cl = ax.contour(x_re, x_im, sigma_grid, levels=levels, colors='k',
                    linewidths=0.5, alpha=0.6)
    ax.clabel(cl, inline=True, fontsize=7, fmt='%d')

    # Colorbar
    cbar = fig.colorbar(cf, ax=ax, shrink=0.85, pad=0.02)
    cbar.set_label(r'$\log_{10}\,\sigma_{\min}(C - z\,I)$', fontsize=11)
    cbar.set_ticks(levels)
    cbar.set_ticklabels([str(lev) for lev in levels])

    # Overlay eigenvalues (only those within the plot region)
    in_region = ((eigenvalues.real >= x_re[0]) & (eigenvalues.real <= x_re[-1]) &
                 (eigenvalues.imag >= x_im[0]) & (eigenvalues.imag <= x_im[-1]))
    eigs_plot = eigenvalues[in_region]

    ax.plot(eigs_plot.real, eigs_plot.imag, 'o', color='white',
            markersize=4, markeredgecolor=NAVY, markeredgewidth=0.8,
            zorder=10, label='Eigenvalues')

    # Stability boundary: Im(lambda) = 0
    ax.axhline(0, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)

    ax.set_xlabel(r'$\mathrm{Re}(\lambda)$')
    ax.set_ylabel(r'$\mathrm{Im}(\lambda)$')
    ax.set_title(
        r'$\varepsilon$-Pseudospectra of Orr--Sommerfeld: '
        r'$R = %d$, $\alpha = %g$, $N = %d$' % (R, alpha, N),
        fontsize=12
    )
    ax.legend(loc='upper left', fontsize=9, framealpha=0.9)

    ax.set_xlim(x_re[0], x_re[-1])
    ax.set_ylim(x_im[0], x_im[-1])

    return fig


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def main():
    """Compute and plot pseudospectra of the Orr-Sommerfeld operator."""

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Physical parameters
    R = 5772           # Critical Reynolds number for plane Poiseuille flow
    alpha = 1.0        # Streamwise wavenumber
    N = 100            # Number of Chebyshev intervals

    print("Pseudospectra of the Orr-Sommerfeld Operator")
    print("=" * 60)
    print(f"  Reynolds number:      R = {R}")
    print(f"  Wavenumber:           alpha = {alpha}")
    print(f"  Chebyshev points:     N = {N}")
    print(f"  Interior points:      {N - 1}")
    print(f"  Matrix C size:        {N - 3} x {N - 3} (after constraint elimination)")
    print()

    # ---- Build reduced matrix C ----
    print("  Building Orr-Sommerfeld matrix C = B_red^{-1} A_red ...")
    C, eigenvalues = build_orr_sommerfeld_matrix(N, R, alpha)
    print(f"  Matrix C shape: {C.shape}")
    print(f"  Number of eigenvalues: {len(eigenvalues)}")

    if len(eigenvalues) > 0:
        lam_crit = eigenvalues[0]
        print(f"  Most unstable eigenvalue:")
        print(f"    lambda = {lam_crit.real:.10f} + {lam_crit.imag:.10f} i")
        if lam_crit.imag < 0:
            print(f"    (stable: Im(lambda) < 0)")
        else:
            print(f"    (UNSTABLE: Im(lambda) > 0)")

    # Print leading eigenvalues
    n_show = min(10, len(eigenvalues))
    print(f"\n  Leading {n_show} eigenvalues (sorted by Im):")
    print(f"  {'Re(lambda)':>16s}  {'Im(lambda)':>16s}")
    print("  " + "-" * 36)
    for k in range(n_show):
        print(f"  {eigenvalues[k].real:>16.10f}  {eigenvalues[k].imag:>16.10f}")
    print()

    # ---- Define grid for pseudospectra ----
    # The eigenvalues lie in Re ~ [0.2, 1.0], Im ~ [-1.0, 0],
    # so we center the grid on this region with some padding.
    nx, ny = 120, 120
    x_re = np.linspace(-0.1, 1.1, nx)
    x_im = np.linspace(-1.0, 0.1, ny)

    print(f"  Computing pseudospectra on {nx} x {ny} grid ...")
    print(f"  Grid: Re in [{x_re[0]}, {x_re[-1]}], Im in [{x_im[0]}, {x_im[-1]}]")
    print(f"  This requires {nx * ny} SVD computations of {C.shape[0]}x{C.shape[0]} matrices.")
    print()

    t_start = time.time()
    sigma_grid = compute_pseudospectra(C, x_re, x_im)
    t_elapsed = time.time() - t_start
    print(f"\n  Pseudospectra computed in {t_elapsed:.1f} seconds.")
    print()

    # ---- Plot ----
    print("  Generating figure ...")
    fig = plot_pseudospectra(x_re, x_im, sigma_grid, eigenvalues, R, alpha, N)

    plt.tight_layout()

    # Save
    output_file = OUTPUT_DIR / 'pseudospectra.pdf'
    fig.savefig(output_file, bbox_inches='tight', pad_inches=0.05)
    fig.savefig(output_file.with_suffix('.png'), bbox_inches='tight',
                pad_inches=0.05, dpi=300)
    print(f"  Figure saved to: {output_file.resolve()}")
    plt.close(fig)

    print()
    print("=" * 60)


if __name__ == '__main__':
    main()
