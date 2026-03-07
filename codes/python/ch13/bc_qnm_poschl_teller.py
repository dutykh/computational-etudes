#!/usr/bin/env python3
"""
bc_qnm_poschl_teller.py
Chapter 13: Boundary Conditions in Spectral Methods

Quasinormal modes (QNMs) of a black hole modelled by the Poeschl-Teller
potential via Chebyshev spectral collocation.

The Regge-Wheeler-like equation reads:

    psi'' + (omega^2 - V(x)) psi = 0,   x in [-L, L]

with the Poeschl-Teller potential V(x) = V_0 sech^2(x).

Boundary conditions encode purely outgoing/ingoing waves:
    psi'(-L) + i omega psi(-L) = 0   (ingoing at left)
    psi'( L) - i omega psi( L) = 0   (outgoing at right)

The substitution lambda = -i omega converts this into a quadratic
eigenvalue problem (QEP):

    (M_0 + lambda M_1 + lambda^2 M_2) psi = 0

which is linearised via the companion matrix approach and solved with
scipy.linalg.eig.

The Poeschl-Teller potential with V_0 = 2 has an exact fundamental QNM:
    omega_0 = sqrt(7)/2 - i/2

Parameters: V_0 = 2, L = 8, N = 80.

Generates Figures 13.6a-c: QNM spectrum, eigenfunctions, and convergence.

Author: Dr. Denys Dutykh
        Mathematics Department
        Khalifa University of Science and Technology
        Abu Dhabi, UAE

Part of the book "Computational Etudes: A Spectral Approach"
https://github.com/dutykh/computational-etudes
"""

import numpy as np
import scipy.linalg
import matplotlib.pyplot as plt
from matplotlib import rcParams
from pathlib import Path

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
OUTPUT_DIR = SCRIPT_DIR / '..' / '..' / '..' / 'textbook' / 'figures' / 'ch13' / 'python'

# Physical parameters
V0 = 2.0                                      # Potential amplitude
OMEGA_EXACT = np.sqrt(7.0) / 2.0 - 0.5j       # Exact fundamental QNM


# -----------------------------------------------------------------------------
# Potential
# -----------------------------------------------------------------------------
def poschl_teller(x, V0=V0):
    """Poeschl-Teller potential V(x) = V_0 sech^2(x)."""
    return V0 / np.cosh(x) ** 2


# -----------------------------------------------------------------------------
# QNM solver
# -----------------------------------------------------------------------------
def solve_qnm(N, L, V0=V0):
    """
    Compute quasinormal modes via a quadratic eigenvalue problem.

    The wave equation  psi'' + (omega^2 - V) psi = 0  with radiation
    BCs is rewritten with lambda = -i omega as

        (M_0 + lambda M_1 + lambda^2 M_2) psi = 0

    and then linearised using the companion matrix approach.

    Parameters
    ----------
    N : int
        Number of Chebyshev intervals.
    L : float
        Half-length of the computational domain [-L, L].
    V0 : float
        Potential amplitude.

    Returns
    -------
    omega : ndarray, shape (2*N+2,)
        Computed quasinormal frequencies (complex).
    psi_all : ndarray, shape (N+1, 2*N+2)
        Eigenvectors (columns are eigenfunctions on the physical grid).
    x_phys : ndarray, shape (N+1,)
        Physical grid points.
    """
    n = N + 1  # total number of collocation points

    # Chebyshev grid and first derivative on [-1, 1]
    D_cheb, x_cheb = cheb_matrix(N)

    # Scale to [-L, L]
    x_phys = L * x_cheb
    D = D_cheb / L        # first derivative on [-L, L]
    D2 = D @ D             # second derivative on [-L, L]

    # Potential at grid points
    V = poschl_teller(x_phys, V0)

    # -----------------------------------------------------------------
    # Assemble QEP matrices  (M0 + lambda M1 + lambda^2 M2) psi = 0
    # -----------------------------------------------------------------
    M0 = np.zeros((n, n), dtype=complex)
    M1 = np.zeros((n, n), dtype=complex)
    M2 = np.zeros((n, n), dtype=complex)

    # Interior rows (i = 1, ..., N-1):
    #   psi'' - V psi - lambda^2 psi = 0
    #   => M0 = D2 - diag(V),  M1 = 0,  M2 = -I
    for i in range(1, N):
        M0[i, :] = D2[i, :]
        M0[i, i] -= V[i]
        M2[i, i] = -1.0

    # Boundary row 0: x_0 = +L  (outgoing wave)
    #   psi'(L) + lambda psi(L) = 0
    #   => M0[0,:] = D[0,:],  M1[0,0] = 1
    M0[0, :] = D[0, :]
    M1[0, 0] = 1.0

    # Boundary row N: x_N = -L  (ingoing wave)
    #   psi'(-L) - lambda psi(-L) = 0
    #   => M0[N,:] = D[N,:],  M1[N,N] = -1
    M0[N, :] = D[N, :]
    M1[N, N] = -1.0

    # -----------------------------------------------------------------
    # Companion linearisation:
    #   A v = lambda B v
    #
    #   A = [[M0,  M1],    B = [[ 0, -M2],
    #        [ 0,   I ]]        [ I,   0 ]]
    # -----------------------------------------------------------------
    Z = np.zeros((n, n), dtype=complex)
    I = np.eye(n, dtype=complex)

    A = np.block([[M0, M1],
                  [Z,  I]])
    B = np.block([[Z,  -M2],
                  [I,   Z]])

    # Solve generalised eigenvalue problem
    lam, vecs = scipy.linalg.eig(A, B)

    # Convert: omega = i * lambda  (suppress warnings from NaN/Inf eigenvalues)
    with np.errstate(invalid='ignore'):
        omega = 1j * lam

    # Extract the physical part (first n components) of each eigenvector
    psi_all = vecs[:n, :]

    return omega, psi_all, x_phys


# -----------------------------------------------------------------------------
# Mode classification
# -----------------------------------------------------------------------------
def classify_modes(omega, max_im=-0.01, max_abs_im=5.0, max_re=10.0):
    """
    Separate physical QNMs from spurious numerical modes.

    Physical QNMs have Im(omega) < 0 (decaying in time) and are not
    too far from the origin.  Modes with large |Re(omega)| are typically
    boundary artefacts from the domain truncation.

    Parameters
    ----------
    omega : ndarray
        Array of computed complex frequencies.
    max_im : float
        Upper bound on Im(omega) for physical modes.
    max_abs_im : float
        Maximum |Im(omega)| to exclude very large spurious modes.
    max_re : float
        Maximum |Re(omega)| to keep.

    Returns
    -------
    phys_idx : ndarray of int
        Indices of physical modes.
    spur_idx : ndarray of int
        Indices of spurious modes.
    """
    finite_mask = np.isfinite(omega)
    phys_mask = (finite_mask &
                 (omega.imag < max_im) &
                 (np.abs(omega.imag) < max_abs_im) &
                 (np.abs(omega.real) < max_re))

    phys_idx = np.where(phys_mask)[0]
    spur_idx = np.where(~phys_mask & finite_mask)[0]

    return phys_idx, spur_idx


def find_fundamental(omega, phys_idx, omega_exact=OMEGA_EXACT):
    """
    Find the fundamental mode (n=0) closest to the known exact value.

    Among the physical modes, the fundamental is identified as the one
    nearest to the exact QNM frequency.  This avoids confusing it with
    spurious boundary artefacts that may have smaller |Im(omega)|.

    Parameters
    ----------
    omega : ndarray
    phys_idx : ndarray of int
    omega_exact : complex
        Reference value for the fundamental mode.

    Returns
    -------
    idx0 : int
        Index into omega of the fundamental mode.
    """
    dist = np.abs(omega[phys_idx] - omega_exact)
    best = phys_idx[np.argmin(dist)]
    return best


def find_overtone(omega, phys_idx, n=1):
    """
    Find the n-th overtone of the Poeschl-Teller potential.

    For V_0 = 2, the exact QNMs are:
        omega_k = sqrt(V_0 - (k + 1/2)^2) - i(k + 1/2),  k = 0, 1, ...
    The first overtone (k=1) is at omega_1 = -i (purely imaginary).

    Parameters
    ----------
    omega : ndarray
    phys_idx : ndarray of int
    n : int
        Overtone number (n=1 is first overtone).

    Returns
    -------
    idx : int or None
        Index into omega of the n-th overtone.
    """
    # Exact overtone for Poeschl-Teller with V0=2
    k = n
    arg = V0 - (k + 0.5) ** 2
    if arg >= 0:
        omega_ot = np.sqrt(arg) - 1j * (k + 0.5)
    else:
        omega_ot = 1j * np.sqrt(-arg) - 1j * (k + 0.5)

    dist = np.abs(omega[phys_idx] - omega_ot)
    if len(dist) == 0:
        return None
    best = phys_idx[np.argmin(dist)]
    return best


# -----------------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------------
def plot_spectrum(omega, fund_idx, output_dir):
    """Figure 1: QNM spectrum in the complex omega-plane."""

    fig, ax = plt.subplots(1, 1, figsize=(7, 5.5))

    # All eigenvalues (grey dots -- all spurious except fundamental)
    finite_mask = np.isfinite(omega)
    om_all = omega[finite_mask]
    clip = (np.abs(om_all.real) < 12) & (np.abs(om_all.imag) < 8)
    ax.plot(om_all[clip].real, om_all[clip].imag, '.',
            color='#BBBBBB', markersize=3, label='Spurious', zorder=1)

    # Fundamental mode (star)
    om0 = omega[fund_idx]
    ax.plot(om0.real, om0.imag, '*', color=CORAL, markersize=16,
            markeredgecolor='white', markeredgewidth=0.5,
            label=r'Fundamental $\omega_0$', zorder=4)

    # Exact fundamental
    ax.plot(OMEGA_EXACT.real, OMEGA_EXACT.imag, '+', color=TEAL,
            markersize=12, markeredgewidth=2.0,
            label=r'Exact $\omega_0$', zorder=5)

    # Reference line Im(omega) = 0
    ax.axhline(0, color='#888888', linewidth=0.5, linestyle='--', zorder=0)

    ax.set_xlabel(r'Re$(\omega)$')
    ax.set_ylabel(r'Im$(\omega)$')
    ax.set_title(u'QNM Spectrum: P\u00f6schl\u2013Teller Potential ($V_0 = %g$)' % V0,
                 fontsize=12)
    ax.legend(loc='lower right', fontsize=9, framealpha=0.9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(True, alpha=0.2)

    # Annotate fundamental
    ax.annotate(
        r'$\omega_0 = \frac{\sqrt{7}}{2} - \frac{i}{2}$',
        xy=(OMEGA_EXACT.real, OMEGA_EXACT.imag),
        xytext=(4.5, -0.15),
        fontsize=10, color=TEAL,
        arrowprops=dict(arrowstyle='->', color=TEAL, alpha=0.7),
        bbox=dict(boxstyle='round,pad=0.3', fc='white', ec=TEAL, alpha=0.85)
    )

    plt.tight_layout()

    output_file = output_dir / 'qnm_spectrum.pdf'
    fig.savefig(output_file, bbox_inches='tight', pad_inches=0.05)
    fig.savefig(output_file.with_suffix('.png'), bbox_inches='tight',
                pad_inches=0.05, dpi=300)
    print(f'Figure saved to: {output_file.resolve()}')
    plt.close(fig)


def plot_eigenfunctions(omega, psi_all, x_phys, fund_idx, ot1_idx, output_dir):
    """Figure 2: Eigenfunctions of fundamental and first overtone."""

    V = poschl_teller(x_phys)
    V_scaled = V / np.max(V)  # normalise to [0, 1] for background

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    titles = [r'Fundamental mode ($n = 0$)', r'First overtone ($n = 1$)']
    indices = [fund_idx, ot1_idx]

    for ax, idx, title in zip(axes, indices, titles):
        if idx is None:
            ax.set_title(title + ' (not found)')
            continue

        psi = psi_all[:, idx]
        # Normalise so max|psi| = 1
        psi = psi / np.max(np.abs(psi))

        om = omega[idx]

        # Shaded potential background
        ax.fill_between(x_phys, 0, V_scaled, color=SKY, alpha=0.15,
                         label=r'$V(x)$ (scaled)')

        # Real and imaginary parts
        ax.plot(x_phys, psi.real, '-', color=NAVY, linewidth=1.5,
                label=r'Re$(\psi)$')
        ax.plot(x_phys, psi.imag, '--', color=CORAL, linewidth=1.5,
                label=r'Im$(\psi)$')

        ax.set_xlabel(r'$x$')
        ax.set_ylabel(r'$\psi(x)$')
        ax.set_title(title + '\n' +
                     r'$\omega = %.4f %+.4fi$' % (om.real, om.imag),
                     fontsize=11)
        ax.legend(loc='upper right', fontsize=8, framealpha=0.9)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.axhline(0, color='#888888', linewidth=0.5)
        ax.grid(True, alpha=0.2)

    fig.suptitle(u'QNM Eigenfunctions: P\u00f6schl\u2013Teller Potential',
                 fontsize=13, y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.93])

    output_file = output_dir / 'qnm_eigenfunctions.pdf'
    fig.savefig(output_file, bbox_inches='tight', pad_inches=0.05)
    fig.savefig(output_file.with_suffix('.png'), bbox_inches='tight',
                pad_inches=0.05, dpi=300)
    print(f'Figure saved to: {output_file.resolve()}')
    plt.close(fig)


def plot_convergence(output_dir):
    """Figure 3: Convergence of fundamental QNM with N and L."""

    # ---- Left panel: error vs N (fixed L = 8) ----
    L_fixed = 8.0
    N_values = list(range(20, 110, 10))
    errors_N = []

    for N_val in N_values:
        omega, _, _ = solve_qnm(N_val, L_fixed)
        phys_idx, _ = classify_modes(omega)
        if len(phys_idx) == 0:
            errors_N.append(np.nan)
            continue
        fund = find_fundamental(omega, phys_idx)
        err = np.abs(omega[fund] - OMEGA_EXACT)
        errors_N.append(max(err, 1e-16))

    # ---- Right panel: error vs L (fixed N = 80) ----
    N_fixed = 80
    L_values = list(range(2, 13))
    errors_L = []

    for L_val in L_values:
        omega, _, _ = solve_qnm(N_fixed, float(L_val))
        phys_idx, _ = classify_modes(omega)
        if len(phys_idx) == 0:
            errors_L.append(np.nan)
            continue
        fund = find_fundamental(omega, phys_idx)
        err = np.abs(omega[fund] - OMEGA_EXACT)
        errors_L.append(max(err, 1e-16))

    # ---- Plot ----
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    ax1 = axes[0]
    ax1.semilogy(N_values, errors_N, 'o-', color=NAVY, linewidth=1.5,
                 markersize=6, markeredgecolor='white', markeredgewidth=0.5)
    ax1.set_xlabel(r'$N$')
    ax1.set_ylabel(r'$|\omega_0^{\mathrm{num}} - \omega_0^{\mathrm{exact}}|$')
    ax1.set_title(r'Convergence vs $N$ (fixed $L = %g$)' % L_fixed, fontsize=11)
    ax1.grid(True, which='both', alpha=0.3)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    ax2 = axes[1]
    ax2.semilogy(L_values, errors_L, 's-', color=CORAL, linewidth=1.5,
                 markersize=6, markeredgecolor='white', markeredgewidth=0.5)
    ax2.set_xlabel(r'$L$')
    ax2.set_ylabel(r'$|\omega_0^{\mathrm{num}} - \omega_0^{\mathrm{exact}}|$')
    ax2.set_title(r'Convergence vs $L$ (fixed $N = %d$)' % N_fixed, fontsize=11)
    ax2.grid(True, which='both', alpha=0.3)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    fig.suptitle(r'Convergence of Fundamental QNM $\omega_0$',
                 fontsize=13, y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.93])

    output_file = output_dir / 'qnm_convergence.pdf'
    fig.savefig(output_file, bbox_inches='tight', pad_inches=0.05)
    fig.savefig(output_file.with_suffix('.png'), bbox_inches='tight',
                pad_inches=0.05, dpi=300)
    print(f'Figure saved to: {output_file.resolve()}')
    plt.close(fig)

    return N_values, errors_N, L_values, errors_L


def main():
    """Compute QNMs and generate all three figures."""

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # ==================================================================
    # Solve the QNM problem
    # ==================================================================
    N = 80
    L = 8.0

    print("Quasinormal Modes of a Black Hole (Poeschl-Teller Potential)")
    print("=" * 60)
    print(f"  V_0 = {V0},  L = {L},  N = {N}")
    print(f"  Exact fundamental: omega_0 = {OMEGA_EXACT}")
    print()

    omega, psi_all, x_phys = solve_qnm(N, L)
    phys_idx, spur_idx = classify_modes(omega)

    print(f"  Total eigenvalues:    {len(omega)}")
    print(f"  Modes with Im(ω)<0:   {len(phys_idx)}")
    print(f"  Modes with Im(ω)≥0:   {len(spur_idx)}")
    print()

    # Fundamental mode
    fund_idx = find_fundamental(omega, phys_idx)
    om0 = omega[fund_idx]
    err0 = np.abs(om0 - OMEGA_EXACT)
    print(f"  Numerical omega_0 =  {om0:.10f}")
    print(f"  |omega_0 - exact|  = {err0:.4e}")
    print()

    # First overtone
    ot1_idx = find_overtone(omega, phys_idx, n=1)
    if ot1_idx is not None:
        om1 = omega[ot1_idx]
        print(f"  First overtone:      omega_1 = {om1:.6f}")
    print()

    # Print a few physical modes sorted by |Im(omega)|
    im_sort = np.argsort(np.abs(omega[phys_idx].imag))
    print("  Physical QNMs (sorted by |Im(omega)|):")
    print("  " + "-" * 50)
    print(f"  {'n':>4}  {'Re(omega)':>14}  {'Im(omega)':>14}")
    print("  " + "-" * 50)
    for k, idx in enumerate(phys_idx[im_sort[:10]]):
        print(f"  {k:>4d}  {omega[idx].real:>14.8f}  {omega[idx].imag:>14.8f}")
    print("  " + "-" * 50)
    print()

    # ==================================================================
    # Generate figures
    # ==================================================================
    print("Generating figures...")
    plot_spectrum(omega, fund_idx, OUTPUT_DIR)
    plot_eigenfunctions(omega, psi_all, x_phys, fund_idx, ot1_idx, OUTPUT_DIR)
    N_vals, err_N, L_vals, err_L = plot_convergence(OUTPUT_DIR)

    # Print convergence tables
    print("\nConvergence Table: omega_0 error vs N (L = 8)")
    print("-" * 30)
    print(f"{'N':>6}  {'Error':>12}")
    print("-" * 30)
    for nv, ev in zip(N_vals, err_N):
        print(f"{nv:>6d}  {ev:>12.4e}")
    print("-" * 30)

    print("\nConvergence Table: omega_0 error vs L (N = 80)")
    print("-" * 30)
    print(f"{'L':>6}  {'Error':>12}")
    print("-" * 30)
    for lv, ev in zip(L_vals, err_L):
        print(f"{lv:>6d}  {ev:>12.4e}")
    print("-" * 30)


if __name__ == '__main__':
    main()
