#!/usr/bin/env python3
"""
ho_kuramoto_sivashinsky.py
Chapter 14: Higher-Order Spectral Methods

Kuramoto-Sivashinsky equation (Computational Etude 14.9):

    u_t + u * u_x + u_xx + u_xxxx = 0

on [0, 32 pi] with periodic boundary conditions. Solved using Fourier
pseudospectral methods in space and the ETDRK4 scheme (Cox & Matthews
2002, with Kassam-Trefethen 2005 contour integral improvement) in time.

The KS equation models laminar flame front instabilities and thin film
flows. In Fourier space, the equation becomes:

    du_hat/dt = L * u_hat + N_hat(u)

where the linear operator L_k = k^2 - k^4 (energy production at low k,
dissipation at high k) and the nonlinear term N_hat = -(i k / 2) FFT(u^2)
uses the conservation form u * u_x = (u^2 / 2)_x.

The ETDRK4 scheme evaluates the phi-functions

    phi_j(z) = integral formulas involving exp(z)

via Kassam & Trefethen's contour integral trick on a circle in the complex
plane, avoiding catastrophic cancellation near z = 0.

The solution exhibits spatiotemporal chaos: cell merging, splitting, and
irregular oscillations characteristic of fourth-order dissipative systems.

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
# ETDRK4 coefficients via Kassam-Trefethen contour integral trick
# -----------------------------------------------------------------------------
def compute_etdrk4_coefficients(Lop, dt, M=32):
    """
    Compute ETDRK4 coefficients using the Kassam-Trefethen (2005) contour
    integral method for numerical stability.

    The standard ETDRK4 formulas involve ratios of the form
    (exp(z) - 1 - z - z^2/2) / z^3 that suffer from catastrophic
    cancellation when z = L_k * dt is near zero. The contour integral
    trick evaluates these on a circle of radius 1 in the complex plane
    around each L_k * dt, averaging over M roots of unity.

    Parameters
    ----------
    Lop : ndarray, shape (N,)
        Linear operator diagonal (real-valued): L_k = k^2 - k^4.
    dt : float
        Time step size.
    M : int, optional
        Number of contour quadrature points (default 32).

    Returns
    -------
    E : ndarray, shape (N,)
        Full-step exponential propagator exp(L * dt).
    E2 : ndarray, shape (N,)
        Half-step exponential propagator exp(L * dt / 2).
    Q : ndarray, shape (N,)
        ETDRK4 coefficient for stage evaluations.
    f1 : ndarray, shape (N,)
        ETDRK4 coefficient for final update (weight on N_v).
    f2 : ndarray, shape (N,)
        ETDRK4 coefficient for final update (weight on N_a + N_b).
    f3 : ndarray, shape (N,)
        ETDRK4 coefficient for final update (weight on N_c).
    """
    N = len(Lop)

    # Full-step and half-step propagators
    E = np.exp(dt * Lop)
    E2 = np.exp(dt * Lop / 2.0)

    # Roots of unity for the contour integral
    r = np.exp(2j * np.pi * np.arange(1, M + 1) / M)  # shape (M,)

    # LR[i, m] = dt * L_i + r_m, evaluated on a circle of radius 1
    # around each dt * L_i in the complex plane
    LR = dt * Lop[:, None] + r[None, :]  # shape (N, M)

    # Contour integral evaluations of the phi-functions:
    #
    # Q  = dt * phi_1(L*dt/2) = dt * (exp(L*dt/2) - 1) / (L*dt)
    #    ~ evaluated at half-step for the intermediate stages
    #
    # f1 = dt * (-4 - LR + exp(LR)*(4 - 3*LR + LR^2)) / LR^3
    # f2 = dt * ( 2 + LR + exp(LR)*(-2 + LR))           / LR^3
    # f3 = dt * (-4 - 3*LR - LR^2 + exp(LR)*(4 - LR))  / LR^3
    #
    # These are averaged over M points on the contour to avoid
    # the 0/0 indeterminate forms at LR = 0.

    Q = dt * np.real(np.mean((np.exp(LR / 2.0) - 1.0) / LR, axis=1))

    f1 = dt * np.real(np.mean(
        (-4.0 - LR + np.exp(LR) * (4.0 - 3.0 * LR + LR**2)) / LR**3,
        axis=1))

    f2 = dt * np.real(np.mean(
        (2.0 + LR + np.exp(LR) * (-2.0 + LR)) / LR**3,
        axis=1))

    f3 = dt * np.real(np.mean(
        (-4.0 - 3.0 * LR - LR**2 + np.exp(LR) * (4.0 - LR)) / LR**3,
        axis=1))

    return E, E2, Q, f1, f2, f3


# -----------------------------------------------------------------------------
# Kuramoto-Sivashinsky solver with ETDRK4
# -----------------------------------------------------------------------------
def solve_ks_etdrk4(N=256, L=32.0 * np.pi, t_final=150.0, dt=0.5,
                    n_snapshots=500):
    """
    Solve the Kuramoto-Sivashinsky equation using Fourier pseudospectral
    spatial discretization and ETDRK4 time integration.

    The equation u_t + u*u_x + u_xx + u_xxxx = 0 is rewritten as:

        du_hat/dt = L_k * u_hat + N_hat(u)

    where L_k = k^2 - k^4 and N_hat = -(ik/2) * FFT(u^2).

    Parameters
    ----------
    N : int
        Number of Fourier modes (grid points).
    L : float
        Domain length [0, L) with periodic BCs.
    t_final : float
        Final integration time.
    dt : float
        Time step size.
    n_snapshots : int
        Number of solution snapshots to store for plotting.

    Returns
    -------
    x : ndarray, shape (N,)
        Spatial grid.
    t_save : ndarray, shape (n_snapshots+1,)
        Times at which snapshots are stored.
    U_save : ndarray, shape (n_snapshots+1, N)
        Solution snapshots u(x, t_j).
    """
    # Spatial grid
    x = L * np.arange(N) / N

    # Wavenumbers: k = (2*pi/L) * [0, 1, 2, ..., N/2-1, -N/2, ..., -1]
    # np.fft.fftfreq(N, d=1/N) gives [0, 1, 2, ..., N/2-1, -N/2, ..., -1]
    k = (2.0 * np.pi / L) * np.fft.fftfreq(N, d=1.0 / N)

    # Linear operator: L_k = k^2 - k^4
    # From u_t = -u*u_x - u_xx - u_xxxx:
    #   -u_xx  -> -(-k^2) u_hat = k^2 u_hat
    #   -u_xxxx -> -(k^4) u_hat = -k^4 u_hat
    # So L_k = k^2 - k^4
    Lop = k**2 - k**4

    # 2/3 dealiasing mask for the quadratic nonlinearity u^2
    # Zero out modes with |k| > N/3 * (2*pi/L)
    dealias = np.ones(N, dtype=bool)
    cutoff = N // 3
    if cutoff > 0:
        dealias[cutoff + 1: N - cutoff] = False

    # Time stepping setup
    n_steps = int(np.ceil(t_final / dt))
    dt = t_final / n_steps  # adjust dt to hit t_final exactly
    save_every = max(1, n_steps // n_snapshots)

    # Compute ETDRK4 coefficients
    E, E2, Q, f1, f2, f3 = compute_etdrk4_coefficients(Lop, dt)

    # Initial condition: u(x, 0) = cos(x/16) * (1 + sin(x/16))
    # This is the initial condition used by Kassam & Trefethen (2005)
    u = np.cos(x / 16.0) * (1.0 + np.sin(x / 16.0))
    v = np.fft.fft(u)

    # Nonlinear term in Fourier space: N_hat = -(ik/2) * FFT(u^2)
    # Uses the conservation form u*u_x = (u^2/2)_x
    def Nhat(v_hat):
        """Evaluate dealiased nonlinear term in Fourier space."""
        u_phys = np.fft.ifft(v_hat).real
        nl = -0.5j * k * np.fft.fft(u_phys**2)
        nl[~dealias] = 0.0  # dealiasing
        return nl

    # Storage for snapshots
    snapshots = [u.copy()]
    snap_times = [0.0]

    # ETDRK4 time integration loop
    # Scheme (Cox & Matthews 2002, Kassam & Trefethen 2005):
    #
    #   Nv = N(v_n)
    #   a  = E2 * v_n + Q * Nv
    #   Na = N(a)
    #   b  = E2 * v_n + Q * Na
    #   Nb = N(b)
    #   c  = E2 * a + Q * (2*Nb - Nv)
    #   Nc = N(c)
    #   v_{n+1} = E * v_n + f1*Nv + 2*f2*(Na + Nb) + f3*Nc

    t = 0.0
    for n in range(n_steps):
        Nv = Nhat(v)
        a = E2 * v + Q * Nv
        Na = Nhat(a)
        b = E2 * v + Q * Na
        Nb = Nhat(b)
        c = E2 * a + Q * (2.0 * Nb - Nv)
        Nc = Nhat(c)
        v = E * v + f1 * Nv + 2.0 * f2 * (Na + Nb) + f3 * Nc

        t += dt

        if (n + 1) % save_every == 0:
            u = np.fft.ifft(v).real
            snapshots.append(u.copy())
            snap_times.append(t)

    return x, np.array(snap_times), np.array(snapshots)


# -----------------------------------------------------------------------------
# Main routine
# -----------------------------------------------------------------------------
def main():
    """Solve the KS equation with ETDRK4 and generate space-time plot."""

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # =========================================================================
    # Solve the Kuramoto-Sivashinsky equation
    # =========================================================================
    N = 256                  # number of Fourier modes
    L = 32.0 * np.pi        # domain length
    t_final = 150.0          # final time
    dt = 0.5                 # time step (ETDRK4 allows larger steps)

    print("Kuramoto-Sivashinsky equation with ETDRK4")
    print(f"  N = {N}, L = {L/np.pi:.0f}*pi, T = {t_final}, dt = {dt}")
    print("  Initial condition: u(x,0) = cos(x/16) * (1 + sin(x/16))")
    print("  Computing...")

    x, t_save, U_save = solve_ks_etdrk4(
        N=N, L=L, t_final=t_final, dt=dt, n_snapshots=500)

    print(f"  {len(t_save)} snapshots saved, t in [{t_save[0]:.1f}, {t_save[-1]:.1f}]")
    print(f"  max |u| = {np.max(np.abs(U_save)):.4f}")

    # =========================================================================
    # Figure: Space-time contour plot
    # =========================================================================
    fig, ax = plt.subplots(figsize=(8, 6))

    X, T = np.meshgrid(x, t_save)

    # Use a symmetric color range based on the 98th percentile
    # (after the initial transient has passed)
    n_transient = len(t_save) // 5
    vmax = np.percentile(np.abs(U_save[n_transient:]), 98)

    pcm = ax.pcolormesh(X / np.pi, T, U_save,
                        cmap='RdBu_r', shading='gouraud',
                        vmin=-vmax, vmax=vmax, rasterized=True)

    ax.set_xlabel(r'$x / \pi$')
    ax.set_ylabel(r'$t$')
    ax.set_title(
        r'Kuramoto--Sivashinsky Equation: ETDRK4 '
        r'($L = 32\pi$, $N = 256$)',
        fontsize=11)

    cbar = fig.colorbar(pcm, ax=ax, shrink=0.85, pad=0.02)
    cbar.set_label(r'$u(x, t)$', fontsize=10)
    cbar.ax.tick_params(labelsize=9)

    # Axis formatting
    ax.set_xlim(0, L / np.pi)
    ax.set_ylim(0, t_final)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    fig.suptitle(
        r'Computational Etude 14.9: Kuramoto--Sivashinsky Equation',
        fontsize=13, y=0.98)

    plt.tight_layout(rect=[0, 0, 1, 0.94])

    # ---- Save ----
    output_file = OUTPUT_DIR / 'kuramoto_sivashinsky.pdf'
    fig.savefig(output_file, bbox_inches='tight', pad_inches=0.05)
    fig.savefig(output_file.with_suffix('.png'), bbox_inches='tight',
                pad_inches=0.05, dpi=300)

    print(f"\nFigure saved to: {output_file.resolve()}")
    plt.close(fig)


if __name__ == '__main__':
    main()
