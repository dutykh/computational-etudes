#!/usr/bin/env python3
"""
navier_stokes_2d.py
Chapter 11: Fourier Pseudospectral Methods for Periodic PDEs

Fourier pseudospectral solver for the 2D incompressible Navier-Stokes
equations in vorticity-streamfunction form on [0, 2*pi)^2, periodic:

    omega_t + J(psi, omega) = nu * Delta omega
    Delta psi = -omega

where J(psi, omega) = psi_x omega_y - psi_y omega_x is the Jacobian.

Uses RK4 with 2/3 dealiasing. The Poisson equation for the streamfunction
is solved trivially in Fourier space: psi_hat = omega_hat / (kx^2 + ky^2).

Generates two figures:
  - ns2d_vorticity.pdf: vorticity snapshots at several times
  - ns2d_energy.pdf: energy and enstrophy decay curves

Author: Dr. Denys Dutykh
Date: February 2026
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

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
NAVY = '#142D6E'
SKY = '#7896D2'
CORAL = '#E74C3C'
TEAL = '#16A085'

SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / '..' / '..' / '..' / 'textbook' / 'figures' / 'ch11' / 'python'


def ns2d_fourier(N=128, nu=1e-3, t_final=20.0, dt=0.01):
    """
    Solve 2D Navier-Stokes in vorticity-streamfunction form:
        omega_t + J(psi, omega) = nu * laplacian(omega)
        laplacian(psi) = -omega
    on [0, 2*pi)^2 with periodic BCs and RK4.

    Parameters:
        N : int - grid points per direction
        nu : float - kinematic viscosity
        t_final : float - final time
        dt : float - time step

    Returns:
        x, y : 1D arrays of grid points
        t_save : snapshot times
        omega_save : list of 2D vorticity fields
        energy_save, enstrophy_save : diagnostic arrays
    """
    L = 2.0 * np.pi

    # 1D grid and wavenumbers
    x = L * np.arange(N) / N
    y = x.copy()
    kx_1d = np.fft.fftfreq(N, d=1.0 / N).astype(float)
    ky_1d = kx_1d.copy()

    # 2D wavenumber grids
    KX, KY = np.meshgrid(kx_1d, ky_1d, indexing='ij')
    K2 = KX**2 + KY**2
    K2[0, 0] = 1.0  # Avoid division by zero for mean mode

    # 2/3 dealiasing mask (2D)
    cutoff = N // 3
    mask = np.ones((N, N), dtype=bool)
    if cutoff > 0:
        mask_1d = np.ones(N, dtype=bool)
        mask_1d[cutoff + 1:N - cutoff] = False
        MX, MY = np.meshgrid(mask_1d, mask_1d, indexing='ij')
        mask = MX & MY

    # Initial condition: random vorticity with zero mean
    np.random.seed(42)
    omega_hat = np.zeros((N, N), dtype=complex)
    # Populate modes with random phases and amplitude ~ k^(-1)
    for i in range(N):
        for j in range(N):
            kk = np.sqrt(KX[i, j]**2 + KY[i, j]**2)
            if 1.0 < kk < 10.0:
                amp = 1.0 / (1.0 + kk)
                phase = 2.0 * np.pi * np.random.rand()
                omega_hat[i, j] = amp * np.exp(1j * phase)
    # Enforce Hermitian symmetry for real-valued omega
    omega = np.fft.ifft2(omega_hat).real
    omega -= np.mean(omega)  # Ensure zero mean
    omega /= np.std(omega)   # Normalize to unit RMS vorticity

    # Diffusive symbol
    diff_symbol = -nu * K2

    def rhs(omega_phys):
        """RHS: -J(psi, omega) + nu * laplacian(omega)."""
        omega_hat = np.fft.fft2(omega_phys)

        # Poisson solve: psi_hat = omega_hat / K2
        psi_hat = omega_hat / K2
        psi_hat[0, 0] = 0.0  # Zero mean streamfunction

        # Velocity: u = psi_y, v = -psi_x
        u = np.fft.ifft2(1j * KY * psi_hat).real
        v = np.fft.ifft2(-1j * KX * psi_hat).real

        # Vorticity gradients
        omega_x = np.fft.ifft2(1j * KX * omega_hat).real
        omega_y = np.fft.ifft2(1j * KY * omega_hat).real

        # Jacobian (advection) in physical space
        jac = u * omega_x + v * omega_y

        # Transform, dealiase, add diffusion
        jac_hat = np.fft.fft2(jac)
        jac_hat[~mask] = 0.0

        f_hat = -jac_hat + diff_symbol * omega_hat
        return np.fft.ifft2(f_hat).real

    # Time stepping
    n_steps = int(np.ceil(t_final / dt))
    dt = t_final / n_steps

    n_save = 200
    save_every = max(1, n_steps // n_save)
    omega_save = [omega.copy()]
    t_save = [0.0]
    energy_save = []
    enstrophy_save = []

    def compute_diagnostics(omega_phys):
        omega_hat = np.fft.fft2(omega_phys)
        psi_hat = omega_hat / K2
        psi_hat[0, 0] = 0.0
        u_hat = 1j * KY * psi_hat
        v_hat = -1j * KX * psi_hat
        u = np.fft.ifft2(u_hat).real
        v = np.fft.ifft2(v_hat).real
        energy = 0.5 * np.mean(u**2 + v**2)
        enstrophy = 0.5 * np.mean(omega_phys**2)
        return energy, enstrophy

    E, Z = compute_diagnostics(omega)
    energy_save.append(E)
    enstrophy_save.append(Z)

    t = 0.0
    for n in range(n_steps):
        k1 = rhs(omega)
        k2 = rhs(omega + 0.5 * dt * k1)
        k3 = rhs(omega + 0.5 * dt * k2)
        k4 = rhs(omega + dt * k3)
        omega = omega + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        t += dt

        if (n + 1) % save_every == 0:
            omega_save.append(omega.copy())
            t_save.append(t)
            E, Z = compute_diagnostics(omega)
            energy_save.append(E)
            enstrophy_save.append(Z)

    return (x, y, np.array(t_save), omega_save,
            np.array(energy_save), np.array(enstrophy_save))


def main():
    N = 128
    nu = 1e-3
    t_final = 20.0

    x, y, t_save, omega_save, energy, enstrophy = ns2d_fourier(
        N=N, nu=nu, t_final=t_final, dt=0.005
    )

    X, Y = np.meshgrid(x, y, indexing='ij')

    # -------------------------------------------------------------------------
    # Figure 1: Vorticity snapshots
    # -------------------------------------------------------------------------
    fig1, axes = plt.subplots(2, 2, figsize=(10, 10))
    t_targets = [0.0, 5.0, 10.0, t_final]
    titles = [f'$t = {t:.0f}$' for t in t_targets]

    vmax = np.percentile(np.abs(omega_save[0]), 95)

    for ax, t_target, title in zip(axes.flat, t_targets, titles):
        idx = np.argmin(np.abs(t_save - t_target))
        omega_snap = omega_save[idx]
        pcm = ax.pcolormesh(X, Y, omega_snap, cmap='RdBu_r',
                            shading='gouraud', vmin=-vmax, vmax=vmax)
        ax.set_title(title, fontsize=12)
        ax.set_xlabel(r'$x$')
        ax.set_ylabel(r'$y$')
        ax.set_aspect('equal')
        fig1.colorbar(pcm, ax=ax, shrink=0.8)

    fig1.suptitle(
        f'2D Navier--Stokes: Decaying Turbulence ($N = {N}$, $\\nu = {nu}$)',
        fontsize=13, y=1.01
    )
    fig1.tight_layout()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig1.savefig(OUTPUT_DIR / 'ns2d_vorticity.pdf',
                 bbox_inches='tight', pad_inches=0.05)
    print(f'Figure saved: {(OUTPUT_DIR / "ns2d_vorticity.pdf").resolve()}')
    plt.close(fig1)

    # -------------------------------------------------------------------------
    # Figure 2: Energy and enstrophy decay
    # -------------------------------------------------------------------------
    fig2, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))

    ax1.plot(t_save[:len(energy)], energy, '-', color=NAVY, linewidth=1.5)
    ax1.set_xlabel(r'$t$')
    ax1.set_ylabel(r'$E(t) = \frac{1}{2}\langle |\mathbf{u}|^2 \rangle$')
    ax1.set_title('Kinetic Energy Decay')
    ax1.grid(True, alpha=0.3)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    ax2.plot(t_save[:len(enstrophy)], enstrophy, '-', color=CORAL, linewidth=1.5)
    ax2.set_xlabel(r'$t$')
    ax2.set_ylabel(r'$Z(t) = \frac{1}{2}\langle \omega^2 \rangle$')
    ax2.set_title('Enstrophy Decay')
    ax2.grid(True, alpha=0.3)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    fig2.tight_layout()
    fig2.savefig(OUTPUT_DIR / 'ns2d_energy.pdf',
                 bbox_inches='tight', pad_inches=0.05)
    print(f'Figure saved: {(OUTPUT_DIR / "ns2d_energy.pdf").resolve()}')
    plt.close(fig2)


if __name__ == '__main__':
    main()
