#!/usr/bin/env python3
"""
fourier_cfl_variable.py
Chapter 17: Stability of Spectral Methods

Computational Etude 17.3 (variable-coefficient extension):
Frozen-coefficient CFL prediction for variable-speed advection.

Two-panel figure:
  (a) Speed profile c(x) = 1 + 0.5*sin(x) on [0, 2*pi) with c_max
      and c_mean annotated.
  (b) Log-log plot of critical dt vs N for RK4 applied to Fourier
      advection u_t + c(x)*u_x = 0.  Experimental binary search
      compared with frozen-coefficient prediction 2.83/(c_max*N/2).

Generates Figure: fourier_cfl_variable.pdf

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

SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / '..' / '..' / '..' / 'textbook' / 'figures' / 'ch17' / 'python'


def variable_rhs(u, c, k):
    """RHS of u_t + c(x)*u_x = 0: returns -c(x) * F^{-1}(ik * F(u))."""
    u_hat = np.fft.fft(u)
    du = np.real(np.fft.ifft(1j * k * u_hat))
    return -c * du


def is_stable(dt, u0, c, k, n_steps=500):
    """Run n_steps of RK4; return True if solution stays bounded."""
    u = u0.copy()
    for _ in range(n_steps):
        k1 = variable_rhs(u, c, k)
        k2 = variable_rhs(u + 0.5 * dt * k1, c, k)
        k3 = variable_rhs(u + 0.5 * dt * k2, c, k)
        k4 = variable_rhs(u + dt * k3, c, k)
        u = u + dt / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        if np.max(np.abs(u)) > 1e6:
            return False
    return True


def find_dt_crit(N, c_func, tol=1e-6):
    """Binary search for the critical dt of variable-coefficient advection."""
    x = 2.0 * np.pi * np.arange(N) / N
    c = c_func(x)
    c_max = np.max(np.abs(c))
    u0 = np.exp(np.cos(x))  # smooth, excites all modes
    k = np.fft.fftfreq(N, d=1.0 / N)

    # Initial bracket
    dt_lo = 0.1 / (c_max * N)
    dt_hi = 10.0 / N

    # Ensure dt_hi is unstable
    while is_stable(dt_hi, u0, c, k):
        dt_hi *= 2.0
        if dt_hi > 100:
            return dt_hi

    # Binary search
    for _ in range(60):
        dt_mid = 0.5 * (dt_lo + dt_hi)
        if is_stable(dt_mid, u0, c, k):
            dt_lo = dt_mid
        else:
            dt_hi = dt_mid
        if dt_hi - dt_lo < tol * dt_lo:
            break

    return 0.5 * (dt_lo + dt_hi)


def main():
    # Speed profile
    c_func = lambda x: 1.0 + 0.5 * np.sin(x)
    c_max = 1.5
    c_mean = 1.0

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

    # --- Panel (a): speed profile ---
    x = np.linspace(0, 2.0 * np.pi, 500)
    c_vals = c_func(x)

    ax1.plot(x, c_vals, color=NAVY, lw=2.0, label=r'$c(x) = 1 + 0.5\,\sin(x)$')
    ax1.axhline(c_max, color=CORAL, ls='--', lw=1.2,
                label=r'$c_{\max} = 1.5$')
    ax1.axhline(c_mean, color=SKY, ls=':', lw=1.2,
                label=r'$\langle c \rangle = 1.0$')
    ax1.set_xlabel(r'$x$')
    ax1.set_ylabel(r'$c(x)$')
    ax1.set_title(r'(a) Variable speed profile')
    ax1.set_xlim(0, 2.0 * np.pi)
    ax1.set_ylim(0.3, 1.7)
    ax1.legend(loc='lower right', fontsize=9)
    ax1.grid(True, alpha=0.3, linewidth=0.5)

    # --- Panel (b): dt_crit vs N ---
    N_values = [16, 32, 64, 128, 256]
    dt_crits = []

    for Nv in N_values:
        dtc = find_dt_crit(Nv, c_func)
        dt_crits.append(dtc)
        print(f"N = {Nv:4d}, dt_crit = {dtc:.6f}, "
              f"dt_crit * c_max * (N/2) = {dtc * c_max * Nv / 2:.4f}")

    dt_crits = np.array(dt_crits)
    N_arr = np.array(N_values, dtype=float)

    ax2.loglog(N_arr, dt_crits, 'o-', color=NAVY, markersize=7, lw=1.5,
               label=r'$\Delta t_{\rm crit}$ (experiment)')

    # Frozen-coefficient prediction: 2.83 / (c_max * N/2)
    N_ref = np.linspace(12, 300, 100)
    ax2.loglog(N_ref, 2.83 / (c_max * N_ref / 2), '--', color=CORAL, lw=1.2,
               label=r'$2.83\, /\, (c_{\max} \cdot N/2)$')

    # Constant-coefficient reference: 2.83 / (N/2)
    ax2.loglog(N_ref, 2.83 / (N_ref / 2), ':', color=SKY, lw=1.2,
               label=r'$2.83\, /\, (N/2)$  [const. coeff.]')

    ax2.set_xlabel(r'$N$')
    ax2.set_ylabel(r'$\Delta t_{\rm crit}$')
    ax2.set_title(r'(b) Frozen-coefficient prediction vs experiment')
    ax2.legend(loc='upper right', fontsize=9)
    ax2.grid(True, which='both', alpha=0.3, linewidth=0.5)

    # Slope annotation
    ax2.text(80, 0.04, r'slope $= -1$', fontsize=11, color=CORAL,
             rotation=-28, ha='center')

    plt.tight_layout()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_DIR / 'fourier_cfl_variable.pdf', bbox_inches='tight')
    fig.savefig(OUTPUT_DIR / 'fourier_cfl_variable.png', bbox_inches='tight')
    print(f"\nFigure saved to {OUTPUT_DIR / 'fourier_cfl_variable.pdf'}")

    plt.show()


if __name__ == '__main__':
    main()
