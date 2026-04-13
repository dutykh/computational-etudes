#!/usr/bin/env python3
"""
fair_comparison.py
Chapter 17: Stability of Spectral Methods

Computational Etude 17.6: Fair Comparison of Time-Stepping Methods.

Two-panel figure comparing Forward Euler, Backward Euler, Crank-Nicolson,
and Integrating Factor (exact exponential) on the 1D Fourier heat equation
u_t = nu * u_xx on [0, 2*pi) with nu = 0.1.

  (a) Error vs dt (log-log).
  (b) Error vs wall-clock time (log-log).

Initial condition: u0 = sin(x) + 0.5*cos(3x).
Exact solution: each Fourier mode k decays as exp(-nu*k^2*T).

Generates Figure 17.6: fair_comparison.pdf

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
import time

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

SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / '..' / '..' / '..' / 'textbook' / 'figures' / 'ch17' / 'python'


def exact_solution_fourier(N, nu, T):
    """Exact solution in Fourier space: each mode decays exponentially."""
    x = np.linspace(0, 2 * np.pi, N, endpoint=False)
    u0 = np.sin(x) + 0.5 * np.cos(3 * x)
    u0_hat = np.fft.fft(u0)

    k = np.fft.fftfreq(N, d=1.0 / N)
    decay = np.exp(-nu * k ** 2 * T)
    u_exact_hat = u0_hat * decay
    return np.real(np.fft.ifft(u_exact_hat))


def forward_euler(N, nu, T, dt):
    """Forward Euler in Fourier space: u_hat^{n+1} = (1 - nu*k^2*dt)*u_hat^n."""
    x = np.linspace(0, 2 * np.pi, N, endpoint=False)
    u0 = np.sin(x) + 0.5 * np.cos(3 * x)
    u_hat = np.fft.fft(u0)

    k = np.fft.fftfreq(N, d=1.0 / N)
    lam = -nu * k ** 2  # eigenvalues of nu * d^2/dx^2

    n_steps = int(round(T / dt))
    dt_actual = T / n_steps

    t0 = time.perf_counter()
    for _ in range(n_steps):
        u_hat = u_hat + dt_actual * lam * u_hat
        if np.max(np.abs(u_hat)) > 1e10:
            return None, time.perf_counter() - t0  # blow-up
    elapsed = time.perf_counter() - t0

    return np.real(np.fft.ifft(u_hat)), elapsed


def backward_euler(N, nu, T, dt):
    """Backward Euler: u_hat^{n+1} = u_hat^n / (1 - dt*lam)."""
    x = np.linspace(0, 2 * np.pi, N, endpoint=False)
    u0 = np.sin(x) + 0.5 * np.cos(3 * x)
    u_hat = np.fft.fft(u0)

    k = np.fft.fftfreq(N, d=1.0 / N)
    lam = -nu * k ** 2

    n_steps = int(round(T / dt))
    dt_actual = T / n_steps
    factor = 1.0 / (1.0 - dt_actual * lam)

    t0 = time.perf_counter()
    for _ in range(n_steps):
        u_hat = u_hat * factor
    elapsed = time.perf_counter() - t0

    return np.real(np.fft.ifft(u_hat)), elapsed


def crank_nicolson(N, nu, T, dt):
    """Crank-Nicolson: u_hat^{n+1} = (1 + dt*lam/2)/(1 - dt*lam/2) * u_hat^n."""
    x = np.linspace(0, 2 * np.pi, N, endpoint=False)
    u0 = np.sin(x) + 0.5 * np.cos(3 * x)
    u_hat = np.fft.fft(u0)

    k = np.fft.fftfreq(N, d=1.0 / N)
    lam = -nu * k ** 2

    n_steps = int(round(T / dt))
    dt_actual = T / n_steps
    factor = (1.0 + 0.5 * dt_actual * lam) / (1.0 - 0.5 * dt_actual * lam)

    t0 = time.perf_counter()
    for _ in range(n_steps):
        u_hat = u_hat * factor
    elapsed = time.perf_counter() - t0

    return np.real(np.fft.ifft(u_hat)), elapsed


def integrating_factor(N, nu, T, dt):
    """Integrating Factor (exact exponential for linear problems).
    u_hat^{n+1} = exp(lam * dt) * u_hat^n."""
    x = np.linspace(0, 2 * np.pi, N, endpoint=False)
    u0 = np.sin(x) + 0.5 * np.cos(3 * x)
    u_hat = np.fft.fft(u0)

    k = np.fft.fftfreq(N, d=1.0 / N)
    lam = -nu * k ** 2

    n_steps = int(round(T / dt))
    dt_actual = T / n_steps
    factor = np.exp(lam * dt_actual)

    t0 = time.perf_counter()
    for _ in range(n_steps):
        u_hat = u_hat * factor
    elapsed = time.perf_counter() - t0

    return np.real(np.fft.ifft(u_hat)), elapsed


def main():
    N = 64
    nu = 0.1
    T = 1.0

    u_exact = exact_solution_fourier(N, nu, T)

    # CFL for forward Euler: dt <= 2 / (nu * (N/2)^2)
    dt_cfl = 2.0 / (nu * (N / 2) ** 2)
    print(f"Forward Euler CFL limit: dt <= {dt_cfl:.6e}")

    methods = {
        'Forward Euler': (forward_euler, NAVY, 'o', '-'),
        'Backward Euler': (backward_euler, TEAL, 's', '-'),
        'Crank-Nicolson': (crank_nicolson, CORAL, '^', '-'),
        'Integrating Factor': (integrating_factor, AMBER, 'D', '-'),
    }

    # dt values: for FE we must stay below CFL, others can use larger dt
    dt_values_all = np.logspace(-5, -1, 30)
    dt_values_fe = dt_values_all[dt_values_all < 0.9 * dt_cfl]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    for name, (method, color, marker, ls) in methods.items():
        dt_vals = dt_values_fe if name == 'Forward Euler' else dt_values_all
        errors = []
        times = []
        dts_used = []

        for dt in dt_vals:
            u_num, elapsed = method(N, nu, T, dt)
            if u_num is None:
                continue
            err = np.max(np.abs(u_num - u_exact))
            if err < 1e-14:
                err = 1e-14
            errors.append(err)
            times.append(elapsed)
            dts_used.append(dt)

        if len(errors) == 0:
            continue

        errors = np.array(errors)
        times = np.array(times)
        dts_used = np.array(dts_used)

        ax1.loglog(dts_used, errors, marker=marker, color=color,
                   markersize=5, lw=1.2, linestyle=ls,
                   label=name, markevery=3)
        ax2.loglog(times, errors, marker=marker, color=color,
                   markersize=5, lw=1.2, linestyle=ls,
                   label=name, markevery=3)

    # Reference slopes
    dt_ref = np.logspace(-5, -1, 50)
    ax1.loglog(dt_ref, 0.5 * dt_ref, '--', color='gray', lw=0.8, alpha=0.7)
    ax1.text(2e-3, 2e-3, r'slope 1', fontsize=8, color='gray', rotation=35)
    ax1.loglog(dt_ref, 2.0 * dt_ref ** 2, ':', color='gray', lw=0.8,
               alpha=0.7)
    ax1.text(5e-3, 8e-5, r'slope 2', fontsize=8, color='gray', rotation=42)

    ax1.set_xlabel(r'$\Delta t$')
    ax1.set_ylabel(r'$\|u_{\rm num} - u_{\rm exact}\|_\infty$')
    ax1.set_title(r'(a) Error vs time step')
    ax1.legend(loc='upper left', fontsize=8)
    ax1.grid(True, which='both', alpha=0.3, linewidth=0.5)
    ax1.set_ylim(1e-15, 1e1)

    ax2.set_xlabel('Wall-clock time (s)')
    ax2.set_ylabel(r'$\|u_{\rm num} - u_{\rm exact}\|_\infty$')
    ax2.set_title(r'(b) Error vs computational cost')
    ax2.legend(loc='upper right', fontsize=8)
    ax2.grid(True, which='both', alpha=0.3, linewidth=0.5)
    ax2.set_ylim(1e-15, 1e1)

    plt.tight_layout()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_DIR / 'fair_comparison.pdf', bbox_inches='tight')
    fig.savefig(OUTPUT_DIR / 'fair_comparison.png', bbox_inches='tight')
    print(f"\nFigure saved to {OUTPUT_DIR / 'fair_comparison.pdf'}")

    plt.show()


if __name__ == '__main__':
    main()
