#!/usr/bin/env python3
"""
kdv_rk_comparison.py
Chapter 17: Stability of Spectral Methods

Computational Etude 17.7: Plain RK4 vs IF-RK4 on the KdV Equation.

Two-panel figure:
  (a) Maximum stable dt vs N (log-log) for plain RK4 and IF-RK4 on
      KdV u_t + 6*u*u_x + u_xxx = 0.  For plain RK4, dt_crit ~ N^{-3}.
      For IF-RK4, dt_crit is much larger (limited by nonlinear term only).
  (b) Solutions at t = 0.5 for N = 128 with both methods (should agree).

Initial condition: u(x,0) = cos(pi*x) on [0, 2) (Zabusky-Kruskal style).

For plain RK4, the linear stability limit is set by the dispersive term
u_xxx whose eigenvalues are i*k^3.  The RK4 stability boundary on the
imaginary axis is |y| <= 2*sqrt(2), giving dt_crit = 2*sqrt(2) / k_max^3.
Since k_max = pi*N/L = pi*N/2, we get dt_crit ~ N^{-3}.

For IF-RK4, the linear part is handled exactly by the integrating factor,
so the stability limit is set by the nonlinear term only (much milder).

Generates Figure 17.7: kdv_rk_comparison.pdf

Author: Dr. Denys Dutykh
        Mathematics Department
        Khalifa University of Science and Technology
        Abu Dhabi, UAE

Part of the book "Computational Etudes: A Spectral Approach"
https://github.com/dutykh/computational-etudes
"""

import numpy as np
from numpy.fft import fft, ifft, fftfreq
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

SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / '..' / '..' / '..' / 'textbook' / 'figures' / 'ch17' / 'python'


def rk4_amplification(z):
    """RK4 amplification factor."""
    return 1.0 + z + z ** 2 / 2.0 + z ** 3 / 6.0 + z ** 4 / 24.0


def rk4_dt_crit_linear(N):
    """Analytical critical dt for RK4 on the linear part of KdV (u_xxx).
    Eigenvalues are lambda_k = -i*k^3, so z = -i*k^3*dt lies on the
    imaginary axis.  RK4 stability on imaginary axis: |y| <= 2*sqrt(2).
    Therefore dt_crit = 2*sqrt(2) / max|k^3| = 2*sqrt(2) / k_max^3."""
    L = 2.0
    k = 2.0 * np.pi * fftfreq(N, d=L / N)
    k_max = np.max(np.abs(k))
    return 2.0 * np.sqrt(2.0) / k_max ** 3


def setup_kdv(N):
    """Set up wavenumbers and dealiasing for KdV on [0, 2)."""
    L = 2.0
    x = np.linspace(0, L, N, endpoint=False)
    k = 2.0 * np.pi * fftfreq(N, d=L / N)
    dealias = np.ones(N)
    kmax = N // 3
    dealias[kmax: N - kmax + 1] = 0.0
    return x, k, dealias


def kdv_ifrk4(N, dt, T):
    """IF-RK4 for KdV: integrating factor removes the linear u_xxx term."""
    x, k, dealias = setup_kdv(N)
    u_hat = fft(np.cos(np.pi * x))

    ik = 1j * k
    lam = -1j * k ** 3

    n_steps = max(1, int(round(T / dt)))
    dt_actual = T / n_steps

    E = np.exp(lam * dt_actual / 2.0)
    E2 = E ** 2

    def NL(v_hat):
        v = np.real(ifft(v_hat))
        vx = np.real(ifft(ik * v_hat))
        return fft(-6.0 * v * vx) * dealias

    for _ in range(n_steps):
        Na = NL(u_hat)
        a = E * u_hat + 0.5 * dt_actual * E * Na
        Nb = NL(a)
        b = E * u_hat + 0.5 * dt_actual * E * Nb
        Nc = NL(b)
        c = E2 * u_hat + dt_actual * E * Nc
        Nd = NL(c)
        u_hat = (E2 * u_hat
                 + dt_actual / 6.0 * (E2 * Na + 2 * E * (Nb + Nc) + Nd))

        if np.max(np.abs(u_hat)) > 1e8:
            return None, x

    u = np.real(ifft(u_hat))
    if np.max(np.abs(u)) > 100:
        return None, x
    return u, x


def kdv_rk4(N, dt, T):
    """Plain RK4 for KdV in Fourier space."""
    x, k, dealias = setup_kdv(N)
    u_hat = fft(np.cos(np.pi * x))

    ik = 1j * k
    ik3 = 1j * k ** 3

    def rhs(v_hat):
        v = np.real(ifft(v_hat))
        vx = np.real(ifft(ik * v_hat))
        return fft(-6.0 * v * vx) * dealias - ik3 * v_hat

    n_steps = max(1, int(round(T / dt)))
    dt_actual = T / n_steps

    for _ in range(n_steps):
        k1 = rhs(u_hat)
        k2 = rhs(u_hat + 0.5 * dt_actual * k1)
        k3 = rhs(u_hat + 0.5 * dt_actual * k2)
        k4 = rhs(u_hat + dt_actual * k3)
        u_hat = u_hat + dt_actual / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)

        if np.max(np.abs(u_hat)) > 1e8:
            return None, x

    u = np.real(ifft(u_hat))
    if np.max(np.abs(u)) > 100:
        return None, x
    return u, x


def find_ifrk4_dt_crit(N, T=0.05, n_bisect=40):
    """Binary search for IF-RK4 critical dt (limited by nonlinear term)."""
    dt_lo = 1e-6
    dt_hi = 0.1

    # Find unstable upper bound
    for _ in range(20):
        u, _ = kdv_ifrk4(N, dt_hi, T)
        if u is None:
            break
        dt_hi *= 2
        if dt_hi > 10:
            return dt_hi

    # Verify dt_lo is stable
    u, _ = kdv_ifrk4(N, dt_lo, T)
    if u is None:
        dt_lo = 1e-8

    for _ in range(n_bisect):
        dt_mid = np.sqrt(dt_lo * dt_hi)
        u, _ = kdv_ifrk4(N, dt_mid, T)
        if u is not None:
            dt_lo = dt_mid
        else:
            dt_hi = dt_mid
        if dt_hi / dt_lo < 1.02:
            break

    return np.sqrt(dt_lo * dt_hi)


def main():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # --- Panel (a): dt_crit vs N ---
    N_values = [32, 64, 128, 256]
    dt_rk4 = []
    dt_ifrk4 = []

    for Nv in N_values:
        # RK4: analytical linear stability limit
        dtc_rk4 = rk4_dt_crit_linear(Nv)
        print(f"N = {Nv:4d}: RK4 dt_crit = {dtc_rk4:.6e}")
        dt_rk4.append(dtc_rk4)

        # IF-RK4: binary search (fast because dt is large)
        dtc_if = find_ifrk4_dt_crit(Nv, T=0.05)
        print(f"         IF-RK4 dt_crit = {dtc_if:.6e}")
        dt_ifrk4.append(dtc_if)

    dt_rk4 = np.array(dt_rk4)
    dt_ifrk4 = np.array(dt_ifrk4)
    N_arr = np.array(N_values, dtype=float)

    ax1.loglog(N_arr, dt_rk4, 'o-', color=NAVY, markersize=7, lw=1.5,
               label=r'Plain RK4')
    ax1.loglog(N_arr, dt_ifrk4, 's-', color=CORAL, markersize=7, lw=1.5,
               label=r'IF-RK4')

    # Reference slopes
    N_ref = np.linspace(20, 300, 100)

    slope_rk4 = np.polyfit(np.log(N_arr), np.log(dt_rk4), 1)[0]
    c_rk4 = np.exp(np.polyfit(np.log(N_arr), np.log(dt_rk4), 1)[1])
    print(f"\nRK4 slope: {slope_rk4:.2f}")
    ax1.loglog(N_ref, c_rk4 * N_ref ** slope_rk4, '--', color=SKY,
               lw=1.0, alpha=0.7,
               label=fr'$\sim N^{{{slope_rk4:.1f}}}$')

    slope_if = np.polyfit(np.log(N_arr), np.log(dt_ifrk4), 1)[0]
    c_if = np.exp(np.polyfit(np.log(N_arr), np.log(dt_ifrk4), 1)[1])
    print(f"IF-RK4 slope: {slope_if:.2f}")
    ax1.loglog(N_ref, c_if * N_ref ** slope_if, ':', color=AMBER,
               lw=1.0, alpha=0.7,
               label=fr'$\sim N^{{{slope_if:.1f}}}$')

    ax1.set_xlabel(r'$N$')
    ax1.set_ylabel(r'$\Delta t_{\rm crit}$')
    ax1.set_title(r'(a) Maximum stable $\Delta t$ vs $N$')
    ax1.legend(loc='center right', fontsize=8)
    ax1.grid(True, which='both', alpha=0.3, linewidth=0.5)

    # --- Panel (b): Solutions at t=0.5 for N=128 ---
    N_sol = 128
    T_sol = 0.5

    dt_rk4_safe = 0.8 * rk4_dt_crit_linear(N_sol)
    print(f"\nComputing RK4 solution at t={T_sol} with N={N_sol}, "
          f"dt={dt_rk4_safe:.6e}...")
    u_rk4, x = kdv_rk4(N_sol, dt_rk4_safe, T_sol)

    # Use a moderate dt for IF-RK4: much larger than the RK4 CFL limit
    # but small enough for the nonlinear term to be resolved accurately
    dt_if_use = 1e-4
    print(f"Computing IF-RK4 solution at t={T_sol} with N={N_sol}, "
          f"dt={dt_if_use:.6e}...")
    u_if, x = kdv_ifrk4(N_sol, dt_if_use, T_sol)

    if u_rk4 is not None:
        ax2.plot(x, u_rk4, '-', color=NAVY, lw=1.8, label='RK4')
    else:
        print("  WARNING: RK4 solution blew up!")
    if u_if is not None:
        ax2.plot(x, u_if, '--', color=CORAL, lw=1.5, label='IF-RK4')
    else:
        print("  WARNING: IF-RK4 solution blew up!")

    ax2.set_xlabel(r'$x$')
    ax2.set_ylabel(r'$u(x, t)$')
    ax2.set_title(fr'(b) KdV solution at $t = {T_sol}$, $N = {N_sol}$')
    ax2.legend(loc='upper right', fontsize=9)
    ax2.grid(True, alpha=0.3, linewidth=0.5)

    plt.tight_layout()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_DIR / 'kdv_rk_comparison.pdf', bbox_inches='tight')
    fig.savefig(OUTPUT_DIR / 'kdv_rk_comparison.png', bbox_inches='tight')
    print(f"\nFigure saved to {OUTPUT_DIR / 'kdv_rk_comparison.pdf'}")

    plt.show()


if __name__ == '__main__':
    main()
