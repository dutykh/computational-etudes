#!/usr/bin/env python3
"""
periodic_pulse_two_grids.py

Chapter 19: Coordinate Transformations
Computational Etude 19.1: A periodic pulse on two grids.

We approximate the sharply localised periodic function

    f(y) = exp(-kappa * (1 - cos(y)))                 y in [-pi, pi],

with kappa large (here kappa = 80) so that the pulse has essentially
vanishing support over most of the period. Two collocation grids compete:

    1. a uniform Fourier grid on y in [-pi, pi];
    2. the arctan/tan mapped grid  y = arctan(L * tan(x)),  L = 0.3,
       which clusters physical points toward y = 0.

For each N in a sweep we compute the L^infty error of the trigonometric
interpolant in physical y-coordinates, sample the function on both grids,
and display coefficient decay. The etude is the chapter's opening shock:
the computational coordinate is part of the numerical method.

Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
Part of "Computational Etudes: A Spectral Approach"
https://github.com/dutykh/computational-etudes
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from map_common import CORAL, GOLD, NAVY, TEAL, output_dir_for, save_fig, setup_matplotlib

SCRIPT_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = output_dir_for(SCRIPT_DIR)


KAPPA = 80.0          # concentration parameter of the periodic pulse
L_MAP = 0.3           # arctan/tan map parameter (L < 1 clusters near y = 0)


def target(y):
    """Periodic, strongly localised test function on [-pi, pi]."""
    return np.exp(-KAPPA * (1.0 - np.cos(y)))


def uniform_fourier_grid(N):
    """Equispaced grid on [-pi, pi), N points (even)."""
    return -np.pi + 2.0 * np.pi * np.arange(N) / N


def arctan_tan_map(x, L):
    """Physical coordinate y corresponding to computational coordinate x.

    We use the 2-pi-periodic variant
        y(x) = 2 * arctan(L * tan(x / 2)),     x, y in [-pi, pi],
    which smoothly bijects the computational interval [-pi, pi] to the
    physical interval [-pi, pi].  For L < 1 it clusters grid points near
    y = 0; for L > 1 it clusters them near y = +-pi.
    """
    return 2.0 * np.arctan(L * np.tan(x / 2.0))


def arctan_tan_inverse(y, L):
    """Inverse of arctan_tan_map: x = 2 * arctan(tan(y/2) / L)."""
    return 2.0 * np.arctan(np.tan(y / 2.0) / L)


def mapped_grid(N, L):
    """Uniform x-grid on the computational period [-pi, pi) mapped to y."""
    x = -np.pi + 2.0 * np.pi * np.arange(N) / N
    return arctan_tan_map(x, L), x


def fourier_interp(y_nodes, f_nodes, y_eval):
    """Trigonometric interpolation of samples f_nodes taken on the grid
    y_nodes = -pi + 2*pi*n/N (n=0..N-1) and extended 2pi-periodically.
    """
    N = len(y_nodes)
    coeffs = np.fft.fft(f_nodes) / N        # DFT of samples
    k = np.fft.fftfreq(N, d=1.0 / N)        # wave-numbers 0, 1, ..., -1
    phase = np.exp(1j * k[:, None] * (y_eval[None, :] + np.pi))
    return (coeffs[:, None] * phase).sum(axis=0).real


def mapped_interp(x_nodes, f_nodes, y_eval, L):
    """Evaluate the mapped Fourier interpolant at arbitrary y_eval."""
    x_eval = arctan_tan_inverse(y_eval, L)
    return fourier_interp(x_nodes, f_nodes, x_eval)


def run_sweep():
    Ns = np.array([8, 12, 16, 24, 32, 48, 64, 96, 128])
    y_eval = np.linspace(-np.pi, np.pi, 4097)
    y_eval[0] += 1e-9
    y_eval[-1] -= 1e-9
    truth = target(y_eval)

    err_uniform, err_mapped = [], []
    for N in Ns:
        # --- uniform Fourier grid
        yU = uniform_fourier_grid(N)
        fU = target(yU)
        approx_U = fourier_interp(yU, fU, y_eval)
        err_uniform.append(np.max(np.abs(approx_U - truth)))

        # --- arctan/tan mapped grid (uniform in x)
        yM, xM = mapped_grid(N, L_MAP)
        fM = target(yM)
        approx_M = mapped_interp(xM, fM, y_eval, L_MAP)
        err_mapped.append(np.max(np.abs(approx_M - truth)))

    return Ns, np.array(err_uniform), np.array(err_mapped)


def coefficient_decay(N):
    """Magnitude |c_k| of Fourier coefficients on the two grids at one N."""
    yU = uniform_fourier_grid(N)
    cU = np.abs(np.fft.fftshift(np.fft.fft(target(yU)) / N))
    yM, xM = mapped_grid(N, L_MAP)
    cM = np.abs(np.fft.fftshift(np.fft.fft(target(yM)) / N))
    k = np.fft.fftshift(np.fft.fftfreq(N, d=1.0 / N))
    return k, cU, cM


def make_figure():
    setup_matplotlib()

    Ns, err_U, err_M = run_sweep()
    Nshow = 32
    yU = uniform_fourier_grid(Nshow)
    yM, _ = mapped_grid(Nshow, L_MAP)
    y_eval = np.linspace(-np.pi, np.pi, 2001)

    fig, axes = plt.subplots(1, 3, figsize=(12.4, 3.6))

    # (a) the pulse and the two grids
    ax = axes[0]
    ax.plot(y_eval, target(y_eval), color=NAVY, lw=1.2, label=r"$f(y)$")
    ax.scatter(yU, np.full_like(yU, -0.08), s=24, marker="x",
               color=CORAL, label=f"uniform, $N={Nshow}$")
    ax.scatter(yM, np.full_like(yM, -0.18), s=24, marker="o",
               color=TEAL, facecolors="none", label=f"arctan/tan, $L={L_MAP}$")
    ax.set_xlim(-np.pi, np.pi)
    ax.set_ylim(-0.28, 1.1)
    ax.set_xlabel(r"$y$")
    ax.set_ylabel(r"$f(y)$, grid points")
    ax.set_title("(a) Pulse + two grids")
    ax.legend(loc="upper right", frameon=False, fontsize=9)

    # (b) error versus N
    ax = axes[1]
    ax.semilogy(Ns, err_U, "-o", color=CORAL, lw=1.1, label="uniform Fourier")
    ax.semilogy(Ns, err_M, "-s", color=TEAL, lw=1.1,
                label=f"arctan/tan, $L={L_MAP}$")
    ax.set_xlabel(r"$N$")
    ax.set_ylabel(r"$\|f - f_N\|_\infty$")
    ax.set_title("(b) Convergence")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(frameon=False, fontsize=9)

    # (c) coefficient decay
    ax = axes[2]
    k, cU, cM = coefficient_decay(96)
    mask = (k >= 0)
    ax.semilogy(k[mask], cU[mask] + 1e-18, "x", color=CORAL, ms=4,
                label="uniform Fourier")
    ax.semilogy(k[mask], cM[mask] + 1e-18, "o", color=TEAL, ms=3,
                mfc="none", label="arctan/tan")
    ax.set_xlabel(r"wave-number $k$")
    ax.set_ylabel(r"$|c_k|$")
    ax.set_title("(c) Coefficient decay, $N=96$")
    ax.grid(True, which="both", alpha=0.3)
    ax.set_ylim(1e-16, 1)
    ax.legend(frameon=False, fontsize=9)

    fig.tight_layout()
    save_fig(fig, OUTPUT_DIR, "periodic_pulse_two_grids")
    plt.close(fig)

    print(f"[19.1] saved figure to {OUTPUT_DIR / 'periodic_pulse_two_grids.pdf'}")
    print("  N, err_uniform, err_mapped:")
    for N, eU, eM in zip(Ns, err_U, err_M):
        print(f"    N={N:4d}   uniform={eU:.3e}   mapped={eM:.3e}")


if __name__ == "__main__":
    make_figure()
