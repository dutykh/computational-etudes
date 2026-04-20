#!/usr/bin/env python3
"""
hermite_width_mismatch.py

Chapter 20: Spectral Methods on Unbounded Intervals.
Computational Etude 20.4: The cost of width mismatch for Hermite functions.

The (normalised) Hermite functions

    psi_n(y) = (pi^0.5  2^n  n!)^(-0.5)  H_n(y)  exp(-y^2/2),

are the eigenfunctions of the harmonic oscillator and hence the
natural basis for Gaussian-like decay.  Their width is fixed at
order 1.  For a target  f(y) = exp(-A * y^2)  of width 1/sqrt(A),
the correct Hermite basis is  psi_n(alpha * y)  with scale alpha
chosen to match f's width; Hermite functions with alpha = 1 will
perform badly if A >> 1 (f too narrow) or A << 1 (f too wide).

We reproduce Boyd (2000) Rule-of-Thumb 14 concretely: match the
built-in length scale of the basis to that of the target.  We sweep
A over three orders of magnitude and compare:

  (A) unscaled Hermite series  sum_n a_n psi_n(y);
  (B) scaled Hermite series   sum_n a_n psi_n(alpha * y)   with
      alpha = sqrt(2 A).

The etude ends with a "quantum harmonic oscillator" demonstration:
for f_true(y) = psi_0(y) = pi^(-1/4) exp(-y^2/2), the Hermite expansion
is trivially exact at N = 0.

Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from numpy.polynomial.hermite_e import hermegauss  # probabilists' version
from numpy.polynomial.hermite import hermgauss      # physicists' version

from unbounded_common import (CORAL, NAVY, TEAL, PURPLE, output_dir_for,
                              save_fig, setup_matplotlib)

SCRIPT_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = output_dir_for(SCRIPT_DIR)


def hermite_psi(n, y):
    """Normalised physicists' Hermite function psi_n(y) = H_n(y) exp(-y^2/2)
    divided by the L^2 normalisation sqrt(2^n n! sqrt(pi)).  Stable
    three-term recurrence in the NORMALISED variant.
    """
    y = np.atleast_1d(y).astype(float)
    psi = np.zeros((n + 1, y.size))
    psi[0] = np.pi ** -0.25 * np.exp(-0.5 * y ** 2)
    if n >= 1:
        psi[1] = np.sqrt(2.0) * y * psi[0]
    for k in range(1, n):
        psi[k + 1] = np.sqrt(2.0 / (k + 1)) * y * psi[k] \
                     - np.sqrt(k / (k + 1)) * psi[k - 1]
    return psi


def hermite_expand(f, N, alpha):
    """Expansion coefficients a_n of  f(y) = sum_n a_n phi_n(y),
    with phi_n(y) = sqrt(alpha) psi_n(alpha y) orthonormal on L^2(R).

    Derivation:  a_n = int f(y) phi_n(y) dy = (1/sqrt(alpha))
    int f(x/alpha) psi_n(x) dx, computed by physicists' Gauss-Hermite
    using int g(x) dx = sum_i w_i exp(x_i^2) g(x_i).
    """
    x, w = hermgauss(N + 32)
    y_nodes = x / alpha
    fv = f(y_nodes)
    psi = hermite_psi(N, x)                           # psi_n(x_i)
    integrand = fv * psi * np.exp(x ** 2)             # shape (N+1, K)
    return (integrand * w).sum(axis=1) / np.sqrt(alpha)


def hermite_eval(coeffs, y, alpha):
    """Evaluate  sum_n a_n sqrt(alpha) psi_n(alpha * y)."""
    psi = hermite_psi(len(coeffs) - 1, alpha * y)
    return np.sqrt(alpha) * (coeffs[:, None] * psi).sum(axis=0)


def max_error(f_exact, coeffs, alpha):
    y = np.linspace(-15.0, 15.0, 4001)
    return np.max(np.abs(hermite_eval(coeffs, y, alpha) - f_exact(y)))


def make_figure():
    setup_matplotlib()

    Ns = np.array([8, 12, 16, 24, 32, 48, 64])
    A_list = [0.1, 0.5, 2.0, 8.0]

    err_alpha1 = {A: [] for A in A_list}
    err_alpha_matched = {A: [] for A in A_list}
    for A in A_list:
        f_exact = lambda y, A=A: np.exp(-A * y ** 2)
        alpha_matched = np.sqrt(2.0 * A)
        for N in Ns:
            c1 = hermite_expand(f_exact, N, 1.0)
            err_alpha1[A].append(max_error(f_exact, c1, 1.0))
            c2 = hermite_expand(f_exact, N, alpha_matched)
            err_alpha_matched[A].append(max_error(f_exact, c2, alpha_matched))

    # quantum oscillator: ground state psi_0
    def f_qho(y):
        return np.pi ** -0.25 * np.exp(-0.5 * y ** 2)
    Ns_q = np.array([1, 2, 4, 8, 16, 32])
    err_qho = []
    for N in Ns_q:
        c = hermite_expand(f_qho, N, 1.0)
        err_qho.append(max_error(f_qho, c, 1.0))

    fig, axes = plt.subplots(1, 3, figsize=(13.4, 3.8))

    ax = axes[0]
    colours = [CORAL, NAVY, TEAL, PURPLE]
    for A, c in zip(A_list, colours):
        ax.semilogy(Ns, err_alpha1[A], "-o", color=c, lw=1.1,
                    label=f"$A={A}$")
    ax.set_xlabel(r"$N$"); ax.set_ylabel(r"$\|f - f_N\|_\infty$")
    ax.set_title(r"(a) $\alpha = 1$ (unscaled)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(frameon=False, fontsize=9, title=r"$e^{-A y^2}$")

    ax = axes[1]
    for A, c in zip(A_list, colours):
        ax.semilogy(Ns, err_alpha_matched[A], "-s", color=c, lw=1.1, mfc="none",
                    label=f"$A={A}$")
    ax.set_xlabel(r"$N$"); ax.set_ylabel(r"$\|f - f_N\|_\infty$")
    ax.set_title(r"(b) $\alpha = \sqrt{2A}$ (matched)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(frameon=False, fontsize=9)

    ax = axes[2]
    ax.semilogy(Ns_q, err_qho, "-D", color=NAVY, lw=1.1,
                label=r"$f = \psi_0$")
    ax.set_xlabel(r"$N$")
    ax.set_ylabel(r"$\|\psi_0 - f_N\|_\infty$")
    ax.set_title("(c) Quantum oscillator: basis matches physics")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(frameon=False, fontsize=9)

    fig.tight_layout()
    save_fig(fig, OUTPUT_DIR, "hermite_width_mismatch")
    plt.close(fig)

    print(f"[20.4] saved figure to {OUTPUT_DIR / 'hermite_width_mismatch.pdf'}")
    for A in A_list:
        print(f"  A={A:.1f}  alpha=1  {err_alpha1[A]}")
        print(f"         matched   {err_alpha_matched[A]}")


if __name__ == "__main__":
    make_figure()
