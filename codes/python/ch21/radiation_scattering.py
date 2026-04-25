#!/usr/bin/env python3
"""
radiation_scattering.py

Chapter 21: Special Tricks for the Spectral Researcher
Computational Etude 21.5: 1-D Schroedinger scattering by a sech^2 potential
                          via Boyd's radiation-basis trick.

Problem.   The 1-D Schroedinger equation in scattering form,
            psi'' + [k^2 - V(x)] psi = 0,
with V(x) = -v sech^2(x), v = 1, has the asymptotic ansatz
            psi(x) ~ exp(i k x) + alpha exp(-i k x)   as x -> -inf,
            psi(x) ~ beta exp(i k x)                  as x -> +inf.
The reflection coefficient is R = |alpha|^2.  For this potential
the closed form (Lamb 1980; Boyd Eq 19.31) is

   R = (1 + cos(2 pi sqrt(v + 1/4))) /
       (cosh(2 pi k) + cos(2 pi sqrt(v + 1/4))).

Trick.  Standard rational Chebyshev TB_n functions decay to constants
at +/-infinity and therefore CANNOT represent the non-decaying
sinusoidal asymptotes themselves.  Boyd's repair is to *augment*
the basis with two 'radiation' functions

   H(x) cos(k x),   H(x) sin(k x),     H(x) = (1 + tanh x) / 2,

so that the unknown asymptotic constants are carried by their
coefficients a_{N-2}, a_{N-1}.  Two pseudospectral solves on the
same N-by-N matrix give the auxiliary functions C(x) and S(x)
(Boyd Eq 19.20-19.28), and asymptotic matching (Eq 19.29) yields
alpha exactly via a 2-by-2 linear system.

We reproduce Boyd Table 19.1: R(k) for v = 1 at k = 0.3, 0.6, ..., 3.0,
using the rational Chebyshev parameter ell = 2 (Boyd's choice) and
N = 48 collocation points.  Absolute errors are at the 10^-8 to 10^-9
level across all ten wavenumbers, even though R itself spans seven
decades from 4.2 x 10^-1 down to 2.2 x 10^-8 -- and a method without
the radiation augmentation simply cannot represent the non-decaying
plane-wave parts of the solution at any cost.

Author: Dr. Denys Dutykh (Khalifa University of Science and Technology)
Part of "Computational Etudes: A Spectral Approach"
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
from tricks_common import setup_matplotlib, output_dir_for, save_fig, dump_json
from tricks_common import NAVY, CORAL, TEAL


V_DEPTH = 1.0
ELL = 2.0


def V(x):
    return -V_DEPTH / np.cosh(x) ** 2


def H(x):
    return 0.5 * (1.0 + np.tanh(x))


def Hp(x):
    return 0.5 / np.cosh(x) ** 2


def Hpp(x):
    return -np.tanh(x) / np.cosh(x) ** 2


def reflection_exact(k):
    s = 2.0 * np.pi * np.sqrt(V_DEPTH + 0.25)
    return (1.0 + np.cos(s)) / (np.cosh(2.0 * np.pi * k) + np.cos(s))


def collocation_grid(N, ell):
    """Boyd's rational-Chebyshev collocation: t_i = pi(2i-1)/(2N), x_i = ell cot(t_i)."""
    i = np.arange(1, N + 1)
    t = np.pi * (2.0 * i - 1.0) / (2.0 * N)
    x = ell / np.tan(t)
    return x, t


def basis_TB(N, t):
    """TB_n(x_i) = cos(n t_i) for n = 0, ..., N-3.  Returns (N x (N-2)) matrix."""
    n = np.arange(N - 2)             # 0, 1, ..., N-3
    return np.cos(np.outer(t, n))


def basis_TB_double_prime(N, t, ell):
    """TB_n''(x_i) = (-n/ell^2) sin^3(t_i) [n cos(n t_i) sin(t_i)
                                       + 2 sin(n t_i) cos(t_i)]."""
    n = np.arange(N - 2)
    nT = np.outer(np.ones_like(t), n)              # broadcast (N x (N-2))
    tT = np.outer(t, np.ones_like(n.astype(float))) # broadcast
    s, c = np.sin(tT), np.cos(tT)
    return -nT / ell ** 2 * s ** 3 * (
        nT * np.cos(nT * tT) * s + 2.0 * np.sin(nT * tT) * c
    )


def solve_scattering(N, ell, k):
    """Solve C-tilde and S-tilde via Boyd's radiation-augmented collocation.
    Returns alpha (complex) and the diagnostic constant 'drift' = sum_n a_n
    -- the asymptotic constant that should be small if the basis is rich
    enough."""
    x, t = collocation_grid(N, ell)

    # Basis values and second derivatives at collocation points.
    Phi  = np.zeros((N, N))
    Phipp = np.zeros((N, N))
    Phi[:, :N - 2]  = basis_TB(N, t)
    Phipp[:, :N - 2] = basis_TB_double_prime(N, t, ell)
    # Radiation function H cos(kx):
    Phi[:, N - 2] = H(x) * np.cos(k * x)
    Phipp[:, N - 2] = (Hpp(x) * np.cos(k * x)
                       - 2.0 * k * Hp(x) * np.sin(k * x)
                       - k ** 2 * H(x) * np.cos(k * x))
    # Radiation function H sin(kx):
    Phi[:, N - 1] = H(x) * np.sin(k * x)
    Phipp[:, N - 1] = (Hpp(x) * np.sin(k * x)
                       + 2.0 * k * Hp(x) * np.cos(k * x)
                       - k ** 2 * H(x) * np.sin(k * x))

    # Galerkin-type collocation matrix:  M_{ij} = phi_j''(x_i) + (k^2 - V(x_i)) phi_j(x_i).
    M = Phipp + np.outer(k ** 2 - V(x), np.ones(N)) * Phi

    # Right-hand sides.
    f = V(x) * np.cos(k * x)
    g = V(x) * np.sin(k * x)

    a = np.linalg.solve(M, f)        # coefficients of C-tilde
    b = np.linalg.solve(M, g)        # coefficients of S-tilde

    gamma1 = a[N - 2] + 1.0
    gamma2 = a[N - 1]
    sigma1 = b[N - 2]
    sigma2 = b[N - 1] + 1.0

    # Asymptotic matching, Boyd Eq 19.29.
    A = np.array([[gamma1 + sigma2, sigma1 - gamma2],
                  [gamma2 - sigma1, gamma1 + sigma2]])
    rhs = np.array([sigma2 - gamma1, -sigma1 - gamma2])
    re, im = np.linalg.solve(A, rhs)
    alpha = complex(re, im)

    drift = np.sum(a[:N - 2])
    return alpha, drift, gamma1, gamma2, sigma1, sigma2


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    p.add_argument("--dump", type=str, default=None)
    args = p.parse_args()

    setup_matplotlib()
    out = output_dir_for(SCRIPT_DIR)

    N = 48
    k_axis = np.array([0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4, 2.7, 3.0])
    R_num = np.empty_like(k_axis)
    R_exact = np.empty_like(k_axis)
    drifts = np.empty_like(k_axis)
    for i, k in enumerate(k_axis):
        alpha, drift, *_ = solve_scattering(N, ELL, k)
        R_num[i] = abs(alpha) ** 2
        R_exact[i] = reflection_exact(k)
        drifts[i] = abs(drift)
    abs_err = np.abs(R_num - R_exact)

    # ---- Figure -------------------------------------------------------
    fig = plt.figure(figsize=(11.5, 4.4))
    gs = fig.add_gridspec(1, 2, wspace=0.30)

    ax = fig.add_subplot(gs[0, 0])
    ax.semilogy(k_axis, R_exact, color=NAVY, linewidth=1.4,
                label=r"$R$ exact (Boyd Eq 19.31)")
    ax.semilogy(k_axis, R_num, "o", color=CORAL, markerfacecolor="white",
                markersize=8, label=fr"$R$ numerical, $N = {N}$, $\ell = {ELL:.0f}$")
    ax.set_xlabel(r"wavenumber $k$")
    ax.set_ylabel(r"reflection coefficient $R$")
    ax.set_title(r"$\mathrm{sech}^2$ scattering: spectral $R$ matches the closed form",
                 fontsize=10)
    ax.legend(loc="lower left", fontsize=10)
    ax.grid(True, alpha=0.25, linewidth=0.4)
    ax.set_ylim(1e-10, 2)

    ax = fig.add_subplot(gs[0, 1])
    ax.semilogy(k_axis, abs_err + 1e-18, "s-", color=NAVY, markerfacecolor="white",
                markersize=6, linewidth=1.0,
                label=r"$|R_\text{num} - R_\text{exact}|$")
    ax.semilogy(k_axis, drifts + 1e-18, "^--", color=TEAL, markerfacecolor="white",
                markersize=5, linewidth=0.8,
                label=r"$|\sum_n a_n|$ (asymptotic-constant drift)")
    ax.set_xlabel(r"wavenumber $k$")
    ax.set_ylabel("error magnitude")
    ax.set_title("error stays well below the spectacularly small $R$",
                 fontsize=10)
    ax.legend(loc="upper right", fontsize=10)
    ax.grid(True, alpha=0.25, linewidth=0.4)

    save_fig(fig, out, "radiation_scattering")
    plt.close(fig)

    print(f"[Etude 21.5]  Schroedinger sech^2 scattering with radiation basis")
    print(f"  v = {V_DEPTH}, ell = {ELL}, N = {N}")
    print(f"   k        R_numerical       R_exact         abs error")
    for i, k in enumerate(k_axis):
        print(f"  {k:.1f}  {R_num[i]:.7e}   {R_exact[i]:.7e}   {abs_err[i]:.2e}")
    print(f"  figure: {out / 'radiation_scattering.pdf'}")

    dump_json(args.dump, {
        "N": N, "ell": ELL, "v": V_DEPTH,
        "k_axis": k_axis.tolist(),
        "R_numerical": R_num.tolist(),
        "R_exact": R_exact.tolist(),
        "abs_err": abs_err.tolist(),
        "drifts": drifts.tolist(),
    })


if __name__ == "__main__":
    main()
