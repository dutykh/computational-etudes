#!/usr/bin/env python3
"""
laguerre_vs_tln.py

Chapter 20: Spectral Methods on Unbounded Intervals.
Computational Etude 20.5: Laguerre functions versus  TL_n  on [0, infty).

Boyd (2000) warns that Laguerre functions  phi_n(y) = exp(-y/2) L_n(y)
perform well ONLY when the target decays exponentially.  The
rational-Chebyshev  TL_n(y) = T_n((y - L)/(y + L))  basis is far more
flexible: it handles exponential decay, algebraic decay, and
asymptotes-to-a-constant at y = infty with equal ease.

We compare both bases on two model functions:

  (i)   f_1(y) = exp(-y)           (pure exponential)
  (ii)  f_2(y) = 1 / (1 + y)       (pure algebraic decay, 1/y at infty)

The etude's lesson: basis choice is asymptotic modelling.  Laguerre
wins on (i) by a small factor; TL_n wins on (ii) spectacularly.

Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from numpy.polynomial.laguerre import laggauss

from unbounded_common import (CORAL, NAVY, TEAL, output_dir_for, save_fig,
                              setup_matplotlib, cheb_eval, dct1_coeffs,
                              tl_map_forward, tl_map_inverse)

SCRIPT_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = output_dir_for(SCRIPT_DIR)


def laguerre_poly(n, y):
    """Laguerre polynomials L_0, L_1, ..., L_n evaluated at y, via the
    three-term recurrence  (n+1) L_{n+1} = (2n+1 - y) L_n - n L_{n-1}.
    """
    y = np.atleast_1d(y).astype(float)
    L = np.zeros((n + 1, y.size))
    L[0] = 1.0
    if n >= 1:
        L[1] = 1.0 - y
    for k in range(1, n):
        L[k + 1] = ((2 * k + 1 - y) * L[k] - k * L[k - 1]) / (k + 1)
    return L


def laguerre_phi(n, y):
    """Laguerre functions phi_n(y) = exp(-y/2) L_n(y) (un-normalised)."""
    return np.exp(-y / 2.0) * laguerre_poly(n, y)


def laguerre_expand(f, N):
    """Expansion coefficients of f in the non-normalised basis
    phi_n(y) = exp(-y/2) L_n(y).  Uses Gauss-Laguerre quadrature.

    Orthogonality:  int_0^infty exp(-y) L_m(y) L_n(y) dy = delta_mn.
    Therefore <f, phi_n>_{L^2(0,infty)} = int_0^infty f(y) exp(-y/2)
    L_n(y) dy, and the expansion is
        f(y) = sum_n a_n phi_n(y),
        a_n = int_0^infty f(y) exp(-y/2) L_n(y) dy / <phi_n, phi_n>
            = 2 int_0^infty f(y) exp(-y/2) L_n(y) dy
    (because  int exp(-y) L_n^2 = 1  <=>  int phi_n^2 = int exp(-y) L_n^2 exp(0) = 1/2).
    We verify this numerically below.
    """
    x, w = laggauss(N + 32)
    # laggauss: int_0^infty e^{-x} g(x) dx ~ sum w_i g(x_i)
    L = laguerre_poly(N, x)
    # a_n = <f, phi_n> / <phi_n, phi_n>
    # <f, phi_n> = int f(y) exp(-y/2) L_n(y) dy
    #            = int e^{-x} [e^{x/2} f(x) L_n(x)] dx  (taking y = x so e^{-x} provided by Laguerre weight)
    # But we want int f * exp(-y/2) L_n WITHOUT the full e^{-y}, so we insert e^{x/2}:
    fv = f(x)
    # int f(y) exp(-y/2) L_n(y) dy = sum_i w_i exp(x_i) * f(x_i) exp(-x_i/2) L_n(x_i)
    #                              = sum_i w_i exp(x_i/2) f(x_i) L_n(x_i)
    inner = (w * np.exp(x / 2.0) * fv) * L
    numer = inner.sum(axis=1)
    # <phi_n, phi_n> = int_0^infty e^{-y} L_n^2 dy = 1
    # But we are using NON-normalised phi_n = e^{-y/2} L_n, so
    # <phi_n, phi_n> = int e^{-y} L_n^2 dy = 1.  (The e^{-y/2} squared gives e^{-y}.)
    return numer


def laguerre_eval(coeffs, y):
    return np.exp(-y / 2.0) * (coeffs[:, None] * laguerre_poly(len(coeffs) - 1, y)).sum(axis=0)


def laguerre_error(f, N):
    coeffs = laguerre_expand(f, N)
    y_fine = np.linspace(0.001, 60.0, 4001)
    truth = f(y_fine)
    approx = laguerre_eval(coeffs, y_fine)
    return np.max(np.abs(approx - truth))


def tln_expand(f, N, L):
    """Expansion of f on [0, infty) in the TL_n basis.  We interpolate
    f on the semi-infinite grid induced by the map y = L (1+x)/(1-x),
    where x are Chebyshev-Gauss-Lobatto interior nodes.
    """
    from cheb_matrix import cheb_matrix  # noqa: E402
    import sys
    sys.path.insert(0, str(SCRIPT_DIR.parent / "ch07"))
    from cheb_matrix import cheb_matrix  # noqa: E402
    _, x = cheb_matrix(N)
    # avoid x = 1 (y = infty): drop first node and prepend to re-use DCT-I;
    # instead we use all x including x = 1 but evaluate f there as the limit
    y = np.where(x < 1.0 - 1e-12, tl_map_forward(x, L), 1e16)
    fv = np.where(x < 1.0 - 1e-12, f(y), 0.0)    # assume f(infty) = 0
    return dct1_coeffs(fv)


def tln_error(f, N, L):
    a = tln_expand(f, N, L)
    y_fine = np.linspace(0.001, 60.0, 4001)
    x_fine = tl_map_inverse(y_fine, L)
    approx = cheb_eval(a, x_fine, N)
    return np.max(np.abs(approx - f(y_fine)))


def make_figure():
    import sys
    sys.path.insert(0, str(SCRIPT_DIR.parent / "ch07"))
    setup_matplotlib()

    Ns = np.array([8, 12, 16, 24, 32, 48, 64, 96])

    err_lag_exp = [laguerre_error(lambda y: np.exp(-y), N) for N in Ns]
    err_lag_alg = [laguerre_error(lambda y: 1.0 / (1.0 + y), N) for N in Ns]

    # pick a good L for TL_n by trial
    L_exp = 2.0
    L_alg = 5.0
    err_tln_exp = [tln_error(lambda y: np.exp(-y), N, L_exp) for N in Ns]
    err_tln_alg = [tln_error(lambda y: 1.0 / (1.0 + y), N, L_alg) for N in Ns]

    fig, axes = plt.subplots(1, 2, figsize=(10.0, 3.8))

    ax = axes[0]
    ax.semilogy(Ns, np.array(err_lag_exp) + 1e-18, "-o", color=CORAL,
                lw=1.1, label="Laguerre")
    ax.semilogy(Ns, np.array(err_tln_exp) + 1e-18, "-s", color=TEAL,
                lw=1.1, mfc="none", label=fr"$TL_n$, $L={L_exp}$")
    ax.set_xlabel(r"$N$"); ax.set_ylabel(r"$\|f - f_N\|_\infty$")
    ax.set_title(r"(a) $f(y) = e^{-y}$ (exponential decay)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(frameon=False, fontsize=9)

    ax = axes[1]
    ax.semilogy(Ns, np.array(err_lag_alg) + 1e-18, "-o", color=CORAL,
                lw=1.1, label="Laguerre")
    ax.semilogy(Ns, np.array(err_tln_alg) + 1e-18, "-s", color=TEAL,
                lw=1.1, mfc="none", label=fr"$TL_n$, $L={L_alg}$")
    ax.set_xlabel(r"$N$"); ax.set_ylabel(r"$\|f - f_N\|_\infty$")
    ax.set_title(r"(b) $f(y) = 1/(1+y)$ (algebraic decay)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(frameon=False, fontsize=9)

    fig.tight_layout()
    save_fig(fig, OUTPUT_DIR, "laguerre_vs_tln")
    plt.close(fig)

    print(f"[20.5] saved figure to {OUTPUT_DIR / 'laguerre_vs_tln.pdf'}")
    for N, le, la, te, ta in zip(Ns, err_lag_exp, err_lag_alg, err_tln_exp, err_tln_alg):
        print(f"  N={N:3d}  lag_exp={le:.2e}  lag_alg={la:.2e}  "
              f"tln_exp={te:.2e}  tln_alg={ta:.2e}")


if __name__ == "__main__":
    make_figure()
