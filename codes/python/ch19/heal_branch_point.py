#!/usr/bin/env python3
"""
heal_branch_point.py

Chapter 19: Coordinate Transformations
Computational Etude 19.5: Healing the square-root branch point.

We approximate the bounded-but-non-analytic function

    g(X) = sqrt(1 - X^2),       X in [-1, 1],

which has square-root branch points at X = +/- 1.  In the native
Chebyshev expansion the coefficients decay only algebraically as
O(1 / n^2) (Boyd Eq. (16.18)); the max-norm error is therefore
O(1 / N).

Under the tanh map  X = tanh(y),  y in (-infty, +infty),  the
transformed function becomes

    g(tanh y) = sech(y),

which is entire in y and decays exponentially fast at +/- infinity.
A rational-Chebyshev (TB_n) approximation on y in (-infty, +infty)
therefore recovers geometric convergence (subgeometric-but-exponential
in Boyd's technical vocabulary).

This etude compares:

 (A) direct Chebyshev approximation of g on [-1, 1];
 (B) Chebyshev collocation in  y  on a truncated window [-Y, Y]  of
     sech(y),  which encapsulates the "implicit tanh-mapped" expansion;
     the window Y is chosen generously (8-12).

Both max-norm errors and cofficient decay are plotted.  The chapter's
lesson: a weak singularity in one coordinate can become irrelevant in
another, and spectral accuracy is a property of the function AS SEEN
BY the basis, not of the function alone.

Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
Part of "Computational Etudes: A Spectral Approach"
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR.parent / "ch07"))
from cheb_matrix import cheb_matrix  # noqa: E402

from map_common import CORAL, NAVY, TEAL, output_dir_for, save_fig, setup_matplotlib

OUTPUT_DIR = output_dir_for(SCRIPT_DIR)


def g_of_X(X):
    return np.sqrt(np.maximum(1.0 - X ** 2, 0.0))


def sech(y):
    return 1.0 / np.cosh(y)


# -------------------------------------------------------------- direct Chebyshev
def chebyshev_coeffs(v):
    """DCT-I coefficients of samples v_k at Chebyshev-GL nodes.

    The normalisation is such that v(x) = sum_n a_n T_n(x) exactly in
    interpolation sense on the N+1 grid points.
    """
    N = len(v) - 1
    # use the symmetric DFT trick of Trefethen p. 74
    V = np.concatenate([v, v[N - 1:0:-1]])
    A = np.real(np.fft.fft(V)) / N
    A[0] *= 0.5
    A[N] *= 0.5
    return A[:N + 1]


def chebyshev_error(N):
    _, x = cheb_matrix(N)
    v = g_of_X(x)
    a = chebyshev_coeffs(v)
    # error on a fine mesh, reconstructing via Clenshaw recurrence
    x_fine = np.linspace(-1 + 1e-10, 1 - 1e-10, 2001)
    # evaluate sum a_n T_n(x) using the recurrence
    T_km2 = np.ones_like(x_fine)
    T_km1 = x_fine.copy()
    val = a[0] * T_km2 + (a[1] if N >= 1 else 0.0) * T_km1
    for n in range(2, N + 1):
        T_k = 2.0 * x_fine * T_km1 - T_km2
        val += a[n] * T_k
        T_km2, T_km1 = T_km1, T_k
    return np.max(np.abs(val - g_of_X(x_fine))), np.abs(a)


# ------------------------------------------- mapped: Chebyshev on y in [-Y, Y]
def mapped_error(N, Y):
    """Chebyshev interpolation of sech(y) on y in [-Y, Y] via affine map."""
    D, xi = cheb_matrix(N)           # xi in [-1, 1]
    y = Y * xi                       # y in [-Y, Y]
    f = sech(y)
    a = chebyshev_coeffs(f)
    # reconstruct on fine mesh
    xi_fine = np.linspace(-1, 1, 2001)
    y_fine = Y * xi_fine
    T_km2 = np.ones_like(xi_fine)
    T_km1 = xi_fine.copy()
    val = a[0] * T_km2 + (a[1] if N >= 1 else 0.0) * T_km1
    for n in range(2, N + 1):
        T_k = 2.0 * xi_fine * T_km1 - T_km2
        val += a[n] * T_k
        T_km2, T_km1 = T_km1, T_k
    return np.max(np.abs(val - sech(y_fine))), np.abs(a)


# ----------------------------------------- pointwise error at fixed N (panel b)
def chebyshev_pointwise_error(N, X_fine):
    """Pointwise error of the direct Chebyshev expansion of sqrt(1 - X^2)."""
    _, x = cheb_matrix(N)
    a = chebyshev_coeffs(g_of_X(x))
    T_km2 = np.ones_like(X_fine)
    T_km1 = X_fine.copy()
    val = a[0] * T_km2 + (a[1] if N >= 1 else 0.0) * T_km1
    for n in range(2, N + 1):
        T_k = 2.0 * X_fine * T_km1 - T_km2
        val += a[n] * T_k
        T_km2, T_km1 = T_km1, T_k
    return np.abs(val - g_of_X(X_fine))


def mapped_pointwise_error(N, Y, X_fine):
    """Pointwise error of the tanh-mapped expansion sech(y), evaluated at
    physical points X = tanh(y).  We invert: given X_fine in (-1, 1), the
    matching y is arctanh(X_fine).  The error in the physical variable
    equals the error of the y-expansion at that y, since
    g(tanh y) = sech(y) exactly."""
    D, xi = cheb_matrix(N)
    y_nodes = Y * xi
    a = chebyshev_coeffs(sech(y_nodes))
    # X_fine in (-1, 1), corresponding y_fine in (-arctanh(1), arctanh(1))
    Xc = np.clip(X_fine, -1.0 + 1e-12, 1.0 - 1e-12)
    y_fine = np.arctanh(Xc)
    # only y_fine in [-Y, Y] are inside the truncation window; outside we
    # have nothing meaningful to report -- mask them.
    inside = np.abs(y_fine) <= Y
    xi_fine = np.where(inside, y_fine / Y, 0.0)
    T_km2 = np.ones_like(xi_fine)
    T_km1 = xi_fine.copy()
    val = a[0] * T_km2 + (a[1] if N >= 1 else 0.0) * T_km1
    for n in range(2, N + 1):
        T_k = 2.0 * xi_fine * T_km1 - T_km2
        val += a[n] * T_k
        T_km2, T_km1 = T_km1, T_k
    err = np.abs(val - sech(y_fine))
    err[~inside] = np.nan        # outside the truncation window
    return err


def make_figure():
    setup_matplotlib()

    Ns = np.array([8, 16, 24, 32, 48, 64, 96, 128])
    err_direct, err_mapped = [], []
    for N in Ns:
        e1, _ = chebyshev_error(N)
        err_direct.append(e1)
        e2, _ = mapped_error(N, Y=10.0)
        err_mapped.append(e2)
    err_direct = np.array(err_direct)
    err_mapped = np.array(err_mapped)

    # coefficient magnitudes at N = 64 for each method
    _, a_direct = chebyshev_error(64)
    _, a_mapped = mapped_error(64, Y=10.0)

    fig, axes = plt.subplots(2, 2, figsize=(11.0, 7.6))

    # (a) the two target functions
    ax = axes[0, 0]
    Xp = np.linspace(-1, 1, 401)
    yp = np.linspace(-6, 6, 401)
    ax.plot(Xp, g_of_X(Xp), color=CORAL, lw=1.4, label=r"$\sqrt{1-X^2}$")
    ax.plot(yp, sech(yp), color=TEAL, lw=1.4, ls="--",
            label=r"$\mathrm{sech}\,y$")
    ax.set_xlabel("coordinate")
    ax.set_ylabel("value")
    ax.set_title(r"(a) Two views of the same function")
    ax.legend(frameon=False, fontsize=10)
    ax.grid(True, alpha=0.3)

    # (b) NEW: pointwise error at fixed N = 32 -- where the methods fail
    ax = axes[0, 1]
    Nshow = 32
    X_fine = np.linspace(-1.0 + 1e-6, 1.0 - 1e-6, 1001)
    err_direct_pw = chebyshev_pointwise_error(Nshow, X_fine)
    err_mapped_pw = mapped_pointwise_error(Nshow, 10.0, X_fine)
    ax.semilogy(X_fine, err_direct_pw + 1e-18, "-", color=CORAL, lw=1.4,
                label=r"direct, $\sqrt{1-X^2}$")
    ax.semilogy(X_fine, err_mapped_pw + 1e-18, "-", color=TEAL, lw=1.4,
                label=r"tanh-mapped, $\mathrm{sech}\,y$")
    ax.set_xlabel(r"$X$")
    ax.set_ylabel(r"$|g(X) - g_N(X)|$")
    ax.set_title(fr"(b) Pointwise error at $N = {Nshow}$")
    ax.set_xlim(-1, 1)
    ax.set_ylim(1e-16, 1e0)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="lower center", frameon=False, fontsize=9)

    # (c) convergence rates
    ax = axes[1, 0]
    ax.loglog(Ns, err_direct, "-o", color=CORAL, lw=1.1,
              label=r"direct, $\sqrt{1-X^2}$")
    ax.loglog(Ns, err_mapped + 1e-18, "-s", color=TEAL, lw=1.1, mfc="none",
              label=r"tanh-mapped, $\mathrm{sech}\,y$")
    ax.loglog(Ns, 1.0 / Ns, ":", color=NAVY, lw=0.8, label=r"$1/N$")
    ax.set_xlabel(r"$N$")
    ax.set_ylabel(r"$\|g - g_N\|_\infty$")
    ax.set_title("(c) Convergence: algebraic vs geometric")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(frameon=False, fontsize=9)

    # (d) coefficient decay at N = 64
    ax = axes[1, 1]
    n_direct = np.arange(len(a_direct))
    n_mapped = np.arange(len(a_mapped))
    ax.semilogy(n_direct, a_direct + 1e-18, "x", color=CORAL, ms=4,
                label="direct Chebyshev")
    ax.semilogy(n_mapped, a_mapped + 1e-18, "o", color=TEAL, ms=3,
                mfc="none", label="tanh-mapped")
    ax.set_xlabel(r"Chebyshev index $n$")
    ax.set_ylabel(r"$|a_n|$")
    ax.set_title("(d) Coefficient decay, $N = 64$")
    ax.grid(True, which="both", alpha=0.3)
    ax.set_ylim(1e-17, 1e1)
    ax.legend(frameon=False, fontsize=9)

    fig.tight_layout()
    save_fig(fig, OUTPUT_DIR, "heal_branch_point")
    plt.close(fig)

    print(f"[19.5] saved figure to {OUTPUT_DIR / 'heal_branch_point.pdf'}")
    for N, ed, em in zip(Ns, err_direct, err_mapped):
        print(f"  N={N:4d}  direct={ed:.3e}  mapped={em:.3e}")


if __name__ == "__main__":
    make_figure()
