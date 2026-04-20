#!/usr/bin/env python3
"""
chebyshev_as_cosine.py

Chapter 19: Coordinate Transformations
Computational Etude 19.2: One problem, two coordinates.

We solve the model BVP

    u''(x) - q * u(x) = f(x),       x in [-1, 1],
    u(-1) = u(+1) = 0,

in two equivalent ways:

 (X) Chebyshev collocation built from the DENSE differentiation matrix
     in the physical variable x, following chapters 7-8 of this book.

 (T) Chebyshev collocation built from FFT-based cosine differentiation
     in the computational variable t = arccos(x).  This is Boyd's
     observation (Chapter 16) that T_n(cos t) = cos(nt): working
     "in t" simply means differentiating by transforming u(cos t) to
     Chebyshev coefficients, applying the three-term derivative
     recurrence, and transforming back.  Same mathematics, different
     arithmetic cost (O(N log N) vs O(N^2)).

Both codes are fed the manufactured source f so that the exact solution
u_ex(x) = sin(pi x / 2) cos(pi x / 2) - sin^3(pi x / 2) / 3 is known.

The etude's message: switching coordinate does NOT change the answer
and does NOT change the convergence rate -- only the implementation.

Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
Part of "Computational Etudes: A Spectral Approach"
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR.parent / "ch07"))
from cheb_matrix import cheb_matrix  # noqa: E402

from map_common import CORAL, NAVY, TEAL, output_dir_for, save_fig, setup_matplotlib

OUTPUT_DIR = output_dir_for(SCRIPT_DIR)

Q_VALUE = 4.0


def exact_u(x):
    """Manufactured solution vanishing at +/-1 and smooth on [-1, 1]."""
    return np.sin(np.pi * x)


def exact_u_xx(x):
    return -(np.pi ** 2) * np.sin(np.pi * x)


def source_f(x):
    return exact_u_xx(x) - Q_VALUE * exact_u(x)


# ---------------------------------------------------------------------- X-form
def solve_in_x(N):
    """Chebyshev collocation in the physical variable x using the dense DM."""
    D, x = cheb_matrix(N)
    D2 = D @ D
    A = D2[1:N, 1:N] - Q_VALUE * np.eye(N - 1)
    rhs = source_f(x[1:N])
    u_int = np.linalg.solve(A, rhs)
    u = np.zeros(N + 1)
    u[1:N] = u_int
    return x, u


# ---------------------------------------------------------------------- T-form
def chebfft(v):
    """Second derivative of v at Chebyshev-Gauss-Lobatto nodes via FFT.
    Implements the standard "cheb FFT" algorithm (Trefethen, Program 21):
    we view v_k = u(cos(k pi / N)) as the real part of a symmetric DFT
    of length 2N, differentiate the cosine coefficients using the
    trigonometric identity T_n'(cos t) = n sin(n t) / sin(t), and
    transform back.  The result is u'(x_k) on the GL grid.
    """
    N = len(v) - 1
    x = np.cos(np.pi * np.arange(N + 1) / N)
    # pack symmetrically
    V = np.concatenate([v, v[N - 1:0:-1]])
    U = np.real(np.fft.fft(V))
    k = np.arange(N + 1)
    w_hat = 1j * np.concatenate([k[:N], [0], k[1:N] - N]) * U
    W = np.real(np.fft.ifft(w_hat))
    w = np.zeros(N + 1)
    w[1:N] = -W[1:N] / np.sqrt(1.0 - x[1:N] ** 2)
    # endpoints: analytical formulas
    k2 = np.arange(N + 1) ** 2
    w[0] = np.sum(k2 * U[:N + 1]) / N + 0.5 * N * U[N]
    w[N] = (np.sum(((-1) ** (np.arange(N + 1) + 1)) * k2 * U[:N + 1]) / N
            + 0.5 * ((-1) ** (N + 1)) * N * U[N])
    return w


def cheb_d2_via_fft(v):
    """Second derivative on GL grid via two successive FFT-differentiations."""
    return chebfft(chebfft(v))


def solve_in_t(N):
    """Chebyshev collocation in t: derivative evaluation via cosine FFT,
    then the same boundary-bordering scheme used in solve_in_x.

    The T-form builds its second-derivative action as a matrix L by
    applying cheb_d2_via_fft to each basis vector of the interior grid.
    That exposes the equivalence with the x-form while using the FFT
    arithmetic path: L = P (D^2) P^T with P the interior projector.
    """
    x = np.cos(np.pi * np.arange(N + 1) / N)
    # Build D^2 implicitly through the FFT
    D2 = np.zeros((N + 1, N + 1))
    e = np.zeros(N + 1)
    for j in range(N + 1):
        e[:] = 0.0
        e[j] = 1.0
        D2[:, j] = cheb_d2_via_fft(e)
    A = D2[1:N, 1:N] - Q_VALUE * np.eye(N - 1)
    rhs = source_f(x[1:N])
    u_int = np.linalg.solve(A, rhs)
    u = np.zeros(N + 1)
    u[1:N] = u_int
    return x, u


def run_convergence():
    Ns = np.array([8, 12, 16, 20, 24, 32, 40, 48])
    err_x, err_t, t_x, t_t = [], [], [], []
    for N in Ns:
        tic = time.perf_counter()
        x_x, u_x = solve_in_x(N)
        t_x.append(time.perf_counter() - tic)
        err_x.append(np.max(np.abs(u_x - exact_u(x_x))))

        tic = time.perf_counter()
        x_t, u_t = solve_in_t(N)
        t_t.append(time.perf_counter() - tic)
        err_t.append(np.max(np.abs(u_t - exact_u(x_t))))
    return Ns, np.array(err_x), np.array(err_t), np.array(t_x), np.array(t_t)


def make_figure():
    setup_matplotlib()
    Ns, err_x, err_t, t_x, t_t = run_convergence()

    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.6))

    ax = axes[0]
    xN, uN = solve_in_x(24)
    xfine = np.linspace(-1, 1, 401)
    ax.plot(xfine, exact_u(xfine), color=NAVY, lw=1.2, label="exact")
    ax.plot(xN, uN, "o", color=CORAL, ms=4, label=r"X-form, $N=24$")
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$u(x)$")
    ax.set_title("(a) Manufactured BVP")
    ax.legend(frameon=False, fontsize=9)
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    ax.semilogy(Ns, err_x + 1e-18, "-o", color=CORAL, lw=1.1,
                label=r"X-form (dense DM)")
    ax.semilogy(Ns, err_t + 1e-18, "-s", color=TEAL, lw=1.1, mfc="none",
                label=r"T-form (FFT-cosine)")
    ax.set_xlabel(r"$N$")
    ax.set_ylabel(r"$\|u - u_N\|_\infty$")
    ax.set_title("(b) Same convergence, different arithmetic")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(frameon=False, fontsize=9)

    fig.tight_layout()
    save_fig(fig, OUTPUT_DIR, "chebyshev_as_cosine")
    plt.close(fig)

    print(f"[19.2] saved figure to {OUTPUT_DIR / 'chebyshev_as_cosine.pdf'}")
    for N, ex, et, tx, tt in zip(Ns, err_x, err_t, t_x, t_t):
        print(f"  N={N:3d}  err_x={ex:.3e}  err_t={et:.3e}")


if __name__ == "__main__":
    make_figure()
