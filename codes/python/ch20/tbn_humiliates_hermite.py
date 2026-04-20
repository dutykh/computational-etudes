#!/usr/bin/env python3
"""
tbn_humiliates_hermite.py

Chapter 20: Spectral Methods on Unbounded Intervals.
Computational Etude 20.6: The function that humiliates Hermite.

Boyd's emblematic example.  The function

    f(y) = 1 / (1 + y^2),       y in (-infty, +infty),

is innocuous -- it is bounded, smooth, and symmetric -- but it has
simple poles at y = +/- i and therefore decays only ALGEBRAICALLY as
1/y^2 at infinity.  The Hermite series

    f ~ sum a_n psi_n(y)                (y native scale)

converges only algebraically because psi_n(y) ~ exp(-y^2/2) decays far
too fast to resolve an algebraically decaying function.  The sinc
series is similarly trapped.

By contrast, the rational-Chebyshev basis TB_n(y) is, for L = 1,
precisely  TB_0(y) - TB_2(y) = 2 * f(y), so f is represented by TWO
basis functions exactly.  This is no coincidence: the TB_n map
y = L cot(t) converts algebraic decay into smooth boundedness in t,
which a Fourier cosine series resolves geometrically.

The etude is the chapter's central shock: basis choice is asymptotic
modelling.

Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from numpy.polynomial.hermite import hermgauss

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR.parent / "ch07"))

from unbounded_common import (CORAL, NAVY, TEAL, output_dir_for, save_fig,
                              setup_matplotlib, cheb_eval, dct1_coeffs,
                              tb_map_forward, tb_map_inverse)
from hermite_width_mismatch import hermite_psi


def target(y):
    return 1.0 / (1.0 + y ** 2)


# ---------------- Hermite expansion (reused from etude 20.4)
def hermite_expand(N, alpha):
    x, w = hermgauss(N + 40)
    y_nodes = x / alpha
    fv = target(y_nodes)
    psi = hermite_psi(N, x)
    integrand = fv * psi * np.exp(x ** 2)
    return (integrand * w).sum(axis=1) / np.sqrt(alpha)


def hermite_error(N, alpha):
    coeffs = hermite_expand(N, alpha)
    y = np.linspace(-40.0, 40.0, 8001)
    psi = hermite_psi(N, alpha * y)
    approx = np.sqrt(alpha) * (coeffs[:, None] * psi).sum(axis=0)
    return np.max(np.abs(approx - target(y)))


# ---------------- sinc expansion (reused from etude 20.3)
def sinc_error(N, h):
    j = np.arange(-N // 2, N // 2 + 1)
    yj = j * h
    fj = target(yj)
    y = np.linspace(-40.0, 40.0, 8001)
    z = (y[:, None] - yj[None, :]) / h
    S = np.sinc(z)
    approx = S @ fj
    return np.max(np.abs(approx - target(y)))


# ---------------- TB_n expansion
def tbn_expand(N, L):
    from cheb_matrix import cheb_matrix
    _, x = cheb_matrix(N)
    y = tb_map_forward(x, L)
    fv = target(y)
    # endpoints x = +/- 1 correspond to y = +/- infty where f = 0
    fv[0] = 0.0
    fv[-1] = 0.0
    return dct1_coeffs(fv)


def tbn_error(N, L):
    a = tbn_expand(N, L)
    y = np.linspace(-40.0, 40.0, 8001)
    x = tb_map_inverse(y, L)
    approx = cheb_eval(a, x, N)
    return np.max(np.abs(approx - target(y)))


def make_figure():
    setup_matplotlib()
    Ns = np.array([8, 12, 16, 24, 32, 48, 64, 96, 128])
    L = 1.0

    err_her = [hermite_error(N, 1.0) for N in Ns]
    err_sinc = [sinc_error(N, np.sqrt(np.pi ** 2 / (2 * N))) for N in Ns]
    err_tbn = [tbn_error(N, L) for N in Ns]

    fig, axes = plt.subplots(1, 2, figsize=(10.0, 3.8))

    ax = axes[0]
    y = np.linspace(-8.0, 8.0, 401)
    ax.plot(y, target(y), color=NAVY, lw=1.2, label=r"$f(y) = 1/(1+y^2)$")
    ax.set_xlabel(r"$y$"); ax.set_ylabel(r"$f$")
    ax.set_title("(a) An innocuous target")
    ax.legend(frameon=False, fontsize=10)
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    ax.loglog(Ns, np.array(err_her) + 1e-18, "-o", color=CORAL, lw=1.1,
              label=r"Hermite, $\alpha = 1$")
    ax.loglog(Ns, np.array(err_sinc) + 1e-18, "-s", color=TEAL,
              lw=1.1, mfc="none", label=r"sinc, $h \sim 1/\sqrt{N}$")
    ax.loglog(Ns, np.array(err_tbn) + 1e-18, "-^", color=NAVY, lw=1.2,
              label=r"$TB_n$, $L=1$")
    ax.set_xlabel(r"$N$"); ax.set_ylabel(r"$\|f - f_N\|_\infty$")
    ax.set_title("(b) Algebraic, subgeometric, geometric")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(frameon=False, fontsize=9)

    fig.tight_layout()
    save_fig(fig, OUTPUT_DIR, "tbn_humiliates_hermite")
    plt.close(fig)

    print(f"[20.6] saved figure to {OUTPUT_DIR / 'tbn_humiliates_hermite.pdf'}")
    for N, h, s, t in zip(Ns, err_her, err_sinc, err_tbn):
        print(f"  N={N:4d}  Her={h:.3e}  sinc={s:.3e}  TB={t:.3e}")


OUTPUT_DIR = output_dir_for(SCRIPT_DIR)

if __name__ == "__main__":
    make_figure()
