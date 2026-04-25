#!/usr/bin/env python3
"""
hilbert_fourier.py

Chapter 21: Special Tricks for the Spectral Researcher
Computational Etude 21.8: Hilbert transform via Fourier multiplier.

The Hilbert transform of a function f is defined by the principal-value
integral

        H{f}(y) = (1/pi) PV \\int  f(x) / (x - y)  dx.

When f is *periodic* with period 2 pi and analytic, its Hilbert
transform on the circle is an exact Fourier multiplier:

        H{f}(y) = -i \\sum_{k != 0} sgn(k) c_k exp(i k y),
        f(x)     =       \\sum_k        c_k exp(i k x).

In code this is three lines: FFT, multiply by -i sgn(k), inverse FFT.
The resulting algorithm inherits all the convergence properties of
the underlying Fourier transform: for analytic periodic f, errors
fall geometrically (in fact super-geometrically here, because the
test function below has factorially-decaying Fourier coefficients).

We test on  f(x) = exp(cos x),  whose Fourier series is

        exp(cos x) = I_0(1) + 2 \\sum_{k=1}^\\infty I_k(1) cos(k x),

so that

        H{f}(y) = 2 \\sum_{k=1}^\\infty I_k(1) sin(k y).

The modified Bessel coefficients  I_k(1)  decay factorially:
        I_k(1) ~ (1/2)^k / k!  as k -> infinity,
which gives the algorithm a particularly clean convergence picture.

Author: Dr. Denys Dutykh (Khalifa University of Science and Technology)
Part of "Computational Etudes: A Spectral Approach"
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.special import iv

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
from tricks_common import setup_matplotlib, output_dir_for, save_fig, dump_json
from tricks_common import NAVY, CORAL, TEAL


def hilbert_via_fft(f_values):
    """Apply the Fourier-multiplier Hilbert transform to a vector of
    samples on a uniform grid spanning one period [-pi, pi)."""
    N = len(f_values)
    F = np.fft.fft(f_values)
    k = np.fft.fftfreq(N, d=1.0 / N)   # integer wavenumbers 0, 1, ..., -1
    multiplier = -1j * np.sign(k)
    multiplier[0] = 0.0
    return np.real(np.fft.ifft(multiplier * F))


def f_test(x):
    return np.exp(np.cos(x))


def hilbert_exact(y, K_max=40):
    """H{exp(cos x)}(y) = 2 sum_{k=1}^infty I_k(1) sin(k y)."""
    val = np.zeros_like(y)
    for k in range(1, K_max + 1):
        val += 2.0 * iv(k, 1.0) * np.sin(k * y)
    return val


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    p.add_argument("--dump", type=str, default=None)
    args = p.parse_args()

    setup_matplotlib()
    out = output_dir_for(SCRIPT_DIR)

    # --- Convergence study --------------------------------------------
    Ns = np.array([4, 6, 8, 10, 12, 14, 16, 20, 24, 32])
    err = np.empty_like(Ns, dtype=float)
    for k, N in enumerate(Ns):
        x = -np.pi + 2.0 * np.pi * np.arange(N) / N
        Hf = hilbert_via_fft(f_test(x))
        Hf_exact = hilbert_exact(x)
        err[k] = np.max(np.abs(Hf - Hf_exact))

    # --- Visual at N = 32 ---------------------------------------------
    N_show = 32
    x_show = -np.pi + 2.0 * np.pi * np.arange(N_show) / N_show
    Hf_show = hilbert_via_fft(f_test(x_show))
    x_dense = np.linspace(-np.pi, np.pi, 401)
    f_dense = f_test(x_dense)
    Hf_dense_exact = hilbert_exact(x_dense)

    # ---- Figure -----------------------------------------------------
    fig = plt.figure(figsize=(11.5, 4.4))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.0], wspace=0.30)

    # Panel A: f and H f
    ax = fig.add_subplot(gs[0, 0])
    ax.plot(x_dense, f_dense, color=NAVY, linewidth=1.4,
            label=r"$f(x) = e^{\cos x}$")
    ax.plot(x_dense, Hf_dense_exact, color=TEAL, linewidth=1.4,
            label=r"$H\{f\}(y)$ exact")
    ax.plot(x_show, Hf_show, "o", color=CORAL, markerfacecolor="white",
            markersize=5, label=fr"$H\{{f\}}_N(y_j)$ via FFT, $N={N_show}$")
    ax.axhline(0.0, color="gray", linewidth=0.4, alpha=0.5)
    ax.set_xlabel(r"$x$ (or $y$)")
    ax.set_ylabel("amplitude")
    ax.set_title(r"$f(x) = e^{\cos x}$ and its Hilbert transform on the circle",
                 fontsize=10)
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(True, alpha=0.25, linewidth=0.4)
    ax.set_xlim(-np.pi - 0.1, np.pi + 0.1)

    # Panel B: convergence
    ax = fig.add_subplot(gs[0, 1])
    ax.semilogy(Ns, err + 1e-18, "o-", color=NAVY, markerfacecolor="white",
                markersize=6, linewidth=1.0,
                label=r"max $|H\{f\}_N - H\{f\}|$ on the grid")
    # Reference: factorial decay of I_k(1).  At high k, I_k(1) ~ (1/2)^k / k!.
    k_ref = np.arange(1, 30)
    bound = 2.0 * iv(Ns / 2, 1.0)
    ax.semilogy(Ns, bound, "--", color=CORAL, linewidth=0.8,
                label=r"$\sim 2\,I_{N/2}(1)$ (Bessel-tail bound)")
    ax.axhline(1e-15, color="gray", linewidth=0.4, alpha=0.5)
    ax.set_xlabel(r"$N$ (Fourier modes per period)")
    ax.set_ylabel("max error")
    ax.set_title("super-geometric convergence on a periodic analytic test",
                 fontsize=10)
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(True, alpha=0.25, linewidth=0.4)
    ax.set_ylim(1e-17, 5)

    save_fig(fig, out, "hilbert_fourier")
    plt.close(fig)

    print(f"[Etude 21.8]  Hilbert transform via Fourier multiplier")
    print(f"  Test: f(x) = exp(cos x), periodic on [-pi, pi]")
    for k, N in enumerate(Ns):
        print(f"  N = {N:3d}: max err = {err[k]:.3e}")
    print(f"  figure: {out / 'hilbert_fourier.pdf'}")

    dump_json(args.dump, {
        "Ns": Ns.tolist(),
        "err": err.tolist(),
    })


if __name__ == "__main__":
    main()
