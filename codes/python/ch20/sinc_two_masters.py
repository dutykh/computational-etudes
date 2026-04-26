#!/usr/bin/env python3
"""
sinc_two_masters.py

Chapter 20: Spectral Methods on Unbounded Intervals.
Computational Etude 20.3: Sinc expansions -- two masters, span and spacing.

The Whittaker cardinal (sinc) expansion on y in (-infty, +infty) is

    f(y) ~ sum_{j=-N/2}^{N/2}  f(j h)  sinc( (y - j h) / h ),
    sinc(z) = sin(pi z) / (pi z).

It is the simplest infinite-interval pseudospectral method.  Unlike
Chebyshev on a finite interval, a sinc series depends on TWO
parameters: the grid spacing h and the truncation N.  The total error
splits into the bandwidth error  E_W(h)  (because the infinite Fourier
integral has been truncated to [-pi/h, pi/h]) and the grid-span error
E_DT(N h)  (because the infinite sum has been cut at |j| > N/2).

Optimally, h is chosen so that the two errors balance; for an
exponentially decaying function like sech(y) this gives h ~ 1/sqrt(N)
and a convergence rate  exp(-q sqrt(N))  (subgeometric).

This etude reproduces Boyd's (2000) subgeometric picture:
 * fix N, vary h: a V-shaped error curve with optimum h^* ~ 1/sqrt(N);
 * follow h = sqrt(pi^2 / (2 N)) (the theoretical optimum for sech):
   the error descends as  exp(-pi sqrt(N/2)).

Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from unbounded_common import (CORAL, NAVY, TEAL, output_dir_for, save_fig,
                              setup_matplotlib)

SCRIPT_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = output_dir_for(SCRIPT_DIR)


def target(y):
    return 1.0 / np.cosh(y)


def sinc_approx(y, N, h):
    """Evaluate the sinc interpolant at points y."""
    j = np.arange(-N // 2, N // 2 + 1)
    yj = j * h
    fj = target(yj)
    # sinc((y - yj) / h)
    z = (y[:, None] - yj[None, :]) / h
    # use numpy sinc which is sin(pi z)/(pi z)
    S = np.sinc(z)
    return S @ fj


def sinc_error(N, h):
    y_fine = np.linspace(-20.0, 20.0, 4001)
    approx = sinc_approx(y_fine, N, h)
    return np.max(np.abs(approx - target(y_fine)))


def make_figure():
    setup_matplotlib()

    # Panel (a): fix N, sweep h
    N_fix = 48
    hs = np.linspace(0.15, 1.8, 80)
    err_sweep = np.array([sinc_error(N_fix, h) for h in hs])

    # Panel (b): follow h_opt(N) ~ sqrt(pi^2 / (2 N)) for sech(y)
    Ns = np.array([8, 12, 16, 24, 32, 48, 64, 96, 128, 192])
    h_opts = np.sqrt(np.pi ** 2 / (2.0 * Ns))
    err_opt = np.array([sinc_error(N, h) for N, h in zip(Ns, h_opts)])

    # 2x2 layout: V-shape | two-master decomposition
    #             convergence | grid cartoon
    fig, axes = plt.subplots(2, 2, figsize=(11.0, 8.0))

    # ---- (a) V-shape: empirical error vs h at fixed N -----------------
    ax = axes[0, 0]
    ax.semilogy(hs, err_sweep, color=CORAL, lw=1.2)
    j_star = int(np.argmin(err_sweep))
    ax.axvline(hs[j_star], color=NAVY, ls=":", lw=1.0,
               label=fr"$h^* = {hs[j_star]:.2f}$ empirical")
    h_theory = np.sqrt(np.pi ** 2 / (2.0 * N_fix))
    ax.axvline(h_theory, color=TEAL, ls="--", lw=1.0,
               label=fr"$\sqrt{{\pi^2/(2N)}} = {h_theory:.2f}$")
    ax.set_xlabel(r"$h$")
    ax.set_ylabel(r"$\|f - f_N\|_\infty$")
    ax.set_title(fr"(a) Fix $N = {N_fix}$, vary $h$")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(frameon=False, fontsize=9)

    # ---- (b) NEW: literal "two masters" decomposition -----------------
    # bandwidth error E_W(h) ~ (4/pi) exp(-pi^2 / (2 h)) -- Stenger asymptotic.
    # grid-span error E_DT(N h) ~ 4 exp(-N h / 2) -- sech tail at edge.
    # The empirical V-shape is the envelope of these two competing
    # exponentials; their crossover sits exactly at the empirical h*.
    ax = axes[0, 1]
    hs_an = np.linspace(hs.min(), hs.max(), 400)
    E_W  = (4.0 / np.pi) * np.exp(-np.pi ** 2 / (2.0 * hs_an))
    E_DT = 4.0 * np.exp(-N_fix * hs_an / 2.0)
    ax.semilogy(hs_an, E_W, color=NAVY, ls="--", lw=1.0,
                label=r"$E_W(h)\sim e^{-\pi^2/(2h)}$ (bandwidth)")
    ax.semilogy(hs_an, E_DT, color=CORAL, ls="--", lw=1.0,
                label=r"$E_{\rm DT}(Nh)\sim e^{-Nh/2}$ (grid span)")
    ax.semilogy(hs_an, E_W + E_DT, color=TEAL, lw=1.4,
                label="sum: $E_W + E_{\\rm DT}$")
    ax.semilogy(hs, err_sweep, color="gray", lw=0.7, alpha=0.7,
                label="empirical V")
    ax.axvline(hs[j_star], color=NAVY, ls=":", lw=0.8)
    ax.set_xlabel(r"$h$")
    ax.set_ylabel("error")
    ax.set_title("(b) Two masters: bandwidth vs grid-span")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(frameon=False, fontsize=8, loc="lower left")
    ax.set_ylim(1e-12, 1e1)

    # ---- (c) Subgeometric convergence at h = h_opt(N) -----------------
    ax = axes[1, 0]
    ax.semilogy(np.sqrt(Ns), err_opt, "-o", color=TEAL, lw=1.1, mfc="none",
                label=r"sinc at $h = \sqrt{\pi^2/(2N)}$")
    ax.semilogy(np.sqrt(Ns), np.exp(-np.pi * np.sqrt(Ns / 2.0)),
                ls=":", color=NAVY, lw=1.0,
                label=r"$\exp(-\pi\sqrt{N/2})$")
    ax.set_xlabel(r"$\sqrt{N}$")
    ax.set_ylabel(r"$\|f - f_N\|_\infty$")
    ax.set_title("(c) Subgeometric convergence at $h = h^*$")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(frameon=False, fontsize=9)

    # ---- (d) Grid cartoon: span and spacing grow together -------------
    ax = axes[1, 1]
    for k, (Nc, offset) in enumerate([(8, 0.7), (32, 0.3)]):
        hk = np.sqrt(np.pi ** 2 / (2.0 * Nc))
        j = np.arange(-Nc // 2, Nc // 2 + 1)
        yj = j * hk
        ax.scatter(yj, np.full_like(yj, offset), s=14,
                   color=(CORAL if k == 0 else TEAL), marker="x",
                   label=fr"$N={Nc}$, $h={hk:.2f}$")
    y_plot = np.linspace(-10.0, 10.0, 401)
    ax.plot(y_plot, 0.05 + 0.6 * target(y_plot), color=NAVY, lw=1.0)
    ax.set_xlim(-10, 10); ax.set_ylim(0, 1.2)
    ax.set_xlabel(r"$y$")
    ax.set_title(r"(d) Grid span and spacing grow together")
    ax.legend(frameon=False, fontsize=9, loc="upper right")

    fig.tight_layout()
    save_fig(fig, OUTPUT_DIR, "sinc_two_masters")
    plt.close(fig)

    print(f"[20.3] saved figure to {OUTPUT_DIR / 'sinc_two_masters.pdf'}")
    for N, h, e in zip(Ns, h_opts, err_opt):
        print(f"  N={N:4d}  h*={h:.3f}  err={e:.3e}")


if __name__ == "__main__":
    make_figure()
