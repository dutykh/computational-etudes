#!/usr/bin/env python3
"""
mathieu_sideband.py

Chapter 21: Special Tricks for the Spectral Researcher
Computational Etude 21.2: Mathieu's equation and sideband truncation.

Mathieu's equation
        u_xx + [lambda - 2 q cos(2 x)] u = 0
on the periodic interval [0, 2 pi] has eigenfunctions ce_n, se_n with
parities determined by n's parity.  The big lesson of Boyd Chapter 19
is that for large n, the dominant Fourier coefficient sits near mode
n -- not near mode 0.  The right truncation is therefore a *band*
of modes around the carrier n, not the usual low-mode truncation.

Reproducing Boyd Figs 19.1 and 19.2 for ce_15:

  Panel A: |a_n| versus n at q = 10, computed from a high-N Galerkin
           solve.  Carrier at n = 15, sidebands at n = 13, 17 (and tiny
           at n = 11, 19).

  Panel B: eigenvalue correction delta(q) = lambda(q) - n^2 for n = 15,
           computed three ways:
              - delta_full   from the high-N Galerkin
              - delta_5      from the 5x5 sideband matrix
              - delta_3      from the inner 3x3 sideband matrix
           For q < 25 the 5x5 result is indistinguishable from the
           reference at the resolution of the plot.

  Panel C: a "breakdown" panel showing what happens at small carrier
           n where q/n^2 is large.  We pick n = 3 and re-run the
           sideband approximations against a high-N reference: the
           5x5 truncation is no longer enough.

The basis is restricted to the cos(odd m x) class to which ce_15
belongs (Boyd's "Even-Odd" parity).

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
from tricks_common import NAVY, CORAL, TEAL, PURPLE


def odd_modes(M):
    """Return the odd integers 1, 3, 5, ..., (2M - 1)."""
    return 2 * np.arange(1, M + 1) - 1


def galerkin_full(N_modes, q):
    """Build the (N_modes x N_modes) Galerkin matrix on the basis
    {cos(1 x), cos(3 x), cos(5 x), ...} for Mathieu's equation.

    The eigenvalue problem is  M a = lambda a, with
        M_{ii} = (2 i - 1)^2,
        M_{i,j} = q  iff |mode_i - mode_j| = 2.
    """
    modes = odd_modes(N_modes)
    M = np.diag(modes.astype(float) ** 2)
    for i in range(N_modes):
        for j in range(N_modes):
            if abs(modes[i] - modes[j]) == 2:
                M[i, j] = q
    return M, modes


def reference_eigenpair(N_modes, q, n_carrier):
    """Diagonalise the full (N_modes x N_modes) Galerkin matrix and
    return (lambda, eigenvector) for the eigenpair whose eigenvector
    has its largest component at mode n_carrier."""
    M, modes = galerkin_full(N_modes, q)
    lam, V = np.linalg.eigh(M)
    # locate carrier index
    idx_carrier = np.where(modes == n_carrier)[0][0]
    # for each eigenvector, find the dominant mode index
    dominant = np.argmax(np.abs(V), axis=0)
    candidates = np.where(dominant == idx_carrier)[0]
    if len(candidates) == 0:
        # Fallback: pick eigenvalue closest to n^2.
        k = np.argmin(np.abs(lam - n_carrier ** 2))
    else:
        # Among candidates, pick the one with the largest |a_carrier|.
        comps = np.abs(V[idx_carrier, candidates])
        k = candidates[np.argmax(comps)]
    return lam[k], V[:, k], modes


def sideband_eigenvalue(n_carrier, q, half_width):
    """Return the sideband-truncation estimate of the eigenvalue near
    n^2.  half_width = 1 gives a 3x3, half_width = 2 gives a 5x5."""
    # mode indices: n - 2 hw, n - 2(hw-1), ..., n, ..., n + 2 hw
    modes = np.array([n_carrier + 2 * k for k in range(-half_width, half_width + 1)],
                     dtype=int)
    M = np.diag(modes.astype(float) ** 2)
    for i in range(len(modes)):
        for j in range(len(modes)):
            if abs(modes[i] - modes[j]) == 2:
                M[i, j] = q
    lam = np.linalg.eigvalsh(M)
    # Pick the eigenvalue closest to n^2.
    return lam[np.argmin(np.abs(lam - n_carrier ** 2))]


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    p.add_argument("--dump", type=str, default=None)
    args = p.parse_args()

    setup_matplotlib()
    out = output_dir_for(SCRIPT_DIR)

    n_carrier = 15
    N_modes = 64                # carries modes 1, 3, ..., 127
    q_demo = 10.0
    q_axis = np.linspace(0.0, 25.0, 200)

    # --- Panel A: coefficient bar chart at q = 10 ---------------------
    lam_demo, vec_demo, modes = reference_eigenpair(N_modes, q_demo, n_carrier)
    # Normalise: enforce a_carrier > 0.
    idx_carrier = np.where(modes == n_carrier)[0][0]
    vec_demo = vec_demo * np.sign(vec_demo[idx_carrier])
    coeff_abs = np.abs(vec_demo)

    # --- Panel B: delta(q) for ce_15 ----------------------------------
    delta_full = np.empty_like(q_axis)
    delta_3 = np.empty_like(q_axis)
    delta_5 = np.empty_like(q_axis)
    for k, q in enumerate(q_axis):
        lam, _, _ = reference_eigenpair(N_modes, q, n_carrier)
        delta_full[k] = lam - n_carrier ** 2
        delta_3[k] = sideband_eigenvalue(n_carrier, q, 1) - n_carrier ** 2
        delta_5[k] = sideband_eigenvalue(n_carrier, q, 2) - n_carrier ** 2

    # --- Panel C: breakdown at small carrier --------------------------
    n_small = 3
    q_small_axis = np.linspace(0.0, 10.0, 200)
    delta_full_small = np.empty_like(q_small_axis)
    delta_3_small = np.empty_like(q_small_axis)
    delta_5_small = np.empty_like(q_small_axis)
    for k, q in enumerate(q_small_axis):
        lam, _, _ = reference_eigenpair(N_modes, q, n_small)
        delta_full_small[k] = lam - n_small ** 2
        delta_3_small[k] = sideband_eigenvalue(n_small, q, 1) - n_small ** 2
        delta_5_small[k] = sideband_eigenvalue(n_small, q, 2) - n_small ** 2

    # --- Figure -------------------------------------------------------
    fig = plt.figure(figsize=(13.5, 4.4))
    gs = fig.add_gridspec(1, 3, wspace=0.32)

    # Panel A: bar chart
    ax = fig.add_subplot(gs[0, 0])
    ax.bar(modes, coeff_abs, width=1.4, color=NAVY, edgecolor=NAVY)
    ax.set_xlim(0, 31)
    ax.set_xlabel(r"Fourier degree $n$  (cos basis, odd $n$)")
    ax.set_ylabel(r"$|a_n|$")
    ax.set_title(rf"Panel A.  $|a_n|$ for $\mathrm{{ce}}_{{15}}$ at $q={q_demo:.0f}$",
                 fontsize=10)
    ax.axvline(n_carrier, color=CORAL, linewidth=0.8, alpha=0.7,
               linestyle="--", label="carrier $n = 15$")
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(True, alpha=0.25, linewidth=0.4, axis="y")

    # Panel B: delta(q) for n=15
    ax = fig.add_subplot(gs[0, 1])
    ax.plot(q_axis, delta_full, color=NAVY, linewidth=1.4,
            label=r"$\delta_\text{full}$ (high-$N$ reference)")
    ax.plot(q_axis, delta_5, color=TEAL, linewidth=1.0, linestyle="--",
            label=r"$\delta_5$ (5$\times$5 sideband)")
    ax.plot(q_axis, delta_3, color=CORAL, linewidth=1.0, linestyle=":",
            label=r"$\delta_3$ (3$\times$3 sideband)")
    ax.set_xlabel(r"$q$")
    ax.set_ylabel(r"$\delta(q) = \lambda(q) - 15^2$")
    ax.set_title(r"Panel B.  $\mathrm{ce}_{15}$: 5$\times$5 already exact at this scale",
                 fontsize=10)
    ax.legend(loc="upper left", fontsize=9)
    ax.grid(True, alpha=0.25, linewidth=0.4)

    # Panel C: breakdown at n=3
    ax = fig.add_subplot(gs[0, 2])
    ax.plot(q_small_axis, delta_full_small, color=NAVY, linewidth=1.4,
            label=r"$\delta_\text{full}$")
    ax.plot(q_small_axis, delta_5_small, color=TEAL, linewidth=1.0, linestyle="--",
            label=r"$\delta_5$")
    ax.plot(q_small_axis, delta_3_small, color=CORAL, linewidth=1.0, linestyle=":",
            label=r"$\delta_3$")
    ax.set_xlabel(r"$q$")
    ax.set_ylabel(r"$\delta(q) = \lambda(q) - 3^2$")
    ax.set_title(r"Panel C.  $\mathrm{ce}_3$: $q/n^2$ large, sideband breaks",
                 fontsize=10)
    ax.legend(loc="lower left", fontsize=9)
    ax.grid(True, alpha=0.25, linewidth=0.4)

    save_fig(fig, out, "mathieu_sideband")
    plt.close(fig)

    print(f"[Etude 21.2]  Mathieu sideband truncation")
    print(f"  ce_{n_carrier}, q = {q_demo}")
    five_largest = np.argsort(coeff_abs)[::-1][:5]
    for j in five_largest:
        print(f"    n = {modes[j]:2d}  |a_n| = {coeff_abs[j]:.4e}")
    print(f"  delta_full at q={q_demo}: {delta_full[np.argmin(np.abs(q_axis - q_demo))]:+.6f}")
    print(f"  delta_5    at q={q_demo}: {delta_5[np.argmin(np.abs(q_axis - q_demo))]:+.6f}")
    print(f"  delta_3    at q={q_demo}: {delta_3[np.argmin(np.abs(q_axis - q_demo))]:+.6f}")
    print(f"  ce_3, q = 5: full = {delta_full_small[np.argmin(np.abs(q_small_axis - 5.0))]:+.4f}, "
          f"5x5 = {delta_5_small[np.argmin(np.abs(q_small_axis - 5.0))]:+.4f}, "
          f"3x3 = {delta_3_small[np.argmin(np.abs(q_small_axis - 5.0))]:+.4f}")
    print(f"  figure: {out / 'mathieu_sideband.pdf'}")

    dump_json(args.dump, {
        "n_carrier": n_carrier, "N_modes": N_modes, "q_demo": q_demo,
        "modes_demo": modes.tolist(),
        "coeff_abs_demo": coeff_abs.tolist(),
        "q_axis": q_axis.tolist(),
        "delta_full": delta_full.tolist(),
        "delta_5": delta_5.tolist(),
        "delta_3": delta_3.tolist(),
        "n_small": n_small,
        "q_small_axis": q_small_axis.tolist(),
        "delta_full_small": delta_full_small.tolist(),
        "delta_5_small": delta_5_small.tolist(),
        "delta_3_small": delta_3_small.tolist(),
    })


if __name__ == "__main__":
    main()
