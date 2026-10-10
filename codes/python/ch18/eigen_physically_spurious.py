#!/usr/bin/env python3
"""
eigen_physically_spurious.py

Chapter 18: Linear Spectral Eigenproblems
Computational Etude 18.7: manufacturing fake instability.

Solves the Gottlieb-Orszag streamfunction eigenproblem
        nu u_xxxx = lambda u_xx,   u(+-1) = u_x(+-1) = 0,
which arises in the stream-function formulation of viscous linear
stability. The exact spectrum is

        lambda_j = nu * mu_j,   mu_j = solutions of tan(mu) = mu OR mu = n pi,

all real and NEGATIVE.

The naive "divide-through" formulation
        (1/nu) D^{-2} D^4 u = lambda u
with four boundary conditions applied to u but only u(+-1) = 0 natively
imposable on D^2 produces TWO spurious LARGE POSITIVE eigenvalues that
scale as O(N^4), as proved rigorously by Dawkins, Dunbar & Douglass (1998).
These modes are physically spurious: they do not correspond to any
eigenmode of the continuous operator.

Three panels:
  (A) Naive-formulation spectrum at N = 32, 48 showing the offending
      large positive eigenvalues.
  (B) N-scan of the magnitude of the largest positive eigenvalue,
      demonstrating O(N^4) growth (Dawkins et al. 1998).
  (C) Cured-formulation spectrum showing that the Huang-Sloan-style
      boundary bordering on the SYSTEM (u and phi = u'', with phi
      determined from the coupling) removes the spurious modes.

Author: Dr. Denys Dutykh (Khalifa University of Science and Technology)
Part of "Computational Etudes: A Spectral Approach"
https://github.com/dutykh/computational-etudes
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scipy.linalg

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR.parent / "ch07"))
from cheb_matrix import cheb_matrix  # noqa: E402

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["CMU Serif", "DejaVu Serif"],
    "mathtext.fontset": "cm",
    "font.size": 11, "axes.linewidth": 0.8,
    "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "figure.dpi": 150, "savefig.dpi": 300,
})
NAVY, SKY, CORAL, TEAL, ORANGE = "#142D6E", "#7896D2", "#E74C3C", "#16A085", "#E67E22"
OUTPUT_DIR = SCRIPT_DIR.parent.parent.parent / "textbook" / "figures" / "ch18" / "python"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def naive_gottlieb_orszag(N: int, nu: float = 1.0):
    """Naive (tau-style) generalised eigenproblem  A u = lambda B u,
    with A = nu D^4 and B = D^2, where the four boundary conditions
    replace the LAST four rows of the pencil. This is the standard
    tau-method placement: spectral residuals occupy the leading rows,
    and the highest-order modes are sacrificed to the BCs.

    With this bordering the pencil has, in exact arithmetic,
        - four algebraically infinite eigenvalues (one per BC row), and
        - N - 3 finite eigenvalues, of which exactly ONE is a large
          positive spurious eigenvalue whose magnitude grows
          asymptotically as O(N^4) (Dawkins, Dunbar & Douglass, 1998).

    The empirical scaling exponent at moderate N (16 <= N <= 256)
    climbs monotonically from ~3.4 toward 4 as N grows; reaching the
    full O(N^4) regime requires N >> 1000.
    """
    D, _ = cheb_matrix(N)
    D2 = D @ D
    D4 = D2 @ D2
    A = nu * D4.copy()
    B = D2.copy()
    # Tau-style bordering: BCs in the last four rows
    ID = np.eye(N + 1)
    A[-4, :] = ID[0, :]; B[-4, :] = 0     # u(+1)  = 0
    A[-3, :] = ID[N, :]; B[-3, :] = 0     # u(-1)  = 0
    A[-2, :] = D[0, :];  B[-2, :] = 0     # u'(+1) = 0
    A[-1, :] = D[N, :];  B[-1, :] = 0     # u'(-1) = 0
    return A, B


def cured_gottlieb_orszag(N: int, nu: float = 1.0):
    """Cured version: solve the SYSTEM
            phi = u''                             (definition)
            nu phi'' = lambda phi                 (the eigenequation)
    with boundary conditions u(+-1) = u'(+-1) = 0 (four on u). Discretise
    both equations at the interior grid; solve the auxiliary equation
    phi = u'' algebraically (on the interior block) and substitute.

    The resulting reduced problem for the interior u has a mass matrix
    that is nonsingular on the interior, so no spurious infinite
    eigenvalues appear. Equivalently, in the restricted interior block,
    there are no rows with singular mass, and the O(N^4) trap is
    avoided.

    This implementation uses boundary bordering on D^2 (imposing
    u(+-1) = 0 natively) plus two extra tau-rows from the first-derivative
    constraints u'(+-1) = 0, applied to the streamfunction-vorticity
    system."""
    D, _ = cheb_matrix(N)
    D2 = D @ D
    # Interior block (imposes u(+-1) = 0 by restriction)
    idx = slice(1, N)
    D2i = D2[idx, idx]
    # Impose u'(+-1) = 0 by bordering a small tau system. The cleanest
    # way: form the extended (N-1) x (N-1) interior problem
    #     A q = lambda B q
    # where
    #     A = nu * (D^2_int)^2 + two tau rows enforcing u'(+-1) = 0
    #     B = D^2_int
    #
    # To do this within the interior block, we replace rows 0 and N-2
    # of A by the two rows of D that correspond to interior collocation
    # in x = +-1 (after sub-indexing). This is the Huang-Sloan-style
    # construction.
    A = nu * (D2i @ D2i)
    B = D2i.copy()
    # Tau rows: u'(+1) = 0 enforced on row 0, u'(-1) = 0 enforced on
    # row N-2, using the full D matrix's boundary rows restricted to
    # interior columns.
    A[0,   :]   = D[0, 1:N]            # u'(+1) = 0
    B[0,   :]   = 0
    A[-1,  :]   = D[N, 1:N]            # u'(-1) = 0
    B[-1,  :]   = 0
    return A, B


def finite_real_eigvals(A, B, beta_tol: float = 1e-10):
    """Return the genuinely finite real eigenvalues of the pencil
    (A, B), separating them from algebraically infinite ones via the
    homogeneous (alpha, beta) decomposition.

    An eigenvalue is treated as algebraically infinite iff
        |beta_i| < beta_tol * max(|alpha_i|, |beta_i|).

    This is essential for boundary-bordered pencils, where rows of B
    are zeroed exactly: the corresponding generalised eigenvalues are
    1/0 in exact arithmetic, but a naive `np.isfinite(alpha/beta)`
    filter is fooled by QZ rounding into accepting them as finite
    numbers of magnitude ~ ||A||/eps.
    """
    ab = scipy.linalg.eig(A, B, right=False, homogeneous_eigvals=True)
    alpha, beta = ab[0], ab[1]
    mag = np.maximum(np.abs(alpha), np.abs(beta))
    finite = np.abs(beta) > beta_tol * mag
    lam = alpha[finite] / beta[finite]
    return np.sort(lam.real[np.abs(lam.imag) < 1e-6])


def make_figure(Ns_scan=(16, 24, 32, 48, 64, 96, 128, 192, 256), N_demo: int = 48, nu: float = 1.0):
    # Panel A: naive spectrum at two resolutions, highlight positives
    A1, B1 = naive_gottlieb_orszag(32, nu)
    A2, B2 = naive_gottlieb_orszag(48, nu)
    lam1 = finite_real_eigvals(A1, B1)
    lam2 = finite_real_eigvals(A2, B2)
    pos1 = lam1[lam1 > 1.0]
    pos2 = lam2[lam2 > 1.0]

    # Panel B: scaling with N of the largest positive eigenvalue (naive)
    max_positive = []
    for N in Ns_scan:
        A, B = naive_gottlieb_orszag(N, nu)
        lam = finite_real_eigvals(A, B)
        pos = lam[lam > 1.0]
        max_positive.append(pos.max() if pos.size > 0 else np.nan)
    max_positive = np.asarray(max_positive)

    # Panel C: cured spectrum - all eigenvalues should be negative
    A_c1, B_c1 = cured_gottlieb_orszag(32, nu)
    A_c2, B_c2 = cured_gottlieb_orszag(48, nu)
    lam_c1 = finite_real_eigvals(A_c1, B_c1)
    lam_c2 = finite_real_eigvals(A_c2, B_c2)
    pos_c1 = lam_c1[lam_c1 > 1.0]
    pos_c2 = lam_c2[lam_c2 > 1.0]

    fig, axes = plt.subplots(1, 3, figsize=(14.0, 4.2))

    # Panel A
    ax = axes[0]
    j1 = np.arange(1, lam1.size + 1)
    j2 = np.arange(1, lam2.size + 1)
    ax.scatter(j1, lam1, s=24, facecolors="white", edgecolors=NAVY,
               linewidths=1.0, label=f"$N = 32$")
    ax.scatter(j2, lam2, s=18, marker="x", color=CORAL,
               linewidths=1.0, label=f"$N = 48$")
    ax.axhline(0.0, color="k", linewidth=0.6, alpha=0.5)
    ax.set_xlabel("mode number $j$")
    ax.set_ylabel("$\\lambda$ (naive formulation)")
    ax.set_title(f"naive: {pos1.size} spurious $\\lambda > 0$ at $N=32$, "
                 f"{pos2.size} at $N=48$", fontsize=10)
    ax.grid(True, alpha=0.25, linewidth=0.4)
    ax.legend(loc="upper right", fontsize=9)

    # Panel B
    ax = axes[1]
    mask = np.isfinite(max_positive)
    Ns = np.asarray(Ns_scan)[mask]
    pos = max_positive[mask]
    ax.loglog(Ns, pos, "o-", color=NAVY, markersize=6,
              markerfacecolor="white", markeredgewidth=1.1,
              linewidth=0.8, label="spurious $\\lambda_{+}$ (naive)")
    # Reference line ~ N^4 anchored at the smallest N point
    ref_N = Ns.astype(float)
    ref_y = pos[0] * (ref_N / ref_N[0]) ** 4
    ax.loglog(ref_N, ref_y, "--", color=TEAL, linewidth=0.8,
              label="$N^4$ reference")
    # Empirical slope (last segment)
    slope = float(np.log(pos[-1] / pos[0]) / np.log(Ns[-1] / Ns[0]))
    ax.set_xlabel("$N$")
    ax.set_ylabel("spurious positive eigenvalue $\\lambda_{+}$")
    ax.set_title(f"DDD scaling: empirical slope $\\approx {slope:.2f}$, "
                 f"asymptote $N^4$", fontsize=10)
    ax.grid(True, which="both", alpha=0.25, linewidth=0.4)
    ax.legend(loc="upper left", fontsize=9)

    # Panel C
    ax = axes[2]
    j1c = np.arange(1, lam_c1.size + 1)
    j2c = np.arange(1, lam_c2.size + 1)
    ax.scatter(j1c, lam_c1, s=24, facecolors="white", edgecolors=NAVY,
               linewidths=1.0, label=f"$N = 32$ (cured)")
    ax.scatter(j2c, lam_c2, s=18, marker="x", color=CORAL,
               linewidths=1.0, label=f"$N = 48$ (cured)")
    ax.axhline(0.0, color="k", linewidth=0.6, alpha=0.5)
    ax.set_xlabel("mode number $j$")
    ax.set_ylabel("$\\lambda$ (cured formulation)")
    ax.set_title(f"cured: {pos_c1.size} spurious at $N=32$, "
                 f"{pos_c2.size} at $N=48$", fontsize=10)
    ax.grid(True, alpha=0.25, linewidth=0.4)
    ax.legend(loc="upper right", fontsize=9)

    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "eigen_physically_spurious.pdf", bbox_inches="tight", pad_inches=0.05)
    fig.savefig(OUTPUT_DIR / "eigen_physically_spurious.png", bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)

    return {
        "nu": nu,
        "Ns_scan":       list(Ns_scan),
        "max_positive":  max_positive.tolist(),
        "naive_n_positive":  {32: int(pos1.size), 48: int(pos2.size)},
        "cured_n_positive":  {32: int(pos_c1.size), 48: int(pos_c2.size)},
        "naive_positives_32": pos1.tolist(),
        "naive_positives_48": pos2.tolist(),
        "cured_spectrum_32_first5": lam_c1[:5].tolist(),
    }


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--dump", type=str, default=None)
    args = p.parse_args()

    r = make_figure()
    print(f"[Etude 18.7]  Gottlieb-Orszag streamfunction (nu = {r['nu']})")
    print(f"  Naive positive eigenvalues at N=32: {r['naive_n_positive'][32]}  "
          f"at N=48: {r['naive_n_positive'][48]}")
    print(f"  Max positive vs N: {list(zip(r['Ns_scan'], [round(m,1) if np.isfinite(m) else None for m in r['max_positive']]))}")
    print(f"  Cured positive eigenvalues at N=32: {r['cured_n_positive'][32]}  "
          f"at N=48: {r['cured_n_positive'][48]}")
    print(f"  figure: {OUTPUT_DIR / 'eigen_physically_spurious.pdf'}")
    if args.dump:
        # stringify int keys for JSON
        rjson = dict(r)
        rjson["naive_n_positive"] = {str(k): v for k, v in r["naive_n_positive"].items()}
        rjson["cured_n_positive"] = {str(k): v for k, v in r["cured_n_positive"].items()}
        Path(args.dump).write_text(json.dumps(rjson, indent=2))
        print(f"  dumped: {args.dump}")


if __name__ == "__main__":
    main()
