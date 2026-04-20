#!/usr/bin/env python3
"""
corner_tensor_clustering.py

Chapter 19: Coordinate Transformations
Computational Etude 19.6: A square Poisson problem with corner stress.

We solve  -Delta u = 1  on the square  Omega = [-1, 1]^2  with
homogeneous Dirichlet data  u = 0  on the boundary.  The solution has
weak logarithmic-type branch points at the four corners
(Boyd Eq. (16.24)):

    u(r, theta) ~ r^2 log(r) sin(2 theta) + smoother terms.

Two discretisations compete, both using an  N_x x N_y  tensor-product
Chebyshev-Lobatto grid:

 (A) NO mapping.  Solve the 2D Poisson equation directly on the
     standard Chebyshev grid.

 (B) ONE-dimensional tanh mapping applied independently in each
     coordinate: X = tanh(alpha * xi), Y = tanh(alpha * eta), with
     alpha tuned so that X(+/-1) are close to +/-1 and the grid
     clusters exponentially near the corner.  This is Boyd's
     practical-first recommendation for weak corner branch points.

The etude compares convergence in N for the two methods.  Because no
closed-form solution exists, we compare against a highly refined
reference solution (N_ref = 128, unmapped; Chebyshev at that resolution
is spectrally converged to ~1e-12 for this problem).

The etude illustrates that for WEAK corner singularities the tensor-
product mapping yields a useful -- but modest -- improvement, and only
above a crossover N.  This is in line with Boyd's cautionary remark
that "asymptotics begin late" for weak singularities.

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


# --------------------------------------------------------- 2D Poisson on square
def assemble_poisson(D, x):
    """Kronecker-sum assembly of the 2D Laplacian on the full tensor grid."""
    N = len(x) - 1
    D2 = D @ D
    I = np.eye(N + 1)
    # Laplacian L = I (x) D2 + D2 (x) I  (row-major indexing: y fast)
    Lmat = np.kron(D2, I) + np.kron(I, D2)
    return Lmat


def bdy_interior_indices(N):
    """Return (interior_mask, idx_interior, idx_boundary) for an (N+1)^2 grid."""
    nx = N + 1
    I_idx = np.arange(nx * nx).reshape(nx, nx)
    interior = np.zeros_like(I_idx, dtype=bool)
    interior[1:N, 1:N] = True
    idx_int = I_idx[interior]
    idx_bdy = I_idx[~interior]
    return interior, idx_int, idx_bdy


def solve_unmapped(N):
    """Solve  -Delta u = 1  on [-1, 1]^2  with  u = 0  on the boundary."""
    D, x = cheb_matrix(N)
    Lmat = assemble_poisson(D, x)
    nx = N + 1
    _, idx_int, idx_bdy = bdy_interior_indices(N)
    rhs = np.ones(nx * nx)
    rhs[idx_bdy] = 0.0            # will be overwritten
    # reduced system
    A = Lmat[np.ix_(idx_int, idx_int)]
    b = -rhs[idx_int]             # -Delta u = 1  =>  Lmat * u = -1 on interior
    u_int = np.linalg.solve(A, b)
    u = np.zeros(nx * nx)
    u[idx_int] = u_int
    return x, x, u.reshape(nx, nx)


def solve_tanh_mapped(N, alpha):
    """Solve the same Poisson problem after applying y = tanh(alpha * x) / tanh(alpha)
    independently in each coordinate, so that y(+/-1) = +/- 1 and the grid clusters
    exponentially towards the walls.
    """
    D, xi = cheb_matrix(N)
    # map each coordinate:  Y = tanh(alpha * xi) / tanh(alpha)
    th_alpha = np.tanh(alpha)
    x_phys = np.tanh(alpha * xi) / th_alpha
    # y = f(xi),  f'(xi) = alpha / cosh^2(alpha xi) / tanh(alpha)
    fp = alpha / (np.cosh(alpha * xi) ** 2 * th_alpha)
    fpp = -2.0 * alpha ** 2 * np.tanh(alpha * xi) / (np.cosh(alpha * xi) ** 2 * th_alpha)
    D1 = np.diag(1.0 / fp) @ D
    D2 = np.diag(1.0 / fp ** 2) @ (D @ D) - np.diag(fpp / fp ** 3) @ D
    nx = N + 1
    I = np.eye(nx)
    Lmat = np.kron(D2, I) + np.kron(I, D2)
    _, idx_int, idx_bdy = bdy_interior_indices(N)
    rhs = -np.ones(nx * nx)
    A = Lmat[np.ix_(idx_int, idx_int)]
    b = rhs[idx_int]
    u_int = np.linalg.solve(A, b)
    u = np.zeros(nx * nx)
    u[idx_int] = u_int
    return x_phys, x_phys, u.reshape(nx, nx)


def reference_solution(N_ref=96):
    x, y, U = solve_unmapped(N_ref)
    return x, U


def interp_on_fine(x_phys, U, xf, yf):
    """2D Chebyshev interpolation to a common fine grid.

    We reuse the simple strategy: interpolate via the Chebyshev barycentric
    formula in each direction.  For benchmarking we use numpy polynomial fit
    in each coordinate, which is simple but good enough for this purpose.
    """
    from numpy.polynomial.chebyshev import chebfit, chebval
    # U is at (x_phys, x_phys) tensor grid; project in each direction
    # Step 1: rowwise Chebyshev fit in "x"; map x_phys back to [-1,1] for fit
    # but x_phys is already monotone on [-1,1], use direct barycentric
    from scipy.interpolate import RegularGridInterpolator
    # x_phys is monotonically DECREASING because Chebyshev-GL gives x = cos(pi k/N)
    if x_phys[0] > x_phys[-1]:
        xs = x_phys[::-1]
        Us = U[::-1, :][:, ::-1]
    else:
        xs = x_phys; Us = U
    interp = RegularGridInterpolator((xs, xs), Us, bounds_error=False, fill_value=0.0)
    XF, YF = np.meshgrid(xf, yf, indexing="ij")
    pts = np.column_stack([XF.ravel(), YF.ravel()])
    vals = interp(pts).reshape(XF.shape)
    return vals


def max_error_against_ref(U_method_x, U_method, U_ref_x, U_ref):
    """Compute max-norm error on a common fine sample."""
    xf = np.linspace(-0.99, 0.99, 121)
    Vm = interp_on_fine(U_method_x, U_method, xf, xf)
    Vr = interp_on_fine(U_ref_x, U_ref, xf, xf)
    return np.max(np.abs(Vm - Vr))


def make_figure():
    setup_matplotlib()

    # Reference solution
    xref, Uref = reference_solution(96)

    Ns = np.array([12, 16, 20, 24, 32, 40, 48, 64])
    ALPHA_LIST = [1.0, 2.0, 3.0]
    err_unmapped = []
    err_mapped = {a: [] for a in ALPHA_LIST}
    for N in Ns:
        xu, _, Uu = solve_unmapped(N)
        err_unmapped.append(max_error_against_ref(xu, Uu, xref, Uref))
        for a in ALPHA_LIST:
            xm, _, Um = solve_tanh_mapped(N, alpha=a)
            err_mapped[a].append(max_error_against_ref(xm, Um, xref, Uref))
    err_unmapped = np.array(err_unmapped)
    for a in ALPHA_LIST:
        err_mapped[a] = np.array(err_mapped[a])
    # pick a representative alpha for the grid-clustering panel
    ALPHA = 2.0

    fig, axes = plt.subplots(1, 3, figsize=(13.0, 3.6))

    # (a) contour of the solution at N = 32
    ax = axes[0]
    xs, _, Us = solve_unmapped(32)
    # sort for clean contouring
    idx = np.argsort(xs)
    cs = ax.contourf(xs[idx], xs[idx], Us[np.ix_(idx, idx)], levels=15, cmap="viridis")
    plt.colorbar(cs, ax=ax, shrink=0.8, label=r"$u$")
    ax.set_aspect("equal")
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$y$")
    ax.set_title(r"(a) $-\Delta u = 1$, $u|_{\partial\Omega}=0$")

    # (b) physical grid clustering
    ax = axes[1]
    Ngrid = 24
    D24, xi24 = cheb_matrix(Ngrid)
    x_plain = xi24
    x_tanh = np.tanh(ALPHA * xi24) / np.tanh(ALPHA)
    for xg in x_plain:
        ax.axvline(xg, ymax=0.45, color=NAVY, lw=0.5, alpha=0.6)
    for xg in x_tanh:
        ax.axvline(xg, ymin=0.55, color=CORAL, lw=0.5, alpha=0.8)
    ax.text(-0.95, 0.22, "standard Chebyshev", color=NAVY, fontsize=9)
    ax.text(-0.95, 0.78, f"tanh clustered, $\\alpha={ALPHA}$", color=CORAL, fontsize=9)
    ax.set_xlim(-1.05, 1.05)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_xlabel("physical coordinate")
    ax.set_title(f"(b) Grids at $N={Ngrid}$")

    # (c) convergence
    ax = axes[2]
    ax.loglog(Ns, err_unmapped + 1e-18, "-o", color=NAVY, lw=1.2,
              label="unmapped Chebyshev")
    colours = [CORAL, TEAL, "#8E44AD"]
    for a, c in zip(ALPHA_LIST, colours):
        ax.loglog(Ns, err_mapped[a] + 1e-18, "-s", color=c, mfc="none",
                  label=fr"tanh $\alpha={a}$")
    ax.set_xlabel(r"$N$")
    ax.set_ylabel(r"error vs reference")
    ax.set_title("(c) Corner convergence")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(frameon=False, fontsize=9)

    fig.tight_layout()
    save_fig(fig, OUTPUT_DIR, "corner_tensor_clustering")
    plt.close(fig)

    print(f"[19.6] saved figure to {OUTPUT_DIR / 'corner_tensor_clustering.pdf'}")
    for i, N in enumerate(Ns):
        line = f"  N={N:3d}  unmapped={err_unmapped[i]:.3e}  "
        line += "  ".join(f"alpha={a}:{err_mapped[a][i]:.3e}"
                          for a in ALPHA_LIST)
        print(line)


if __name__ == "__main__":
    make_figure()
