#!/usr/bin/env python3
"""
symbolic_boundary_layer.py

Chapter 21: Special Tricks for the Spectral Researcher
Computational Etude 21.9: a boundary-layer BVP solved symbolically by a
                          four-term Galerkin method.

Carrier and Pearson's textbook example (Boyd Eq 20.3):

        epsilon^2 u_xx - u = -1,    u(-1) = u(1) = 0,

with exact solution u(x) = 1 - cosh(x/epsilon) / cosh(1/epsilon).
For small epsilon the solution has a steep boundary layer of width
~ epsilon at each endpoint.

Galerkin precepts for symbolic work (Boyd Chapter 20):
   * use Galerkin (not collocation) at small N: the residual is
     orthogonalised against the test functions, so the resulting
     algebra is more compact;
   * incorporate the boundary conditions into the basis by writing
     u = (1 - x^2) p(x) -- 'basis recombination' instead of boundary
     bordering;
   * exploit the x -> -x symmetry of the problem by using only even
     powers in p(x);
   * apply UNWEIGHTED moments (test functions {1, x^2, x^4, ...})
     because we know the action is in the boundary layers and weighting
     would dilute their contribution.

Following Boyd Eq 20.5, take the four-term symmetric ansatz
        u_4(x) = (1 - x^2) (a_0 + a_2 x^2 + a_4 x^4 + a_6 x^6).
Substituting into the ODE gives the residual R_4(x) (Eq 20.6).  Demand
        integral_{-1}^{1} x^{2j} R_4 dx = 0,    j = 0, 1, 2, 3,
giving four linear equations in (a_0, a_2, a_4, a_6).  Solve symbolically
in epsilon.  The coefficients turn out to be RATIONAL functions of
epsilon (Boyd Eq 20.8 with a common denominator D(epsilon) of degree 8).

We then evaluate the resulting analytic approximation pointwise on a
fine grid and compare with the exact cosh solution for epsilon ranging
from 1 down to 1/20, reproducing Boyd Table 20.3.

Author: Dr. Denys Dutykh (Khalifa University of Science and Technology)
Part of "Computational Etudes: A Spectral Approach"
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import sympy as sp

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
from tricks_common import setup_matplotlib, output_dir_for, save_fig, dump_json
from tricks_common import NAVY, CORAL, TEAL


def galerkin_4_term():
    """Solve Boyd Eq 20.5-20.7 symbolically.  Returns
    (u4_expr_in_x_eps, dict_a_in_eps, denom_D_eps)."""
    x, eps = sp.symbols("x epsilon", positive=True)
    a = sp.symbols("a_0 a_2 a_4 a_6")

    # Trial solution u_4 with built-in boundary conditions.
    p = a[0] + a[1] * x ** 2 + a[2] * x ** 4 + a[3] * x ** 6
    u4 = (1 - x ** 2) * p

    # Residual of the original equation.
    R4 = eps ** 2 * sp.diff(u4, x, 2) - u4 + 1

    # Symmetric Galerkin moments: integrate against x^{2j} for j=0,1,2,3.
    eqs = [sp.integrate(x ** (2 * j) * R4, (x, -1, 1)) for j in range(4)]
    sol = sp.solve(eqs, a, dict=True)[0]

    # Common-denominator form.
    sol_simpl = {ai: sp.together(sp.simplify(sol[ai])) for ai in a}
    denom = sp.together(sol_simpl[a[0]]).as_numer_denom()[1]
    return u4.subs(sol_simpl), sol_simpl, denom, eps, x


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    p.add_argument("--dump", type=str, default=None)
    args = p.parse_args()

    setup_matplotlib()
    out = output_dir_for(SCRIPT_DIR)

    u4_expr, a_dict, D_expr, eps, x = galerkin_4_term()
    print("Symbolic four-term Galerkin solution:")
    for ai, val in a_dict.items():
        print(f"  {ai} = {val}")
    print(f"  common denominator D(eps) = {sp.expand(D_expr)}")

    u_exact_expr = 1 - sp.cosh(x / eps) / sp.cosh(1 / eps)
    u4_func = sp.lambdify((x, eps), u4_expr, "numpy")
    u_exact_func = sp.lambdify((x, eps), u_exact_expr, "numpy")

    eps_table = [sp.Rational(1, 20), sp.Rational(3, 40),
                 sp.Rational(1, 10), sp.Rational(3, 20),
                 sp.Rational(1, 5),  sp.Rational(1, 4),
                 sp.Rational(3, 10), sp.Rational(2, 5),
                 sp.Rational(1, 2),  sp.Rational(3, 4),
                 sp.Integer(1)]
    eps_floats = np.array([float(e) for e in eps_table])

    x_eval = np.linspace(-1.0, 1.0, 4001)
    err_table = []
    for e in eps_floats:
        u4_v = u4_func(x_eval, e)
        ue_v = u_exact_func(x_eval, e)
        err_table.append(float(np.max(np.abs(u4_v - ue_v))))
    err_arr = np.array(err_table)

    # Visual: solutions at eps = 1/10 and eps = 1/20
    e_show = [0.1, 0.05]
    sols = []
    exact = []
    for e in e_show:
        sols.append(u4_func(x_eval, e))
        exact.append(u_exact_func(x_eval, e))

    fig = plt.figure(figsize=(11.5, 4.4))
    gs = fig.add_gridspec(1, 2, wspace=0.30)

    # Panel A: solutions
    ax = fig.add_subplot(gs[0, 0])
    for i, (e, u4_v, ue_v) in enumerate(zip(e_show, sols, exact)):
        ax.plot(x_eval, ue_v, color=NAVY,  linewidth=1.2,
                linestyle=("-" if i == 0 else "--"),
                label=fr"$u_\text{{exact}},\ \epsilon = 1/{int(round(1/e))}$")
        ax.plot(x_eval[::80], u4_v[::80], "o", color=CORAL,
                markersize=5, markerfacecolor="white",
                label=fr"$u_4,\ \epsilon = 1/{int(round(1/e))}$" if i == 0 else None)
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$u(x)$")
    ax.set_title(r"$\epsilon^2 u'' - u = -1$, $u(\pm 1) = 0$:"
                 r" four-term Galerkin vs exact",
                 fontsize=10)
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(True, alpha=0.25, linewidth=0.4)

    # Panel B: max error vs eps
    ax = fig.add_subplot(gs[0, 1])
    ax.loglog(eps_floats, err_arr, "o-", color=NAVY, markerfacecolor="white",
              markersize=6, linewidth=1.0,
              label=r"max $|u_4 - u_\mathrm{exact}|$")
    ax.set_xlabel(r"$\epsilon$")
    ax.set_ylabel(r"$L^\infty$ error")
    ax.set_title("Boyd Table 20.3: four-term symbolic accuracy "
                 "across two decades of $\\epsilon$", fontsize=10)
    ax.legend(loc="upper left", fontsize=10)
    ax.grid(True, alpha=0.25, linewidth=0.4)

    save_fig(fig, out, "symbolic_boundary_layer")
    plt.close(fig)

    print(f"\n[Etude 21.9]  symbolic boundary-layer BVP, four-term Galerkin")
    print(f"   epsilon       L_inf error")
    for e, err in zip(eps_floats, err_arr):
        print(f"   {e:8.5f}    {err:.5e}")
    print(f"  figure: {out / 'symbolic_boundary_layer.pdf'}")

    dump_json(args.dump, {
        "epsilon": eps_floats.tolist(),
        "Linf_error": err_arr.tolist(),
        "a0_expr": str(a_dict[list(a_dict.keys())[0]]),
        "denominator_D": str(sp.expand(D_expr)),
    })


if __name__ == "__main__":
    main()
