# symbolic_boundary_layer.jl
# Chapter 21, Etude 21.9: Carrier-Pearson boundary-layer BVP solved
# symbolically by a four-term Galerkin method.  Reproduces Boyd Table 20.3.
#
# Trial: u(x) = (1 - x^2)(a_0 + a_2 x^2 + a_4 x^4 + a_6 x^6).
# Test functions: {1, x^2, x^4, x^6} (unweighted moments).
# Solve the 4-by-4 symbolic linear system; coefficients come out as
# rational functions of epsilon, agreeing with Boyd Eq 20.8.
#
# Author: Dr. Denys Dutykh (Khalifa University of Science and Technology)

using SymPy
include(joinpath(@__DIR__, "tricks_common.jl"))
apply_theme!()
const OUT = output_dir(@__DIR__)


function galerkin_4_term()
    @syms x::real epsilon::real
    a = [SymPy.symbols("a_$(2k)") for k in 0:3]
    p = sum(a[k+1] * x ^ (2k) for k in 0:3)
    u = (1 - x^2) * p

    R = epsilon^2 * diff(u, x, 2) - u + 1
    eqs = [integrate(x^(2j) * R, (x, -1, 1)) for j in 0:3]
    sol = solve(eqs, a)             # Dict-like

    return u, sol, x, epsilon, a
end


function main_run(dump_path)
    u_expr, sol_dict, x, epsilon, a_syms = galerkin_4_term()
    println("Symbolic four-term Galerkin solution:")
    for (k, ak) in enumerate(a_syms)
        val = simplify(sol_dict[ak])
        println("  $ak = $val")
    end

    u4 = simplify(u_expr.subs([(ak, sol_dict[ak]) for ak in a_syms]))
    u_exact = 1 - cosh(x / epsilon) / cosh(1 / epsilon)
    u4_func = lambdify(u4, [x, epsilon])
    ue_func = lambdify(u_exact, [x, epsilon])

    eps_table = [1//20, 3//40, 1//10, 3//20, 1//5,
                 1//4, 3//10, 2//5, 1//2, 3//4, 1]
    eps_floats = Float64.(eps_table)

    x_eval = collect(range(-1.0, 1.0; length = 4001))
    err_table = Float64[]
    for e in eps_floats
        u4_v = [u4_func(xx, e) for xx in x_eval]
        ue_v = [ue_func(xx, e) for xx in x_eval]
        push!(err_table, maximum(abs.(u4_v .- ue_v)))
    end

    e_show = [0.1, 0.05]
    sols = [[u4_func(xx, e) for xx in x_eval] for e in e_show]
    exact = [[ue_func(xx, e) for xx in x_eval] for e in e_show]

    fig = Figure(size = (1150, 440))

    ax1 = Axis(fig[1, 1]; xlabel = "x", ylabel = "u(x)",
               title = "epsilon^2 u'' - u = -1, u(+/-1) = 0:  4-term Galerkin vs exact")
    for i in 1:length(e_show)
        e = e_show[i]
        lines!(ax1, x_eval, exact[i]; color = NAVY, linewidth = 1.2,
               linestyle = (i == 1 ? :solid : :dash),
               label = "u_exact, eps = 1/$(round(Int, 1/e))")
        scatter!(ax1, x_eval[1:80:end], sols[i][1:80:end];
                 color = :white, strokecolor = CORAL, strokewidth = 1.0,
                 markersize = 5,
                 label = (i == 1 ? "u_4 (4-term Galerkin)" : nothing))
    end
    axislegend(ax1, position = :rt, labelsize = 9)

    ax2 = Axis(fig[1, 2]; xlabel = "epsilon", ylabel = "L_inf error",
               yscale = log10, xscale = log10,
               title = "Boyd Table 20.3: 4-term symbolic accuracy across epsilon")
    scatterlines!(ax2, eps_floats, err_table; color = NAVY,
        markercolor = :white, strokecolor = NAVY, strokewidth = 1.0,
        markersize = 6, linewidth = 1.0,
        label = "max |u_4 - u_exact|")
    axislegend(ax2, position = :lt, labelsize = 10)

    save(joinpath(OUT, "symbolic_boundary_layer.pdf"), fig)
    save(joinpath(OUT, "symbolic_boundary_layer.png"), fig, px_per_unit = 4)

    println("\n[Etude 21.9]  symbolic boundary-layer BVP, 4-term Galerkin")
    for (e, err) in zip(eps_floats, err_table)
        println("   eps = $e   L_inf err = $err")
    end
    println("  figure: ", joinpath(OUT, "symbolic_boundary_layer.pdf"))

    if dump_path !== nothing
        write_json(dump_path, Dict("epsilon" => eps_floats,
            "Linf_error" => err_table))
    end
end


function parse_args()
    dp = nothing; i = 1
    while i <= length(ARGS)
        if ARGS[i] == "--dump" && i < length(ARGS)
            dp = ARGS[i+1]; i += 2
        else; i += 1; end
    end
    return dp
end


if abspath(PROGRAM_FILE) == @__FILE__
    main_run(parse_args())
end
