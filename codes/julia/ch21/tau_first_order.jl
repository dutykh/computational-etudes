# tau_first_order.jl
# Chapter 21, Etude 21.11: tau method on u' + u = 0, u(-1) = 1.
#
# Original: u'(x) + u(x) = 0, exact solution u(x) = exp(-(x+1)).
# Modified: v'(x) + v(x) = tau T_N(x), v(-1) = 1.
# Solve the (N+2) x (N+2) augmented system [D + I | -T_N(x)] [v; tau] = b.
#
# Author: Dr. Denys Dutykh (Khalifa University of Science and Technology)

include(joinpath(@__DIR__, "tricks_common.jl"))
include(joinpath(@__DIR__, "..", "ch07", "cheb_matrix.jl"))
apply_theme!()
const OUT = output_dir(@__DIR__)


function tau_solve(N::Int)
    D, x = cheb_matrix(N)
    L = D + I                              # L v = 0  (original RHS)
    TN = cos.(N .* acos.(x))                # T_N at CGL points
    A = zeros(N + 2, N + 2)
    b = zeros(N + 2)
    A[1:N+1, 1:N+1] = L
    A[1:N+1, N+2] = -TN
    A[N+2, :] .= 0.0
    A[N+2, N+1] = 1.0    # v at x = -1 (last CGL node)
    b[N+2] = 1.0
    sol = A \ b
    return x, sol[1:N+1], sol[N+2]
end


function standard_pseudospectral(N::Int)
    D, x = cheb_matrix(N)
    L = D + I
    L[N+1, :] .= 0.0
    L[N+1, N+1] = 1.0
    rhs = zeros(N + 1)
    rhs[N+1] = 1.0
    return x, L \ rhs
end


function make_figure(; dump_path = nothing)
    Ns = [4, 6, 8, 12, 16, 24, 32]
    tau_vals = Float64[]
    err_tau  = Float64[]
    err_pseu = Float64[]
    for N in Ns
        x, v, tau = tau_solve(N)
        push!(tau_vals, abs(tau))
        u_exact = exp.(-(x .+ 1.0))
        push!(err_tau, maximum(abs.(v .- u_exact)))
        x2, u_ps = standard_pseudospectral(N)
        push!(err_pseu, maximum(abs.(u_ps .- exp.(-(x2 .+ 1.0)))))
    end

    N_show = 16
    x_show, v_show, tau_show = tau_solve(N_show)
    x_dense = collect(range(-1, 1; length = 401))
    u_exact_dense = exp.(-(x_dense .+ 1.0))

    fig = Figure(size = (1150, 440))

    ax1 = Axis(fig[1, 1]; xlabel = "N", ylabel = "|tau| and pointwise error",
               yscale = log10,
               title = "geometric decay of tau, geometric convergence",
               limits = ((3, 33), (1e-16, 5.0)))
    scatterlines!(ax1, Ns, tau_vals; color = NAVY, markercolor = :white,
        strokecolor = NAVY, strokewidth = 1.0, markersize = 6, linewidth = 1.0,
        label = "|tau| (Lanczos perturbation)")
    scatterlines!(ax1, Ns, err_tau; color = TEAL, marker = :rect,
        markercolor = :white, strokecolor = TEAL, strokewidth = 1.0,
        markersize = 5, linewidth = 0.9, label = "max |v_N - u_exact|")
    scatterlines!(ax1, Ns, err_pseu; color = CORAL, marker = :utriangle,
        markercolor = :white, strokecolor = CORAL, strokewidth = 1.0,
        markersize = 5, linewidth = 0.9, label = "max |u_N^p-s - u_exact|")
    axislegend(ax1, position = :rt, labelsize = 9)

    ax2 = Axis(fig[1, 2]; xlabel = "x", ylabel = "u(x)",
               title = "N = $N_show:  tau ~ $(round(tau_show; sigdigits=3))")
    lines!(ax2, x_dense, u_exact_dense; color = NAVY, linewidth = 1.4,
           label = "u_exact = exp(-(x+1))")
    scatter!(ax2, x_show, v_show; color = :white, strokecolor = TEAL,
             strokewidth = 1.0, markersize = 6,
             label = "v_N at CGL nodes")
    axislegend(ax2, position = :rt, labelsize = 9)

    save(joinpath(OUT, "tau_first_order.pdf"), fig)
    save(joinpath(OUT, "tau_first_order.png"), fig, px_per_unit = 4)

    println("[Etude 21.11]  tau method on u' + u = 0, u(-1) = 1")
    for (k, N) in enumerate(Ns)
        println("  N = $N:  |tau| = $(tau_vals[k]), err_tau = $(err_tau[k]), "
                * "err_p-s = $(err_pseu[k])")
    end
    println("  figure: ", joinpath(OUT, "tau_first_order.pdf"))

    if dump_path !== nothing
        write_json(dump_path, Dict("Ns" => Ns,
            "tau" => tau_vals, "err_tau" => err_tau,
            "err_pseudospectral" => err_pseu))
    end
end


using LinearAlgebra

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
    make_figure(; dump_path = parse_args())
end
