# semi_infinite_compare.jl
# Chapter 19: Coordinate Transformations
# Computational Etude 19.4: Infinity without truncation.
#
# Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)

using CairoMakie
using LinearAlgebra
using Printf

include(joinpath(@__DIR__, "..", "ch07", "cheb_matrix.jl"))

function solve_truncation(N, L)
    D, xi = cheb_matrix(N)
    y = 0.5 * L .* (xi .+ 1)
    Dy = (2 / L) .* D
    Dy2 = Dy * Dy
    A = Dy2[2:N, 2:N] - I(N - 1)
    rhs = -Dy2[2:N, 1] .* 0.0 .- Dy2[2:N, N + 1] .* 1.0
    u_int = A \ rhs
    u = zeros(N + 1); u[1] = 0.0; u[N + 1] = 1.0; u[2:N] = u_int
    return y, u
end

function solve_algebraic(N, ell)
    D, x = cheb_matrix(N)
    fp  = @. 2ell / (1 - x)^2
    fpp = @. 4ell / (1 - x)^3
    Dy = Diagonal(1 ./ fp) * D
    Dy2 = Diagonal(1 ./ fp .^ 2) * (D * D) - Diagonal(fpp ./ fp .^ 3) * D
    y = @. ell * (1 + x) / (1 - x)
    A = Dy2[2:N, 2:N] - I(N - 1)
    rhs = -Dy2[2:N, 1] .* 0.0 .- Dy2[2:N, N + 1] .* 1.0
    u_int = A \ rhs
    u = zeros(N + 1); u[1] = 0.0; u[N + 1] = 1.0; u[2:N] = u_int
    return y, u
end

function max_err(y, u)
    mask = isfinite.(y) .& (y .< 1e6)
    return maximum(abs.(u[mask] .- exp.(-y[mask])))
end

function run()
    outdir = joinpath(@__DIR__, "..", "..", "..",
                      "textbook", "figures", "ch19", "julia")
    mkpath(outdir)
    Ns = [12, 16, 20, 24, 32, 40, 48, 64]
    L_trunc = [10, 20, 40]   # truncation half-widths
    ell_map = [1, 2, 4, 8]   # map parameters
    E_trunc = zeros(length(Ns), length(L_trunc))
    E_map   = zeros(length(Ns), length(ell_map))
    for (i, N) in enumerate(Ns)
        for (j, L) in enumerate(L_trunc)
            y, u = solve_truncation(N, L); E_trunc[i, j] = maximum(abs.(u .- exp.(-y)))
        end
        for (j, ell) in enumerate(ell_map)
            y, u = solve_algebraic(N, ell); E_map[i, j] = max_err(y, u)
        end
    end

    NAVY = colorant"#142D6E"; CORAL = colorant"#E74C3C"; TEAL = colorant"#16A085"
    ORANGE = colorant"#E67E22"; GOLD = colorant"#D4A017"; PURPLE = colorant"#8E44AD"
    fig = Figure(size=(1300, 340))

    ax1 = Axis(fig[1, 1], xlabel="y", ylabel="u(y)", title="(a) Solution at N=24")
    y_t, u_t = solve_truncation(24, 20); y_m, u_m = solve_algebraic(24, 2)
    yplot = collect(range(0, 12, length=401))
    lines!(ax1, yplot, exp.(-yplot); color=NAVY, linewidth=1.2)
    scatter!(ax1, y_t, u_t; color=CORAL, markersize=5)
    y_m_clip = clamp.(y_m, 0, 12)
    scatter!(ax1, y_m_clip, u_m; color=TEAL, markersize=5)

    ax2 = Axis(fig[1, 2], xlabel="N", ylabel="max error", yscale=log10,
               title="(b) Truncation")
    colours_t = [CORAL, ORANGE, GOLD]
    for (j, (L, c)) in enumerate(zip(L_trunc, colours_t))
        scatterlines!(ax2, Ns, E_trunc[:, j]; color=c, label="L=$L")
    end
    axislegend(ax2; position=:rt)

    ax3 = Axis(fig[1, 3], xlabel="N", ylabel="max error", yscale=log10,
               title="(c) Algebraic map")
    colours_m = [CORAL, TEAL, PURPLE, NAVY]
    for (j, (ell, c)) in enumerate(zip(ell_map, colours_m))
        scatterlines!(ax3, Ns, max.(E_map[:, j], 1e-18); color=c, label="ℓ=$ell")
    end
    axislegend(ax3; position=:rt)

    save(joinpath(outdir, "semi_infinite_compare.pdf"), fig)
    save(joinpath(outdir, "semi_infinite_compare.png"), fig)
    @printf("[19.4-julia] saved figure\n")
end

if abspath(PROGRAM_FILE) == @__FILE__
    run()
end
