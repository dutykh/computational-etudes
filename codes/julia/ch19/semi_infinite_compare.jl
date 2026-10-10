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

    Y_LIM = 12.0
    Ngrid = 24
    _, x_demo = cheb_matrix(Ngrid)
    x_line = collect(range(-1.0, 1.0 - 1e-12, length=401))

    fig = Figure(size=(1100, 760))

    # (a) Solution at N=24
    ax1 = Axis(fig[1, 1]; xlabel="y", ylabel="u(y)",
               title="(a) Solution at N = 24")
    y_t, u_t = solve_truncation(24, 20); y_m, u_m = solve_algebraic(24, 2)
    yplot = collect(range(0, 12, length=401))
    h_ex = lines!(ax1, yplot, exp.(-yplot); color=NAVY, linewidth=1.4)
    h_tr = scatter!(ax1, y_t, u_t; color=CORAL, marker=:circle, markersize=7)
    y_m_clip = clamp.(y_m, 0, 12)
    h_al = scatter!(ax1, y_m_clip, u_m; color=TEAL, marker=:rect, markersize=7)
    axislegend(ax1, [h_ex, h_tr, h_al],
               ["exact", "truncation L = 20", "algebraic ell = 2"];
               position=:rt, framevisible=false)

    # (b) NEW: algebraic map shape y(x) for the four ell values
    colours_m = [CORAL, TEAL, PURPLE, NAVY]
    ax2 = Axis(fig[1, 2];
               xlabel="computational coordinate x",
               ylabel="physical y",
               title="(b) Algebraic map  y = ell (1+x)/(1-x)",
               limits=((-1.10, 1.10), (-1.6, Y_LIM + 1.0)))
    legendH = []
    for (j, ell) in enumerate(ell_map)
        y_curve = ell .* (1 .+ x_line) ./ (1 .- x_line)
        visible = y_curve .<= Y_LIM
        push!(legendH,
              lines!(ax2, x_line[visible], y_curve[visible];
                     color=colours_m[j], linewidth=1.6,
                     label="ell = $ell"))
    end
    scatter!(ax2, x_demo, fill(-0.7, length(x_demo));
             color=NAVY, marker=:vline, markersize=12)
    text!(ax2, -1.0, -1.1; text="CGL nodes x_j",
          color=NAVY, fontsize=10)
    hlines!(ax2, [0.0]; color=(:grey, 0.4), linewidth=0.4)
    vlines!(ax2, [0.0]; color=(:grey, 0.4), linewidth=0.4)
    axislegend(ax2, legendH,
               ["ell = $e" for e in ell_map];
               position=:lt, framevisible=false)

    # (c) Truncation: error vs N
    ax3 = Axis(fig[2, 1]; xlabel="N", ylabel="max error", yscale=log10,
               title="(c) Truncation: error vs N")
    colours_t = [CORAL, ORANGE, GOLD]
    for (j, (L, c)) in enumerate(zip(L_trunc, colours_t))
        scatterlines!(ax3, Ns, E_trunc[:, j]; color=c, label="trunc. L = $L")
    end
    axislegend(ax3; position=:rt, framevisible=false)

    # (d) Algebraic map: error vs N
    ax4 = Axis(fig[2, 2]; xlabel="N", ylabel="max error", yscale=log10,
               title="(d) Algebraic map: error vs N")
    for (j, (ell, c)) in enumerate(zip(ell_map, colours_m))
        scatterlines!(ax4, Ns, max.(E_map[:, j], 1e-18);
                      color=c, label="algebraic ell = $ell")
    end
    axislegend(ax4; position=:rt, framevisible=false)

    save(joinpath(outdir, "semi_infinite_compare.pdf"), fig)
    save(joinpath(outdir, "semi_infinite_compare.png"), fig)
    @printf("[19.4-julia] saved figure\n")
end

if abspath(PROGRAM_FILE) == @__FILE__
    run()
end
