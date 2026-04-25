# map1d_toolkit.jl
# Chapter 19: Coordinate Transformations
# Computational Etude 19.3: Build a reusable one-dimensional map toolkit.
#
# Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)

using CairoMakie
using LinearAlgebra
using Printf

include(joinpath(@__DIR__, "..", "ch07", "cheb_matrix.jl"))

struct Map1D
    forward::Function
    inverse::Function
    fprime::Function
    fdoubleprime::Function
end

function derivative_matrices(m::Map1D, Dx, x)
    fp  = m.fprime.(x)
    fpp = m.fdoubleprime.(x)
    Dx2 = Dx * Dx
    Dy  = Diagonal(1.0 ./ fp) * Dx
    Dy2 = Diagonal(1.0 ./ fp .^ 2) * Dx2 - Diagonal(fpp ./ fp .^ 3) * Dx
    return Dy, Dy2
end

function algebraic_semi_infinite(ell)
    Map1D(x -> ell*(1+x)/(1-x), y -> (y-ell)/(y+ell),
          x -> 2ell/(1-x)^2, x -> 4ell/(1-x)^3)
end
function tanh_map()
    Map1D(x -> tanh(x), y -> atanh(y),
          x -> 1 / cosh(x)^2, x -> -2*tanh(x)/cosh(x)^2)
end

function convergence_alg(N, ell)
    Dx, x = cheb_matrix(N)
    m = algebraic_semi_infinite(ell)
    Dy, Dy2 = derivative_matrices(m, Dx, x)
    y = m.forward.(x)
    u = exp.(-y); du_ex = -exp.(-y); d2u_ex = exp.(-y)
    mask = isfinite.(y) .& (abs.(y) .< 1e6)
    e1 = maximum(abs.((Dy * u)[mask] .- du_ex[mask]))
    e2 = maximum(abs.((Dy2 * u)[mask] .- d2u_ex[mask]))
    return e1, e2
end
function convergence_tanh(N)
    Dx, x = cheb_matrix(N)
    m = tanh_map()
    Dy, Dy2 = derivative_matrices(m, Dx, x)
    y = m.forward.(x)
    u = 1.0 ./ (1 .+ y .^ 2)
    du_ex = -2 .* y ./ (1 .+ y .^ 2) .^ 2
    d2u_ex = (6 .* y .^ 2 .- 2) ./ (1 .+ y .^ 2) .^ 3
    e1 = maximum(abs.(Dy * u .- du_ex))
    e2 = maximum(abs.(Dy2 * u .- d2u_ex))
    return e1, e2
end

function run()
    outdir = joinpath(@__DIR__, "..", "..", "..",
                      "textbook", "figures", "ch19", "julia")
    mkpath(outdir)

    Ns = [8, 12, 16, 20, 24, 32, 40, 48, 64]
    alg1 = zeros(length(Ns)); alg2 = zeros(length(Ns))
    tnh1 = zeros(length(Ns)); tnh2 = zeros(length(Ns))
    for (i, N) in enumerate(Ns)
        alg1[i], alg2[i] = convergence_alg(N, 2.0)
        tnh1[i], tnh2[i] = convergence_tanh(N)
    end

    NAVY = colorant"#142D6E"; CORAL = colorant"#E74C3C"; TEAL = colorant"#16A085"

    Ngrid = 24
    ELL = 2.0
    Y_LIM_ALG = 12.0
    _, x_demo = cheb_matrix(Ngrid)
    x_line = collect(range(-1.0, 1.0 - 1e-12, length=401))

    fig = Figure(size=(1100, 760))

    # (a) tanh map
    ax1 = Axis(fig[1, 1];
               xlabel="computational coordinate x",
               ylabel="physical y",
               title="(a) tanh map: x_j (bottom) -> y_j (right)",
               limits=((-1.20, 1.20), (-1.30, 1.20)))
    y_tanh_line = tanh.(x_line)
    y_tanh_grid = tanh.(x_demo)
    h_t = lines!(ax1, x_line, y_tanh_line; color=TEAL, linewidth=1.6,
                 label="y = tanh(x)")
    scatter!(ax1, x_demo, fill(-1.18, length(x_demo));
             color=NAVY, marker=:vline, markersize=12)
    scatter!(ax1, fill(1.10, length(y_tanh_grid)), y_tanh_grid;
             color=TEAL, marker=:hline, markersize=12)
    for j in 1:3:length(x_demo)
        lines!(ax1, [x_demo[j], x_demo[j], 1.10],
                    [-1.18, y_tanh_grid[j], y_tanh_grid[j]];
                    color=TEAL, linewidth=0.4, linestyle=:dot)
    end
    hlines!(ax1, [0.0]; color=(:grey, 0.4), linewidth=0.4)
    vlines!(ax1, [0.0]; color=(:grey, 0.4), linewidth=0.4)
    axislegend(ax1, [h_t], ["y = tanh(x)"];
               position=:lt, framevisible=false)

    # (b) algebraic semi-infinite map
    ax2 = Axis(fig[1, 2];
               xlabel="computational coordinate x",
               ylabel="physical y",
               title="(b) algebraic map (ell = 2): clusters near y = 0",
               limits=((-1.20, 1.20), (-1.6, Y_LIM_ALG + 1.0)))
    y_alg_line = ELL .* (1 .+ x_line) ./ (1 .- x_line)
    visible_line = y_alg_line .<= Y_LIM_ALG
    h_a = lines!(ax2, x_line[visible_line], y_alg_line[visible_line];
                 color=CORAL, linewidth=1.6, label="y = ell(1+x)/(1-x)")
    x_finite = x_demo[1:end-1]
    y_alg_grid = ELL .* (1 .+ x_finite) ./ (1 .- x_finite)
    visible_grid = y_alg_grid .<= Y_LIM_ALG
    scatter!(ax2, x_demo, fill(-0.7, length(x_demo));
             color=NAVY, marker=:vline, markersize=12)
    scatter!(ax2, fill(1.10, sum(visible_grid)), y_alg_grid[visible_grid];
             color=CORAL, marker=:hline, markersize=12)
    for j in 1:length(x_finite)
        if y_alg_grid[j] <= Y_LIM_ALG && j % 2 == 0
            lines!(ax2, [x_finite[j], x_finite[j], 1.10],
                        [-0.7, y_alg_grid[j], y_alg_grid[j]];
                        color=CORAL, linewidth=0.4, linestyle=:dot)
        end
    end
    n_off = sum(.!visible_grid) + 1   # +1 for x=1 endpoint at infinity
    text!(ax2, 0.10, Y_LIM_ALG - 1.6;
          text="y -> infinity ($n_off ticks beyond view)",
          color=CORAL, fontsize=10)
    hlines!(ax2, [0.0]; color=(:grey, 0.4), linewidth=0.4)
    vlines!(ax2, [0.0]; color=(:grey, 0.4), linewidth=0.4)
    axislegend(ax2, [h_a], ["y = ell(1+x)/(1-x)"];
               position=:lt, framevisible=false)

    # (c) First-derivative convergence
    ax3 = Axis(fig[2, 1];
               xlabel="N", ylabel="max error of u'",
               title="(c) First derivative: D_y u",
               yscale=log10)
    scatterlines!(ax3, Ns, alg1 .+ 1e-18;
                  color=CORAL, marker=:circle, label="algebraic")
    scatterlines!(ax3, Ns, tnh1 .+ 1e-18;
                  color=TEAL, marker=:rect, label="tanh")
    axislegend(ax3; position=:rt, framevisible=false)

    # (d) Second-derivative convergence -- the N^2-lag panel
    ax4 = Axis(fig[2, 2];
               xlabel="N", ylabel="max error of u''",
               title="(d) Second derivative: D_y^2 u (N^2-lag)",
               yscale=log10)
    scatterlines!(ax4, Ns, alg2 .+ 1e-18;
                  color=CORAL, marker=:circle, linestyle=:dash, label="algebraic")
    scatterlines!(ax4, Ns, tnh2 .+ 1e-18;
                  color=TEAL, marker=:rect, linestyle=:dash, label="tanh")
    axislegend(ax4; position=:rt, framevisible=false)

    save(joinpath(outdir, "map1d_toolkit.pdf"), fig)
    save(joinpath(outdir, "map1d_toolkit.png"), fig)
    @printf("[19.3-julia] saved figure\n")
end

if abspath(PROGRAM_FILE) == @__FILE__
    run()
end
