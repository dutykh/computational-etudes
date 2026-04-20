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

function algebraic_semi_infinite(L)
    Map1D(x -> L*(1+x)/(1-x), y -> (y-L)/(y+L),
          x -> 2L/(1-x)^2, x -> 4L/(1-x)^3)
end
function tanh_map()
    Map1D(x -> tanh(x), y -> atanh(y),
          x -> 1 / cosh(x)^2, x -> -2*tanh(x)/cosh(x)^2)
end

function convergence_alg(N, L)
    Dx, x = cheb_matrix(N)
    m = algebraic_semi_infinite(L)
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
    fig = Figure(size=(940, 340))

    ax1 = Axis(fig[1, 1], xlabel="coordinate value", title="(a) Grid clustering, N=24",
               limits=((-1.2, 12.2), (-1.5, 0.6)))
    _, x_demo = cheb_matrix(24)
    y_alg = [min(2.0*(1+xi)/(1-xi), 12.0) for xi in x_demo]
    y_tanh = tanh.(x_demo)
    for xg in x_demo;   vlines!(ax1, xg; color=NAVY, linewidth=1, ymin=0.5, ymax=0.6); end
    for xg in y_alg;    vlines!(ax1, xg; color=CORAL, linewidth=1, ymin=0.30, ymax=0.40); end
    for xg in y_tanh;   vlines!(ax1, xg; color=TEAL, linewidth=1, ymin=0.10, ymax=0.20); end

    ax2 = Axis(fig[1, 2], xlabel="N", ylabel="mapped-deriv error",
               yscale=log10, title="(b) Validation")
    scatterlines!(ax2, Ns, alg1 .+ 1e-18; color=CORAL, label="alg D1")
    scatterlines!(ax2, Ns, alg2 .+ 1e-18; color=CORAL, linestyle=:dash, label="alg D2")
    scatterlines!(ax2, Ns, tnh1 .+ 1e-18; color=TEAL, label="tanh D1")
    scatterlines!(ax2, Ns, tnh2 .+ 1e-18; color=TEAL, linestyle=:dash, label="tanh D2")
    axislegend(ax2; position=:rt)

    save(joinpath(outdir, "map1d_toolkit.pdf"), fig)
    save(joinpath(outdir, "map1d_toolkit.png"), fig)
    @printf("[19.3-julia] saved figure\n")
end

if abspath(PROGRAM_FILE) == @__FILE__
    run()
end
