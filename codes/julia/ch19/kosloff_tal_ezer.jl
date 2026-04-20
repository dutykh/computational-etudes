# kosloff_tal_ezer.jl
# Chapter 19: Coordinate Transformations
# Computational Etude 19.8: Accuracy versus timestep.
#
# Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)

using CairoMakie
using FFTW
using LinearAlgebra
using Printf

include(joinpath(@__DIR__, "..", "ch07", "cheb_matrix.jl"))

function kte_grid(N, beta)
    D, xi = cheb_matrix(N)
    denom = asin(1 - beta)
    y = denom == 0 ? xi : asin.((1 - beta) .* xi) ./ denom
    return y, xi, D
end

function kte_derivative(N, beta)
    D, xi = cheb_matrix(N)
    denom = asin(1 - beta)
    y = asin.((1 - beta) .* xi) ./ denom
    fp = (1 - beta) ./ (sqrt.(1 .- (1 - beta)^2 .* xi .^ 2) .* denom)
    Dy = Diagonal(1 ./ fp) * D
    return y, xi, D, Dy
end

function cheb_coeff_last(v)
    N = length(v) - 1
    V = vcat(v, reverse(v[2:N]))
    A = real.(fft(V)) ./ N
    A[1] *= 0.5; A[N + 1] *= 0.5
    return A[N + 1]
end

function run()
    outdir = joinpath(@__DIR__, "..", "..", "..",
                      "textbook", "figures", "ch19", "julia")
    mkpath(outdir)

    Ns = collect(8:2:96)
    min_std = zeros(length(Ns)); min_opt = zeros(length(Ns)); min_cons = zeros(length(Ns))
    stiff_std = zeros(length(Ns)); stiff_opt = zeros(length(Ns)); stiff_cons = zeros(length(Ns))
    for (i, N) in enumerate(Ns)
        y, _, D = kte_grid(N, 0.0)
        min_std[i] = minimum(diff(sort(y)))
        stiff_std[i] = maximum(abs.(eigvals(D)))

        beta_opt = 1 - cos(1 / N)
        y, _, _, Dy = kte_derivative(N, beta_opt)
        min_opt[i] = minimum(diff(sort(y)))
        stiff_opt[i] = maximum(abs.(eigvals(Dy)))

        beta_cons = 1 - cos(0.5)
        y, _, _, Dy = kte_derivative(N, beta_cons)
        min_cons[i] = minimum(diff(sort(y)))
        stiff_cons[i] = maximum(abs.(eigvals(Dy)))
    end

    Ns_c = collect(3:2:63)
    aN_agg = zeros(length(Ns_c)); aN_cons = zeros(length(Ns_c))
    for (i, N) in enumerate(Ns_c)
        y_a, _, _ = kte_grid(N, 1 - cos(1 / N)); aN_agg[i]  = abs(cheb_coeff_last(y_a))
        y_c, _, _ = kte_grid(N, 1 - cos(0.5));  aN_cons[i] = abs(cheb_coeff_last(y_c))
    end

    NAVY = colorant"#142D6E"; CORAL = colorant"#E74C3C"; TEAL = colorant"#16A085"
    fig = Figure(size=(1340, 360))

    ax1 = Axis(fig[1, 1], xlabel="N", ylabel="min spacing",
               xscale=log10, yscale=log10, title="(a) Minimum grid spacing")
    scatterlines!(ax1, Ns, min_std; color=NAVY, label="standard")
    scatterlines!(ax1, Ns, min_opt; color=CORAL, label="KTE 1/N²")
    scatterlines!(ax1, Ns, min_cons; color=TEAL, label="KTE 1-cos(1/2)")
    lines!(ax1, Ns, 2.0 ./ Ns; color=:gray, linestyle=:dot, label="2/N")
    axislegend(ax1; position=:lb)

    ax2 = Axis(fig[1, 2], xlabel="N", ylabel="ρ(D)",
               xscale=log10, yscale=log10, title="(b) Stiffness")
    scatterlines!(ax2, Ns, stiff_std; color=NAVY, label="standard")
    scatterlines!(ax2, Ns, stiff_opt; color=CORAL, label="KTE 1/N²")
    scatterlines!(ax2, Ns, stiff_cons; color=TEAL, label="KTE 1-cos(1/2)")
    axislegend(ax2; position=:lt)

    ax3 = Axis(fig[1, 3], xlabel="N", ylabel="|a_N| of f(y)=y",
               xscale=log10, yscale=log10, title="(c) Accuracy destroyed")
    scatter!(ax3, Ns_c, max.(aN_agg, 1e-18); color=CORAL, label="aggressive")
    scatter!(ax3, Ns_c, max.(aN_cons, 1e-18); color=TEAL, label="conservative")
    lines!(ax3, Ns_c, 0.488 ./ (Ns_c .^ 2); color=NAVY, linestyle=:dot, label="0.488/N²")
    axislegend(ax3; position=:lb)

    save(joinpath(outdir, "kosloff_tal_ezer.pdf"), fig)
    save(joinpath(outdir, "kosloff_tal_ezer.png"), fig)
    @printf("[19.8-julia] saved figure\n")
end

if abspath(PROGRAM_FILE) == @__FILE__
    run()
end
