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
    fig = Figure(size=(1100, 760))

    # (a) NEW: KTE map shape y(xi) for the three regimes
    Ndemo = 32
    _, xi_demo = cheb_matrix(Ndemo)
    x_line = collect(range(-1, 1, length=401))
    ax1 = Axis(fig[1, 1];
               xlabel="computational coordinate xi",
               ylabel="physical coordinate y",
               title="(a) KTE map  y = arcsin((1-beta) xi)/arcsin(1-beta)",
               limits=((-1.10, 1.10), (-1.40, 1.10)),
               aspect=1)
    h_std = lines!(ax1, x_line, x_line; color=NAVY, linewidth=1.4)
    beta_a = 1 - cos(1 / Ndemo)
    y_a = asin.((1 - beta_a) .* x_line) ./ asin(1 - beta_a)
    h_agg = lines!(ax1, x_line, y_a; color=CORAL, linewidth=1.6)
    beta_c = 1 - cos(0.5)
    y_c = asin.((1 - beta_c) .* x_line) ./ asin(1 - beta_c)
    h_con = lines!(ax1, x_line, y_c; color=TEAL, linewidth=1.6)
    scatter!(ax1, xi_demo, fill(-1.12, length(xi_demo));
             color=NAVY, marker=:vline, markersize=10)
    text!(ax1, -1.0, -1.30; text="CGL nodes xi_j",
          color=NAVY, fontsize=10)
    hlines!(ax1, [0.0]; color=(:grey, 0.4), linewidth=0.4)
    vlines!(ax1, [0.0]; color=(:grey, 0.4), linewidth=0.4)
    axislegend(ax1, [h_std, h_agg, h_con],
               ["standard (beta = 0): identity",
                "aggressive beta ~ 1/N^2 (N = $Ndemo)",
                "conservative beta = 1 - cos(1/2)"];
               position=:lt, framevisible=false)

    # (b) Min spacing vs N
    ax2 = Axis(fig[1, 2]; xlabel="N", ylabel="min spacing",
               xscale=log10, yscale=log10,
               title="(b) Minimum spacing vs N")
    scatterlines!(ax2, Ns, min_std; color=NAVY, label="standard")
    scatterlines!(ax2, Ns, min_opt; color=CORAL, label="KTE 1/N^2")
    scatterlines!(ax2, Ns, min_cons; color=TEAL, label="KTE 1-cos(1/2)")
    lines!(ax2, Ns, 2.0 ./ Ns; color=:gray, linestyle=:dot, label="2/N")
    axislegend(ax2; position=:lb, framevisible=false)

    # (c) Stiffness vs N
    ax3 = Axis(fig[2, 1]; xlabel="N", ylabel="rho(D)",
               xscale=log10, yscale=log10, title="(c) Stiffness rho(D) vs N")
    scatterlines!(ax3, Ns, stiff_std; color=NAVY, label="standard")
    scatterlines!(ax3, Ns, stiff_opt; color=CORAL, label="KTE 1/N^2")
    scatterlines!(ax3, Ns, stiff_cons; color=TEAL, label="KTE 1-cos(1/2)")
    axislegend(ax3; position=:lt, framevisible=false)

    # (d) Coefficient |a_N| of f(y) = y
    ax4 = Axis(fig[2, 2]; xlabel="N", ylabel="|a_N| of f(y) = y",
               xscale=log10, yscale=log10,
               title="(d) Accuracy destroyed by aggressive scaling")
    scatter!(ax4, Ns_c, max.(aN_agg, 1e-18); color=CORAL,
             marker=:circle, label="aggressive")
    scatter!(ax4, Ns_c, max.(aN_cons, 1e-18); color=TEAL,
             marker=:rect, label="conservative")
    lines!(ax4, Ns_c, 0.488 ./ (Ns_c .^ 2); color=NAVY,
           linestyle=:dot, label="0.488/N^2")
    axislegend(ax4; position=:lb, framevisible=false)

    save(joinpath(outdir, "kosloff_tal_ezer.pdf"), fig)
    save(joinpath(outdir, "kosloff_tal_ezer.png"), fig)
    @printf("[19.8-julia] saved figure\n")
end

if abspath(PROGRAM_FILE) == @__FILE__
    run()
end
