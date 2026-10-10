#!/usr/bin/env julia
#=
chebyshev_stiffness.jl
Chapter 17: Stability of Spectral Methods

Computational Etude 17.4: Chebyshev Stiffness and Eigenvalue Structure.

Three figures:
  (1) cheb_eigenvalues.pdf -- Eigenvalues of Chebyshev D^2 (Dirichlet BCs)
      for N=32 plotted on the real axis.
  (2) cheb_eigenvectors.pdf -- Physical eigenvector (smooth) and outlier
      eigenvector (boundary-localised).
  (3) cheb_cfl_scaling.pdf -- Log-log plot of spectral radius rho(D2) vs N,
      showing N^4 scaling.

Generates Figures 17.4a,b,c

Author: Dr. Denys Dutykh
        Mathematics Department
        Khalifa University of Science and Technology
        Abu Dhabi, UAE

Part of the book "Computational Etudes: A Spectral Approach"
https://github.com/dutykh/computational-etudes
=#

using LinearAlgebra
using Plots
using LaTeXStrings
using Printf

gr()

const NAVY  = RGB(20/255, 45/255, 110/255)
const SKY   = RGB(120/255, 150/255, 210/255)
const CORAL = RGB(231/255, 76/255, 60/255)
const TEAL  = RGB(22/255, 160/255, 133/255)
const AMBER = RGB(230/255, 126/255, 34/255)

const SCRIPT_DIR = @__DIR__
const OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch17", "julia")
mkpath(OUTPUT_DIR)

"""Chebyshev differentiation matrix of size (N+1) x (N+1)."""
function cheb(N)
    if N == 0; return [0.0;;], [1.0]; end
    x = [cos(pi * j / N) for j in 0:N]
    c = [j == 0 || j == N ? 2.0 : 1.0 for j in 0:N] .* [(-1.0)^j for j in 0:N]
    X = repeat(x, 1, N+1)
    dX = X - X' + I(N+1)
    D = (c * (1.0 ./ c)') ./ dX
    D .-= Diagonal(vec(sum(D, dims=2)))
    return D, x
end

function figure_eigenvalues()
    N = 32
    D, x = cheb(N)
    D2 = D * D
    D2_int = D2[2:N, 2:N]

    eigs = sort(real.(eigvals(D2_int)))

    # Physical eigenvalue approximations: -(n*pi/2)^2 for domain [-1,1]
    n_phys = 1:N-1
    phys_approx = -(n_phys .* pi ./ 2.0) .^ 2

    # Separate physical from outliers
    threshold = -(N * pi / 2.0)^2 * 0.5
    phys_mask = eigs .> threshold
    out_mask  = .!phys_mask

    p = plot(xlabel=L"\lambda\ \mathrm{(eigenvalue)}", ylabel="",
             title=L"Eigenvalues of Chebyshev $D^2$ with Dirichlet BCs ($N = 32$)",
             yticks=[], legend=:topleft, legendfontsize=9,
             grid=true, gridalpha=0.3, size=(1000, 350))
    hline!(p, [0], color=:gray, lw=0.5, label="")

    # Physical modes
    scatter!(p, eigs[phys_mask], zeros(sum(phys_mask)),
             color=NAVY, ms=6, markerstrokecolor=:white, markerstrokewidth=0.5,
             label="physical modes")
    # Outlier modes
    if any(out_mask)
        scatter!(p, eigs[out_mask], zeros(sum(out_mask)),
                 color=CORAL, ms=6, markershape=:square,
                 markerstrokecolor=:white, markerstrokewidth=0.5,
                 label="outlier modes")
    end

    # Mark physical approximations
    nmark = min(7, N-1)
    scatter!(p, phys_approx[1:nmark], fill(0.05, nmark),
             color=TEAL, ms=4, markershape=:dtriangle, alpha=0.7,
             label=L"-n^2\pi^2/4\ \mathrm{(predicted)}")

    savefig(p, joinpath(OUTPUT_DIR, "cheb_eigenvalues.pdf"))
    savefig(p, joinpath(OUTPUT_DIR, "cheb_eigenvalues.png"))
    println("Figure saved to $(joinpath(OUTPUT_DIR, "cheb_eigenvalues.pdf"))")
    return eigs
end

function figure_eigenvectors()
    N = 32
    D, x = cheb(N)
    D2 = D * D
    D2_int = D2[2:N, 2:N]
    x_int = x[2:N]

    F = eigen(D2_int)
    idx = sortperm(real.(F.values))
    eigs_sorted = real.(F.values[idx])
    vecs_sorted = F.vectors[:, idx]

    # Physical eigenvector: mode N/4 = 8 (0-indexed mode 7 in sorted order)
    mode_idx = N ÷ 4   # 1-based index for mode 8
    v_phys = real.(vecs_sorted[:, mode_idx])
    v_phys ./= maximum(abs.(v_phys))
    if v_phys[length(v_phys) ÷ 4] < 0
        v_phys .*= -1
    end

    # Full vector including BCs
    v_full = zeros(N + 1)
    v_full[2:N] .= v_phys

    p1 = plot(x, v_full, color=NAVY, lw=1.8, label="",
              xlabel=L"$x$", ylabel=L"$v(x)$",
              title="(a) Physical eigenvector (mode $(N÷4), " *
                    latexstring("\\lambda = $(round(eigs_sorted[mode_idx], digits=1))") * ")",
              xlims=(-1.05, 1.05), grid=true, gridalpha=0.3)
    scatter!(p1, x_int, v_phys, color=NAVY, ms=3,
             markerfacecolor=:white, markerstrokewidth=1.0, label="")

    # Outlier eigenvector: most negative eigenvalue (index 1)
    v_out = real.(vecs_sorted[:, 1])
    v_out ./= maximum(abs.(v_out))
    v_full_out = zeros(N + 1)
    v_full_out[2:N] .= v_out

    p2 = plot(x, v_full_out, color=CORAL, lw=1.8, label="",
              xlabel=L"$x$", ylabel=L"$v(x)$",
              title="(b) Outlier eigenvector (" *
                    latexstring("\\lambda = $(round(Int, eigs_sorted[1]))") * ")",
              xlims=(-1.05, 1.05), grid=true, gridalpha=0.3)
    scatter!(p2, x_int, v_out, color=CORAL, ms=3,
             markerfacecolor=:white, markerstrokewidth=1.0, label="")

    fig = plot(p1, p2, layout=(1, 2), size=(1100, 450))

    savefig(fig, joinpath(OUTPUT_DIR, "cheb_eigenvectors.pdf"))
    savefig(fig, joinpath(OUTPUT_DIR, "cheb_eigenvectors.png"))
    println("Figure saved to $(joinpath(OUTPUT_DIR, "cheb_eigenvectors.pdf"))")
end

function figure_cfl_scaling()
    N_values = [8, 12, 16, 20, 24, 32, 48, 64]
    rho_values = Float64[]

    for N in N_values
        D, _ = cheb(N)
        D2 = D * D
        D2_int = D2[2:N, 2:N]
        eigs = eigvals(D2_int)
        rho = maximum(abs.(eigs))
        push!(rho_values, rho)
        @printf("N = %3d, rho(D2) = %.2e, rho/N^4 = %.4f\n", N, rho, rho / N^4)
    end

    N_arr = Float64.(N_values)

    p = plot(N_arr, rho_values, xscale=:log10, yscale=:log10,
             marker=:circle, ms=6, lw=1.5, color=NAVY,
             label=L"\rho(\tilde{D}_2)",
             xlabel=L"$N$", ylabel=L"\rho(\tilde{D}_2)",
             title=L"Spectral radius of Chebyshev $D^2$ (Dirichlet BCs)",
             legend=:topleft, legendfontsize=10,
             grid=true, gridalpha=0.3, minorgrid=true, size=(700, 500))

    # Reference line: 0.048 * N^4
    N_ref = range(6, 80, length=100)
    plot!(p, N_ref, 0.048 .* N_ref .^ 4, ls=:dash, lw=1.2, color=CORAL,
          label=L"0.048\, N^4")

    # Slope annotation
    annotate!(p, 30, 1e5, text("slope = 4", 10, CORAL, rotation=40))

    savefig(p, joinpath(OUTPUT_DIR, "cheb_cfl_scaling.pdf"))
    savefig(p, joinpath(OUTPUT_DIR, "cheb_cfl_scaling.png"))
    println("Figure saved to $(joinpath(OUTPUT_DIR, "cheb_cfl_scaling.pdf"))")
end

function main()
    figure_eigenvalues()
    figure_eigenvectors()
    figure_cfl_scaling()
end

main()
