#!/usr/bin/env julia
#=
pseudospectra_demo.jl
Chapter 17: Stability of Spectral Methods

Computational Etude 17.5: Pseudospectra of Normal vs Nonnormal Matrices.

Two-panel figure comparing pseudospectra:
  (a) Normal diagonal matrix with eigenvalues -1, -2, ..., -16.
  (b) Nonnormal upper triangular matrix with same diagonal but A[i,i+1] = 10.

Pseudospectra computed via resolvent norm:
  ||(zI - A)^{-1}|| = 1 / sigma_min(zI - A)

Generates Figure 17.5: pseudospectra_comparison.pdf

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

gr()

const NAVY = RGB(20/255, 45/255, 110/255)
const SKY  = RGB(120/255, 150/255, 210/255)
const CORAL = RGB(231/255, 76/255, 60/255)
const DARK_NAVY = RGB(10/255, 26/255, 64/255)

const SCRIPT_DIR = @__DIR__
const OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch17", "julia")
mkpath(OUTPUT_DIR)

"""
Compute log10(resolvent norm) on a grid.
resolvent_norm(z) = ||(zI - A)^{-1}|| = 1 / sigma_min(zI - A)
"""
function compute_pseudospectra(A, x_range, y_range; n_grid=200)
    xr = range(x_range[1], x_range[2], length=n_grid)
    yr = range(y_range[1], y_range[2], length=n_grid)

    n = size(A, 1)
    log_resolvent = zeros(n_grid, n_grid)

    for j in 1:n_grid
        for i in 1:n_grid
            z = complex(xr[j], yr[i])
            M = z * I(n) - A
            smin = svdvals(M)[end]   # smallest singular value
            if smin > 0
                log_resolvent[i, j] = log10(1.0 / smin)
            else
                log_resolvent[i, j] = 15.0  # effectively infinity
            end
        end
    end

    return collect(xr), collect(yr), log_resolvent
end

function main()
    N = 16

    # Normal matrix: diagonal with eigenvalues -1, -2, ..., -N
    A_normal = diagm(-collect(1.0:N))

    # Nonnormal matrix: same diagonal, upper diagonal entries = 10
    A_nonnormal = copy(A_normal)
    for i in 1:N-1
        A_nonnormal[i, i+1] = 10.0
    end

    eigenvalues = -collect(1.0:N)

    println("Computing pseudospectra of normal matrix...")
    x_range_n  = (-18, 4)
    y_range_n  = (-8, 8)
    X_n, Y_n, log_res_n = compute_pseudospectra(A_normal, x_range_n, y_range_n, n_grid=250)

    println("Computing pseudospectra of nonnormal matrix...")
    x_range_nn = (-18, 10)
    y_range_nn = (-12, 12)
    X_nn, Y_nn, log_res_nn = compute_pseudospectra(A_nonnormal, x_range_nn, y_range_nn, n_grid=250)

    # Contour levels: log10(resolvent_norm) = 1, 2, 3
    levels = [1.0, 2.0, 3.0]

    # --- Panel (a): Normal matrix ---
    p1 = plot(xlabel=L"\mathrm{Re}(z)", ylabel=L"\mathrm{Im}(z)",
              title=L"(a) Normal ($A = \mathrm{diag}(-1,\ldots,-16)$)",
              aspect_ratio=:equal, grid=true, gridalpha=0.15)
    contour!(p1, X_n, Y_n, log_res_n, levels=levels,
             c=cgrad([SKY, NAVY, DARK_NAVY]), lw=[1.0, 1.2, 1.5],
             colorbar=false)
    scatter!(p1, eigenvalues, zeros(N), color=CORAL, ms=4,
             markerstrokecolor=:white, markerstrokewidth=0.5,
             label="eigenvalues")
    hline!(p1, [0], color=:gray, lw=0.5, label="")
    vline!(p1, [0], color=:gray, lw=0.5, label="")

    # --- Panel (b): Nonnormal matrix ---
    p2 = plot(xlabel=L"\mathrm{Re}(z)", ylabel=L"\mathrm{Im}(z)",
              title=L"(b) Nonnormal ($A_{i,i+1} = 10$, same diagonal)",
              aspect_ratio=:equal, grid=true, gridalpha=0.15)
    contour!(p2, X_nn, Y_nn, log_res_nn, levels=levels,
             c=cgrad([SKY, NAVY, DARK_NAVY]), lw=[1.0, 1.2, 1.5],
             colorbar=false)
    scatter!(p2, eigenvalues, zeros(N), color=CORAL, ms=4,
             markerstrokecolor=:white, markerstrokewidth=0.5,
             label="eigenvalues")
    hline!(p2, [0], color=:gray, lw=0.5, label="")
    vline!(p2, [0], color=:gray, lw=0.5, label="")

    # Highlight pseudospectral bulge
    annotate!(p2, 2, 3, text("pseudospectral\nbulge into RHP", 8, CORAL, :center))

    fig = plot(p1, p2, layout=(1, 2), size=(1300, 500))

    savefig(fig, joinpath(OUTPUT_DIR, "pseudospectra_comparison.pdf"))
    savefig(fig, joinpath(OUTPUT_DIR, "pseudospectra_comparison.png"))
    println("Figure saved to $(joinpath(OUTPUT_DIR, "pseudospectra_comparison.pdf"))")

    display(fig)
end

main()
