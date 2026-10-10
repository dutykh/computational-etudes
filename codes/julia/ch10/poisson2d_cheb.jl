# poisson2d_cheb.jl
#
# Chapter 10: Spectral PDE Solvers with Chebyshev and Fourier Grids
#
# 2D Poisson equation solver on Chebyshev tensor grid.
#
# PDE: -Lap(u) = f, -1 < x, y < 1
# BCs: u = 0 on boundary (Dirichlet)
#
# Uses manufactured solution:
#     u(x, y) = (1 - x^2)(1 - y^2) cos(pi*x/2) cos(pi*y/2)
#
# which automatically satisfies homogeneous Dirichlet BCs.
#
# Generates Figure 10.7: Exact, numerical, and error comparison.
#
# Author: Dr. Denys Dutykh
#         Mathematics Department
#         Khalifa University of Science and Technology
#         Abu Dhabi, UAE
#
# Part of the book "Computational Etudes: A Spectral Approach"
# https://github.com/dutykh/computational-etudes

using CairoMakie
using LaTeXStrings
using Colors
using LinearAlgebra
using Printf

# Import cheb_matrix
include(joinpath(@__DIR__, "chebfft.jl"))
include(joinpath(@__DIR__, "..", "ch07", "cheb_matrix.jl"))

# -----------------------------------------------------------------------------
# Publication-quality CairoMakie configuration
# -----------------------------------------------------------------------------
set_theme!(Theme(
    fontsize = 11,
    fonts = (regular = "CMU Serif", bold = "CMU Serif Bold", italic = "CMU Serif Italic"),
    Axis = (xlabelsize = 11, ylabelsize = 11, titlesize = 12,
            xticklabelsize = 10, yticklabelsize = 10,
            spinewidth = 0.8, xtickwidth = 0.8, ytickwidth = 0.8),
    Legend = (labelsize = 10,),
))

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
const NAVY = colorant"#142D6E"

# Output path
SCRIPT_DIR = @__DIR__
OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch10", "julia")
mkpath(OUTPUT_DIR)

# -----------------------------------------------------------------------------
# Manufactured solution
# -----------------------------------------------------------------------------
"""
    exact_solution(xx, yy)

Manufactured exact solution that satisfies u = 0 on the boundary.
"""
exact_solution(xx, yy) = (1 .- xx.^2) .* (1 .- yy.^2) .* cos.(pi .* xx ./ 2) .* cos.(pi .* yy ./ 2)

"""
    exact_laplacian(xx, yy)

Analytical Laplacian of the manufactured solution.

`u = g(x)*h(y)` with `g(x) = (1-x^2)*cos(pi*x/2)`, `h(y) = (1-y^2)*cos(pi*y/2)`.
`Lap(u) = g''(x)*h(y) + g(x)*h''(y)`.
"""
function exact_laplacian(xx, yy)
    p = π
    cx = cos.(p .* xx ./ 2); cy = cos.(p .* yy ./ 2)
    sx = sin.(p .* xx ./ 2); sy = sin.(p .* yy ./ 2)
    gx  = (1 .- xx.^2) .* cx
    hy  = (1 .- yy.^2) .* cy
    gpp = -(2 .+ p^2/4 .* (1 .- xx.^2)) .* cx .+ 2p .* xx .* sx
    hpp = -(2 .+ p^2/4 .* (1 .- yy.^2)) .* cy .+ 2p .* yy .* sy
    return gpp .* hy .+ gx .* hpp
end

# -----------------------------------------------------------------------------
# Poisson solver
# -----------------------------------------------------------------------------
"""
    poisson2d_cheb(; N=32)

Solve 2D Poisson equation with manufactured solution.

The discrete system is:
    A * u_vec = f_vec
where A is the negative Kronecker sum of D2_int with itself.

Returns `(xx, yy, U, U_exact, error)`.
"""
function poisson2d_cheb(; N=32)
    # Build differentiation matrices
    D, x = cheb_matrix(N)
    D2 = D * D
    xx = [x[i] for i in 1:N+1, j in 1:N+1]
    yy = [x[j] for i in 1:N+1, j in 1:N+1]

    # Exact solution
    U_exact = exact_solution(xx, yy)

    # Interior indices
    ii = 2:N
    D2_int = D2[ii, ii]
    I_int = Matrix{Float64}(I, N - 1, N - 1)

    # Kronecker sum Laplacian for interior (represents Delta)
    L = kron(D2_int, I_int) + kron(I_int, D2_int)

    # For -Delta u = f, we have -L * u = f, so use A = -L
    A = -L

    # RHS from analytical Laplacian (not discrete, to test convergence)
    f_full = -exact_laplacian(xx, yy)
    f_int = f_full[ii, ii]

    # Solve A * u = f where A = -Delta
    u_vec = A \ vec(f_int)

    # Reconstruct full solution
    U = zeros(N + 1, N + 1)
    U[ii, ii] = reshape(u_vec, N - 1, N - 1)

    # Compute error
    error = maximum(abs.(U .- U_exact))

    return xx, yy, U, U_exact, error
end

# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
function main()
    # Solve Poisson equation
    N = 32
    xx, yy, U, U_exact, error = poisson2d_cheb(N=N)

    @printf("N = %d, Max error = %.2e\n", N, error)

    # Convergence study upfront so panel (d) can plot it
    println("\nConvergence Study:")
    println("-" ^ 40)
    @printf("%6s  %14s  %10s\n", "N", "Max Error", "Ratio")
    println("-" ^ 40)

    Ns_conv = [8, 12, 16, 20, 24, 28, 32, 40, 48]
    errs_conv = Float64[]
    prev_error = nothing
    for Nv in Ns_conv
        _, _, _, _, err_val = poisson2d_cheb(N=Nv)
        push!(errs_conv, err_val)
        if prev_error !== nothing && err_val > 1e-15
            ratio = prev_error / err_val
        else
            ratio = NaN
        end
        @printf("%6d  %14.6e  %10.2f\n", Nv, err_val, ratio)
        prev_error = err_val
    end

    # Create 2x2 figure
    fig = Figure(size = (1100, 900))

    vmin, vmax = minimum(U_exact), maximum(U_exact)
    n_levels = 30

    # ---- (a) Exact solution -----------------------------------------
    ax1 = Axis(fig[1, 1],
               xlabel = L"x", ylabel = L"y",
               title  = "(a) Exact solution u_exact(x, y)",
               aspect = DataAspect())
    contourf!(ax1, xx[:, 1], yy[1, :], U_exact,
              levels = range(vmin, vmax, length=n_levels), colormap = :RdBu)
    contour!(ax1, xx[:, 1], yy[1, :], U_exact,
             levels = 10, color = (:black, 0.5), linewidth = 0.3)

    # ---- (b) Numerical solution -------------------------------------
    ax2 = Axis(fig[1, 2],
               xlabel = L"x", ylabel = L"y",
               title  = "(b) Numerical solution u_N(x, y) at N = $N",
               aspect = DataAspect())
    contourf!(ax2, xx[:, 1], yy[1, :], U,
              levels = range(vmin, vmax, length=n_levels), colormap = :RdBu)
    contour!(ax2, xx[:, 1], yy[1, :], U,
             levels = 10, color = (:black, 0.5), linewidth = 0.3)

    # ---- (c) Pointwise log10 error ----------------------------------
    ax3 = Axis(fig[2, 1],
               xlabel = L"x", ylabel = L"y",
               title  = "(c) log₁₀ |u_N - u_exact|  (max = $(@sprintf("%.2e", error)))",
               aspect = DataAspect())
    err = abs.(U .- U_exact)
    err_plot = log10.(max.(err, 1e-16))
    contourf!(ax3, xx[:, 1], yy[1, :], err_plot,
              levels = 20, colormap = :viridis)

    # ---- (d) NEW: spectral convergence ------------------------------
    ax4 = Axis(fig[2, 2],
               xlabel = "N (Chebyshev degree, each direction)",
               ylabel = "max-norm error",
               yscale = log10,
               title  = "(d) spectral convergence: error vs N")
    scatterlines!(ax4, Ns_conv, errs_conv .+ 1e-18; color = NAVY,
                  markercolor = :white, strokecolor = NAVY, strokewidth = 1.0,
                  markersize = 6, linewidth = 1.2,
                  label = "max |u_N - u_exact|")
    hlines!(ax4, [1e-15]; color = :gray, linestyle = :dot, linewidth = 0.6,
            label = "machine epsilon")
    axislegend(ax4; position = :rt, labelsize = 9)

    Label(fig[0, :], "2D Poisson Equation: Spectral Solution", fontsize = 14)

    # Save figure
    save(joinpath(OUTPUT_DIR, "poisson2d_solution.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "poisson2d_solution.png"), fig, px_per_unit = 4)
    println("Figure saved to: $(joinpath(OUTPUT_DIR, "poisson2d_solution.pdf"))")
end

main()
