# ho_biharmonic_poisson.jl
#
# Chapter 14: Higher-Order Spectral Methods
#
# Biharmonic Poisson BVP on a clamped square plate with manufactured solution:
#
#     Delta^2 u = f(x,y),   (x, y) in (-1, 1)^2
#     u = 0,  du/dn = 0      on the boundary
#
# Manufactured exact solution:
#     u(x,y) = sin^2(pi*x) * sin^2(pi*y)
#
# Right-hand side:
#     f(x,y) = 4*pi^4 * [4*cos(2*pi*x)*cos(2*pi*y)
#                         - cos(2*pi*x) - cos(2*pi*y)]
#
# Generates Figure: biharmonic_poisson.pdf
#   Left panel:  contourf of numerical solution at N = 24
#   Right panel: convergence (max error vs N, semilogy)
#
# Author: Dr. Denys Dutykh
#         Mathematics Department
#         Khalifa University of Science and Technology
#         Abu Dhabi, UAE
#
# Part of the book "Computational Etudes: A Spectral Approach"
# https://github.com/dutykh/computational-etudes

using CairoMakie
using Colors
using LinearAlgebra
using Printf

# Include Chebyshev matrix from Chapter 7
const SCRIPT_DIR = @__DIR__
include(joinpath(SCRIPT_DIR, "..", "ch07", "cheb_matrix.jl"))

# Output directory
const OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch14", "julia")

# Book color scheme
const NAVY   = colorant"#142D6E"
const SKY    = colorant"#7896D2"
const CORAL  = colorant"#E74C3C"
const TEAL   = colorant"#16A085"
const PURPLE = colorant"#8E44AD"
const ORANGE = colorant"#E67E22"

# Publication-quality theme
set_theme!(Theme(
    fontsize = 11,
    fonts = (;
        regular = "CMU Serif",
        bold    = "CMU Serif Bold",
        italic  = "CMU Serif Italic"
    ),
    Axis = (;
        xlabelsize = 11, ylabelsize = 11,
        xticklabelsize = 10, yticklabelsize = 10,
        titlesize = 12,
        spinewidth = 0.8,
        xtickwidth = 0.8, ytickwidth = 0.8
    )
))

# =============================================================================
# 1D clamped biharmonic operator via polynomial trick
# =============================================================================
function build_clamped_biharmonic(N)
    D, x = cheb_matrix(N)

    D2 = D * D
    D3 = D2 * D
    D4full = D3 * D

    x2 = x .^ 2
    s = zeros(N + 1)
    s[2:N] = 1.0 ./ (1.0 .- x2[2:N])
    S = Diagonal(s)

    D4mat = (Diagonal(1.0 .- x2) * D4full
             - 8.0 * Diagonal(x) * D3
             - 12.0 * D2) * S

    D4 = D4mat[2:N, 2:N]
    D2int = D2[2:N, 2:N]

    return D4, D2int, x
end

# =============================================================================
# Exact solution and RHS
# =============================================================================
u_exact(x, y) = sin(π * x)^2 * sin(π * y)^2

function f_rhs(x, y)
    c2x = cos(2π * x)
    c2y = cos(2π * y)
    return 4π^4 * (4c2x * c2y - c2x - c2y)
end

# =============================================================================
# Solver
# =============================================================================
function solve_biharmonic_poisson(N)
    D4, D2int, x = build_clamped_biharmonic(N)

    M = N - 1
    I_M = Matrix(1.0I, M, M)

    # 2D biharmonic via Kronecker products
    L = kron(I_M, D4) + kron(D4, I_M) + 2kron(D2int, D2int)

    # Interior grid
    xi = x[2:N]
    XX = [xi[j] for i in 1:M, j in 1:M]
    YY = [xi[i] for i in 1:M, j in 1:M]

    # RHS and solve
    f_vec = vec(f_rhs.(XX, YY))
    u_vec = L \ f_vec

    # Error
    u_ex_vec = vec(u_exact.(XX, YY))
    err = maximum(abs.(u_vec - u_ex_vec))

    # Reshape to full grid
    U = zeros(N + 1, N + 1)
    U[2:N, 2:N] = reshape(u_vec, (M, M))

    return U, x, err
end

# =============================================================================
# Main
# =============================================================================
function main()
    mkpath(OUTPUT_DIR)

    # ----- Convergence study -----
    N_values = collect(6:2:24)
    errors = Float64[]

    println("Biharmonic Poisson BVP with Manufactured Solution")
    println("=" ^ 55)
    println("  u(x,y) = sin²(πx) · sin²(πy)")
    println("  f(x,y) = 4π⁴(4cos(2πx)cos(2πy) - cos(2πx) - cos(2πy))")
    println()
    println("  Convergence study:")
    println("  " * "-"^40)
    @printf("  %4s   %12s   %12s\n", "N", "Interior pts", "Max error")
    println("  " * "-"^40)

    U_plot = nothing
    x_plot = nothing

    for N in N_values
        U, x, err = solve_biharmonic_poisson(N)
        push!(errors, err)
        @printf("  %4d   %12d   %12.4e\n", N, (N - 1)^2, err)

        if N == N_values[end]
            U_plot = U
            x_plot = x
        end
    end

    println("  " * "-"^40)
    println()

    # ----- Figure: two panels -----
    fig = Figure(size = (900, 400))

    # Left panel: contourf of solution
    ax1 = Axis(fig[1, 1],
               aspect = DataAspect(),
               title = @sprintf("Numerical solution (N = %d)", N_values[end]),
               xlabel = L"x", ylabel = L"y",
               xticks = [-1, 0, 1], yticks = [-1, 0, 1])

    max_val = maximum(abs.(U_plot))
    lvls = range(-max_val, max_val, length = 21)
    cf = contourf!(ax1, x_plot, x_plot, U_plot', levels = lvls, colormap = :RdBu)
    contour!(ax1, x_plot, x_plot, U_plot', levels = 8,
             color = NAVY, linewidth = 0.6)
    lines!(ax1, [-1, 1, 1, -1, -1], [-1, -1, 1, 1, -1],
           color = :black, linewidth = 1.0)
    xlims!(ax1, -1.05, 1.05)
    ylims!(ax1, -1.05, 1.05)
    Colorbar(fig[1, 1][1, 2], cf, width = 10)

    # Right panel: convergence
    ax2 = Axis(fig[1, 2],
               yscale = log10,
               title = "Spectral convergence",
               xlabel = L"N",
               ylabel = "Max error",
               xticks = N_values)

    scatterlines!(ax2, N_values, errors,
                  color = CORAL, markersize = 7, linewidth = 1.5,
                  label = "Max error")
    hlines!(ax2, [2.2e-16], color = :gray, linestyle = :dot,
            linewidth = 0.8, label = "Machine ε")
    axislegend(ax2, position = :rt, framevisible = true)

    # Save
    save(joinpath(OUTPUT_DIR, "biharmonic_poisson.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "biharmonic_poisson.png"), fig, px_per_unit = 4)
    @printf("  Figure saved to: %s\n", joinpath(OUTPUT_DIR, "biharmonic_poisson.pdf"))

    println("=" ^ 55)
end

main()
