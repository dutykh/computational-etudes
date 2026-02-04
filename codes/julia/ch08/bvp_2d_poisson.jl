# bvp_2d_poisson.jl
#
# Solves the 2D Poisson equation using tensor product Chebyshev methods:
#
#     u_xx + u_yy = f(x,y),  (x,y) in (-1,1)^2,  u = 0 on boundary
#
# Test problem: f(x,y) = -2*pi^2 * sin(pi*x) * sin(pi*y)
# Exact solution: u(x,y) = sin(pi*x) * sin(pi*y)
#
# This demonstrates tensor product grids and Kronecker product operators.
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

# Import Chebyshev matrix function from Chapter 7
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
# Book color scheme
# -----------------------------------------------------------------------------
const NAVY   = colorant"#142D6E"
const SKY    = colorant"#7896D2"
const CORAL  = colorant"#E74C3C"
const TEAL   = colorant"#16A085"
const PURPLE = colorant"#8E44AD"
const ORANGE = colorant"#E67E22"

# Output paths
const SCRIPT_DIR = @__DIR__
const OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch08", "julia")

# -----------------------------------------------------------------------------
# Problem definition
# -----------------------------------------------------------------------------
"""f(x,y) = -2*pi^2 * sin(pi*x) * sin(pi*y)."""
rhs_function(x, y) = -2π^2 * sin(π * x) * sin(π * y)

"""u(x,y) = sin(pi*x) * sin(pi*y)."""
exact_solution(x, y) = sin(π * x) * sin(π * y)

# -----------------------------------------------------------------------------
# Solver
# -----------------------------------------------------------------------------
"""
    solve_2d_poisson(N)

Solve 2D Poisson equation u_xx + u_yy = f with u = 0 on boundary.

Uses Kronecker product formulation: L = I ⊗ D^2 + D^2 ⊗ I

# Returns
- `X, Y`: 2D meshgrid of Chebyshev points
- `U`: solution on the grid
- `L_mat`: 2D Laplacian operator matrix
"""
function solve_2d_poisson(N)
    D2, _, x = cheb_second_derivative_matrix(N)

    # Interior points only
    n_int = N - 1
    D2_int = D2[2:N, 2:N]
    x_int = x[2:N]

    # Build 2D Laplacian using Kronecker products
    # L = I ⊗ D^2 + D^2 ⊗ I
    Imat = I(n_int)
    L_mat = kron(Imat, D2_int) + kron(D2_int, Imat)

    # Create meshgrid for interior points
    X_int = [x_int[j] for _ in 1:n_int, j in 1:n_int]
    Y_int = [x_int[i] for i in 1:n_int, _ in 1:n_int]

    # Right-hand side (flattened in column-major order, Julia default)
    F = [rhs_function(X_int[i, j], Y_int[i, j]) for i in 1:n_int, j in 1:n_int]
    f_vec = vec(F)    # Column-major (Julia default)

    # Solve the linear system
    u_vec = L_mat \ f_vec

    # Reshape solution
    U_int = reshape(u_vec, n_int, n_int)

    # Embed in full grid with boundary conditions
    U = zeros(N + 1, N + 1)
    U[2:N, 2:N] = U_int

    # Create full meshgrid
    X = [x[j] for _ in 1:N+1, j in 1:N+1]
    Y = [x[i] for i in 1:N+1, _ in 1:N+1]

    return X, Y, U, L_mat
end

function main()
    mkpath(OUTPUT_DIR)

    N = 16

    # Solve problem
    X, Y, U, L_mat = solve_2d_poisson(N)
    U_exact = [exact_solution(X[i, j], Y[i, j]) for i in 1:N+1, j in 1:N+1]
    error = maximum(abs.(U .- U_exact))

    # =========================================================================
    # Figure 1: Tensor product grid
    # =========================================================================
    fig1 = Figure(size=(520, 520))

    ax1 = Axis(fig1[1, 1],
        xlabel = L"x",
        ylabel = L"y",
        title = @sprintf("Chebyshev Tensor Product Grid (N = %d)", N),
        titlesize = 12,
        aspect = DataAspect())

    # Add grid lines
    for i in 1:N+1
        lines!(ax1, X[i, :], Y[i, :], color=SKY, linewidth=0.5, alpha=0.6)
        lines!(ax1, X[:, i], Y[:, i], color=SKY, linewidth=0.5, alpha=0.6)
    end

    # Plot grid points
    scatter!(ax1, vec(X), vec(Y), color=TEAL, markersize=6,
             strokecolor=:white, strokewidth=0.5)

    xlims!(ax1, -1.1, 1.1)
    ylims!(ax1, -1.1, 1.1)
    hidespines!(ax1, :t, :r)

    # Add annotation about clustering
    text!(ax1, 0.3, 0.5,
        text="Grid clusters\nnear boundaries",
        fontsize=9, color=RGBf(0.4, 0.4, 0.4))

    save(joinpath(OUTPUT_DIR, "tensor_grid.pdf"), fig1)
    save(joinpath(OUTPUT_DIR, "tensor_grid.png"), fig1, px_per_unit=4)
    println("Figure saved to: ", joinpath(OUTPUT_DIR, "tensor_grid.pdf"))

    # =========================================================================
    # Figure 2: 2D Poisson solution
    # =========================================================================
    fig2 = Figure(size=(950, 430))

    # 3D surface plot
    ax2a = Axis3(fig2[1, 1],
        xlabel = L"x",
        ylabel = L"y",
        zlabel = L"u(x,y)",
        title = @sprintf("Solution (max error: %.2e)", error),
        titlesize = 11,
        elevation = 25.0 * π / 180,
        azimuth = 45.0 * π / 180)

    # Create finer grid for plotting
    x_fine = range(-1, 1, length=100)
    U_exact_fine = [exact_solution(xi, yj) for xi in x_fine, yj in x_fine]

    surface!(ax2a, collect(x_fine), collect(x_fine), U_exact_fine,
             colormap=:viridis, alpha=0.8)
    scatter!(ax2a, vec(X), vec(Y), vec(U), color=CORAL, markersize=4)

    # Contour plot
    ax2b = Axis(fig2[1, 2],
        xlabel = L"x",
        ylabel = L"y",
        title = "Contour Plot with Grid Points",
        titlesize = 11,
        aspect = DataAspect())

    levels = range(-1, 1, length=21)
    contourf!(ax2b, collect(x_fine), collect(x_fine), U_exact_fine,
              levels=collect(levels), colormap=:viridis)
    contour!(ax2b, collect(x_fine), collect(x_fine), U_exact_fine,
             levels=collect(levels), color=(:white, 0.5), linewidth=0.5)

    # Mark numerical points
    scatter!(ax2b, vec(X), vec(Y), color=(CORAL, 0.6), markersize=4,
             strokewidth=0)

    Colorbar(fig2[1, 3], limits=(-1, 1), colormap=:viridis, label=L"u(x,y)")

    Label(fig2[0, :],
        L"2D Poisson: $\nabla^2 u = -2\pi^2 \sin(\pi x)\sin(\pi y)$",
        fontsize=13)

    save(joinpath(OUTPUT_DIR, "poisson_2d.pdf"), fig2)
    save(joinpath(OUTPUT_DIR, "poisson_2d.png"), fig2, px_per_unit=4)
    println("Figure saved to: ", joinpath(OUTPUT_DIR, "poisson_2d.pdf"))

    # =========================================================================
    # Figure 3: Laplacian sparsity pattern
    # =========================================================================
    fig3 = Figure(size=(520, 520))

    n_int = N - 1
    ax3 = Axis(fig3[1, 1],
        title = @sprintf("2D Laplacian Structure (N = %d, size %d×%d)", N, n_int^2, n_int^2),
        titlesize = 11,
        xlabel = "Column index",
        ylabel = "Row index",
        yreversed = true,
        aspect = DataAspect())

    # Show sparsity pattern (find nonzero entries)
    rows_nz = Int[]
    cols_nz = Int[]
    for j in 1:n_int^2
        for i in 1:n_int^2
            if abs(L_mat[i, j]) > 1e-14
                push!(rows_nz, i)
                push!(cols_nz, j)
            end
        end
    end

    scatter!(ax3, cols_nz, rows_nz, color=NAVY, markersize=1.5)

    # Add annotation
    text!(ax3, n_int^2 / 2, n_int^2 + n_int^2 * 0.08,
        text=L"L = I \otimes D^2 + D^2 \otimes I",
        fontsize=11, align=(:center, :top))

    save(joinpath(OUTPUT_DIR, "laplacian_sparsity.pdf"), fig3)
    save(joinpath(OUTPUT_DIR, "laplacian_sparsity.png"), fig3, px_per_unit=4)
    println("Figure saved to: ", joinpath(OUTPUT_DIR, "laplacian_sparsity.pdf"))

    # Print convergence table
    println("\nConvergence Table: 2D Poisson Equation")
    println("-" ^ 50)
    @printf("%6s %12s %14s\n", "N", "Grid Size", "Max Error")
    println("-" ^ 50)

    for N_val in [4, 8, 12, 16, 20, 24, 32]
        X_n, Y_n, U_n, _ = solve_2d_poisson(N_val)
        U_exact_n = [exact_solution(X_n[i, j], Y_n[i, j]) for i in 1:N_val+1, j in 1:N_val+1]
        err = maximum(abs.(U_n .- U_exact_n))
        grid_size = (N_val + 1)^2
        @printf("%6d %12d %14.2e\n", N_val, grid_size, err)
    end

    println("-" ^ 50)
end

main()
