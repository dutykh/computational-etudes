# bc_laplace_2d.jl
#
# Chapter 13: Boundary Conditions in Spectral Methods
#
# Laplace equation with piecewise boundary data on [-1,1]^2:
#
#     u_xx + u_yy = 0,  (x,y) in (-1,1)^2
#
# Boundary conditions:
#     u(x, 1)  = sin^4(pi*x)  for -1 < x < 0  (top, left half)
#     u(1, y)  = (1/5)*sin(3*pi*y)              (right boundary)
#     u = 0    everywhere else on the boundary
#
# Solved via boundary bordering (Method II):
#     1. Build full 2D Laplacian L = I x D^2 + D^2 x I
#     2. Replace boundary rows with identity rows and set boundary data
#     3. Solve the global linear system
#
# Generates Figure 13.5: 2D Laplace with piecewise BCs.
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

# Publication-quality CairoMakie configuration
set_theme!(Theme(
    fontsize = 11,
    fonts = (regular = "CMU Serif", bold = "CMU Serif Bold", italic = "CMU Serif Italic"),
    Axis = (xlabelsize = 11, ylabelsize = 11, titlesize = 12,
            xticklabelsize = 10, yticklabelsize = 10,
            spinewidth = 0.8, xtickwidth = 0.8, ytickwidth = 0.8),
    Legend = (labelsize = 10,),
))

# Book color scheme
const NAVY   = colorant"#142D6E"
const SKY    = colorant"#7896D2"
const CORAL  = colorant"#E74C3C"
const TEAL   = colorant"#16A085"
const PURPLE = colorant"#8E44AD"
const ORANGE = colorant"#E67E22"

# Output paths
const SCRIPT_DIR = @__DIR__
const OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch13", "julia")

# -----------------------------------------------------------------------------
# Boundary data functions
# -----------------------------------------------------------------------------

"""Top boundary u(x, 1) = sin^4(pi*x) for x < 0, else 0."""
function bc_top(x)
    val = zeros(length(x))
    for i in eachindex(x)
        if x[i] < 0
            val[i] = sin(pi * x[i])^4
        end
    end
    return val
end

"""Right boundary u(1, y) = (1/5)*sin(3*pi*y)."""
bc_right(y) = (1.0 / 5.0) .* sin.(3.0 .* pi .* y)

# -----------------------------------------------------------------------------
# Solver
# -----------------------------------------------------------------------------
"""
    solve_laplace_2d(N)

Solve the 2D Laplace equation with piecewise boundary data.

Uses boundary bordering (Method II): build the full (N+1)^2 system,
then replace each boundary row with an identity row and the
corresponding boundary value.
"""
function solve_laplace_2d(N)
    D2, D, x = cheb_second_derivative_matrix(N)
    n = N + 1  # number of points per direction

    # Build 2D Laplacian: L = I x D2 + D2 x I
    I_n = Matrix(1.0I, n, n)
    L = kron(I_n, D2) + kron(D2, I_n)

    # Create meshgrid (column-major: x varies fastest)
    # Grid point (i, j) -> index (j-1)*n + i
    # x-coordinate: x[i], y-coordinate: x[j]
    xx = zeros(n * n)
    yy = zeros(n * n)
    for j in 1:n
        for i in 1:n
            k = (j - 1) * n + i
            xx[k] = x[i]
            yy[k] = x[j]
        end
    end

    # Right-hand side (zero for Laplace equation)
    rhs = zeros(n * n)

    # Identify boundary nodes and apply boundary conditions
    for k in 1:n*n
        xk, yk = xx[k], yy[k]
        is_boundary = (abs(abs(xk) - 1.0) < 1e-14 || abs(abs(yk) - 1.0) < 1e-14)

        if is_boundary
            # Replace row with identity row
            L[k, :] .= 0.0
            L[k, k] = 1.0

            # Set boundary value
            if abs(yk - 1.0) < 1e-14
                # Top boundary: y = 1
                rhs[k] = xk < 0 ? sin(pi * xk)^4 : 0.0
            elseif abs(xk - 1.0) < 1e-14
                # Right boundary: x = 1
                rhs[k] = (1.0 / 5.0) * sin(3.0 * pi * yk)
            else
                # All other boundaries: u = 0
                rhs[k] = 0.0
            end
        end
    end

    # Solve the linear system
    u_vec = L \ rhs

    # Reshape to 2D: U[i, j] = u at (x[i], x[j])
    U = reshape(u_vec, n, n)

    # Create meshgrid arrays for plotting
    XX = [x[i] for i in 1:n, j in 1:n]
    YY = [x[j] for i in 1:n, j in 1:n]

    return XX, YY, U
end

# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
function main()
    mkpath(OUTPUT_DIR)

    N = 24
    XX, YY, U = solve_laplace_2d(N)

    # =========================================================================
    # Figure: 3D surface plot
    # =========================================================================
    fig = Figure(size = (650, 500))

    ax = Axis3(fig[1, 1],
               xlabel = L"x",
               ylabel = L"y",
               zlabel = L"u(x,y)",
               title = L"Laplace Equation: $\nabla^2 u = 0$ with Piecewise BCs ($N = %$(N)$)",
               azimuth = -55 * pi / 180,
               elevation = 25 * pi / 180)

    # Custom colormap using book colors
    book_cmap = cgrad([NAVY, TEAL, SKY, colorant"#FFFFFF"])

    surface!(ax, XX, YY, U, colormap = book_cmap, alpha = 0.9)

    # Save figure
    save(joinpath(OUTPUT_DIR, "laplace_2d_mixed.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "laplace_2d_mixed.png"), fig, px_per_unit = 4)
    println("Figure saved to: $(joinpath(OUTPUT_DIR, "laplace_2d_mixed.pdf"))")

    # =========================================================================
    # Print summary
    # =========================================================================
    println("\n2D Laplace Equation with Piecewise Boundary Data")
    println("="^55)
    @printf("  Grid size:          N = %d  (%d unknowns)\n", N, (N+1)^2)
    println("  Domain:             [-1, 1]^2")
    @printf("  Solution range:     [%.6f, %.6f]\n", minimum(U), maximum(U))
    println("  BCs applied:")
    println("    Top (y=1, x<0):   sin^4(pi*x)")
    println("    Right (x=1):      (1/5)*sin(3*pi*y)")
    println("    Elsewhere:        u = 0")
    println("="^55)

    # Verify boundary conditions
    n = N + 1
    x = [cos(pi * j / N) for j in 0:N]

    # Top boundary check (j=1, y=x[1]=1)
    u_top = U[:, 1]
    bc_top_exact = bc_top(x)
    top_err = maximum(abs.(u_top .- bc_top_exact))

    # Right boundary check (i=1, x=x[1]=1)
    u_right = U[1, :]
    bc_right_exact = bc_right(x)
    right_err = maximum(abs.(u_right .- bc_right_exact))

    # Bottom boundary check (j=n, y=x[N+1]=-1)
    u_bottom = U[:, n]
    bottom_err = maximum(abs.(u_bottom))

    # Left boundary check (i=n, x=x[N+1]=-1)
    u_left = U[n, :]
    left_err = maximum(abs.(u_left))

    println("\n  BC verification:")
    @printf("    Top boundary error:    %.2e\n", top_err)
    @printf("    Right boundary error:  %.2e\n", right_err)
    @printf("    Bottom boundary error: %.2e\n", bottom_err)
    @printf("    Left boundary error:   %.2e\n", left_err)
end

main()
