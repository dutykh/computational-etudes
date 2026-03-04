# bvp_helmholtz.jl
#
# Solves the Helmholtz equation:
#
#     u_xx + u_yy + k^2 * u = f(x,y),  (x,y) in (-1,1)^2,  u = 0 on boundary
#
# with localized Gaussian forcing:
#     f(x,y) = exp(-20*[(x-0.3)^2 + (y+0.4)^2])
#
# Uses k = 7 to demonstrate near-resonance behavior with the (2,4) eigenmode.
# (Eigenvalue for mode (m,n) is k^2 = (pi/2)^2*(m^2 + n^2), so (2,4) gives
# k ~ 7.02)
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
using Interpolations
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
"""
Localized Gaussian forcing:
f(x,y) = exp(-sigma*[(x-x0)^2 + (y-y0)^2])
"""
function forcing_function(x, y; x0=0.3, y0=-0.4, sigma=20.0)
    return exp(-sigma * ((x - x0)^2 + (y - y0)^2))
end

# -----------------------------------------------------------------------------
# Solver
# -----------------------------------------------------------------------------
"""
    solve_helmholtz(N, k)

Solve Helmholtz equation u_xx + u_yy + k^2*u = f with u = 0 on boundary.

# Arguments
- `N::Int`: number of intervals in each direction
- `k::Float64`: wave number

# Returns
- `X, Y`: 2D meshgrid of Chebyshev points
- `U`: solution on the grid
"""
function solve_helmholtz(N, k)
    D2, _, x = cheb_second_derivative_matrix(N)

    # Interior points only
    n_int = N - 1
    D2_int = D2[2:N, 2:N]
    x_int = x[2:N]

    # Build 2D Laplacian + k^2*I using Kronecker products
    # L = I ⊗ D^2 + D^2 ⊗ I + k^2*(I ⊗ I)
    Imat = I(n_int)
    L = kron(Imat, D2_int) + kron(D2_int, Imat) + k^2 * I(n_int^2)

    # Create meshgrid for interior points
    X_int = [x_int[j] for _ in 1:n_int, j in 1:n_int]
    Y_int = [x_int[i] for i in 1:n_int, _ in 1:n_int]

    # Right-hand side (flattened in column-major order)
    F = [forcing_function(X_int[i, j], Y_int[i, j]) for i in 1:n_int, j in 1:n_int]
    f_vec = vec(F)

    # Solve the linear system
    u_vec = L \ f_vec

    # Reshape solution
    U_int = reshape(u_vec, n_int, n_int)

    # Embed in full grid with boundary conditions
    U = zeros(N + 1, N + 1)
    U[2:N, 2:N] = U_int

    # Create full meshgrid
    X = [x[j] for _ in 1:N+1, j in 1:N+1]
    Y = [x[i] for i in 1:N+1, _ in 1:N+1]

    return X, Y, U
end

function main()
    mkpath(OUTPUT_DIR)

    N = 32
    k = 7.0    # Near resonance with (2,4) mode

    # Theoretical resonance wavenumbers
    println("Resonance wavenumbers for first few modes:")
    println("-" ^ 50)
    @printf("%15s %12s\n", "Mode (m,n)", "k_mn")
    println("-" ^ 50)
    for m in 1:5
        for n in 1:5
            k_mn = (π / 2) * sqrt(m^2 + n^2)
            if k_mn < 10
                @printf("%15s %12.4f\n", "($m,$n)", k_mn)
            end
        end
    end
    println("-" ^ 50)

    # Solve Helmholtz equation
    X, Y, U = solve_helmholtz(N, k)

    # Create figure
    fig = Figure(size=(1000, 430))

    # Panel 1: 3D surface
    ax1 = Axis3(fig[1, 1],
        xlabel = L"x",
        ylabel = L"y",
        zlabel = L"u(x,y)",
        title = @sprintf("Solution (k = %.0f)", k),
        titlesize = 11,
        elevation = 25.0 * π / 180,
        azimuth = 45.0 * π / 180)

    # Get grid vectors for surface plotting
    x_grid = X[1, :]
    y_grid = Y[:, 1]

    surface!(ax1, x_grid, y_grid, U, colormap=:RdBu, alpha=0.9)

    # Panel 2: Contour plot
    ax2 = Axis(fig[1, 2],
        xlabel = L"x",
        ylabel = L"y",
        title = "Contour Plot",
        titlesize = 11,
        aspect = DataAspect())

    # Fine interpolation for smoother contours
    x_sorted = sort(X[1, :])
    y_sorted = sort(Y[:, 1])
    sort_x_idx = sortperm(X[1, :])
    sort_y_idx = sortperm(Y[:, 1])
    U_sorted = U[sort_y_idx, sort_x_idx]

    x_fine = range(-1, 1, length=200)
    y_fine = range(-1, 1, length=200)

    # Interpolate using Interpolations.jl
    itp = interpolate((collect(x_sorted), collect(y_sorted)), U_sorted', Gridded(Linear()))
    etp = extrapolate(itp, Line())
    U_fine = [etp(xi, yi) for xi in x_fine, yi in y_fine]

    vmax = maximum(abs.(U_fine))
    levels_vals = range(-vmax, vmax, length=31)

    contourf!(ax2, collect(x_fine), collect(y_fine), U_fine',
              levels=collect(levels_vals), colormap=:RdBu)
    contour!(ax2, collect(x_fine), collect(y_fine), U_fine',
             levels=collect(levels_vals[1:2:end]), color=(:black, 0.5), linewidth=0.3)

    # Mark forcing location
    scatter!(ax2, [0.3], [-0.4], color=ORANGE, markersize=15, marker=:star5,
             strokecolor=:white, strokewidth=1, label="Forcing center")

    axislegend(ax2, position=:rt, labelsize=9)

    Colorbar(fig[1, 3], limits=(-vmax, vmax), colormap=:RdBu, label=L"u(x,y)")

    # Main title
    Label(fig[0, :],
        @sprintf("Helmholtz Equation: ∇²u + k²u = f(x,y) (k = %.0f, near (2,4) resonance at k = 7.02)", k),
        fontsize=12)

    # Save figure
    save(joinpath(OUTPUT_DIR, "helmholtz.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "helmholtz.png"), fig, px_per_unit=4)

    println("\nFigure saved to: ", joinpath(OUTPUT_DIR, "helmholtz.pdf"))

    # Compare solutions at different k values
    println("\nSolution amplitude at forcing location for different k:")
    println("-" ^ 50)
    @printf("%8s %16s %15s\n", "k", "|u| at center", "Near mode")
    println("-" ^ 50)

    # Find forcing point in grid
    idx_x = argmin(abs.(X[1, :] .- 0.3))
    idx_y = argmin(abs.(Y[:, 1] .- (-0.4)))

    for k_test in [3.0, 5.0, 7.0, 7.02, 8.0, 10.0]
        _, _, U_test = solve_helmholtz(N, k_test)
        u_at_force = abs(U_test[idx_y, idx_x])

        # Find nearest eigenvalue
        nearest_mode = ""
        min_dist = Inf
        for m in 1:9
            for n in 1:9
                k_mn = (π / 2) * sqrt(m^2 + n^2)
                if abs(k_mn - k_test) < min_dist
                    min_dist = abs(k_mn - k_test)
                    nearest_mode = "($m,$n)"
                end
            end
        end

        @printf("%8.2f %16.6f %15s\n", k_test, u_at_force, nearest_mode)
    end

    println("-" ^ 50)
    println("\nNote: Amplitude increases dramatically near resonance values.")
end

main()
