# bvp_nonlinear.jl
#
# Solves the Bratu equation (nonlinear BVP) using Newton iteration:
#
#     u_xx + lambda * exp(u) = 0,  x in (-1, 1),  u(+/-1) = 0
#
# with lambda = 0.5 (subcritical case, unique solution exists).
#
# IMPORTANT: The critical value for [-1, 1] domain is lambda_c ~ 0.878.
# For lambda > lambda_c, no solution exists! We use lambda = 0.5 for a
# stable solution.
#
# This demonstrates that nonlinear BVPs are easily handled with spectral
# methods using Newton iteration.
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
# Solvers
# -----------------------------------------------------------------------------
"""
    solve_bratu_newton(N; lam=0.5, tol=1e-12, max_iter=30)

Solve the Bratu equation u_xx + lambda*exp(u) = 0 with u(+/-1) = 0 using Newton.

Newton iteration for F(u) = D^2*u + lambda*exp(u) = 0:
    J*delta_u = -F(u)
    J = D^2 + lambda*diag(exp(u))
    u_new = u + alpha*delta_u  (with line search for robustness)
"""
function solve_bratu_newton(N; lam=0.5, tol=1e-12, max_iter=30)
    D2, _, x = cheb_second_derivative_matrix(N)

    # Extract interior system
    D2_int = D2[2:N, 2:N]

    # Initial guess: u = 0
    u = zeros(N + 1)

    residuals = Float64[]

    for it in 1:max_iter
        # Residual: F(u) = D^2*u + lambda*exp(u) (interior points)
        exp_u = exp.(u[2:N])
        F = D2_int * u[2:N] + lam * exp_u

        res_norm = maximum(abs.(F))
        push!(residuals, res_norm)

        if res_norm < tol
            return x, u, residuals, it
        end

        # Jacobian: J = D^2 + lambda*diag(exp(u))
        J = D2_int + lam * Diagonal(exp_u)

        # Newton step
        delta_u = J \ (-F)

        # Backtracking line search to ensure residual decreases
        alpha = 1.0
        u_trial = copy(u)
        for _ in 1:10    # Max 10 backtracking steps
            u_trial[2:N] = u[2:N] + alpha * delta_u
            exp_u_trial = exp.(u_trial[2:N])
            F_trial = D2_int * u_trial[2:N] + lam * exp_u_trial
            res_trial = maximum(abs.(F_trial))

            if res_trial < res_norm
                break
            end
            alpha *= 0.5
        end

        # Update interior points with damped step
        u[2:N] = u[2:N] + alpha * delta_u
    end

    return x, u, residuals, max_iter
end

"""
    solve_bratu_fixedpoint(N; lam=0.5, tol=1e-10, max_iter=100)

Solve Bratu equation using fixed-point iteration for comparison.
Fixed-point iteration: D^2 u_new = -lambda*exp(u_old)
"""
function solve_bratu_fixedpoint(N; lam=0.5, tol=1e-10, max_iter=100)
    D2, _, x = cheb_second_derivative_matrix(N)
    D2_int = D2[2:N, 2:N]

    # Initial guess
    u = zeros(N + 1)

    residuals = Float64[]

    for it in 1:max_iter
        # RHS
        rhs = -lam * exp.(u[2:N])

        # Solve
        u_new_int = D2_int \ rhs

        # Compute change
        change = maximum(abs.(u_new_int .- u[2:N]))
        push!(residuals, change)

        u[2:N] = u_new_int

        if change < tol
            return x, u, residuals, it
        end
    end

    return x, u, residuals, max_iter
end

function main()
    mkpath(OUTPUT_DIR)

    # Create figure
    fig = Figure(size=(950, 400))

    # Panel 1: Solution for various N
    ax1 = Axis(fig[1, 1],
        xlabel = L"x",
        ylabel = L"u(x)",
        title = L"Bratu Equation: $u_{xx} + \lambda e^u = 0$",
        titlesize = 11)

    N_values = [8, 16, 32]
    colors = [CORAL, TEAL, PURPLE]

    for (N, color) in zip(N_values, colors)
        x, u, residuals, n_iter = solve_bratu_newton(N)
        scatterlines!(ax1, x, u, color=color, linewidth=1.2, markersize=5,
                      strokecolor=:white, strokewidth=0.3,
                      label=@sprintf("N = %d (%d iter)", N, n_iter))
    end

    xlims!(ax1, -1.05, 1.05)
    axislegend(ax1, position=:cb, labelsize=9, nbanks=3)
    hidespines!(ax1, :t, :r)

    # Panel 2: Convergence history
    ax2 = Axis(fig[1, 2],
        xlabel = "Iteration",
        ylabel = "Residual norm",
        title = "Convergence History (N = 16)",
        titlesize = 11,
        yscale = log10)

    # Compare fixed-point vs Newton for N=16
    N = 16
    _, _, res_fp, n_fp = solve_bratu_fixedpoint(N, tol=1e-14, max_iter=50)
    _, _, res_newton, n_newton = solve_bratu_newton(N, tol=1e-14, max_iter=20)

    scatterlines!(ax2, 1:length(res_fp), res_fp, color=CORAL,
                  linewidth=1.5, markersize=7, label="Fixed-point")
    scatterlines!(ax2, 1:length(res_newton), res_newton, color=TEAL,
                  linewidth=1.5, markersize=7, marker=:rect, label="Newton")

    ylims!(ax2, 1e-16, 10)
    axislegend(ax2, position=:rt, labelsize=9)
    hidespines!(ax2, :t, :r)

    # Add tolerance line
    hlines!(ax2, [1e-12], color=:gray, linestyle=:dash, linewidth=1)
    text!(ax2, 2, 3e-12, text="tol", fontsize=9, color=:gray)

    # Main title
    Label(fig[0, :],
        L"Nonlinear BVP: Bratu Equation with $\lambda = 0.5$",
        fontsize=13)

    # Save figure
    save(joinpath(OUTPUT_DIR, "bratu.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "bratu.png"), fig, px_per_unit=4)

    println("Figure saved to: ", joinpath(OUTPUT_DIR, "bratu.pdf"))

    # Generate convergence table
    println("\n" * "=" ^ 60)
    println("Table 7.1: Bratu Equation Convergence (Newton)")
    println("=" ^ 60)
    @printf("%6s %12s %16s\n", "N", "Iterations", "u(0)")
    println("-" ^ 60)

    for N in [4, 8, 16, 32, 64]
        x, u, _, n_iter = solve_bratu_newton(N)
        # Find index closest to x=0
        idx_0 = argmin(abs.(x))
        @printf("%6d %12d %16.10f\n", N, n_iter, u[idx_0])
    end

    println("=" ^ 60)
    println("\nNote: Newton iteration converges quadratically (in 5-8 iterations)")
    println("while fixed-point converges linearly (in ~18 iterations).")
end

main()
