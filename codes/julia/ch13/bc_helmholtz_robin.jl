# bc_helmholtz_robin.jl
#
# Chapter 13: Boundary Conditions in Spectral Methods
#
# Helmholtz equation with Robin boundary conditions (boundary bordering / tau method).
#
#     -u_xx + k^2 u = f(x),  x in (-1, 1)
#     u'(-1) - alpha u(-1) = 0   (Robin at left)
#     u(1) = 0                    (Dirichlet at right)
#
# Parameters: k = 2, f(x) = exp(-x^2)
#
# The full (N+1) x (N+1) system is built and boundary rows are overwritten:
#     Row 1   (x = 1):   identity row for Dirichlet u(1) = 0
#     Row N+1 (x = -1):  D[N+1,:] - alpha * e_{N+1}^T for Robin condition
#
# Generates:
#     Figure 13.2: Helmholtz solutions for various alpha and convergence
#     Figure 13.3: Condition number vs alpha
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
# Problem setup
# -----------------------------------------------------------------------------

"""Right-hand side: f(x) = exp(-x^2)."""
rhs_function(x) = exp.(-x.^2)

"""
    solve_helmholtz_robin(N; k=2.0, alpha=5.0)

Solve -u_xx + k^2 u = f(x) with Robin/Dirichlet BCs using boundary bordering.

Row 1   (x=1, Dirichlet):  L[1,:] = e_1^T          rhs[1] = 0
Row N+1 (x=-1, Robin):     L[N+1,:] = D[N+1,:] - alpha * e_{N+1}^T   rhs[N+1] = 0
"""
function solve_helmholtz_robin(N; k=2.0, alpha=5.0)
    D2, D, x = cheb_second_derivative_matrix(N)
    n = N + 1

    # Build operator: L = -D^2 + k^2 I
    I_mat = Matrix(1.0I, n, n)
    L = -D2 + k^2 * I_mat
    rhs = rhs_function(x)

    # Row 1 (x[1] = 1): Dirichlet u(1) = 0
    L[1, :] .= 0.0
    L[1, 1] = 1.0
    rhs[1] = 0.0

    # Row N+1 (x[N+1] = -1): Robin u'(-1) - alpha * u(-1) = 0
    L[n, :] = D[n, :] - alpha * I_mat[n, :]
    rhs[n] = 0.0

    # Solve
    u = L \ rhs

    return x, u
end

# -----------------------------------------------------------------------------
# Barycentric interpolation helper
# -----------------------------------------------------------------------------
function cheb_interpolate(x_src, u_src, x_target)
    n = length(x_src)
    w = [(-1.0)^j for j in 0:n-1]
    w[1] *= 0.5
    w[end] *= 0.5

    u_target = similar(x_target)
    for (k, xt) in enumerate(x_target)
        idx = findfirst(i -> abs(x_src[i] - xt) < 1e-15, 1:n)
        if idx !== nothing
            u_target[k] = u_src[idx]
        else
            num = sum(w[j] / (xt - x_src[j]) * u_src[j] for j in 1:n)
            den = sum(w[j] / (xt - x_src[j]) for j in 1:n)
            u_target[k] = num / den
        end
    end
    return u_target
end

# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
function main()
    mkpath(OUTPUT_DIR)

    k = 2.0

    # =========================================================================
    # Figure 1: Solutions for various alpha + convergence
    # =========================================================================
    fig1 = Figure(size = (900, 380))

    # ---- Panel 1: Solutions for different alpha ----
    ax1 = Axis(fig1[1, 1],
               xlabel = L"x",
               ylabel = L"u(x)",
               title = L"Solutions ($N = 32$, $k = 2$)")

    alpha_values = [0, 1, 5, 20]
    colors_alpha = [NAVY, TEAL, CORAL, PURPLE]

    N_plot = 32
    for (alpha, col) in zip(alpha_values, colors_alpha)
        x, u = solve_helmholtz_robin(N_plot, k=k, alpha=Float64(alpha))
        scatterlines!(ax1, x, u, color = col, linewidth = 1.5,
                      markersize = 4, markercolor = col,
                      strokecolor = :white, strokewidth = 0.3,
                      label = L"\alpha = %$(alpha)")
    end

    xlims!(ax1, -1.05, 1.05)
    axislegend(ax1, position = :lt, framevisible = true, labelsize = 9)
    hidespines!(ax1, :t, :r)

    # ---- Panel 2: Convergence for alpha = 5 ----
    ax2 = Axis(fig1[1, 2],
               xlabel = L"N",
               ylabel = "Max Error",
               title = L"Convergence ($\alpha = 5$)",
               yscale = log10)

    alpha_conv = 5.0

    # Reference solution with high N
    N_ref = 100
    x_ref, u_ref = solve_helmholtz_robin(N_ref, k=k, alpha=alpha_conv)

    N_values = collect(4:2:40)
    errors = Float64[]

    for N_val in N_values
        x_n, u_n = solve_helmholtz_robin(N_val, k=k, alpha=alpha_conv)
        u_ref_n = cheb_interpolate(x_ref, u_ref, x_n)
        err = max(maximum(abs.(u_n .- u_ref_n)), 1e-16)
        push!(errors, err)
    end

    scatterlines!(ax2, N_values, errors, color = CORAL, linewidth = 1.5,
                  markersize = 6, markercolor = CORAL,
                  strokecolor = :white, strokewidth = 0.5)

    xlims!(ax2, 0, 44)
    ylims!(ax2, 1e-16, 1e1)
    hidespines!(ax2, :t, :r)

    # Machine precision line
    hlines!(ax2, [2.2e-16], color = :gray, linestyle = :dot, linewidth = 1)
    text!(ax2, 42, 4e-16, text = "Machine precision", fontsize = 9,
          color = :gray, align = (:right, :bottom))

    # Main title
    Label(fig1[0, :],
          L"Helmholtz with Robin BC: $-u_{xx} + k^2 u = e^{-x^2}$",
          fontsize = 13)

    # Save figure 1
    save(joinpath(OUTPUT_DIR, "helmholtz_robin.pdf"), fig1)
    save(joinpath(OUTPUT_DIR, "helmholtz_robin.png"), fig1, px_per_unit = 4)
    println("Figure 1 saved to: $(joinpath(OUTPUT_DIR, "helmholtz_robin.pdf"))")

    # =========================================================================
    # Figure 2: Condition number vs alpha
    # =========================================================================
    fig2 = Figure(size = (500, 380))

    ax3 = Axis(fig2[1, 1],
               xlabel = L"\alpha",
               ylabel = L"\mathrm{cond}(L)",
               title = L"Condition Number vs Robin Parameter $\alpha$",
               xscale = log10,
               yscale = log10)

    alpha_range = 10.0 .^ range(-2, 2, length=80)
    N_cond_values = [16, 32, 64]
    colors_cond = [TEAL, CORAL, PURPLE]

    for (N_c, col) in zip(N_cond_values, colors_cond)
        cond_numbers = Float64[]
        for alpha_val in alpha_range
            D2, D, x = cheb_second_derivative_matrix(N_c)
            n = N_c + 1
            I_mat = Matrix(1.0I, n, n)
            L = -D2 + k^2 * I_mat

            # Apply BCs
            L[1, :] .= 0.0
            L[1, 1] = 1.0
            L[n, :] = D[n, :] - alpha_val * I_mat[n, :]

            push!(cond_numbers, cond(L))
        end

        lines!(ax3, alpha_range, cond_numbers, color = col, linewidth = 1.5,
               label = L"N = %$(N_c)")
    end

    axislegend(ax3, position = :rb, framevisible = true, labelsize = 10)
    hidespines!(ax3, :t, :r)

    # Save figure 2
    save(joinpath(OUTPUT_DIR, "robin_conditioning.pdf"), fig2)
    save(joinpath(OUTPUT_DIR, "robin_conditioning.png"), fig2, px_per_unit = 4)
    println("Figure 2 saved to: $(joinpath(OUTPUT_DIR, "robin_conditioning.pdf"))")

    # Print convergence table
    @printf("\nConvergence Table: Helmholtz with Robin BC (alpha = %.0f)\n", alpha_conv)
    println("-"^40)
    @printf("%6s %14s\n", "N", "Max Error")
    println("-"^40)
    for (N_val, err) in zip(N_values, errors)
        @printf("%6d %14.2e\n", N_val, err)
    end
    println("-"^40)
end

main()
