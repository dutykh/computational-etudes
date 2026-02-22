# bc_inhom_lifting.jl
#
# Chapter 13: Boundary Conditions in Spectral Methods
#
# Inhomogeneous Poisson equation solved via the lifting technique.
#
#     u_xx = exp(4x),  x in (-1, 1),  u(-1) = 0,  u(1) = 1
#
# The lifting function w(x) = (x + 1)/2 satisfies the boundary conditions:
#     w(-1) = 0,  w(1) = 1
#
# Since w'' = 0, the homogeneous variable v = u - w satisfies:
#     v_xx = exp(4x),  v(-1) = 0,  v(1) = 0
#
# After solving for v with homogeneous Dirichlet BCs, we reconstruct u = v + w.
#
# Generates Figure 13.1: Inhomogeneous Poisson with lifting.
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

"""Lifting function w(x) = (x + 1)/2 satisfying w(-1)=0, w(1)=1."""
lifting_function(x) = (x .+ 1.0) ./ 2.0

"""Right-hand side: f(x) = exp(4x)."""
rhs_function(x) = exp.(4.0 .* x)

"""
    solve_poisson_lifting(N)

Solve u_xx = exp(4x) with u(-1)=0, u(1)=1 via lifting.

Steps:
    1. Lifting: w(x) = (x+1)/2, so v = u - w satisfies v(+-1) = 0
    2. Since w'' = 0: v_xx = exp(4x)
    3. Solve for v using homogeneous Dirichlet BCs (matrix stripping)
    4. Reconstruct u = v + w
"""
function solve_poisson_lifting(N)
    D2, D, x = cheb_second_derivative_matrix(N)

    # Interior system for homogeneous Dirichlet BCs
    # Julia 1-based: x[1]=1 (right), x[N+1]=-1 (left)
    # Interior points: indices 2:N
    D2_int = D2[2:N, 2:N]
    f_int = rhs_function(x[2:N])

    # Solve v_xx = f(x) with v(+-1) = 0
    v_int = D2_int \ f_int

    # Assemble full v with BCs
    v = zeros(N + 1)
    v[2:N] = v_int

    # Reconstruct u = v + w
    u = v .+ lifting_function(x)

    return x, u
end

"""Compute high-resolution reference solution."""
function compute_reference_solution(N_ref=100)
    return solve_poisson_lifting(N_ref)
end

# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
function main()
    mkpath(OUTPUT_DIR)

    # Compute reference solution on high-resolution grid
    N_ref = 100
    x_ref, u_ref = compute_reference_solution(N_ref)

    # Fine plotting grid via barycentric interpolation
    # Use Chebyshev polynomial fit for interpolation
    x_fine = collect(range(-1, 1, length=500))

    # Barycentric interpolation from Chebyshev grid
    function cheb_interpolate(x_src, u_src, x_target)
        n = length(x_src)
        w = [(-1.0)^j for j in 0:n-1]
        w[1] *= 0.5
        w[end] *= 0.5

        u_target = similar(x_target)
        for (k, xt) in enumerate(x_target)
            # Check if xt coincides with a grid point
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

    u_fine = cheb_interpolate(x_ref, u_ref, x_fine)

    # Create figure with two panels
    fig = Figure(size = (900, 380))

    # ---- Panel 1: Solution for N=16 ----
    N = 16
    x_num, u_num = solve_poisson_lifting(N)

    # Interpolate reference to comparison grid
    u_exact_grid = cheb_interpolate(x_ref, u_ref, x_num)
    error_val = maximum(abs.(u_num .- u_exact_grid))

    ax1 = Axis(fig[1, 1],
               xlabel = L"x",
               ylabel = L"u(x)",
               title = @sprintf("Solution (N = %d, max error: %.2e)", N, error_val))

    lines!(ax1, x_fine, u_fine, color = NAVY, linewidth = 1.5, label = "Reference")
    scatter!(ax1, x_num, u_num, color = TEAL, markersize = 8,
             strokecolor = :white, strokewidth = 0.5, label = "Spectral")

    xlims!(ax1, -1.05, 1.05)
    axislegend(ax1, position = :lt, framevisible = true, labelsize = 9)
    hidespines!(ax1, :t, :r)

    # ---- Panel 2: Convergence study ----
    ax2 = Axis(fig[1, 2],
               xlabel = L"N",
               ylabel = "Max Error",
               title = "Convergence",
               yscale = log10)

    N_values = collect(4:2:40)
    errors = Float64[]

    for N_val in N_values
        x_n, u_n = solve_poisson_lifting(N_val)
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
    Label(fig[0, :],
          L"Inhomogeneous Poisson: $u_{xx} = e^{4x}$, $u(-1)=0$, $u(1)=1$ (Lifting)",
          fontsize = 13)

    # Save figure
    save(joinpath(OUTPUT_DIR, "inhom_lifting.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "inhom_lifting.png"), fig, px_per_unit = 4)
    println("Figure saved to: $(joinpath(OUTPUT_DIR, "inhom_lifting.pdf"))")

    # Print convergence table
    println("\nConvergence Table: Inhomogeneous Poisson with Lifting")
    println("-"^40)
    @printf("%6s %14s\n", "N", "Max Error")
    println("-"^40)
    for (N_val, err) in zip(N_values, errors)
        @printf("%6d %14.2e\n", N_val, err)
    end
    println("-"^40)
end

main()
