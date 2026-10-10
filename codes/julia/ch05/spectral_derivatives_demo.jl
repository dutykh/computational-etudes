# spectral_derivatives_demo.jl
#
# Demonstrates the periodic spectral differentiation matrix by computing
# first and second derivatives of a smooth periodic function.
#
# Test function: u(x) = exp(sin^2(x))
# - First derivative:  u'(x)  = sin(2x) * exp(sin^2(x))
# - Second derivative: u''(x) = [sin^2(2x) + 2*cos(2x)] * exp(sin^2(x))
#
# With just N = 64 grid points, the spectral method achieves near machine precision!
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

# -----------------------------------------------------------------------------
# Publication-quality CairoMakie configuration
# -----------------------------------------------------------------------------
set_theme!(Theme(
    fontsize = 11,
    fonts = (regular = "CMU Serif", bold = "CMU Serif Bold", italic = "CMU Serif Italic"),
    Axis = (xlabelsize = 12, ylabelsize = 12, titlesize = 12,
            xticklabelsize = 10, yticklabelsize = 10,
            spinewidth = 0.8, xtickwidth = 0.8, ytickwidth = 0.8),
    Legend = (labelsize = 10,),
))

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
const NAVY   = colorant"#142D6E"
const SKY    = colorant"#7896D2"
const CORAL  = colorant"#E74C3C"
const TEAL   = colorant"#16A085"
const PURPLE = colorant"#8E44AD"

SCRIPT_DIR = @__DIR__
OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch05", "julia")
mkpath(OUTPUT_DIR)

# -----------------------------------------------------------------------------
# Periodic spectral differentiation matrix
# -----------------------------------------------------------------------------
"""
    spectral_diff_periodic(N)

Construct the periodic spectral differentiation matrix for N equispaced
points on [0, 2pi).
"""
function spectral_diff_periodic(N)
    N < 2 && error("Need at least 2 points")

    h = 2pi / N
    x = h .* (0:N-1)

    D = zeros(N, N)
    for i in 1:N
        for j in 1:N
            if i != j
                diff = i - j
                D[i, j] = 0.5 * ((-1)^diff) / tan(diff * pi / N)
            end
        end
    end

    return D, collect(x)
end

# -----------------------------------------------------------------------------
# Main computation and visualization
# -----------------------------------------------------------------------------
function main()
    # Number of grid points
    N = 64

    # Construct the spectral differentiation matrix
    D, x = spectral_diff_periodic(N)

    # Fine grid for plotting exact solutions
    x_fine = range(0, 2pi, length = 200)

    # Test function: u(x) = exp(sin^2(x))
    u      = exp.(sin.(x) .^ 2)
    u_fine = exp.(sin.(x_fine) .^ 2)

    # Exact first derivative: u'(x) = sin(2x) * exp(sin^2(x))
    u1_exact = sin.(2 .* x) .* exp.(sin.(x) .^ 2)
    u1_fine  = sin.(2 .* x_fine) .* exp.(sin.(x_fine) .^ 2)

    # Exact second derivative: u''(x) = [sin^2(2x) + 2*cos(2x)] * exp(sin^2(x))
    u2_exact = (sin.(2 .* x) .^ 2 .+ 2 .* cos.(2 .* x)) .* exp.(sin.(x) .^ 2)
    u2_fine  = (sin.(2 .* x_fine) .^ 2 .+ 2 .* cos.(2 .* x_fine)) .* exp.(sin.(x_fine) .^ 2)

    # Numerical derivatives using spectral differentiation matrix
    u1_num = D * u           # First derivative: D * u
    u2_num = D * (D * u)     # Second derivative: D^2 * u

    # Compute errors
    error_u1 = maximum(abs.(u1_num .- u1_exact))
    error_u2 = maximum(abs.(u2_num .- u2_exact))

    # -------------------------------------------------------------------------
    # Create figure
    # -------------------------------------------------------------------------
    fig = Figure(size = (880, 380))

    # Pi tick positions and labels
    pi_ticks = ([0, pi/2, pi, 3pi/2, 2pi],
                [L"0", L"\pi/2", L"\pi", L"3\pi/2", L"2\pi"])

    # Left panel: u(x) and u'(x)
    ax1 = Axis(fig[1, 1],
               xlabel = L"x",
               ylabel = L"u(x),\; u'(x)",
               title = "First Derivative (N = $N)",
               xticks = pi_ticks)

    lines!(ax1, collect(x_fine), collect(u_fine),
           color = NAVY, linewidth = 2,
           label = L"u(x) = e^{\sin^2 x}")
    lines!(ax1, collect(x_fine), collect(u1_fine),
           color = CORAL, linewidth = 2,
           label = L"u'(x) exact")
    scatter!(ax1, x, u1_num,
             color = :white, markersize = 7,
             strokecolor = TEAL, strokewidth = 2,
             label = L"u'(x) spectral")

    xlims!(ax1, 0, 2pi)
    axislegend(ax1, position = :rt, labelsize = 9)
    ax1.xgridvisible = true
    ax1.xgridcolor = RGBAf(0, 0, 0, 0.1)
    ax1.ygridvisible = true
    ax1.ygridcolor = RGBAf(0, 0, 0, 0.1)

    # Add error annotation
    text!(ax1, 0.05, 0.05,
          text = @sprintf("Max error: %.2e", error_u1),
          fontsize = 10, color = TEAL,
          space = :relative,
          align = (:left, :bottom))

    # Right panel: u''(x)
    ax2 = Axis(fig[1, 2],
               xlabel = L"x",
               ylabel = L"u''(x)",
               title = "Second Derivative (N = $N)",
               xticks = pi_ticks)

    lines!(ax2, collect(x_fine), collect(u2_fine),
           color = PURPLE, linewidth = 2,
           label = L"u''(x) exact")
    scatter!(ax2, x, u2_num,
             color = :white, markersize = 7, marker = :rect,
             strokecolor = TEAL, strokewidth = 2,
             label = L"u''(x) spectral")

    xlims!(ax2, 0, 2pi)
    axislegend(ax2, position = :rt, labelsize = 9)
    ax2.xgridvisible = true
    ax2.xgridcolor = RGBAf(0, 0, 0, 0.1)
    ax2.ygridvisible = true
    ax2.ygridcolor = RGBAf(0, 0, 0, 0.1)

    # Add error annotation
    text!(ax2, 0.05, 0.05,
          text = @sprintf("Max error: %.2e", error_u2),
          fontsize = 10, color = TEAL,
          space = :relative,
          align = (:left, :bottom))

    # -------------------------------------------------------------------------
    # Save figure
    # -------------------------------------------------------------------------
    save(joinpath(OUTPUT_DIR, "spectral_derivatives_demo.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "spectral_derivatives_demo.png"), fig, px_per_unit = 4)
    println("Figure saved to: $(joinpath(OUTPUT_DIR, "spectral_derivatives_demo.pdf"))")

    # -------------------------------------------------------------------------
    # Print results
    # -------------------------------------------------------------------------
    println("\n" * "=" ^ 60)
    println("Spectral Differentiation Demo")
    println("=" ^ 60)
    println("Test function: u(x) = exp(sin^2(x))")
    println("Number of grid points: N = $N")
    println("-" ^ 60)
    @printf("First derivative  max error: %.4e\n", error_u1)
    @printf("Second derivative max error: %.4e\n", error_u2)
    println("-" ^ 60)
    println("With only 64 points, we achieve near machine precision!")
    println("=" ^ 60)

    # -------------------------------------------------------------------------
    # Convergence study for table
    # -------------------------------------------------------------------------
    println("\nConvergence Study (for table):")
    println("-" ^ 50)
    @printf("%6s  %14s  %14s\n", "N", "Error u'", "Error u''")
    println("-" ^ 50)

    for N_test in [8, 16, 32, 64]
        D_test, x_test = spectral_diff_periodic(N_test)
        u_test = exp.(sin.(x_test) .^ 2)
        u1_exact_test = sin.(2 .* x_test) .* exp.(sin.(x_test) .^ 2)
        u2_exact_test = (sin.(2 .* x_test) .^ 2 .+ 2 .* cos.(2 .* x_test)) .* exp.(sin.(x_test) .^ 2)
        u1_num_test = D_test * u_test
        u2_num_test = D_test * (D_test * u_test)
        err1 = maximum(abs.(u1_num_test .- u1_exact_test))
        err2 = maximum(abs.(u2_num_test .- u2_exact_test))
        @printf("%6d  %14.4e  %14.4e\n", N_test, err1, err2)
    end

    println("-" ^ 50)
end

main()
