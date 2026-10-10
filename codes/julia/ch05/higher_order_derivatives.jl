# higher_order_derivatives.jl
#
# Demonstrates spectral differentiation for higher-order derivatives (up to 4th order)
# and compares D^2 construction methods: matrix squaring vs direct formula.
#
# The test function is u(x) = exp(-sin(2x)) on [0, 2pi), which is smooth, periodic,
# and has higher frequency content than simpler test functions.
#
# This script generates:
# 1. A 2x2 figure showing derivatives from 1st to 4th order
# 2. A comparison figure for D^2 construction methods
# 3. A convergence table showing accuracy vs derivative order
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
    Axis = (xlabelsize = 11, ylabelsize = 11, titlesize = 12,
            xticklabelsize = 10, yticklabelsize = 10,
            spinewidth = 0.8, xtickwidth = 0.8, ytickwidth = 0.8),
    Legend = (labelsize = 10,),
))

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
const N_DEMO = 32   # Grid points for demonstration figure
const N_FINE = 500  # Fine grid for plotting exact solutions

# Book color scheme
const NAVY   = colorant"#142D6E"
const SKY    = colorant"#7896D2"
const CORAL  = colorant"#E74C3C"
const TEAL   = colorant"#16A085"
const PURPLE = colorant"#8E44AD"
const ORANGE = colorant"#E67E22"

# Output paths
SCRIPT_DIR = @__DIR__
OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch05", "julia")
mkpath(OUTPUT_DIR)

# -----------------------------------------------------------------------------
# Test function and its exact derivatives
# -----------------------------------------------------------------------------

"""u(x) = exp(-sin(2x))"""
test_function(x) = exp(-sin(2x))

"""u'(x) = -2*cos(2x)*exp(-sin(2x))"""
exact_derivative_1(x) = -2 * cos(2x) * exp(-sin(2x))

"""u''(x) = 4*(sin(2x) + cos^2(2x))*exp(-sin(2x))"""
function exact_derivative_2(x)
    s = sin(2x)
    c = cos(2x)
    return 4 * (s + c^2) * exp(-s)
end

"""u'''(x) = 8*cos(2x)*(sin^2(2x) - 3*sin(2x))*exp(-sin(2x))"""
function exact_derivative_3(x)
    s = sin(2x)
    c = cos(2x)
    return 8 * c * (s^2 - 3s) * exp(-s)
end

"""u''''(x) computed from the pattern."""
function exact_derivative_4(x)
    s = sin(2x)
    c = cos(2x)
    term = -s^3 + 3s^2 + 5s * c^2 - 3c^2 - s^2 * c^2
    return 16 * term * exp(-s)
end

# -----------------------------------------------------------------------------
# Spectral differentiation matrix (periodic)
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
# Figure 1: 2x2 figure showing derivatives from 1st to 4th order
# -----------------------------------------------------------------------------
function create_derivatives_figure()
    # Construct differentiation matrix
    D, x = spectral_diff_periodic(N_DEMO)

    # Sample test function at grid points
    u = test_function.(x)

    # Compute numerical derivatives via matrix powers
    u1_num = D * u
    u2_num = D * (D * u)
    u3_num = D * (D * (D * u))
    u4_num = D * (D * (D * (D * u)))

    # Fine grid for exact solutions
    x_fine = collect(range(0, 2pi, length = N_FINE))
    u1_exact_fine = exact_derivative_1.(x_fine)
    u2_exact_fine = exact_derivative_2.(x_fine)
    u3_exact_fine = exact_derivative_3.(x_fine)
    u4_exact_fine = exact_derivative_4.(x_fine)

    # Exact values at grid points for error computation
    u1_exact_grid = exact_derivative_1.(x)
    u2_exact_grid = exact_derivative_2.(x)
    u3_exact_grid = exact_derivative_3.(x)
    u4_exact_grid = exact_derivative_4.(x)

    # Compute errors
    err1 = maximum(abs.(u1_num .- u1_exact_grid))
    err2 = maximum(abs.(u2_num .- u2_exact_grid))
    err3 = maximum(abs.(u3_num .- u3_exact_grid))
    err4 = maximum(abs.(u4_num .- u4_exact_grid))

    # Create figure
    fig = Figure(size = (820, 620))

    # Supertitle
    Label(fig[0, 1:2],
          L"Higher-Order Spectral Derivatives of $u(x) = e^{-\sin(2x)}$ ($N = %$(N_DEMO)$)",
          fontsize = 13)

    # Pi tick positions and labels
    pi_ticks = ([0, pi/2, pi, 3pi/2, 2pi],
                [L"0", L"\pi/2", L"\pi", L"3\pi/2", L"2\pi"])

    # Derivative data for each panel
    derivatives = [
        (u1_exact_fine, u1_num, err1, L"u'(x)",    "First derivative"),
        (u2_exact_fine, u2_num, err2, L"u''(x)",   "Second derivative"),
        (u3_exact_fine, u3_num, err3, L"u'''(x)",  "Third derivative"),
        (u4_exact_fine, u4_num, err4, L"u''''(x)", "Fourth derivative"),
    ]

    colors_num = [CORAL, TEAL, PURPLE, ORANGE]

    for (idx, (exact_fine, num, err, ylabel, title)) in enumerate(derivatives)
        row = (idx - 1) ÷ 2 + 1
        col = (idx - 1) % 2 + 1

        ax = Axis(fig[row, col],
                  xlabel = L"x",
                  ylabel = ylabel,
                  title = title,
                  xticks = pi_ticks)

        # Plot exact solution
        lines!(ax, x_fine, exact_fine,
               color = NAVY, linewidth = 1.5, label = "Exact")

        # Plot numerical solution
        scatter!(ax, x, num,
                 color = colors_num[idx], markersize = 5, label = "Spectral")

        xlims!(ax, 0, 2pi)

        axislegend(ax, position = :rt, labelsize = 9)

        # Add error annotation
        text!(ax, 0.95, 0.05,
              text = @sprintf("Max error: %.2e", err),
              fontsize = 9, color = :black,
              space = :relative,
              align = (:right, :bottom),
              )

        # Clean styling
        hidespines!(ax, :t, :r)
        ax.leftspinecolor   = RGBf(0.267, 0.267, 0.267)
        ax.bottomspinecolor = RGBf(0.267, 0.267, 0.267)
        ax.xtickcolor = RGBf(0.267, 0.267, 0.267)
        ax.ytickcolor = RGBf(0.267, 0.267, 0.267)
        ax.ygridvisible = true
        ax.ygridcolor   = RGBAf(0.533, 0.533, 0.533, 0.2)
        ax.ygridstyle   = :solid
        ax.xgridvisible = true
        ax.xgridcolor   = RGBAf(0.533, 0.533, 0.533, 0.2)
        ax.xgridstyle   = :solid
    end

    # Save figure
    output_file = joinpath(OUTPUT_DIR, "higher_order_derivatives.pdf")
    save(output_file, fig)
    save(joinpath(OUTPUT_DIR, "higher_order_derivatives.png"), fig, px_per_unit = 4)
    println("Figure saved to: $output_file")

    return err1, err2, err3, err4
end

# -----------------------------------------------------------------------------
# Figure 2: D^2 matrix squaring demonstration (3-panel)
# -----------------------------------------------------------------------------
function create_d2_comparison_figure()
    N = 16  # Grid size for visualization

    # Construct D and D^2
    D, x = spectral_diff_periodic(N)
    D2 = D * D

    # Test second derivative accuracy
    u = test_function.(x)
    u2_num = D2 * u
    u2_exact_grid = exact_derivative_2.(x)
    error_d2 = maximum(abs.(u2_num .- u2_exact_grid))

    # Fine grid for plotting
    x_fine = collect(range(0, 2pi, length = N_FINE))

    # Create figure
    fig = Figure(size = (960, 340))

    # -------------------------------------------------------------------------
    # Panel 1: Structure of D^2 matrix
    # -------------------------------------------------------------------------
    ax1 = Axis(fig[1, 1],
               xlabel = L"Column $j$",
               ylabel = L"Row $i$",
               title = L"$D^2 = D \cdot D$ (matrix structure)",
               yreversed = true,
               aspect = DataAspect())

    vmax = maximum(abs.(D2))
    hm = heatmap!(ax1, 0:N-1, 0:N-1, D2',
                  colormap = :RdBu,
                  colorrange = (-vmax, vmax))

    Colorbar(fig[1, 1], hm, width = 10, tellwidth = true,
             vertical = true, halign = :right)

    ax1.xticks = [0, N ÷ 2, N - 1]
    ax1.yticks = [0, N ÷ 2, N - 1]

    # -------------------------------------------------------------------------
    # Panel 2: Eigenvalue spectrum of D^2
    # -------------------------------------------------------------------------
    ax2 = Axis(fig[1, 2],
               xlabel = "Index",
               ylabel = "Eigenvalue",
               title = L"Eigenvalues of $D^2$")

    eigvals_d2 = eigvals(D2)
    eigvals_sorted = sort(real.(eigvals_d2))

    # Theoretical eigenvalues for d^2/dx^2: -k^2 for k = -N/2+1, ..., N/2
    k_vals = (-N ÷ 2 + 1):(N ÷ 2)
    theoretical_eigvals = sort(Float64.(-k_vals .^ 2))

    scatter!(ax2, 1:N, eigvals_sorted,
             color = NAVY, markersize = 6, label = "Numerical")
    scatter!(ax2, 1:N, theoretical_eigvals,
             color = CORAL, markersize = 5, marker = :rect,
             strokecolor = CORAL, strokewidth = 1.5,
             label = L"Theoretical $-k^2$")

    axislegend(ax2, position = :lt, labelsize = 9)

    # Clean styling
    hidespines!(ax2, :t, :r)
    ax2.leftspinecolor   = RGBf(0.267, 0.267, 0.267)
    ax2.bottomspinecolor = RGBf(0.267, 0.267, 0.267)
    ax2.ygridvisible = true
    ax2.ygridcolor   = RGBAf(0.533, 0.533, 0.533, 0.3)
    ax2.ygridstyle   = :solid
    ax2.xgridvisible = true
    ax2.xgridcolor   = RGBAf(0.533, 0.533, 0.533, 0.3)
    ax2.xgridstyle   = :solid

    # -------------------------------------------------------------------------
    # Panel 3: Second derivative accuracy
    # -------------------------------------------------------------------------
    pi_ticks = ([0, pi/2, pi, 3pi/2, 2pi],
                [L"0", L"\pi/2", L"\pi", L"3\pi/2", L"2\pi"])

    ax3 = Axis(fig[1, 3],
               xlabel = L"x",
               ylabel = L"u''(x)",
               title = @sprintf("Second derivative (error = %.2e)", error_d2),
               xticks = pi_ticks)

    lines!(ax3, x_fine, exact_derivative_2.(x_fine),
           color = NAVY, linewidth = 1.5, label = "Exact")
    scatter!(ax3, x, u2_num,
             color = TEAL, markersize = 6, label = L"D^2 u")

    xlims!(ax3, 0, 2pi)
    axislegend(ax3, position = :rt, labelsize = 9)

    # Clean styling
    hidespines!(ax3, :t, :r)
    ax3.leftspinecolor   = RGBf(0.267, 0.267, 0.267)
    ax3.bottomspinecolor = RGBf(0.267, 0.267, 0.267)
    ax3.ygridvisible = true
    ax3.ygridcolor   = RGBAf(0.533, 0.533, 0.533, 0.2)
    ax3.ygridstyle   = :solid
    ax3.xgridvisible = true
    ax3.xgridcolor   = RGBAf(0.533, 0.533, 0.533, 0.2)
    ax3.xgridstyle   = :solid

    # Supertitle
    Label(fig[0, 1:3],
          L"Matrix Squaring for Second Derivatives ($N = %$(N)$)",
          fontsize = 13)

    # Save figure
    output_file = joinpath(OUTPUT_DIR, "d2_comparison.pdf")
    save(output_file, fig)
    save(joinpath(OUTPUT_DIR, "d2_comparison.png"), fig, px_per_unit = 4)
    println("Figure saved to: $output_file")

    return error_d2
end

# -----------------------------------------------------------------------------
# Convergence table
# -----------------------------------------------------------------------------
function generate_convergence_table()
    N_values = [8, 16, 32, 64]
    results = Tuple{Int, Float64, Float64, Float64, Float64}[]

    println("\n" * "=" ^ 70)
    println("Convergence Table: Higher-Order Derivatives")
    println("Test function: u(x) = exp(-sin(2x))")
    println("=" ^ 70)
    @printf("%6s %14s %14s %14s %14s\n", "N", "Error u'", "Error u''", "Error u'''", "Error u''''")
    println("-" ^ 70)

    for N in N_values
        D, x = spectral_diff_periodic(N)
        u = test_function.(x)

        # Numerical derivatives
        u1_num = D * u
        u2_num = D * (D * u)
        u3_num = D * (D * (D * u))
        u4_num = D * (D * (D * (D * u)))

        # Exact derivatives at grid points
        u1_exact = exact_derivative_1.(x)
        u2_exact = exact_derivative_2.(x)
        u3_exact = exact_derivative_3.(x)
        u4_exact = exact_derivative_4.(x)

        # Errors
        err1 = maximum(abs.(u1_num .- u1_exact))
        err2 = maximum(abs.(u2_num .- u2_exact))
        err3 = maximum(abs.(u3_num .- u3_exact))
        err4 = maximum(abs.(u4_num .- u4_exact))

        push!(results, (N, err1, err2, err3, err4))
        @printf("%6d %14.2e %14.2e %14.2e %14.2e\n", N, err1, err2, err3, err4)
    end

    println("=" ^ 70)

    return results
end

# -----------------------------------------------------------------------------
# Main function
# -----------------------------------------------------------------------------
function main()
    println("=" ^ 70)
    println("Higher-Order Derivatives Demonstration")
    println("=" ^ 70)

    # Create main figure
    println("\n1. Creating higher-order derivatives figure...")
    err1, err2, err3, err4 = create_derivatives_figure()
    @printf("   Errors: u'=%.2e, u''=%.2e, u'''=%.2e, u''''=%.2e\n",
            err1, err2, err3, err4)

    # Create D^2 comparison figure
    println("\n2. Creating D^2 matrix squaring demonstration...")
    d2_error = create_d2_comparison_figure()
    @printf("   Second derivative error via D^2: %.2e\n", d2_error)

    # Generate convergence table
    println("\n3. Generating convergence table...")
    results = generate_convergence_table()

    println("\n" * "=" ^ 70)
    println("All figures generated successfully!")
    println("=" ^ 70)
end

main()
