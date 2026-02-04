# cheb_convergence.jl
#
# Demonstrates spectral convergence rates for four functions of increasing
# smoothness. Shows the fundamental relationship between function regularity
# and convergence rate of spectral differentiation.
#
# Test functions (different from Trefethen's standard examples):
# 1. |x|^(5/2) - Third derivative of bounded variation (algebraic: ~N^{-2.5})
# 2. exp(-1/(1-x^2)) - Smooth bump function, C^inf but not analytic (superalgebraic)
# 3. tanh(5x) - Analytic in [-1,1], poles at +-i*pi/10 (exponential)
# 4. x^8 - Polynomial of degree 8 (exact for N >= 8)
#
# This script generates Figure 6.5 and Table 6.2 for Chapter 6.
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
using Printf

# Import Chebyshev matrix function
include(joinpath(@__DIR__, "cheb_matrix.jl"))

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
const PURPLE = colorant"#9B59B6"
const ORANGE = colorant"#E67E22"

# Output paths
SCRIPT_DIR = @__DIR__
OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch07", "julia")
mkpath(OUTPUT_DIR)

# -----------------------------------------------------------------------------
# Test functions and their exact derivatives
# -----------------------------------------------------------------------------

# Function 1: |x|^(5/2) - C^2, 3rd derivative in BV
f1(x) = abs(x)^2.5
f1_prime(x) = 2.5 * abs(x)^1.5 * sign(x)

# Function 2: exp(-1/(1-x^2)) - Smooth bump, C^inf but not analytic at boundaries
function f2(x)
    abs(x) < 1.0 ? exp(-1.0 / (1.0 - x^2)) : 0.0
end

function f2_prime(x)
    if abs(x) < 1.0
        inner = 1.0 - x^2
        return exp(-1.0 / inner) * (-2.0 * x) / inner^2
    else
        return 0.0
    end
end

# Function 3: tanh(5x) - Analytic, poles at x = +-i*pi/10
f3(x) = tanh(5.0 * x)
f3_prime(x) = 5.0 / cosh(5.0 * x)^2

# Function 4: x^8 - Polynomial degree 8
f4(x) = x^8
f4_prime(x) = 8.0 * x^7

# List of test functions
const TEST_FUNCTIONS = [
    (name = L"|x|^{5/2}",         f = f1, fp = f1_prime,
     smoothness = L"C^2, 3rd deriv in BV", rate = L"O(N^{-2.5})",
     color = CORAL),
    (name = L"e^{-1/(1-x^2)}",    f = f2, fp = f2_prime,
     smoothness = L"C^\infty, not analytic", rate = L"faster than any $N^{-k}$",
     color = TEAL),
    (name = L"\tanh(5x)",          f = f3, fp = f3_prime,
     smoothness = "Analytic",       rate = L"O(\rho^{-N})",
     color = PURPLE),
    (name = L"x^8",                f = f4, fp = f4_prime,
     smoothness = "Polynomial",     rate = L"Exact for $N \geq 8$",
     color = ORANGE),
]

# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
function main()
    # Range of N values
    N_values = [4, 6, 8, 10, 12, 16, 20, 24, 32, 40, 48, 64, 80, 96, 128]

    # Create 2x2 figure
    fig = Figure(size = (700, 550))

    # Main title
    Label(fig[0, :], "Spectral Convergence: Four Functions of Increasing Smoothness",
          fontsize = 13, font = "CMU Serif Bold")

    positions = [(1, 1), (1, 2), (2, 1), (2, 2)]

    for (idx, (tf, pos)) in enumerate(zip(TEST_FUNCTIONS, positions))
        ax = Axis(fig[pos[1], pos[2]],
                  xlabel = L"N",
                  ylabel = "Max Error",
                  title = tf.name,
                  yscale = log10)
        hidespines!(ax, :t, :r)
        xlims!(ax, 0, 135)
        ylims!(ax, 1e-16, 1e2)

        # Grid
        ax.xgridvisible = true
        ax.xgridcolor = RGBAf(0.533, 0.533, 0.533, 0.3)
        ax.ygridvisible = true
        ax.ygridcolor = RGBAf(0.533, 0.533, 0.533, 0.3)

        errors = Float64[]

        for Nv in N_values
            D, x = cheb_matrix(Nv)
            v = tf.f.(x)
            w = D * v
            w_exact = tf.fp.(x)

            # Compute max error (avoiding boundary issues for bump function)
            if idx == 2  # Bump function
                interior = findall(xi -> xi > -0.99 && xi < 0.99, x)
                if !isempty(interior)
                    err = maximum(abs.(w[interior] .- w_exact[interior]))
                else
                    err = maximum(abs.(w .- w_exact))
                end
            else
                err = maximum(abs.(w .- w_exact))
            end

            push!(errors, max(err, 1e-16))
        end

        # Plot convergence
        scatterlines!(ax, Float64.(N_values), errors, color = tf.color,
                      linewidth = 1.5, markersize = 6,
                      strokecolor = :white, strokewidth = 0.5)

        # Add reference lines for theoretical rates
        if idx == 1  # |x|^(5/2) - algebraic O(N^-2.5)
            mask = N_values .>= 16
            N_ref = Float64.(N_values[mask])
            err_ref = errors[mask]
            if !isempty(err_ref)
                C = err_ref[1] * N_ref[1]^2.5
                lines!(ax, N_ref, C .* N_ref .^ (-2.5),
                       color = :gray, linewidth = 1, linestyle = :dash, alpha = 0.7,
                       label = L"O(N^{-2.5})")
            end
        elseif idx == 4  # x^8 - exact for N >= 8
            hlines!(ax, 1e-14, color = :gray, linewidth = 1, linestyle = :dash, alpha = 0.7)
        end

        # Add convergence rate annotation
        text!(ax, 125, 10^1.5,
              text = string(tf.rate),
              fontsize = 9, align = (:right, :top))
    end

    rowgap!(fig.layout, 1, 5)
    rowgap!(fig.layout, 2, 15)

    # Save figure
    save(joinpath(OUTPUT_DIR, "convergence_waterfall.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "convergence_waterfall.png"), fig, px_per_unit = 4)
    println("Figure saved to: $(joinpath(OUTPUT_DIR, "convergence_waterfall.pdf"))")

    # Generate convergence table
    println("\n" * "=" ^ 80)
    println("Convergence Table: Chebyshev Spectral Differentiation")
    println("=" ^ 80)

    # Table header
    @printf("%6s  %14s  %14s  %14s  %14s\n",
            "N", "|x|^(5/2)", "bump", "tanh(5x)", "x^8")
    println("-" ^ 80)

    # Table rows
    for Nv in N_values
        D, x = cheb_matrix(Nv)
        errs = Float64[]
        for (idx, tf) in enumerate(TEST_FUNCTIONS)
            v = tf.f.(x)
            w = D * v
            w_exact = tf.fp.(x)

            if idx == 2  # Bump function
                interior = findall(xi -> xi > -0.99 && xi < 0.99, x)
                if !isempty(interior)
                    err = maximum(abs.(w[interior] .- w_exact[interior]))
                else
                    err = maximum(abs.(w .- w_exact))
                end
            else
                err = maximum(abs.(w .- w_exact))
            end
            push!(errs, err)
        end

        @printf("%6d", Nv)
        for err in errs
            if err < 1e-14
                @printf("  %14s", "< 1e-14")
            else
                @printf("  %14.2e", err)
            end
        end
        println()
    end

    println("=" ^ 80)

    # Summary
    println("\nSummary:")
    println("-" ^ 80)
    @printf("%-25s %-25s %-20s\n", "Function", "Smoothness", "Expected Rate")
    println("-" ^ 80)
    summaries = [
        ("|x|^(5/2)",       "C^2, 3rd deriv in BV",    "O(N^(-2.5))"),
        ("exp(-1/(1-x^2))", "C^inf, not analytic",      "faster than any N^(-k)"),
        ("tanh(5x)",        "Analytic",                  "O(rho^(-N))"),
        ("x^8",             "Polynomial",                "Exact for N >= 8"),
    ]
    for (name, smooth, rate) in summaries
        @printf("%-25s %-25s %-20s\n", name, smooth, rate)
    end
    println("-" ^ 80)
end

main()
