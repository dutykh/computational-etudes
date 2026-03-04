# convergence_rates.jl
#
# Demonstration of spectral differentiation convergence rates for functions
# of varying smoothness, illustrating Theorems 3 and 4.
#
# Three test functions on [0, 2pi]:
#     1. |sin(x)|^3       - Finite regularity, algebraic convergence O(N^{-3})
#     2. 1/(1+sin^2(x/2)) - Analytic in strip, geometric convergence O(c^{-N})
#     3. exp(sin(x))       - Entire function, super-geometric convergence
#
# The plot shows differentiation error vs N on a semi-log scale, demonstrating
# how function smoothness determines the convergence rate of spectral methods.
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
using LaTeXStrings

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
const NAVY   = colorant"#142D6E"
const SKY    = colorant"#7896D2"
const CORAL  = colorant"#E74C3C"
const TEAL   = colorant"#16A085"
const PURPLE = colorant"#9B59B6"

SCRIPT_DIR = @__DIR__
OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch06", "julia")
mkpath(OUTPUT_DIR)

# -----------------------------------------------------------------------------
# Test functions and their exact derivatives
# -----------------------------------------------------------------------------
"""Finite regularity: |sin(x)|^3"""
f1(x) = abs(sin(x))^3

"""Exact derivative of |sin(x)|^3: 3 sin(x) |sin(x)| cos(x)"""
f1_deriv(x) = 3 * sin(x) * abs(sin(x)) * cos(x)

"""Analytic in strip: 1/(1 + sin^2(x/2))"""
f2(x) = 1.0 / (1.0 + sin(x / 2.0)^2)

"""Exact derivative of 1/(1 + sin^2(x/2))"""
f2_deriv(x) = -sin(x / 2.0) * cos(x / 2.0) / (1.0 + sin(x / 2.0)^2)^2

"""Entire function: exp(sin(x))"""
f3(x) = exp(sin(x))

"""Exact derivative of exp(sin(x))"""
f3_deriv(x) = cos(x) * exp(sin(x))

# -----------------------------------------------------------------------------
# Spectral differentiation matrix (periodic)
# -----------------------------------------------------------------------------
"""
    spectral_diff_periodic(N)

Construct the periodic spectral differentiation matrix for `N` points.

Returns `(D, x)` where `D` is the N x N differentiation matrix and
`x` contains the grid points on [0, 2pi).
"""
function spectral_diff_periodic(N)
    h = 2pi / N
    x = h * collect(0:N-1)

    # Build differentiation matrix using the cotangent formula
    D = zeros(N, N)
    for i in 1:N
        for j in 1:N
            if i != j
                D[i, j] = 0.5 * ((-1)^(i - j)) / tan((i - j) * pi / N)
            end
        end
    end

    return D, x
end

# -----------------------------------------------------------------------------
# Convergence study
# -----------------------------------------------------------------------------
"""
    compute_errors(N_values)

Compute differentiation errors for each function at various `N`.
Returns three vectors of max-norm errors.
"""
function compute_errors(N_values)
    errors1 = Float64[]
    errors2 = Float64[]
    errors3 = Float64[]

    for N in N_values
        D, x = spectral_diff_periodic(N)

        # Function 1: |sin(x)|^3
        u1 = f1.(x)
        du1_exact = f1_deriv.(x)
        du1_spectral = D * u1
        push!(errors1, maximum(abs.(du1_spectral - du1_exact)))

        # Function 2: 1/(1 + sin^2(x/2))
        u2 = f2.(x)
        du2_exact = f2_deriv.(x)
        du2_spectral = D * u2
        push!(errors2, maximum(abs.(du2_spectral - du2_exact)))

        # Function 3: exp(sin(x))
        u3 = f3.(x)
        du3_exact = f3_deriv.(x)
        du3_spectral = D * u3
        push!(errors3, maximum(abs.(du3_spectral - du3_exact)))
    end

    return errors1, errors2, errors3
end

# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
function main()
    # Grid sizes to test
    N_values = [6, 8, 10, 12, 14, 16, 20, 24, 28, 32, 40, 48, 56, 64]

    # Compute errors
    errors1, errors2, errors3 = compute_errors(N_values)

    # Create figure
    fig = Figure(size = (600, 410))
    ax = Axis(fig[1, 1],
              xlabel = L"Number of grid points $N$",
              ylabel = L"Differentiation error $\Vert f' - Df\Vert_\infty$",
              title  = "Spectral Differentiation Convergence: Smoothness Matters",
              yscale = log10,
              xlabelpadding = 3, ylabelpadding = 3)

    # Plot errors
    scatterlines!(ax, N_values, errors1,
                  color = TEAL, linewidth = 1.5, markersize = 7, marker = :circle,
                  label = L"$|\sin(x)|^3$ (finite regularity)")
    scatterlines!(ax, N_values, errors2,
                  color = CORAL, linewidth = 1.5, markersize = 7, marker = :rect,
                  label = L"$1/(1+\sin^2(x/2))$ (analytic in strip)")
    scatterlines!(ax, N_values, errors3,
                  color = NAVY, linewidth = 1.5, markersize = 7, marker = :diamond,
                  label = L"$\exp(\sin(x))$ (entire)")

    # Add theoretical reference lines
    N_ref = range(8, 64, length=100)

    # O(N^{-3}) reference for finite regularity
    C1 = errors1[4] * N_values[4]^3  # index 4 corresponds to N=12
    lines!(ax, collect(N_ref), C1 ./ Float64.(collect(N_ref)).^3,
           linestyle = :dash, color = (TEAL, 0.5), linewidth = 1)
    text!(ax, 66, C1 / 66^3 * 1.5, text = L"O(N^{-3})", fontsize = 10,
          color = TEAL, align = (:left, :bottom))

    # Geometric decay reference for analytic strip
    # Fit the exponential decay rate
    idx_fit = findall(i -> N_values[i] >= 10 && N_values[i] <= 40 && errors2[i] > 1e-15,
                      1:length(N_values))
    if length(idx_fit) >= 2
        N_fit = Float64.(N_values[idx_fit])
        log_err = log.(errors2[idx_fit])
        # Linear fit: log(error) = a + b*N
        A = hcat(ones(length(N_fit)), N_fit)
        coeffs = A \ log_err
        decay_rate = -coeffs[2]
        C2 = exp(coeffs[1])
        lines!(ax, collect(N_ref), C2 .* exp.(-decay_rate .* collect(N_ref)),
               linestyle = :dash, color = (CORAL, 0.5), linewidth = 1)
        text!(ax, 30, C2 * exp(-decay_rate * 30) * 3,
              text = latexstring("O(e^{-$(round(decay_rate, digits=2)) N})"),
              fontsize = 10, color = CORAL)
    end

    # Machine epsilon line
    hlines!(ax, [2.2e-16], color = :gray, linestyle = :dot, linewidth = 1)
    text!(ax, 66, 4e-16, text = "Machine precision", fontsize = 9, color = :gray,
          align = (:right, :bottom))

    # Add annotation for super-geometric convergence
    # Index 7 corresponds to N=20
    text!(ax, 28, errors3[7],
          text = "Super-geometric\nconvergence",
          fontsize = 10, color = NAVY, align = (:left, :center))

    # Axis settings
    xlims!(ax, 0, 70)
    ylims!(ax, 1e-16, 1e1)

    # Legend
    axislegend(ax, position = :rt, framevisible = true,
               framecolor = :transparent, bgcolor = RGBAf(1, 1, 1, 0.9),
               padding = (8, 8, 4, 4))

    # Clean styling
    hidespines!(ax, :t, :r)
    ax.leftspinecolor   = RGBf(0.267, 0.267, 0.267)
    ax.bottomspinecolor = RGBf(0.267, 0.267, 0.267)
    ax.xtickcolor = RGBf(0.267, 0.267, 0.267)
    ax.ytickcolor = RGBf(0.267, 0.267, 0.267)

    # Grid
    ax.ygridvisible = true
    ax.ygridcolor   = RGBAf(0.533, 0.533, 0.533, 0.2)
    ax.ygridstyle   = :solid
    ax.xgridvisible = true
    ax.xgridcolor   = RGBAf(0.533, 0.533, 0.533, 0.2)
    ax.xgridstyle   = :solid

    # Save figure
    save(joinpath(OUTPUT_DIR, "convergence_rates.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "convergence_rates.png"), fig, px_per_unit = 4)
    println("Figure saved to: $(joinpath(OUTPUT_DIR, "convergence_rates.pdf"))")

    # Print results table
    println("\nConvergence Study Results:")
    println("=" ^ 65)
    println(lpad("N", 4), "  ", lpad("|sin(x)|^3", 14), "  ",
            lpad("1/(1+sin^2)", 14), "  ", lpad("exp(sin)", 14))
    println("-" ^ 65)
    for (i, N) in enumerate(N_values)
        @printf("%4d  %14.6e  %14.6e  %14.6e\n", N, errors1[i], errors2[i], errors3[i])
    end
    println("=" ^ 65)
end

main()
