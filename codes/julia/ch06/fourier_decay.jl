# fourier_decay.jl
#
# Demonstration of the decay hierarchy for Fourier coefficients, illustrating
# Theorem 1: the connection between function smoothness and spectral decay.
#
# Three test functions on [0, 2pi]:
#     1. |sin(x)|^3    - Finite regularity (C^2 but not C^3)
#                        Decay rate: O(k^{-4})
#     2. 1/(1+sin^2(x/2)) - Analytic in a strip in the complex plane
#                        Decay rate: O(c^{-|k|}) geometric
#     3. exp(sin(x))   - Entire function (analytic everywhere)
#                        Decay rate: faster than any exponential
#
# The plot shows |f_hat_k| vs k on a log scale, revealing:
#     - Algebraic decay appears as a straight line on log-log scale
#     - Geometric decay appears as a straight line on semi-log scale
#     - Super-geometric decay curves downward on semi-log scale
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
using FFTW
using LaTeXStrings
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
const NAVY   = colorant"#142D6E"
const SKY    = colorant"#7896D2"
const CORAL  = colorant"#E74C3C"
const TEAL   = colorant"#16A085"
const PURPLE = colorant"#9B59B6"

SCRIPT_DIR = @__DIR__
OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch06", "julia")
mkpath(OUTPUT_DIR)

# -----------------------------------------------------------------------------
# Test functions
# -----------------------------------------------------------------------------
"""
    f_finite_regularity(x)

|sin(x)|^3 - has bounded third derivative, fourth derivative in BV.
Decay rate: O(k^{-4})
"""
f_finite_regularity(x) = abs(sin(x))^3

"""
    f_analytic_strip(x)

1/(1 + sin^2(x/2)) - analytic in a strip |Im(z)| < a for some a > 0.
Has poles where sin^2(z/2) = -1, i.e., at complex values.
Decay rate: O(c^{-|k|}) geometric
"""
f_analytic_strip(x) = 1.0 / (1.0 + sin(x / 2.0)^2)

"""
    f_entire(x)

exp(sin(x)) - entire function (analytic throughout the complex plane).
Decay rate: faster than any exponential (super-geometric)
"""
f_entire(x) = exp(sin(x))

# -----------------------------------------------------------------------------
# Fourier coefficient computation
# -----------------------------------------------------------------------------
"""
    compute_fourier_coefficients(f, N)

Compute the Fourier coefficients of a 2pi-periodic function `f`
sampled at `N` equispaced points.

Returns `(k, f_hat_abs)` where `k` contains wavenumbers 0, 1, ..., N/2
and `f_hat_abs` the absolute values of the corresponding Fourier coefficients.
"""
function compute_fourier_coefficients(f, N)
    x = range(0, 2pi, length=N+1)[1:N]  # endpoint=false
    f_vals = f.(x)

    # Compute FFT and normalize
    f_hat = fft(f_vals) / N

    # Take absolute values and keep only positive frequencies
    f_hat_abs = abs.(f_hat[1:N÷2+1])

    # Wavenumbers
    k = 0:N÷2

    return collect(k), f_hat_abs
end

# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
function main()
    N = 128  # Use enough points to see the decay clearly

    # Compute Fourier coefficients for each function
    k1, f1_hat = compute_fourier_coefficients(f_finite_regularity, N)
    k2, f2_hat = compute_fourier_coefficients(f_analytic_strip, N)
    k3, f3_hat = compute_fourier_coefficients(f_entire, N)

    # Create figure
    fig = Figure(size = (600, 410))
    ax = Axis(fig[1, 1],
              xlabel = L"Wavenumber $k$",
              ylabel = L"Fourier coefficient magnitude $|\hat{f}_k|$",
              title  = "The Decay Hierarchy: Smoothness Determines Spectral Decay",
              yscale = log10,
              xlabelpadding = 3,
              ylabelpadding = 3)

    # Plot coefficients (skip k=0 for cleaner visualization on log scale)
    scatterlines!(ax, k1[2:end], f1_hat[2:end],
                  color = TEAL, linewidth = 1.5, markersize = 5, marker = :circle,
                  label = L"$|\sin(x)|^3$ (finite regularity)")
    scatterlines!(ax, k2[2:end], f2_hat[2:end],
                  color = CORAL, linewidth = 1.5, markersize = 5, marker = :rect,
                  label = L"$1/(1+\sin^2(x/2))$ (analytic in strip)")
    scatterlines!(ax, k3[2:end], f3_hat[2:end],
                  color = NAVY, linewidth = 1.5, markersize = 5, marker = :diamond,
                  label = L"$\exp(\sin(x))$ (entire)")

    # Add theoretical decay rate references

    # O(k^{-4}) reference for finite regularity
    k_ref_alg = 5:59
    C1 = f1_hat[11] * 10^4  # index 11 corresponds to k=10
    lines!(ax, collect(k_ref_alg), C1 ./ Float64.(k_ref_alg).^4,
           linestyle = :dash, color = (TEAL, 0.5), linewidth = 1)
    text!(ax, 55, C1 / 55^4 * 2, text = L"O(k^{-4})", fontsize = 10, color = TEAL)

    # Geometric decay reference for analytic strip
    # Find decay rate by fitting
    k_fit = 10:24
    log_f2 = log.(f2_hat[k_fit .+ 1] .+ 1e-20)  # +1 for 1-based indexing
    # Linear fit: log(f_hat) = a + b*k
    A = hcat(ones(length(k_fit)), Float64.(collect(k_fit)))
    coeffs = A \ log_f2
    decay_rate = -coeffs[2]
    # Only draw reference line in the visible range (above machine precision)
    k_ref_exp = 5:34
    lines!(ax, collect(k_ref_exp),
           exp(coeffs[1]) .* exp.(-decay_rate .* Float64.(collect(k_ref_exp))),
           linestyle = :dash, color = (CORAL, 0.6), linewidth = 1.2)
    # Position label near the data
    text!(ax, 18, exp(coeffs[1]) * exp(-decay_rate * 18) * 2.5,
          text = latexstring("O(e^{-$(round(decay_rate, digits=2)) k})"),
          fontsize = 10, color = CORAL)

    # Machine epsilon line
    hlines!(ax, [1e-16], color = :gray, linestyle = :dot, linewidth = 1)
    text!(ax, 62, 2e-16, text = "Machine precision", fontsize = 9, color = :gray,
          align = (:right, :bottom))

    # Axis settings
    xlims!(ax, 0, 65)
    ylims!(ax, 1e-17, 1e1)

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
    save(joinpath(OUTPUT_DIR, "decay_hierarchy.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "decay_hierarchy.png"), fig, px_per_unit = 4)
    println("Figure saved to: $(joinpath(OUTPUT_DIR, "decay_hierarchy.pdf"))")

    # Print coefficient values for verification
    println("\nFourier Coefficient Magnitudes:")
    println("=" ^ 60)
    println(lpad("k", 4), "  ", lpad("|sin(x)|^3", 14), "  ",
            lpad("1/(1+sin^2)", 14), "  ", lpad("exp(sin)", 14))
    println("-" ^ 60)
    for k in [1, 2, 5, 10, 20, 30, 40, 50]
        if k < length(f1_hat)
            Printf.@printf("%4d  %14.6e  %14.6e  %14.6e\n",
                           k, f1_hat[k+1], f2_hat[k+1], f3_hat[k+1])
        end
    end
    println("=" ^ 60)
end

main()
