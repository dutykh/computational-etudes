# convergence_comparison.jl
#
# Compares the convergence of polynomial interpolation on equispaced vs
# Chebyshev nodes for the Runge function.
#
# The interpolation error bound is:
#     ||f - p_N||_inf <= (1 + Lambda_N) E_N(f)
#
# where:
#     - Lambda_N is the Lebesgue constant
#     - E_N(f) is the best polynomial approximation error
#
# For the Runge function f(x) = 1/(1 + 25x^2) with poles at +/-0.2i:
#     - Chebyshev: error ~ O(rho^{-N}) where rho ~ 1.04
#     - Equispaced: diverges due to 2^N growth of Lambda_N
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
using LaTeXStrings
using Statistics

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
const N_FINE = 2000

# Book color scheme
const NAVY  = colorant"#142D6E"
const SKY   = colorant"#7896D2"
const CORAL = colorant"#E74C3C"
const TEAL  = colorant"#16A085"

# Output path
const SCRIPT_DIR = @__DIR__
const OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch04", "julia")

# -----------------------------------------------------------------------------
# Functions
# -----------------------------------------------------------------------------
runge(x) = 1.0 / (1.0 + 25.0 * x^2)
equispaced_nodes(N) = collect(range(-1, 1, length=N+1))
chebyshev_nodes(N) = [cos(j * π / N) for j in 0:N]

function lagrange_interpolate(x_nodes, f_nodes, x_eval)
    n = length(x_nodes)
    p_eval = zeros(length(x_eval))

    for k in 1:n
        L_k = ones(length(x_eval))
        for j in 1:n
            if j != k
                L_k .*= (x_eval .- x_nodes[j]) ./ (x_nodes[k] - x_nodes[j])
            end
        end
        p_eval .+= f_nodes[k] .* L_k
    end

    return p_eval
end

function max_interpolation_error(node_func, N; n_eval=N_FINE)
    x_eval = collect(range(-1, 1, length=n_eval))
    x_nodes = node_func(N)
    f_nodes = runge.(x_nodes)
    p_interp = lagrange_interpolate(x_nodes, f_nodes, x_eval)
    return maximum(abs.(runge.(x_eval) .- p_interp))
end

# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
function main()
    # Range of polynomial degrees
    N_values = 2:2:50

    # Compute errors
    errors_equi = [max_interpolation_error(equispaced_nodes, N) for N in N_values]
    errors_cheb = [max_interpolation_error(chebyshev_nodes, N) for N in N_values]

    # Create figure
    fig = Figure(size = (800, 550))
    ax = Axis(fig[1, 1],
              xlabel = L"N" * " (polynomial degree)",
              ylabel = L"\Vert f - p_N\Vert_\infty" * " (max interpolation error)",
              title = "Convergence: Equispaced vs Chebyshev Interpolation",
              limits = ((0, 52), (1e-6, 1e4)),
              yscale = log10)

    # Plot errors
    scatterlines!(ax, collect(N_values), errors_equi, color = CORAL,
                  linewidth = 1.5, markersize = 8, marker = :circle,
                  label = "Equispaced nodes")
    scatterlines!(ax, collect(N_values), errors_cheb, color = TEAL,
                  linewidth = 1.5, markersize = 8, marker = :rect,
                  label = "Chebyshev nodes")

    # Add theoretical convergence rate for Chebyshev
    rho = abs(0.2im + sqrt(complex(-1.04)))
    N_theory = collect(range(10, 50, length=100))

    # Fit the constant
    idx_fit = collect(N_values) .>= 20
    if any(idx_fit)
        C_fit = mean(errors_cheb[idx_fit] .* rho .^ collect(N_values)[idx_fit])
        lines!(ax, N_theory, C_fit .* rho .^ (-N_theory),
               linestyle = :dash, color = TEAL, linewidth = 1.5, alpha = 0.5,
               label = latexstring("Theory: O(\\rho^{-N}), \\rho \\approx $(round(rho, digits=2))"))
    end

    # Reference line
    hlines!(ax, 1, color = colorant"#888888", linewidth = 0.5, linestyle = :dot,
            label = "Error = 1 (useless)")

    # Annotations
    text!(ax, 30.0, 1e3, text = "Divergence!", fontsize = 10, color = CORAL)
    text!(ax, 32.0, 1e-4, text = "Convergence", fontsize = 10, color = TEAL)

    # Formula box
    text!(ax, 0.05, 0.05, text = L"f(x) = \frac{1}{1 + 25x^2}",
          fontsize = 12, space = :relative, align = (:left, :bottom))

    axislegend(ax, position = :rt, framevisible = false, bgcolor = (:white, 0.9))

    # Save figure
    mkpath(OUTPUT_DIR)
    save(joinpath(OUTPUT_DIR, "convergence_comparison.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "convergence_comparison.png"), fig, px_per_unit = 4)

    println("Figure saved to: ", joinpath(OUTPUT_DIR, "convergence_comparison.pdf"))

    # Print summary table
    println("\nConvergence Comparison (Runge function):")
    println("-" ^ 55)
    println("$(rpad("N", 4)) $(lpad("Equispaced", 15)) $(lpad("Chebyshev", 15)) $(lpad("Ratio", 12))")
    println("-" ^ 55)
    for i in 1:5:length(N_values)
        N = N_values[i]
        err_e = errors_equi[i]
        err_c = errors_cheb[i]
        ratio = err_c > 0 ? err_e / err_c : Inf
        println("$(lpad(N, 4)) $(lpad(string(round(err_e, sigdigits=2)), 15)) $(lpad(string(round(err_c, sigdigits=2)), 15)) $(lpad(string(round(ratio, sigdigits=2)), 12))")
    end
    println("-" ^ 55)

    # Analysis
    println("\nAnalysis:")
    println("  Runge function poles at z = +/-0.2i")
    println("  Chebyshev convergence parameter: rho = $(round(rho, digits=4))")
    println("  Expected rate: error ~ rho^{-N} = $(round(1/rho, digits=4))^N")

    # Verify convergence rate
    if length(errors_cheb) > 5
        log_errors = log.(errors_cheb[6:end])
        N_fit = collect(N_values)[6:end]
        # Linear fit in log space
        A = hcat(ones(length(N_fit)), Float64.(N_fit))
        coeffs = A \ log_errors
        rho_estimated = exp(-coeffs[2])
        println("  Estimated from data: rho = $(round(rho_estimated, digits=4))")
    end

    return fig
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
