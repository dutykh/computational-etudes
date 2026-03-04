# convergence_zoom.jl
#
# Zoomed view of Chebyshev interpolation convergence for the Runge function,
# without the divergent equispaced curve.
#
# This allows better visualization of the geometric convergence rate O(rho^{-N})
# where rho ~ 1.22 for the Runge function.
#
# The Runge function f(x) = 1/(1 + 25x^2) has poles at +/-0.2i, giving:
#     rho = |0.2i + sqrt((0.2i)^2 - 1)| ~ 1.22
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
    # Range of polynomial degrees (extended for better view)
    N_values = 2:2:70

    # Compute errors for Chebyshev only
    errors_cheb = [max_interpolation_error(chebyshev_nodes, N) for N in N_values]

    # Theoretical convergence parameter
    rho = abs(0.2im + sqrt(complex(-1.04)))

    # Create figure
    fig = Figure(size = (700, 500))
    ax = Axis(fig[1, 1],
              xlabel = L"N" * " (polynomial degree)",
              ylabel = L"\Vert f - p_N\Vert_\infty" * " (max interpolation error)",
              title = "Chebyshev Interpolation Convergence (Runge Function)",
              limits = ((0, 72), (1e-16, 1e1)),
              yscale = log10)

    # Plot Chebyshev errors
    scatterlines!(ax, collect(N_values), errors_cheb, color = TEAL,
                  linewidth = 1.5, markersize = 8, marker = :rect,
                  label = "Chebyshev interpolation")

    # Add theoretical convergence rate
    N_theory = collect(range(10, 70, length=100))
    idx_fit = collect(N_values) .>= 20
    C_fit = mean(errors_cheb[idx_fit] .* rho .^ collect(N_values)[idx_fit])
    lines!(ax, N_theory, C_fit .* rho .^ (-N_theory),
           linestyle = :dash, color = NAVY, linewidth = 2, alpha = 0.7,
           label = latexstring("Theory: C \\cdot \\rho^{-N}, \\rho \\approx $(round(rho, digits=2))"))

    # Machine precision line
    hlines!(ax, 1e-15, color = colorant"#888888", linewidth = 0.8, linestyle = :dot,
            label = "Machine precision")

    # Annotation for convergence rate
    text!(ax, 35.0, 1e-6,
          text = latexstring("O(\\rho^{-N}) \\text{ with } \\rho \\approx $(round(rho, digits=2))"),
          fontsize = 11, color = TEAL)

    # Formula box
    text!(ax, 0.05, 0.05,
          text = L"f(x) = \frac{1}{1 + 25x^2}" * "\nPoles at " * L"z = \pm 0.2i",
          fontsize = 11, space = :relative, align = (:left, :bottom))

    axislegend(ax, position = :rt, framevisible = false, bgcolor = (:white, 0.9))

    # Save figure
    mkpath(OUTPUT_DIR)
    save(joinpath(OUTPUT_DIR, "convergence_zoom.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "convergence_zoom.png"), fig, px_per_unit = 4)

    println("Figure saved to: ", joinpath(OUTPUT_DIR, "convergence_zoom.pdf"))

    # Print summary table
    println("\nChebyshev Convergence for Runge Function:")
    println("-" ^ 50)
    println("$(rpad("N", 4)) $(lpad("Error", 15)) $(lpad("Theoretical", 15))")
    println("-" ^ 50)
    for i in 1:5:length(N_values)
        N = N_values[i]
        err = errors_cheb[i]
        theory = C_fit * rho^(-N)
        println("$(lpad(N, 4)) $(lpad(string(round(err, sigdigits=2)), 15)) $(lpad(string(round(theory, sigdigits=2)), 15))")
    end
    println("-" ^ 50)

    # Analysis
    println("\nAnalysis:")
    println("  Runge function poles at z = +/-0.2i")
    println("  Theoretical rho = $(round(rho, digits=6))")

    # Estimate rho from the data
    log_errors = log.(errors_cheb[6:end])
    N_fit = collect(N_values)[6:end]
    A = hcat(ones(length(N_fit)), Float64.(N_fit))
    coeffs = A \ log_errors
    rho_estimated = exp(-coeffs[2])
    println("  Estimated from data: rho = $(round(rho_estimated, digits=6))")
    println("  Error reduction per degree: factor of $(round(rho_estimated, digits=3))")

    return fig
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
