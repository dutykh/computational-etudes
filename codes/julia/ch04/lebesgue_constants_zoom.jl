# lebesgue_constants_zoom.jl
#
# Zoomed comparison of Lebesgue constant growth for Legendre and Chebyshev nodes.
#
# This figure omits the exponentially growing equispaced curve to allow better
# visualization of the difference between:
#     - Chebyshev: Lambda_N = (2/pi) ln(N) + O(1)
#     - Legendre:  Lambda_N ~ O(sqrt(N))
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
using Polynomials

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
# Book color scheme
const NAVY = colorant"#142D6E"
const SKY  = colorant"#7896D2"
const TEAL = colorant"#16A085"

# Output path
const SCRIPT_DIR = @__DIR__
const OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch04", "julia")

# -----------------------------------------------------------------------------
# Node generation
# -----------------------------------------------------------------------------
chebyshev_nodes(N) = [cos(j * π / N) for j in 0:N]

function legendre_nodes(N)
    N == 0 && return [0.0]
    N == 1 && return [-1.0, 1.0]

    # Build P_N via recurrence
    P_prev = [1.0]
    P_curr = [0.0, 1.0]

    for n in 1:N-1
        P_next = ((2n + 1) .* [0; P_curr] .- n .* [P_prev; zeros(2)]) ./ (n + 1)
        P_prev = P_curr
        P_curr = P_next
    end

    P_N_poly = Polynomial(P_curr)
    P_N_deriv = derivative(P_N_poly)

    interior = roots(P_N_deriv)
    interior = real.(filter(r -> abs(imag(r)) < 1e-10, interior))

    return sort(vcat([-1.0], interior, [1.0]))
end

# -----------------------------------------------------------------------------
# Lebesgue function computation
# -----------------------------------------------------------------------------
function lagrange_basis(x_nodes, k, x_eval)
    n = length(x_nodes)
    L_k = ones(length(x_eval))
    for j in 1:n
        if j != k
            L_k .*= (x_eval .- x_nodes[j]) ./ (x_nodes[k] - x_nodes[j])
        end
    end
    return L_k
end

function lebesgue_function(x_nodes, x_eval)
    N = length(x_nodes) - 1
    Lambda = zeros(length(x_eval))
    for k in 1:N+1
        Lambda .+= abs.(lagrange_basis(x_nodes, k, x_eval))
    end
    return Lambda
end

function lebesgue_constant(x_nodes; n_eval=2000)
    x_eval = collect(range(-1, 1, length=n_eval))
    Lambda = lebesgue_function(x_nodes, x_eval)
    return maximum(Lambda)
end

# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
function main()
    # Create figure (single panel)
    fig = Figure(size = (700, 500))
    ax = Axis(fig[1, 1],
              xlabel = L"N" * " (polynomial degree)",
              ylabel = L"\Lambda_N" * " (Lebesgue constant)",
              title = "Lebesgue Constants: Legendre vs Chebyshev",
              limits = ((0, 52), (0, 8)))

    # Compute Lebesgue constants for range of N
    N_values = 2:50

    Lambda_cheb_values = [lebesgue_constant(chebyshev_nodes(n)) for n in N_values]
    Lambda_leg_values = [lebesgue_constant(legendre_nodes(n)) for n in N_values]

    # Plot computed values
    scatterlines!(ax, collect(N_values), Lambda_leg_values,
                  color = TEAL, linewidth = 1.5, markersize = 8, marker = :rect,
                  label = "Legendre")
    scatterlines!(ax, collect(N_values), Lambda_cheb_values,
                  color = SKY, linewidth = 1.5, markersize = 8, marker = :utriangle,
                  label = "Chebyshev")

    # Plot asymptotic curves
    N_asy = collect(range(3, 50, length=100))

    # Chebyshev: (2/pi) ln(N) + constant
    Lambda_cheb_asy = (2/π) .* log.(N_asy) .+ 0.6
    lines!(ax, N_asy, Lambda_cheb_asy, linestyle = :dash, color = SKY,
           linewidth = 1.5, alpha = 0.6, label = L"(2/\pi)\ln N + 0.6")

    # Legendre: c * sqrt(N) with fitted constant
    c_leg = Lambda_leg_values[end] / sqrt(N_values[end])
    Lambda_leg_asy = c_leg .* sqrt.(N_asy)
    lines!(ax, N_asy, Lambda_leg_asy, linestyle = :dash, color = TEAL,
           linewidth = 1.5, alpha = 0.6,
           label = latexstring("$(round(c_leg, digits=2))\\sqrt{N}"))

    # Growth rate annotations with arrows
    text!(ax, 38.0, 6.5, text = L"O(\sqrt{N})", fontsize = 11, color = TEAL)
    text!(ax, 38.0, 2.0, text = L"O(\ln N)", fontsize = 11, color = SKY)

    axislegend(ax, position = :lt, framevisible = false, labelsize = 10)

    # Save figure
    mkpath(OUTPUT_DIR)
    save(joinpath(OUTPUT_DIR, "lebesgue_constants_zoom.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "lebesgue_constants_zoom.png"), fig, px_per_unit = 4)

    println("Figure saved to: ", joinpath(OUTPUT_DIR, "lebesgue_constants_zoom.pdf"))

    # Print comparison table
    println("\nLebesgue Constants Comparison (Legendre vs Chebyshev):")
    println("-" ^ 50)
    println("$(rpad("N", 4)) $(lpad("Legendre", 15)) $(lpad("Chebyshev", 15)) $(lpad("Ratio", 10))")
    println("-" ^ 50)
    for i in 1:10:length(N_values)
        n = N_values[i]
        ratio = Lambda_leg_values[i] / Lambda_cheb_values[i]
        println("$(lpad(n, 4)) $(lpad(round(Lambda_leg_values[i], digits=4), 15)) $(lpad(round(Lambda_cheb_values[i], digits=4), 15)) $(lpad(round(ratio, digits=2), 10))")
    end
    println("-" ^ 50)

    return fig
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
