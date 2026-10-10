# lebesgue_functions.jl
#
# Visualizes and compares Lebesgue functions for different node distributions.
#
# The Lebesgue function is:
#     Lambda_N(x) = sum_{k=0}^{N} |L_k(x)|
#
# where L_k(x) are the Lagrange basis polynomials. The Lebesgue constant is:
#     Lambda_N = max_{x in [-1,1]} Lambda_N(x)
#
# The Lebesgue constant bounds the interpolation error:
#     ||f - p_N||_inf <= (1 + Lambda_N) E_N(f)
#
# Asymptotic growth rates:
#     - Chebyshev: Lambda_N = (2/pi) ln(N) + O(1)
#     - Legendre:  Lambda_N ~ O(sqrt(N))
#     - Equispaced: Lambda_N ~ 2^N / (N ln N) (exponential!)
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
const N_FINE = 1000

# Book color scheme
const NAVY  = colorant"#142D6E"
const SKY   = colorant"#7896D2"
const CORAL = colorant"#E74C3C"
const TEAL  = colorant"#16A085"

# Output path
const SCRIPT_DIR = @__DIR__
const OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch04", "julia")

# -----------------------------------------------------------------------------
# Node generation
# -----------------------------------------------------------------------------
equispaced_nodes(N) = collect(range(-1, 1, length=N+1))
chebyshev_nodes(N) = [cos(j * π / N) for j in 0:N]

"""
    legendre_nodes(N)

Legendre-Gauss-Lobatto nodes on [-1, 1].
These are +/-1 plus the N-1 zeros of P'_N(x).
"""
function legendre_nodes(N)
    N == 0 && return [0.0]
    N == 1 && return [-1.0, 1.0]

    # Build Legendre polynomial P_N in the Legendre basis
    # Coefficients: all zeros except the N-th is 1
    coeffs = zeros(N + 1)
    coeffs[N + 1] = 1.0

    # Convert to standard polynomial using recurrence relation
    # P_0(x) = 1, P_1(x) = x, (n+1)P_{n+1}(x) = (2n+1)x P_n(x) - n P_{n-1}(x)
    P_prev = [1.0]       # P_0
    P_curr = [0.0, 1.0]  # P_1

    for n in 1:N-1
        P_next = ((2n + 1) .* [0; P_curr] .- n .* [P_prev; zeros(2)]) ./ (n + 1)
        # Pad P_prev for next iteration
        P_prev = P_curr
        P_curr = P_next
    end

    # P_curr is now P_N(x) in monomial basis
    # Compute derivative
    P_N_poly = Polynomial(P_curr)
    P_N_deriv = derivative(P_N_poly)

    # Find roots of P'_N
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
    # Fine grid for plotting
    x_fine = collect(range(-1, 1, length=N_FINE))

    # Create figure with two panels
    fig = Figure(size = (1100, 450))

    # ==========================================================================
    # Left panel: Lebesgue functions for N = 10
    # ==========================================================================
    N_plot = 10

    x_equi = equispaced_nodes(N_plot)
    x_cheb = chebyshev_nodes(N_plot)
    x_leg = legendre_nodes(N_plot)

    Lambda_equi = lebesgue_function(x_equi, x_fine)
    Lambda_cheb = lebesgue_function(x_cheb, x_fine)
    Lambda_leg = lebesgue_function(x_leg, x_fine)

    ax1 = Axis(fig[1, 1],
               xlabel = L"x", ylabel = L"\Lambda_N(x)",
               title = "Lebesgue Functions (N = $N_plot)",
               limits = ((-1, 1), (0, 35)))

    lines!(ax1, x_fine, Lambda_equi, color = CORAL, linewidth = 1.5,
           label = latexstring("Equispaced (\\Lambda_{$N_plot} = $(round(maximum(Lambda_equi), digits=1)))"))
    lines!(ax1, x_fine, Lambda_leg, color = TEAL, linewidth = 1.5,
           label = latexstring("Legendre (\\Lambda_{$N_plot} = $(round(maximum(Lambda_leg), digits=2)))"))
    lines!(ax1, x_fine, Lambda_cheb, color = SKY, linewidth = 1.5,
           label = latexstring("Chebyshev (\\Lambda_{$N_plot} = $(round(maximum(Lambda_cheb), digits=2)))"))

    hlines!(ax1, 1, color = colorant"#888888", linewidth = 0.5, linestyle = :dash)

    # Annotation
    text!(ax1, -0.6, 25.0, text = "Peaks at boundaries\nfor equispaced nodes",
          fontsize = 9, color = :gray)

    axislegend(ax1, position = :ct, framevisible = false, labelsize = 9)

    # ==========================================================================
    # Right panel: Lebesgue constant growth with N
    # ==========================================================================
    N_values = 2:30

    Lambda_equi_values = [lebesgue_constant(equispaced_nodes(n)) for n in N_values]
    Lambda_cheb_values = [lebesgue_constant(chebyshev_nodes(n)) for n in N_values]
    Lambda_leg_values = [lebesgue_constant(legendre_nodes(n)) for n in N_values]

    ax2 = Axis(fig[1, 2],
               xlabel = L"N" * " (polynomial degree)",
               ylabel = L"\Lambda_N" * " (Lebesgue constant)",
               title = "Growth of Lebesgue Constants",
               limits = ((0, 32), nothing),
               yscale = log10)

    scatterlines!(ax2, collect(N_values), Lambda_equi_values,
                  color = CORAL, linewidth = 1.5, markersize = 6, marker = :circle,
                  label = "Equispaced")
    scatterlines!(ax2, collect(N_values), Lambda_leg_values,
                  color = TEAL, linewidth = 1.5, markersize = 6, marker = :rect,
                  label = "Legendre")
    scatterlines!(ax2, collect(N_values), Lambda_cheb_values,
                  color = SKY, linewidth = 1.5, markersize = 6, marker = :utriangle,
                  label = "Chebyshev")

    # Plot asymptotic curve for Chebyshev
    N_asy = range(5, 30, length=100)
    Lambda_cheb_asy = (2/π) .* log.(N_asy) .+ 0.6
    lines!(ax2, collect(N_asy), Lambda_cheb_asy, linestyle = :dash, color = SKY,
           linewidth = 1, alpha = 0.5, label = L"(2/\pi)\ln N")

    # Growth rate annotations
    text!(ax2, 25.0, 1e6, text = L"O(2^N)", fontsize = 10, color = CORAL)
    text!(ax2, 25.0, 4.0, text = L"O(\sqrt{N})", fontsize = 10, color = TEAL)
    text!(ax2, 25.0, 2.5, text = L"O(\ln N)", fontsize = 10, color = SKY)

    axislegend(ax2, position = :lt, framevisible = false, labelsize = 9)

    # Save figure
    mkpath(OUTPUT_DIR)
    save(joinpath(OUTPUT_DIR, "lebesgue_functions.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "lebesgue_functions.png"), fig, px_per_unit = 4)

    println("Figure saved to: ", joinpath(OUTPUT_DIR, "lebesgue_functions.pdf"))

    # Print summary table
    println("\nLebesgue Constants Comparison:")
    println("-" ^ 65)
    println("$(rpad("N", 4)) $(lpad("Equispaced", 15)) $(lpad("Legendre", 15)) $(lpad("Chebyshev", 15))")
    println("-" ^ 65)
    for i in 1:5:length(N_values)
        n = N_values[i]
        println("$(lpad(n, 4)) $(lpad(round(Lambda_equi_values[i], digits=2), 15)) $(lpad(round(Lambda_leg_values[i], digits=4), 15)) $(lpad(round(Lambda_cheb_values[i], digits=4), 15))")
    end
    println("-" ^ 65)

    println("\nAsymptotic Growth Rates:")
    println("  Chebyshev:  Lambda_N ~ (2/pi) ln(N) + O(1)")
    println("  Legendre:   Lambda_N ~ O(sqrt(N))")
    println("  Equispaced: Lambda_N ~ 2^(N+1) / (eN ln N)")

    return fig
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
