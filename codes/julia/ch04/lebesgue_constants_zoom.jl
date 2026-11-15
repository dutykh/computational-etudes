# lebesgue_constants_zoom.jl
#
# Zoomed comparison of Lebesgue constant growth for Legendre-Gauss-Lobatto and
# Chebyshev-Gauss-Lobatto nodes.
#
# This figure omits the exponentially growing equispaced curve, so that the two
# logarithmic growths can be told apart:
#     - Chebyshev-Gauss-Lobatto: Lambda_N = (2/pi)(ln N + gamma + ln(8/pi)) + o(1)
#     - Legendre-Gauss-Lobatto:  Lambda_N = O(ln N), close to (2/pi) ln(N+1) + 0.685,
#       the estimate of Hesthaven drawn from numerical experiments
# The Legendre-Gauss-Lobatto nodes are the Fekete points of [-1, 1] (Fejer),
# and their constant stays about 0.28 below the Chebyshev one. Neither family
# is optimal: every node set has Lambda_N > (2/pi) ln(N+1) + 0.5212 (Brutman,
# 1978), and both curves stay a little above this bound.
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
using LinearAlgebra

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

"""
    legendre_nodes(N)

Legendre-Gauss-Lobatto nodes on [-1, 1]: the endpoints and the N - 1 zeros of
P_N'. These zeros are those of the Gegenbauer polynomial C^(3/2)_(N-1), hence
the eigenvalues of its symmetric tridiagonal Jacobi matrix, with zero diagonal
and off-diagonal entries sqrt(k(k + 2) / ((2k + 1)(2k + 3))), k = 1, ..., N - 2
(Golub-Welsch). Unlike the roots of P_N' in the monomial basis, whose
coefficients grow exponentially with N, this eigenvalue problem is perfectly
conditioned and gives the nodes to machine precision for every N.
"""
function legendre_nodes(N)
    N == 0 && return [0.0]
    N == 1 && return [-1.0, 1.0]

    k = 1:N-2
    offdiag = sqrt.(k .* (k .+ 2) ./ ((2 .* k .+ 1) .* (2 .* k .+ 3)))
    interior = eigvals(SymTridiagonal(zeros(N - 1), collect(Float64, offdiag)))

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
              xlabel = L"N\text{ (polynomial degree)}",
              ylabel = L"\Lambda_N\text{ (Lebesgue constant)}",
              title = "Lebesgue Constants: Legendre vs Chebyshev",
              limits = ((0, 52), (1.0, 4.4)))

    # Compute Lebesgue constants for range of N
    N_values = 2:50

    Lambda_cheb_values = [lebesgue_constant(chebyshev_nodes(n)) for n in N_values]
    Lambda_leg_values = [lebesgue_constant(legendre_nodes(n)) for n in N_values]

    # Plot computed values
    scatterlines!(ax, collect(N_values), Lambda_leg_values,
                  color = TEAL, linewidth = 1.5, markersize = 8, marker = :rect,
                  label = "Legendre–Gauss–Lobatto")
    scatterlines!(ax, collect(N_values), Lambda_cheb_values,
                  color = SKY, linewidth = 1.5, markersize = 8, marker = :utriangle,
                  label = "Chebyshev–Gauss–Lobatto")

    # Logarithmic asymptotics of the two families
    N_asy = collect(range(2, 50, length=200))

    # Chebyshev-Gauss-Lobatto: (2/pi)(ln N + gamma + ln(8/pi)) + o(1)
    Lambda_cheb_asy = (2/π) .* (log.(N_asy) .+ Base.MathConstants.eulergamma .+ log(8/π))
    lines!(ax, N_asy, Lambda_cheb_asy, linestyle = :dash, color = SKY,
           linewidth = 1.5, alpha = 0.6,
           label = L"(2/\pi)(\ln N + \gamma + \ln(8/\pi))")

    # Legendre-Gauss-Lobatto: the estimate (2/pi) ln(N+1) + 0.685 of Hesthaven
    Lambda_leg_asy = (2/π) .* log.(N_asy .+ 1) .+ 0.685
    lines!(ax, N_asy, Lambda_leg_asy, linestyle = :dash, color = TEAL,
           linewidth = 1.5, alpha = 0.6, label = L"(2/\pi)\ln(N+1) + 0.685")

    # Lower bound for every node set, hence for the optimal nodes (Brutman)
    Lambda_low = (2/π) .* log.(N_asy .+ 1) .+ 0.5212
    lines!(ax, N_asy, Lambda_low, linestyle = :dot, color = colorant"#555555",
           linewidth = 1.5,
           label = L"(2/\pi)\ln(N+1) + 0.5212\text{, lower bound for all nodes}")

    axislegend(ax, position = :lt, framevisible = false, labelsize = 10)

    # Save figure
    mkpath(OUTPUT_DIR)
    save(joinpath(OUTPUT_DIR, "lebesgue_constants_zoom.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "lebesgue_constants_zoom.png"), fig, px_per_unit = 4)

    println("Figure saved to: ", joinpath(OUTPUT_DIR, "lebesgue_constants_zoom.pdf"))

    # Print comparison table
    println("\nLebesgue Constants Comparison (Legendre vs Chebyshev):")
    println("-" ^ 66)
    println("$(rpad("N", 4)) $(lpad("Legendre", 12)) $(lpad("Chebyshev", 12)) $(lpad("Difference", 12)) $(lpad("Lower bound", 13))")
    println("-" ^ 66)
    for i in 1:8:length(N_values)
        n = N_values[i]
        low = (2/π) * log(n + 1) + 0.5212
        diff = Lambda_cheb_values[i] - Lambda_leg_values[i]
        println("$(lpad(n, 4)) $(lpad(round(Lambda_leg_values[i], digits=4), 12)) $(lpad(round(Lambda_cheb_values[i], digits=4), 12)) $(lpad(round(diff, digits=4), 12)) $(lpad(round(low, digits=4), 13))")
    end
    println("-" ^ 66)

    return fig
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
