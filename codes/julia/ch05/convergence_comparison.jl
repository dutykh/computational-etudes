# convergence_comparison.jl
#
# The "Computational Etude" of Chapter 5: Comparing the accuracy of finite
# difference methods (orders 2, 4, 6) against the spectral method for
# differentiating a smooth periodic function.
#
# Test Function:
#     u(x) = 1 / (2 + sin(x))   on [0, 2pi)
#     u'(x) = -cos(x) / (2 + sin(x))^2
#
# This function is:
#     - Smooth (analytic) and periodic
#     - Has a non-trivial Fourier spectrum
#     - Complex-plane singularities at x = -pi/2 +/- i*arcsinh(2)
#
# Expected Results:
#     - FD2: Error ~ O(N^{-2})  [algebraic]
#     - FD4: Error ~ O(N^{-4})  [algebraic]
#     - FD6: Error ~ O(N^{-6})  [algebraic]
#     - Spectral: Error ~ O(c^{-N})  [geometric/exponential]
#
# The spectral method reaches machine precision around N ~ 20, while
# finite difference methods require much larger N for comparable accuracy.
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

# Include shared utilities
include(joinpath(@__DIR__, "fdweights.jl"))

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

    return D, x
end

# -----------------------------------------------------------------------------
# Test function and its derivative
# -----------------------------------------------------------------------------
u_func(x)        = 1.0 / (2.0 + sin(x))
u_deriv_exact(x) = -cos(x) / (2.0 + sin(x))^2

# -----------------------------------------------------------------------------
# Finite difference matrix construction
# -----------------------------------------------------------------------------
"""
    fd_diff_periodic(N, order)

Construct a periodic finite difference differentiation matrix.
"""
function fd_diff_periodic(N, order)
    h = 2pi / N

    stencil_half = order ÷ 2
    stencil_nodes = h .* (-stencil_half:stencil_half)

    weights = fdweights(0.0, collect(stencil_nodes), 1)

    D = zeros(N, N)
    for i in 1:N
        for (k, w) in enumerate(weights)
            j = mod(i - 1 + k - 1 - stencil_half, N) + 1
            D[i, j] += w
        end
    end

    return D
end

# -----------------------------------------------------------------------------
# Convergence study
# -----------------------------------------------------------------------------
"""
    compute_errors(N_values, orders)

Compute maximum differentiation errors for various N and FD orders.
"""
function compute_errors(N_values, orders)
    errors_fd = Dict(order => Float64[] for order in orders)
    errors_spectral = Float64[]

    for N in N_values
        h = 2pi / N
        x = h .* (0:N-1)

        u = u_func.(x)
        du_exact = u_deriv_exact.(x)

        # Finite difference methods
        for order in orders
            if order < N
                D_fd = fd_diff_periodic(N, order)
                du_fd = D_fd * u
                err = maximum(abs.(du_fd .- du_exact))
            else
                err = NaN
            end
            push!(errors_fd[order], err)
        end

        # Spectral method
        D_spec, _ = spectral_diff_periodic(N)
        du_spec = D_spec * u
        err_spec = maximum(abs.(du_spec .- du_exact))
        push!(errors_spectral, err_spec)
    end

    return errors_fd, errors_spectral
end

# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
function main()
    # Grid sizes to test
    N_values = [4, 6, 8, 10, 12, 14, 16, 20, 24, 28, 32, 40, 48, 56, 64]
    orders = [2, 4, 6]

    # Compute errors
    errors_fd, errors_spectral = compute_errors(N_values, orders)

    # Create figure
    fig = Figure(size = (640, 440))
    ax = Axis(fig[1, 1],
              xlabel = L"Number of grid points $N$",
              ylabel = L"Maximum error $\|u' - Du\|_\infty$",
              title = "Differentiation Error: Finite Differences vs Spectral Method",
              yscale = log10)

    # Plot finite difference errors
    colors_fd = Dict(2 => CORAL, 4 => TEAL, 6 => SKY)
    markers_fd = Dict(2 => :circle, 4 => :rect, 6 => :utriangle)

    for order in orders
        scatterlines!(ax, N_values, errors_fd[order],
                      marker = markers_fd[order], color = colors_fd[order],
                      linewidth = 1.5, markersize = 6,
                      label = "FD$order (order $order)")
    end

    # Plot spectral errors
    scatterlines!(ax, N_values, errors_spectral,
                  marker = :diamond, color = NAVY, linewidth = 2, markersize = 7,
                  label = "Spectral")

    # Add theoretical reference lines
    N_ref = range(8, 64, length = 100)

    # FD2: O(N^{-2})
    C2 = errors_fd[2][5] * N_values[5]^2
    lines!(ax, collect(N_ref), C2 ./ N_ref .^ 2,
           linestyle = :dash, color = (CORAL, 0.4), linewidth = 1)

    # FD4: O(N^{-4})
    C4 = errors_fd[4][5] * N_values[5]^4
    lines!(ax, collect(N_ref), C4 ./ N_ref .^ 4,
           linestyle = :dash, color = (TEAL, 0.4), linewidth = 1)

    # FD6: O(N^{-6})
    C6 = errors_fd[6][5] * N_values[5]^6
    lines!(ax, collect(N_ref), C6 ./ N_ref .^ 6,
           linestyle = :dash, color = (SKY, 0.4), linewidth = 1)

    # Machine epsilon line
    hlines!(ax, [2.2e-16], color = :gray, linestyle = :dot, linewidth = 1)
    text!(ax, 66, 4e-16,
          text = "Machine precision", fontsize = 9, color = :gray,
          align = (:right, :bottom))

    # Add slope annotations
    text!(ax, 68, C2 / 68^2 * 1.8,
          text = L"O(N^{-2})", fontsize = 10, color = CORAL,
          align = (:right, :bottom))
    text!(ax, 68, C4 / 68^4 * 2.5,
          text = L"O(N^{-4})", fontsize = 10, color = TEAL,
          align = (:right, :bottom))
    text!(ax, 68, C6 / 68^6 * 3.0,
          text = L"O(N^{-6})", fontsize = 10, color = SKY,
          align = (:right, :bottom))

    # Add "spectral accuracy" annotation
    text!(ax, 12, 5e-10,
          text = "Spectral\naccuracy!", fontsize = 10, color = NAVY,
          font = :bold, align = (:center, :center))

    # Add test function info
    text!(ax, 0.98, 0.98,
          text = L"u(x) = \frac{1}{2 + \sin(x)}, \quad [0, 2\pi)",
          fontsize = 9,
          space = :relative,
          align = (:right, :top))

    # Axis settings
    xlims!(ax, 0, 70)
    ylims!(ax, 1e-16, 1e1)
    axislegend(ax, position = (0.55, 0.85), labelsize = 9)

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

    # Save figure
    save(joinpath(OUTPUT_DIR, "convergence_comparison.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "convergence_comparison.png"), fig, px_per_unit = 4)
    println("Figure saved to: $(joinpath(OUTPUT_DIR, "convergence_comparison.pdf"))")

    # Print results table
    println("\nConvergence Study Results:")
    println("=" ^ 70)
    @printf("%4s  %12s  %12s  %12s  %12s\n", "N", "FD2", "FD4", "FD6", "Spectral")
    println("-" ^ 70)
    for (i, N) in enumerate(N_values)
        @printf("%4d  %12.4e  %12.4e  %12.4e  %12.4e\n",
                N, errors_fd[2][i], errors_fd[4][i], errors_fd[6][i], errors_spectral[i])
    end
    println("=" ^ 70)

    # Compute convergence rates
    println("\nConvergence Rates (estimated from last two points):")
    for order in orders
        e1, e2 = errors_fd[order][end-1], errors_fd[order][end]
        n1, n2 = N_values[end-1], N_values[end]
        if e1 > 0 && e2 > 0
            rate = log(e1 / e2) / log(n2 / n1)
            @printf("  FD%d: %.2f (expected: %d)\n", order, rate, order)
        end
    end

    # Spectral rate (geometric)
    e1, e2 = errors_spectral[7], errors_spectral[9]  # N=16, N=24
    if e1 > 1e-15 && e2 > 1e-15
        rate = (log(e1) - log(e2)) / (N_values[9] - N_values[7])
        @printf("  Spectral: geometric rate ~ exp(%.3f*N)\n", -rate)
    end
end

main()
