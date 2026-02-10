# periodic_cardinal_functions.jl
#
# Visualizes periodic cardinal functions (discrete Dirichlet kernel) for equispaced
# nodes on a periodic domain.
#
# The periodic cardinal function is:
#     phi_j(x) = sin(N(x - x_j)/2) / (N sin((x - x_j)/2))
#
# These functions satisfy phi_j(x_k) = delta_{jk} (Kronecker delta), making them
# the "cardinal functions" for trigonometric interpolation on periodic domains.
# The interpolant is:
#     p(x) = sum_{j=0}^{N-1} u_j * phi_j(x)
#
# Unlike Lagrange basis polynomials on non-periodic domains, the periodic cardinal
# functions are bounded and exhibit damped oscillations away from their peak,
# reflecting the sin(x)/x-like structure of the discrete Dirichlet kernel.
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
const N      = 16     # Number of grid points (should be even)
const N_FINE = 1000   # Number of points for smooth plotting

# Book color scheme
const NAVY   = colorant"#142D6E"
const SKY    = colorant"#7896D2"
const CORAL  = colorant"#E74C3C"
const TEAL   = colorant"#16A085"
const PURPLE = colorant"#8E44AD"
const ORANGE = colorant"#E67E22"

SCRIPT_DIR = @__DIR__
OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch05", "julia")
mkpath(OUTPUT_DIR)

# -----------------------------------------------------------------------------
# Periodic cardinal function computation
# -----------------------------------------------------------------------------
"""
    periodic_cardinal(x, x_j, N)

Compute the periodic cardinal function phi_j(x) centered at x_j.

    phi_j(x) = sin(N(x - x_j)/2) / (N sin((x - x_j)/2))

This is the discrete Dirichlet kernel, satisfying phi_j(x_k) = delta_{jk}
for equispaced nodes x_k = 2*pi*k/N on [0, 2*pi).

# Arguments
- `x`: points at which to evaluate the cardinal function
- `x_j`: center point (node location)
- `N`: number of grid points

# Returns
- `phi`: values of phi_j at x
"""
function periodic_cardinal(x::AbstractVector, x_j::Real, N::Int)
    theta = (x .- x_j) ./ 2.0

    phi = zeros(length(x))

    for i in eachindex(x)
        if abs(sin(theta[i])) < 1e-14
            # Singularity: use L'Hopital's rule
            # lim_{theta->0} sin(N*theta) / (N*sin(theta)) = 1
            phi[i] = 1.0
        else
            phi[i] = sin(N * theta[i]) / (N * sin(theta[i]))
        end
    end

    return phi
end

# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
function main()
    # Grid points on [0, 2*pi)
    h = 2pi / N
    x_nodes = h .* (0:N-1) |> collect

    # Fine grid for plotting
    x_fine = collect(range(0, 2pi, length = N_FINE))

    # Create figure
    fig = Figure(size = (900, 450))

    ax = Axis(fig[1, 1],
              xlabel = L"x",
              ylabel = L"\phi_j(x)",
              title = L"Periodic Cardinal Functions ($N = %$N$)",
              xticks = ([0, pi/2, pi, 3pi/2, 2pi],
                        [L"0", L"\pi/2", L"\pi", L"3\pi/2", L"2\pi"]))

    # Indices of cardinal functions to plot (spread across the domain)
    # Note: Python uses 0-based indices [0,4,8,12]; Julia uses 1-based [1,5,9,13]
    indices_to_plot = [1, 5, 9, 13]
    # Labels use 0-based indexing to match mathematical convention
    labels_to_plot  = [0, 4, 8, 12]
    colors = [CORAL, TEAL, SKY, PURPLE]

    # Plot cardinal functions
    for (idx, (j, lbl, col)) in enumerate(zip(indices_to_plot, labels_to_plot, colors))
        phi_j = periodic_cardinal(x_fine, x_nodes[j], N)
        lines!(ax, x_fine, phi_j, color = col, linewidth = 1.5,
               label = L"\phi_{%$lbl}(x)")
    end

    # Mark all nodes on x-axis
    scatter!(ax, x_nodes, zeros(N), color = NAVY, markersize = 5,
             label = "Nodes")

    # Highlight the peak values at each plotted cardinal function's node
    for (j, col) in zip(indices_to_plot, colors)
        scatter!(ax, [x_nodes[j]], [1.0], color = col, markersize = 8,
                 strokecolor = :white, strokewidth = 1.5)
    end

    # Reference lines
    hlines!(ax, [0.0], color = RGBf(0.533, 0.533, 0.533), linewidth = 0.5)
    hlines!(ax, [1.0], color = RGBf(0.533, 0.533, 0.533), linewidth = 0.5,
            linestyle = :dash)

    # Axis limits
    xlims!(ax, 0, 2pi)
    ylims!(ax, -0.25, 1.15)

    # Legend
    axislegend(ax, position = :rt, labelsize = 9, nbanks = 2)

    # Add annotation about cardinal property
    # Arrow from annotation text to the peak of phi_4 (index 5 in Julia)
    arrows!(ax, [x_nodes[5] + 0.8], [0.85], [-0.75], [0.12],
            color = :gray, linewidth = 0.8)
    text!(ax, x_nodes[5] + 0.85, 0.85,
          text = L"\phi_j(x_k) = \delta_{jk}",
          fontsize = 10, color = :gray,
          align = (:left, :center))

    # Clean styling
    hidespines!(ax, :t, :r)
    ax.leftspinecolor   = RGBf(0.267, 0.267, 0.267)
    ax.bottomspinecolor = RGBf(0.267, 0.267, 0.267)
    ax.xtickcolor = RGBf(0.267, 0.267, 0.267)
    ax.ytickcolor = RGBf(0.267, 0.267, 0.267)
    ax.xgridvisible = true
    ax.xgridcolor   = RGBAf(0.533, 0.533, 0.533, 0.2)
    ax.xgridstyle   = :solid
    ax.ygridvisible = true
    ax.ygridcolor   = RGBAf(0.533, 0.533, 0.533, 0.2)
    ax.ygridstyle   = :solid

    # Save figure
    save(joinpath(OUTPUT_DIR, "periodic_cardinal_functions.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "periodic_cardinal_functions.png"), fig, px_per_unit = 4)
    println("Figure saved to: $(joinpath(OUTPUT_DIR, "periodic_cardinal_functions.pdf"))")

    # -------------------------------------------------------------------------
    # Verify cardinal property
    # -------------------------------------------------------------------------
    println("\nVerification of cardinal property for N = $N:")
    println("-" ^ 50)
    @printf("%4s %12s %25s\n", "j", "phi_j(x_j)", "max |phi_j(x_k)|, k!=j")
    println("-" ^ 50)
    for j in 1:N
        phi_at_nodes = periodic_cardinal(x_nodes, x_nodes[j], N)
        phi_at_own = phi_at_nodes[j]
        phi_at_others = [phi_at_nodes[k] for k in 1:N if k != j]
        max_at_others = maximum(abs.(phi_at_others))
        @printf("%4d %12.6f %25.2e\n", j - 1, phi_at_own, max_at_others)
    end
    println("-" ^ 50)

    # Print maximum oscillation amplitude
    println("\nMaximum |phi_j(x)| away from x_j:")
    for (j, lbl) in zip(indices_to_plot, labels_to_plot)
        phi_j = periodic_cardinal(x_fine, x_nodes[j], N)
        # Find points away from the node
        away_from_node = abs.(x_fine .- x_nodes[j]) .> 0.1
        if any(away_from_node)
            max_away = maximum(abs.(phi_j[away_from_node]))
            @printf("  phi_%d: max |phi_%d(x)| = %.4f\n", lbl, lbl, max_away)
        end
    end
end

main()
