# heat_equation_evolution.jl
#
# Visualizes the time evolution of the heat equation solution using the
# truncated Fourier series derived in Chapter 2. The initial condition is
# a triangle wave, which demonstrates how higher frequency modes decay
# faster than lower ones.
#
# Heat equation (periodic BCs on [0, 2pi]):
#     u_t = u_xx
#
# Analytical solution:
#     u(x,t) = a_0 + Sum (a_n cos(nx) + b_n sin(nx)) exp(-n^2 t)
#
# Triangle wave initial condition:
#     f(x) = pi - |x - pi|  on [0, 2pi]
#
# Fourier coefficients (computed analytically):
#     a_0 = pi/2
#     a_n = (4/pi) / n^2  for odd n, 0 for even n
#     b_n = 0             (symmetric function)
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
const N_MODES  = 50      # Number of Fourier modes in truncation
const N_POINTS = 1000    # Spatial grid points for smooth curves
const TIMES    = [0.0, 0.005, 0.02, 0.05, 0.1, 0.3]  # Time snapshots

# Book color scheme
const NAVY = colorant"#142D6E"
const SKY  = colorant"#7896D2"

# Output path (relative to this script's location)
SCRIPT_DIR = @__DIR__
OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch02", "julia")
mkpath(OUTPUT_DIR)

# -----------------------------------------------------------------------------
# Fourier coefficients for the triangle wave
# -----------------------------------------------------------------------------
"""
    triangle_fourier_coefficients(n_modes)

Compute Fourier coefficients for the triangle wave f(x) = pi - |x - pi|.

The Fourier series is:
    f(x) = pi/2 + (4/pi) Sum_{k=1}^inf cos((2k-1)x) / (2k-1)^2

Returns (a0, a_n, b_n) where a0 is the constant term,
a_n[n] is the coefficient of cos(nx), and b_n[n] is the coefficient of sin(nx).
"""
function triangle_fourier_coefficients(n_modes)
    a0 = pi / 2

    a_n = zeros(n_modes + 1)  # a_n[n+1] = coefficient of cos(nx)
    b_n = zeros(n_modes + 1)  # b_n[n+1] = coefficient of sin(nx)

    for n in 1:n_modes
        if n % 2 == 1  # odd n only
            a_n[n+1] = 4.0 / (pi * n^2)
        end
    end

    return a0, a_n, b_n
end

# -----------------------------------------------------------------------------
# Evaluate the truncated Fourier series solution
# -----------------------------------------------------------------------------
"""
    heat_solution(x, t, a0, a_n, b_n)

Evaluate the heat equation solution at positions `x` and time `t`.

    u(x,t) = a_0 + Sum_{n=1}^N (a_n cos(nx) + b_n sin(nx)) exp(-n^2 t)
"""
function heat_solution(x, t, a0, a_n, b_n)
    u = fill(a0, length(x))

    n_modes = length(a_n) - 1
    for n in 1:n_modes
        decay = exp(-n^2 * t)
        @. u += (a_n[n+1] * cos(n * x) + b_n[n+1] * sin(n * x)) * decay
    end

    return u
end

# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
function main()
    # Compute Fourier coefficients
    a0, a_n, b_n = triangle_fourier_coefficients(N_MODES)

    # Spatial grid
    x = range(0, 2pi, length=N_POINTS)

    # Create figure with golden ratio proportions
    fig = Figure(size = (600, 370))
    ax = Axis(fig[1, 1],
              xlabel = L"x",
              ylabel = L"u(x,t)",
              xlabelpadding = 3,
              ylabelpadding = 3)

    # Elegant color palette: navy to sky blue gradient
    n_times = length(TIMES)
    colors = [RGBf(
        (20 + (i-1)/(n_times-1) * 100) / 255,
        (45 + (i-1)/(n_times-1) * 105) / 255,
        (110 + (i-1)/(n_times-1) * 100) / 255
    ) for i in 1:n_times]

    # Varying line widths
    line_widths = range(2.0, 1.2, length=n_times)

    # Plot solution at each time
    for (i, t) in enumerate(TIMES)
        u = heat_solution(collect(x), t, a0, a_n, b_n)

        # Create time label
        if t == 0.0
            label = L"t = 0"
        elseif t < 0.01
            label = latexstring("t = $(Printf.@sprintf("%.3f", t))")
        elseif t < 0.1
            label = latexstring("t = $(Printf.@sprintf("%.2f", t))")
        else
            label = latexstring("t = $(Printf.@sprintf("%.1f", t))")
        end

        lines!(ax, collect(x), u, color = colors[i], linewidth = line_widths[i],
               label = label)
    end

    # Axis limits
    xlims!(ax, 0, 2pi)
    ylims!(ax, -0.1, 3.4)

    # Custom x-ticks with pi notation
    ax.xticks = ([0, pi/2, pi, 3pi/2, 2pi],
                 [L"0", L"\frac{\pi}{2}", L"\pi", L"\frac{3\pi}{2}", L"2\pi"])

    # Hide top and right spines
    hidespines!(ax, :t, :r)
    ax.leftspinecolor  = RGBf(0.267, 0.267, 0.267)
    ax.bottomspinecolor = RGBf(0.267, 0.267, 0.267)
    ax.xtickcolor = RGBf(0.267, 0.267, 0.267)
    ax.ytickcolor = RGBf(0.267, 0.267, 0.267)

    # Subtle horizontal grid
    ax.ygridvisible = true
    ax.ygridcolor   = RGBAf(0.533, 0.533, 0.533, 0.2)
    ax.ygridstyle   = :solid
    ax.xgridvisible = false

    # Legend
    axislegend(ax, position = :ct, nbanks = 3, framevisible = true,
               framecolor = :transparent, bgcolor = RGBAf(1, 1, 1, 0.9),
               padding = (8, 8, 4, 4), rowgap = 0)

    # Save figures
    save(joinpath(OUTPUT_DIR, "heat_evolution.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "heat_evolution.png"), fig, px_per_unit = 4)
    println("Figure saved to: $(joinpath(OUTPUT_DIR, "heat_evolution.pdf"))")
end

main()
