# wave_equation_evolution.jl
#
# Visualizes the time evolution of the wave equation solution using the
# truncated Fourier sine series derived in Chapter 2. The initial condition is
# a plucked string (triangular displacement), which demonstrates standing wave
# oscillations.
#
# Wave equation (Dirichlet BCs on [0, L]):
#     u_tt = c^2 u_xx
#
# Analytical solution:
#     u(x,t) = Sum (a_n cos(w_n t) + b_n sin(w_n t)) sin(n pi x / L)
#
# where w_n = c n pi / L
#
# Plucked string initial condition (triangle):
#     f(x) = (2h/L)x           for 0 <= x <= L/2
#     f(x) = 2h(1 - x/L)      for L/2 <= x <= L
#     g(x) = 0                 (initially at rest)
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
const L_DOMAIN = Float64(pi)  # Domain length [0, L]
const C_SPEED  = 1.0     # Wave speed
const H_HEIGHT = 1.0     # Height of plucked string at center

# Period of oscillation: T = 2L/c
const T_PERIOD = 2 * L_DOMAIN / C_SPEED

# Time snapshots: show half a period of oscillation
const TIMES = [0.0, T_PERIOD/8, T_PERIOD/4, 3*T_PERIOD/8, T_PERIOD/2]

# Book color scheme
const NAVY = colorant"#142D6E"
const SKY  = colorant"#7896D2"

# Output path (relative to this script's location)
SCRIPT_DIR = @__DIR__
OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch02", "julia")
mkpath(OUTPUT_DIR)

# -----------------------------------------------------------------------------
# Fourier coefficients for the plucked string (triangle)
# -----------------------------------------------------------------------------
"""
    plucked_string_coefficients(n_modes, L, h)

Compute Fourier sine coefficients for the plucked string initial condition.

    f(x) = (2h/L)x           for 0 <= x <= L/2
    f(x) = 2h(1 - x/L)      for L/2 <= x <= L

The Fourier sine series is:
    f(x) = Sum a_n sin(n pi x / L)

where a_n = (8h/(n^2 pi^2)) sin(n pi/2)

Returns (a_n, b_n) where a_n[n+1] is the coefficient of sin(n pi x / L)
and b_n[n+1] = 0 since the string is initially at rest.
"""
function plucked_string_coefficients(n_modes, L, h)
    a_n = zeros(n_modes + 1)  # a_n[n+1] = coefficient of sin(n*pi*x/L)
    b_n = zeros(n_modes + 1)  # b_n[n+1] = 0 since initially at rest

    for n in 1:n_modes
        if n % 2 == 1  # odd n only
            # sin(n*pi/2) = 1 for n = 1,5,9,... and -1 for n = 3,7,11,...
            sign_val = (n % 4 == 1) ? 1 : -1
            a_n[n+1] = sign_val * 8.0 * h / (n^2 * pi^2)
        end
    end

    return a_n, b_n
end

# -----------------------------------------------------------------------------
# Evaluate the truncated Fourier series solution
# -----------------------------------------------------------------------------
"""
    wave_solution(x, t, a_n, b_n, L, c)

Evaluate the wave equation solution at positions `x` and time `t`.

    u(x,t) = Sum_{n=1}^N (a_n cos(w_n t) + b_n sin(w_n t)) sin(n pi x / L)

where w_n = c n pi / L
"""
function wave_solution(x, t, a_n, b_n, L, c)
    u = zeros(length(x))

    n_modes = length(a_n) - 1
    for n in 1:n_modes
        omega_n = c * n * pi / L
        spatial  = @. sin(n * pi * x / L)
        temporal = a_n[n+1] * cos(omega_n * t) + b_n[n+1] * sin(omega_n * t)
        @. u += temporal * spatial
    end

    return u
end

# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
function main()
    # Compute Fourier coefficients
    a_n, b_n = plucked_string_coefficients(N_MODES, L_DOMAIN, H_HEIGHT)

    # Spatial grid
    x = range(0, L_DOMAIN, length=N_POINTS)

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
        u = wave_solution(collect(x), t, a_n, b_n, L_DOMAIN, C_SPEED)

        # Create time label as fraction of period
        if t == 0.0
            label = L"t = 0"
        else
            frac = t / T_PERIOD
            if abs(frac - 1/8) < 0.001
                label = L"t = T/8"
            elseif abs(frac - 1/4) < 0.001
                label = L"t = T/4"
            elseif abs(frac - 3/8) < 0.001
                label = L"t = 3T/8"
            elseif abs(frac - 1/2) < 0.001
                label = L"t = T/2"
            else
                label = latexstring("t = $(@sprintf("%.3f", t))")
            end
        end

        lines!(ax, collect(x), u, color = colors[i], linewidth = line_widths[i],
               label = label)
    end

    # Axis limits
    xlims!(ax, 0, L_DOMAIN)
    ylims!(ax, -1.2, 1.2)

    # Custom x-ticks with pi notation
    ax.xticks = ([0, L_DOMAIN/4, L_DOMAIN/2, 3*L_DOMAIN/4, L_DOMAIN],
                 [L"0", L"\frac{\pi}{4}", L"\frac{\pi}{2}", L"\frac{3\pi}{4}", L"\pi"])

    # Hide top and right spines
    hidespines!(ax, :t, :r)
    ax.leftspinecolor   = RGBf(0.267, 0.267, 0.267)
    ax.bottomspinecolor = RGBf(0.267, 0.267, 0.267)
    ax.xtickcolor = RGBf(0.267, 0.267, 0.267)
    ax.ytickcolor = RGBf(0.267, 0.267, 0.267)

    # Subtle horizontal grid
    ax.ygridvisible = true
    ax.ygridcolor   = RGBAf(0.533, 0.533, 0.533, 0.2)
    ax.ygridstyle   = :solid
    ax.xgridvisible = false

    # Reference line at y = 0
    hlines!(ax, [0.0], color = RGBAf(0.533, 0.533, 0.533, 0.3), linewidth = 0.5)

    # Legend
    axislegend(ax, position = :ct, nbanks = 3, framevisible = true,
               framecolor = :transparent, bgcolor = RGBAf(1, 1, 1, 0.9),
               padding = (8, 8, 4, 4), rowgap = 0)

    # Save figures
    save(joinpath(OUTPUT_DIR, "wave_evolution.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "wave_evolution.png"), fig, px_per_unit = 4)
    println("Figure saved to: $(joinpath(OUTPUT_DIR, "wave_evolution.pdf"))")
end

main()
