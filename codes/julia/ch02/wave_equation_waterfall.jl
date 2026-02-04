# wave_equation_waterfall.jl
#
# Waterfall (3D surface) plot showing the complete space-time evolution of the
# wave equation solution. This visualization demonstrates the standing wave
# oscillations of a plucked string over time.
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
const N_MODES  = 50    # Number of Fourier modes in truncation
const NX       = 150   # Spatial grid points
const NT       = 100   # Time grid points
const L_DOMAIN = Float64(pi)  # Domain length [0, L]
const C_SPEED  = 1.0   # Wave speed
const H_HEIGHT = 1.0   # Height of plucked string at center

# Period of oscillation: T = 2L/c
const T_PERIOD = 2 * L_DOMAIN / C_SPEED
const T_MAX    = T_PERIOD  # Show one full period

# Book color scheme
const NAVY = colorant"#142D6E"
const SKY  = colorant"#7896D2"

# Output path (relative to this script's location)
SCRIPT_DIR = @__DIR__
OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch02", "julia")
mkpath(OUTPUT_DIR)

# -----------------------------------------------------------------------------
# Fourier coefficients for the plucked string
# -----------------------------------------------------------------------------
"""
    plucked_string_coefficients(n_modes, L, h)

Compute Fourier sine coefficients for the plucked string initial condition.
"""
function plucked_string_coefficients(n_modes, L, h)
    a_n = zeros(n_modes + 1)
    b_n = zeros(n_modes + 1)

    for n in 1:n_modes
        if n % 2 == 1  # odd n only
            sign_val = (n % 4 == 1) ? 1 : -1
            a_n[n+1] = sign_val * 8.0 * h / (n^2 * pi^2)
        end
    end

    return a_n, b_n
end

# -----------------------------------------------------------------------------
# Evaluate the wave equation solution
# -----------------------------------------------------------------------------
"""
    wave_solution(x, t, a_n, b_n, L, c)

Evaluate the wave equation solution at positions `x` and time `t`.

    u(x,t) = Sum_{n=1}^N (a_n cos(w_n t) + b_n sin(w_n t)) sin(n pi x / L)
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

    # Create 2D grid for space-time
    x = range(0, L_DOMAIN, length=NX)
    t = range(0, T_MAX, length=NT)

    # Evaluate solution on grid
    U = zeros(NT, NX)
    for i in 1:NT
        U[i, :] = wave_solution(collect(x), t[i], a_n, b_n, L_DOMAIN, C_SPEED)
    end

    # Create custom colormap (navy to white to sky)
    navy_rgb  = (20/255, 45/255, 110/255)
    white_rgb = (1.0, 1.0, 1.0)
    sky_rgb   = (120/255, 150/255, 210/255)

    n_colors = 256
    colors_list = Vector{RGBf}(undef, n_colors)
    for i in 1:n_colors÷2
        t_val = (i - 1) / (n_colors÷2 - 1)
        colors_list[i] = RGBf(
            navy_rgb[1] * (1 - t_val) + white_rgb[1] * t_val,
            navy_rgb[2] * (1 - t_val) + white_rgb[2] * t_val,
            navy_rgb[3] * (1 - t_val) + white_rgb[3] * t_val
        )
    end
    for i in 1:n_colors÷2
        t_val = (i - 1) / (n_colors÷2 - 1)
        colors_list[n_colors÷2 + i] = RGBf(
            white_rgb[1] * (1 - t_val) + sky_rgb[1] * t_val,
            white_rgb[2] * (1 - t_val) + sky_rgb[2] * t_val,
            white_rgb[3] * (1 - t_val) + sky_rgb[3] * t_val
        )
    end
    custom_cmap = cgrad(colors_list)

    # Create figure with 3D surface
    fig = Figure(size = (700, 500))
    ax = Axis3(fig[1, 1],
               xlabel = L"x", ylabel = L"t", zlabel = L"u(x,t)",
               xlabeloffset = 40, ylabeloffset = 40, zlabeloffset = 50,
               xticklabelsize = 9, yticklabelsize = 9, zticklabelsize = 9,
               azimuth = -60 * pi / 180 + pi,
               elevation = 25 * pi / 180)

    # Plot the surface
    surface!(ax, collect(x), collect(t), U',
             colormap = custom_cmap, alpha = 0.9)

    # Add wireframe for better depth perception
    wireframe!(ax, collect(x), collect(t), U',
               color = RGBAf(0.5, 0.5, 0.5, 0.3), linewidth = 0.2,
               overdraw = false)

    # Custom x-ticks with pi notation
    ax.xticks = ([0, L_DOMAIN/4, L_DOMAIN/2, 3*L_DOMAIN/4, L_DOMAIN],
                 [L"0", L"\frac{\pi}{4}", L"\frac{\pi}{2}", L"\frac{3\pi}{4}", L"\pi"])

    # Custom t-ticks as fractions of period
    ax.yticks = ([0, T_PERIOD/4, T_PERIOD/2, 3*T_PERIOD/4, T_PERIOD],
                 [L"0", L"\frac{T}{4}", L"\frac{T}{2}", L"\frac{3T}{4}", L"T"])

    # Save figures
    save(joinpath(OUTPUT_DIR, "wave_waterfall.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "wave_waterfall.png"), fig, px_per_unit = 4)
    println("Figure saved to: $(joinpath(OUTPUT_DIR, "wave_waterfall.pdf"))")
end

main()
