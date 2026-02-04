# heat_equation_waterfall.jl
#
# Waterfall (3D surface) plot showing the complete space-time evolution of the
# heat equation solution. This visualization demonstrates how the initial
# triangle wave smooths out over time as higher frequency modes decay faster.
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
const N_MODES = 50    # Number of Fourier modes in truncation
const NX      = 150   # Spatial grid points
const NT      = 100   # Time grid points
const T_MAX   = 0.3   # Maximum time

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
"""
function triangle_fourier_coefficients(n_modes)
    a0 = pi / 2
    a_n = zeros(n_modes + 1)
    b_n = zeros(n_modes + 1)

    for n in 1:n_modes
        if n % 2 == 1  # odd n only
            a_n[n+1] = 4.0 / (pi * n^2)
        end
    end

    return a0, a_n, b_n
end

# -----------------------------------------------------------------------------
# Evaluate the heat equation solution
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

    # Create 2D grid for space-time
    x = range(0, 2pi, length=NX)
    t = range(0, T_MAX, length=NT)

    # Evaluate solution on grid
    U = zeros(NT, NX)
    for i in 1:NT
        U[i, :] = heat_solution(collect(x), t[i], a0, a_n, b_n)
    end

    # Create custom colormap (navy to white to sky)
    navy_rgb = (20/255, 45/255, 110/255)
    white_rgb = (1.0, 1.0, 1.0)
    sky_rgb = (120/255, 150/255, 210/255)

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
    ax.xticks = ([0, pi/2, pi, 3pi/2, 2pi],
                 [L"0", L"\frac{\pi}{2}", L"\pi", L"\frac{3\pi}{2}", L"2\pi"])

    # Save figures
    save(joinpath(OUTPUT_DIR, "heat_waterfall.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "heat_waterfall.png"), fig, px_per_unit = 4)
    println("Figure saved to: $(joinpath(OUTPUT_DIR, "heat_waterfall.pdf"))")
end

main()
