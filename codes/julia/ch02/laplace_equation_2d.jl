# laplace_equation_2d.jl
#
# Visualizes the solution of the Laplace equation in a periodic strip using the
# truncated Fourier series derived in Chapter 2. The boundary condition is a
# smooth function demonstrating how different Fourier modes decay toward the
# interior.
#
# Laplace equation in strip D = [0, 2pi] x [0, 1]:
#     u_xx + u_yy = 0
#
# Boundary conditions:
#     u(x, 0) = f(x) = sin(x) + 0.5*sin(3x)   (prescribed at bottom)
#     u(x, 1) = 0                               (zero at top)
#     u periodic in x with period 2pi
#
# Analytical solution:
#     u(x,y) = a_0(1-y) + Sum [a_n cos(nx) + b_n sin(nx)] sinh(n(1-y))/sinh(n)
#
# where a_n, b_n are Fourier coefficients of f(x).
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
const NX      = 300   # Grid points in x direction
const NY      = 150   # Grid points in y direction

# Book color scheme
const NAVY = colorant"#142D6E"
const SKY  = colorant"#7896D2"

# Output path (relative to this script's location)
SCRIPT_DIR = @__DIR__
OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch02", "julia")
mkpath(OUTPUT_DIR)

# -----------------------------------------------------------------------------
# Boundary data and Fourier coefficients
# -----------------------------------------------------------------------------
"""
    boundary_function(x)

Boundary data at y = 0: f(x) = sin(x) + 0.5*sin(3x).
"""
function boundary_function(x)
    return sin.(x) .+ 0.5 .* sin.(3 .* x)
end

"""
    boundary_fourier_coefficients(n_modes)

Compute Fourier coefficients for f(x) = sin(x) + 0.5*sin(3x).

For this specific function:
    a_0 = 0
    a_n = 0 for all n (no cosine terms)
    b_1 = 1, b_3 = 0.5, b_n = 0 otherwise

Returns (a0, a_n, b_n).
"""
function boundary_fourier_coefficients(n_modes)
    a0 = 0.0
    a_n = zeros(n_modes + 1)
    b_n = zeros(n_modes + 1)

    # f(x) = sin(x) + 0.5*sin(3x)
    b_n[2] = 1.0    # b_1 (1-indexed: index 2 corresponds to n=1)
    b_n[4] = 0.5    # b_3 (1-indexed: index 4 corresponds to n=3)

    return a0, a_n, b_n
end

# -----------------------------------------------------------------------------
# Evaluate the Laplace equation solution
# -----------------------------------------------------------------------------
"""
    laplace_solution(X, Y, a0, a_n, b_n)

Evaluate the Laplace equation solution in the strip.

    u(x,y) = a_0(1-y) + Sum_{n=1}^N [a_n cos(nx) + b_n sin(nx)] sinh(n(1-y))/sinh(n)

Parameters:
    X, Y: coordinate arrays (2D)
    a0, a_n, b_n: Fourier coefficients

Returns:
    U: solution values (same shape as X, Y)
"""
function laplace_solution(X, Y, a0, a_n, b_n)
    # Start with constant mode
    U = a0 .* (1.0 .- Y)

    n_modes = length(a_n) - 1
    for n in 1:n_modes
        # Skip if both coefficients are zero
        if abs(a_n[n+1]) < 1e-15 && abs(b_n[n+1]) < 1e-15
            continue
        end

        # y-dependent factor: sinh(n(1-y)) / sinh(n)
        sinh_n = sinh(n)
        if sinh_n > 1e-10
            y_factor = sinh.(n .* (1.0 .- Y)) ./ sinh_n
        else
            y_factor = zeros(size(Y))
        end

        # Add contribution from this mode
        @. U += (a_n[n+1] * cos(n * X) + b_n[n+1] * sin(n * X)) * y_factor
    end

    return U
end

# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
function main()
    # Compute Fourier coefficients
    a0, a_n, b_n = boundary_fourier_coefficients(N_MODES)

    # Create 2D grid
    x = range(0, 2pi, length=NX)
    y = range(0, 1,   length=NY)
    X = [xi for yi in y, xi in x]
    Y = [yi for yi in y, xi in x]

    # Evaluate solution on grid
    U = laplace_solution(X, Y, a0, a_n, b_n)

    # Create custom diverging colormap (navy to white to sky)
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

    # Determine symmetric color limits
    vmax = maximum(abs.(U))

    # Create figure
    fig = Figure(size = (600, 350))
    ax = Axis(fig[1, 1],
              xlabel = L"x",
              ylabel = L"y",
              xlabelpadding = 3,
              ylabelpadding = 3)

    # Filled contour plot (Makie expects z[x_idx, y_idx], so transpose U)
    levels_fill = range(-vmax, vmax, length=25)
    cf = contourf!(ax, collect(x), collect(y), U',
                   levels = levels_fill, colormap = custom_cmap,
                   extendlow = :auto, extendhigh = :auto)

    # Add contour lines for emphasis
    levels_line = range(-vmax, vmax, length=13)
    contour!(ax, collect(x), collect(y), U',
             levels = levels_line, color = RGBAf(0.2, 0.2, 0.2, 0.5),
             linewidth = 0.3)

    # Colorbar
    Colorbar(fig[1, 2], cf, label = L"u(x,y)", labelsize = 11,
             ticklabelsize = 9, width = 12)

    # Custom x-ticks with pi notation
    ax.xticks = ([0, pi/2, pi, 3pi/2, 2pi],
                 [L"0", L"\frac{\pi}{2}", L"\pi", L"\frac{3\pi}{2}", L"2\pi"])

    # y-ticks
    ax.yticks = [0, 0.25, 0.5, 0.75, 1.0]

    # Spine styling (keep all four spines for the strip domain)
    ax.topspinecolor    = RGBf(0.267, 0.267, 0.267)
    ax.rightspinecolor  = RGBf(0.267, 0.267, 0.267)
    ax.leftspinecolor   = RGBf(0.267, 0.267, 0.267)
    ax.bottomspinecolor = RGBf(0.267, 0.267, 0.267)
    ax.topspinevisible  = true
    ax.rightspinevisible = true
    ax.spinewidth = 0.5
    ax.xtickcolor = RGBf(0.267, 0.267, 0.267)
    ax.ytickcolor = RGBf(0.267, 0.267, 0.267)

    # Save figures
    save(joinpath(OUTPUT_DIR, "laplace_solution.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "laplace_solution.png"), fig, px_per_unit = 4)
    println("Figure saved to: $(joinpath(OUTPUT_DIR, "laplace_solution.pdf"))")
end

main()
