# two_views_function.jl
#
# Chapter 9: Physical and Fourier Space on Grids
#
# Demonstrates the two views of a function: physical space and Fourier space.
# Shows a smooth function and its rapidly decaying Fourier transform.
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
using FFTW

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
const NAVY   = colorant"#142D6E"
const SKY    = colorant"#7896D2"
const CORAL  = colorant"#E74C3C"
const TEAL   = colorant"#16A085"
const PURPLE = colorant"#9B59B6"
const ORANGE = colorant"#E67E22"

# Output path
const SCRIPT_DIR = @__DIR__
const OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch09", "julia")

# -----------------------------------------------------------------------------
# Compute spectrum via FFT
# -----------------------------------------------------------------------------
"""
    compute_spectrum(u_func; L=20.0, N=1024)

Approximate the Fourier transform of `u_func` via FFT.

# Arguments
- `u_func`: callable function to transform
- `L`: domain is [-L/2, L/2]
- `N`: number of grid points

# Returns
- `x, u, k, u_hat`: grid, values, wavenumbers, Fourier coefficients
"""
function compute_spectrum(u_func; L=20.0, N=1024)
    x  = range(-L/2, stop=L/2, length=N+1)[1:N]  # endpoint=false equivalent
    dx = x[2] - x[1]
    u  = u_func.(x)

    # Approximate continuous FT via FFT
    k     = 2π .* fftfreq(N, 1/dx)
    u_hat = dx .* fft(u)

    return collect(x), u, k, u_hat
end

# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
function main()
    # Test function: Gaussian-modulated cosine
    u(x) = exp(-x^2) * cos(3x)

    # Compute physical and Fourier representations
    x, u_vals, k, u_hat = compute_spectrum(u; L=20.0, N=2^10)

    # Shift for plotting (center k=0)
    k_shifted     = fftshift(k)
    u_hat_shifted = fftshift(u_hat)

    # Create figure with two panels
    fig = Figure(size = (700, 280))

    # Left: Physical space
    ax1 = Axis(fig[1, 1],
               xlabel = L"x",
               ylabel = L"u(x)",
               title  = "Physical space",
               limits = ((-8, 8), nothing))
    lines!(ax1, x, u_vals, color = NAVY, linewidth = 1.2)
    hlines!(ax1, [0.0], color = :black, linewidth = 0.5)

    # Right: Fourier space (log scale)
    ax2 = Axis(fig[1, 2],
               xlabel = L"k",
               ylabel = L"|\hat{u}(k)|",
               title  = "Fourier space (approximate spectrum)",
               yscale = log10,
               limits = ((-15, 15), (1e-16, 10)))
    lines!(ax2, k_shifted, abs.(u_hat_shifted), color = NAVY, linewidth = 1.2)

    # Save figure
    mkpath(OUTPUT_DIR)
    save(joinpath(OUTPUT_DIR, "two_views_function.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "two_views_function.png"), fig, px_per_unit = 4)

    println("Figure saved to: ", joinpath(OUTPUT_DIR, "two_views_function.pdf"))

    return fig
end

# Run if executed as script
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
