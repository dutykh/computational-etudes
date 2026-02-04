# zero_padding_interpolation.jl
#
# Chapter 9: Physical and Fourier Space on Grids
#
# Demonstrates band-limited interpolation via zero-padding in Fourier space.
# This is Etude 4 from the chapter.
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
# Zero-padding interpolation
# -----------------------------------------------------------------------------
"""
    zero_pad_interpolate(v; q=4)

Interpolate a periodic function by zero-padding in Fourier space.

# Arguments
- `v`: function values on coarse grid (length N)
- `q`: upsampling factor (fine grid has q*N points)

# Returns
- `v_fine`: interpolated values on fine grid
"""
function zero_pad_interpolate(v; q=4)
    N = length(v)
    M = q * N

    # Forward FFT
    v_hat = fft(v)

    # Zero-pad in Fourier space
    v_hat_padded = zeros(ComplexF64, M)

    # Copy low frequencies
    # Positive frequencies: indices 1, 2, ..., N/2  (Julia 1-indexed)
    v_hat_padded[1:N÷2] = v_hat[1:N÷2]
    # Negative frequencies: last N/2-1 entries
    v_hat_padded[end-N÷2+2:end] = v_hat[end-N÷2+2:end]

    # Handle Nyquist mode: split between +N/2 and -N/2
    v_hat_padded[N÷2+1]     = v_hat[N÷2+1] / 2
    v_hat_padded[M-N÷2+1]   = v_hat[N÷2+1] / 2

    # Inverse FFT and scale
    v_fine = real.(ifft(v_hat_padded)) .* q

    return v_fine
end

# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
function main()
    # Coarse grid parameters
    N = 32
    x_coarse = 2π .* collect(0:N-1) ./ N

    # Test function: exp(sin(x))
    v = exp.(sin.(x_coarse))

    # Interpolate via zero-padding
    q = 4
    M = q * N
    x_fine = 2π .* collect(0:M-1) ./ M
    v_fine = zero_pad_interpolate(v; q=q)

    # Exact function for comparison
    v_exact = exp.(sin.(x_fine))

    # Compute interpolation error
    error = maximum(abs.(v_fine .- v_exact))
    println("N = $N coarse points, M = $M fine points")
    println("Maximum interpolation error: $(round(error, sigdigits=2))")

    # Create figure
    fig = Figure(size = (600, 320))
    ax = Axis(fig[1, 1],
              xlabel = L"x",
              ylabel = L"\exp(\sin x)",
              title  = "Periodic band-limited interpolation via FFT zero-padding",
              limits = ((0, 2π), nothing))

    # Plot interpolant
    lines!(ax, x_fine, v_fine, color = NAVY, linewidth = 1.2,
           label = "Band-limited interpolant")

    # Plot coarse samples
    scatter!(ax, x_coarse, v, color = :white, markersize = 8,
             strokecolor = :black, strokewidth = 2,
             label = "Coarse samples (N = $N)")

    # Legend
    axislegend(ax, position = :rt, framevisible = false, bgcolor = (:white, 0.9))

    # Save figure
    mkpath(OUTPUT_DIR)
    save(joinpath(OUTPUT_DIR, "zero_padding_interpolation.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "zero_padding_interpolation.png"), fig, px_per_unit = 4)

    println("Figure saved to: ", joinpath(OUTPUT_DIR, "zero_padding_interpolation.pdf"))

    return fig
end

# Run if executed as script
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
