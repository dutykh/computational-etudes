# fft_aliasing.jl
#
# Chapter 9: Physical and Fourier Space on Grids
#
# Demonstrates aliasing in FFT of high-frequency data.
# sin(17x) sampled at N=32 points shows spike at k=-15 due to aliasing.
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
# Main visualization
# -----------------------------------------------------------------------------
function main()
    # Grid parameters
    N = 32
    x = 2π .* collect(0:N-1) ./ N

    # High-frequency function: k=17 > Nyquist limit N/2=16
    u = sin.(17 .* x)

    # Compute FFT
    u_hat = fft(u) ./ N

    # Wavenumber vector
    k = collect(0:N-1)
    k[k .> N÷2] .-= N   # Convert to -N/2+1, ..., N/2

    # Shift for plotting
    k_shifted     = fftshift(k)
    u_hat_shifted = fftshift(u_hat)

    # Find the aliased location
    println("N = $N, Nyquist limit = $(N÷2)")
    println("True wavenumber: 17")
    println("Alias: 17 = -15 + 32, so appears at k = -15")
    idx_max = argmax(abs.(u_hat_shifted))
    println("Maximum |u_hat_k| at k = $(k_shifted[idx_max])")

    # Create figure
    fig = Figure(size = (600, 320))
    ax = Axis(fig[1, 1],
              xlabel = L"k$ (wavenumber)$",
              ylabel = L"|\hat{u}_k|",
              title  = L"FFT of $\sin(17x)$ with $N = 32$ points: aliasing to $k = -15$",
              limits = ((-N÷2 - 1, N÷2 + 1), (0, 0.6)))

    # Stem plot of spectrum
    stem!(ax, k_shifted, abs.(u_hat_shifted),
          color = NAVY, stemcolor = NAVY, markersize = 6)

    # Highlight the aliased peak
    aliased_k   = -15
    aliased_idx = findfirst(k_shifted .== aliased_k)
    scatter!(ax, [aliased_k], [abs(u_hat_shifted[aliased_idx])],
             color = CORAL, markersize = 10,
             label = "Aliased peak at k = $aliased_k")

    # Add annotation arrow
    text!(ax, aliased_k + 5, 0.4,
          text = L"$k = 17$ aliases to $k = -15$",
          fontsize = 10)
    arrows!(ax, [aliased_k + 4.5], [0.39], [-4.0], [-0.05],
            color = CORAL, linewidth = 1.0)

    # Legend
    axislegend(ax, position = :rt, framevisible = false, bgcolor = (:white, 0.9))

    # Save figure
    mkpath(OUTPUT_DIR)
    save(joinpath(OUTPUT_DIR, "fft_aliasing.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "fft_aliasing.png"), fig, px_per_unit = 4)

    println("Figure saved to: ", joinpath(OUTPUT_DIR, "fft_aliasing.pdf"))

    return fig
end

# Run if executed as script
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
