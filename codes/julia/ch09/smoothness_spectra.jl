# smoothness_spectra.jl
#
# Chapter 9: Physical and Fourier Space on Grids
#
# Demonstrates the relationship between function smoothness and spectral decay.
# Shows spectra of: exp(sin(x)), periodic hat, and square wave.
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
# Compute spectrum of a periodic function
# -----------------------------------------------------------------------------
"""
    compute_spectrum(f_func; N=256)

Compute Fourier spectrum of a periodic function on [0, 2*pi).

# Returns
- `k`: wavenumbers (positive only, excluding k=0)
- `f_hat_mag`: magnitude of Fourier coefficients
"""
function compute_spectrum(f_func; N=256)
    x = 2π .* collect(0:N-1) ./ N
    f = f_func.(x)
    f_hat = fft(f) ./ N

    # Take positive wavenumbers only (excluding k=0)
    k         = collect(1:N÷2-1)
    f_hat_mag = abs.(f_hat[2:N÷2])

    return k, f_hat_mag
end

# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
function main()
    N = 256

    # Define three functions with different smoothness

    # 1. Analytic: exp(sin(x))
    f1(x) = exp(sin(x))

    # 2. C^0 but not C^1: periodic hat function centered at pi
    f2(x) = max(0.0, 1.0 - abs(x - π) / (π/2))

    # 3. Discontinuous: square wave
    f3(x) = sign(sin(x))

    # Compute spectra
    k1, spec1 = compute_spectrum(f1; N=N)
    k2, spec2 = compute_spectrum(f2; N=N)
    k3, spec3 = compute_spectrum(f3; N=N)

    # Create figure with 3 panels
    fig = Figure(size = (840, 280))

    # Panel 1: exp(sin(x)) - semilog plot (exponential decay)
    ax1 = Axis(fig[1, 1],
               xlabel = L"k",
               ylabel = L"|\hat{f}_k|",
               title  = L"$\exp(\sin x)$ (analytic)",
               yscale = log10,
               limits = ((0, N÷2), (1e-16, 10)))
    scatterlines!(ax1, k1, spec1, color = NAVY, markersize = 4, linewidth = 0.8)

    # Panel 2: Hat function - log-log plot (algebraic decay ~k^{-2})
    ax2 = Axis(fig[1, 2],
               xlabel = L"k",
               ylabel = L"|\hat{f}_k|",
               title  = L"Hat function ($C^0$)",
               xscale = log10, yscale = log10,
               limits = ((1, N÷2), nothing))
    scatterlines!(ax2, k2, spec2, color = NAVY, markersize = 4, linewidth = 0.8)
    # Reference line for k^{-2}
    k_ref = [5.0, 100.0]
    lines!(ax2, k_ref, 0.5 .* k_ref .^ (-2), color = CORAL, linewidth = 1.5,
           linestyle = :dash, label = L"k^{-2}")
    axislegend(ax2, position = :rt, framevisible = false,
               bgcolor = (:white, 0.9), labelsize = 9)

    # Panel 3: Square wave - log-log plot (algebraic decay ~k^{-1})
    ax3 = Axis(fig[1, 3],
               xlabel = L"k",
               ylabel = L"|\hat{f}_k|",
               title  = "Square wave (discontinuous)",
               xscale = log10, yscale = log10,
               limits = ((1, N÷2), nothing))
    scatterlines!(ax3, k3, spec3, color = NAVY, markersize = 4, linewidth = 0.8)
    # Reference line for k^{-1}
    k_ref2 = [1.0, 100.0]
    lines!(ax3, k_ref2, 1.5 .* k_ref2 .^ (-1), color = CORAL, linewidth = 1.5,
           linestyle = :dash, label = L"k^{-1}")
    axislegend(ax3, position = :rt, framevisible = false,
               bgcolor = (:white, 0.9), labelsize = 9)

    # Save figure
    mkpath(OUTPUT_DIR)
    save(joinpath(OUTPUT_DIR, "smoothness_spectra.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "smoothness_spectra.png"), fig, px_per_unit = 4)

    println("Figure saved to: ", joinpath(OUTPUT_DIR, "smoothness_spectra.pdf"))

    return fig
end

# Run if executed as script
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
