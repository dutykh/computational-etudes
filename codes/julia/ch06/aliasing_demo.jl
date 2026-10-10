# aliasing_demo.jl
#
# Visualization of the aliasing phenomenon (Theorem 2: Poisson summation).
#
# When a function is sampled at N equispaced points, frequencies that differ
# by multiples of N become indistinguishable. High-frequency content "folds"
# onto low frequencies, corrupting the discrete Fourier transform.
#
# The aliasing formula:
#     f_tilde_k = sum_{j=-infty}^{infty} f_hat_{k + jN}
#
# This demo shows:
#     1. A function with high-frequency content sampled coarsely
#     2. The "folding" diagram showing how frequencies overlap
#     3. Comparison of true spectrum vs aliased spectrum
#
# For smooth functions, aliasing is harmless because high-frequency
# coefficients are negligible. For non-smooth functions, aliasing
# introduces significant errors.
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
const NAVY   = colorant"#142D6E"
const SKY    = colorant"#7896D2"
const CORAL  = colorant"#E74C3C"
const TEAL   = colorant"#16A085"
const PURPLE = colorant"#9B59B6"

SCRIPT_DIR = @__DIR__
OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch06", "julia")
mkpath(OUTPUT_DIR)

# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
function main()
    """Create a three-panel figure demonstrating aliasing."""

    fig = Figure(size = (900, 310))

    # -------------------------------------------------------------------------
    # Panel 1: Physical space - coarse sampling misses oscillations
    # -------------------------------------------------------------------------
    ax1 = Axis(fig[1, 1],
               xlabel = L"x",
               ylabel = L"f(x)",
               title  = "Coarse Sampling Misses High Frequencies",
               xlabelpadding = 3, ylabelpadding = 3)

    # A function with both low and high frequency content
    f(x) = sin(x) + 0.5 * sin(10 * x)

    x_fine = range(0, 2pi, length=500)
    N_coarse = 12  # Coarse sampling
    x_coarse = range(0, 2pi, length=N_coarse+1)[1:N_coarse]  # endpoint=false

    # Plot fine function
    lines!(ax1, collect(x_fine), f.(x_fine),
           color = NAVY, linewidth = 1.5, label = "True function")

    # Plot coarse samples
    scatter!(ax1, collect(x_coarse), f.(x_coarse),
             color = CORAL, markersize = 8, label = L"Sampled ($N=12$)")

    # Interpolant through coarse samples (trigonometric)
    # The aliased version sees sin(10x) as sin(10x - 12x) = sin(-2x)
    # So the perceived function is sin(x) + 0.5*sin(-2x) = sin(x) - 0.5*sin(2x)
    f_aliased(x) = sin(x) - 0.5 * sin(2 * x)

    lines!(ax1, collect(x_fine), f_aliased.(x_fine),
           color = TEAL, linewidth = 1.5, linestyle = :dash,
           label = "Aliased interpolant")

    xlims!(ax1, 0, 2pi)
    ylims!(ax1, -1.8, 1.8)

    axislegend(ax1, position = :rt, framevisible = true,
               framecolor = :transparent, bgcolor = RGBAf(1, 1, 1, 0.9),
               padding = (6, 6, 3, 3), labelsize = 9)

    # -------------------------------------------------------------------------
    # Panel 2: Frequency folding schematic
    # -------------------------------------------------------------------------
    ax2 = Axis(fig[1, 2],
               xlabel = L"Wavenumber $k$",
               title  = L"Frequency Folding ($N=12$)",
               xlabelpadding = 3)

    N = 12  # Number of sample points
    nyquist = N ÷ 2  # Nyquist frequency

    # Draw the frequency axis
    hlines!(ax2, [0.0], color = :gray, linewidth = 0.5)

    # Mark the Nyquist boundaries
    vlines!(ax2, [-nyquist, nyquist], color = :gray, linestyle = :dash, linewidth = 1)

    # Shade the resolved region
    band!(ax2, [-nyquist, nyquist], [-1.0, -1.0], [1.0, 1.0],
          color = (TEAL, 0.2))
    text!(ax2, 0, 0.8, text = "Resolved\nfrequencies", align = (:center, :center),
          fontsize = 10, color = TEAL)

    # Show folding arrows for k=10
    k_original = 10
    k_aliased = k_original - N  # 10 - 12 = -2

    # Arrow from k=10 to k=-2
    arrows!(ax2, [k_original], [0.3], [k_aliased - k_original], [0.0],
            color = CORAL, linewidth = 2, arrowsize = 10)
    scatter!(ax2, [k_original], [0.3], color = CORAL, markersize = 8)
    text!(ax2, k_original, 0.45, text = L"k=10", align = (:center, :bottom),
          fontsize = 10, color = CORAL)
    text!(ax2, (k_original + k_aliased) / 2, 0.5, text = "folds to",
          align = (:center, :bottom), fontsize = 9, color = CORAL)

    # Arrow from k=-10 to k=2
    k_original2 = -10
    k_aliased2 = k_original2 + N  # -10 + 12 = 2

    arrows!(ax2, [k_original2], [-0.3], [k_aliased2 - k_original2], [0.0],
            color = PURPLE, linewidth = 2, arrowsize = 10)
    scatter!(ax2, [k_original2], [-0.3], color = PURPLE, markersize = 8)
    text!(ax2, k_original2, -0.45, text = L"k=-10", align = (:center, :top),
          fontsize = 10, color = PURPLE)

    # Labels for Nyquist boundaries
    text!(ax2, -nyquist, -0.7, text = L"-N/2 = -6", align = (:center, :center),
          fontsize = 10, color = :gray)
    text!(ax2, nyquist, -0.7, text = L"N/2 = 6", align = (:center, :center),
          fontsize = 10, color = :gray)

    xlims!(ax2, -18, 18)
    ylims!(ax2, -1, 1)
    hideydecorations!(ax2)

    # -------------------------------------------------------------------------
    # Panel 3: True vs aliased spectrum
    # -------------------------------------------------------------------------
    ax3 = Axis(fig[1, 3],
               xlabel = L"Wavenumber $k$",
               ylabel = L"|\hat{g}_k|",
               title  = "True vs. Aliased Spectrum",
               yscale = log10,
               xlabelpadding = 3, ylabelpadding = 3)

    # A function with significant high-frequency content
    g(x) = 1.0 / (1.5 + cos(x))

    # Compute "true" spectrum with many points
    N_fine = 128
    x_fine_s = range(0, 2pi, length=N_fine+1)[1:N_fine]
    g_fine = g.(x_fine_s)
    g_hat_fine = fft(g_fine) / N_fine

    # Compute aliased spectrum with few points
    N_coarse_s = 16
    x_coarse_s = range(0, 2pi, length=N_coarse_s+1)[1:N_coarse_s]
    g_coarse = g.(x_coarse_s)
    g_hat_coarse = fft(g_coarse) / N_coarse_s

    # Plot true spectrum (positive frequencies only)
    k_fine_s = 0:N_fine÷2
    scatterlines!(ax3, collect(k_fine_s), abs.(g_hat_fine[1:N_fine÷2+1]),
                  color = NAVY, markersize = 4, marker = :circle, linewidth = 1,
                  label = L"True spectrum ($N=128$)")

    # Plot aliased spectrum
    k_coarse_s = 0:N_coarse_s÷2
    scatterlines!(ax3, collect(k_coarse_s), abs.(g_hat_coarse[1:N_coarse_s÷2+1]),
                  color = CORAL, markersize = 6, marker = :rect, linewidth = 1.5,
                  label = L"Aliased spectrum ($N=16$)")

    # Mark Nyquist frequency
    vlines!(ax3, [N_coarse_s ÷ 2], color = :gray, linestyle = :dash, linewidth = 1)
    text!(ax3, N_coarse_s ÷ 2 + 1, 1e-1,
          text = "Nyquist\n(N/2=8)", fontsize = 9, color = :gray)

    xlims!(ax3, -1, 50)
    ylims!(ax3, 1e-10, 1e1)

    axislegend(ax3, position = :rt, framevisible = true,
               framecolor = :transparent, bgcolor = RGBAf(1, 1, 1, 0.9),
               padding = (6, 6, 3, 3), labelsize = 9)

    # Clean styling for all panels
    for ax in [ax1, ax2, ax3]
        hidespines!(ax, :t, :r)
        ax.leftspinecolor   = RGBf(0.267, 0.267, 0.267)
        ax.bottomspinecolor = RGBf(0.267, 0.267, 0.267)
        ax.xtickcolor = RGBf(0.267, 0.267, 0.267)
        ax.ytickcolor = RGBf(0.267, 0.267, 0.267)
    end

    # Save figure
    save(joinpath(OUTPUT_DIR, "aliasing_visualization.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "aliasing_visualization.png"), fig, px_per_unit = 4)
    println("Figure saved to: $(joinpath(OUTPUT_DIR, "aliasing_visualization.pdf"))")
end

main()
