# chebfft_accuracy.jl
#
# Chapter 10: Spectral PDE Solvers with Chebyshev and Fourier Grids
#
# Creates Figure 10.2: Demonstration of spectral accuracy for Chebyshev
# differentiation via FFT.
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

# Import chebfft function
include(joinpath(@__DIR__, "chebfft.jl"))

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
const NAVY  = colorant"#142D6E"
const SKY   = colorant"#7896D2"
const CORAL = colorant"#E74C3C"
const TEAL  = colorant"#16A085"

# Output path
SCRIPT_DIR = @__DIR__
OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch10", "julia")
mkpath(OUTPUT_DIR)

# -----------------------------------------------------------------------------
# Test functions
# -----------------------------------------------------------------------------
f(x) = exp(x) * cos(4x)
f_prime(x) = exp(x) * (cos(4x) - 4 * sin(4x))

# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
function main()
    fig = Figure(size = (960, 280))

    # Fine grid for plotting
    x_fine = range(-1, 1, length=500)

    # =========================================================================
    # Panel 1: Function with Chebyshev nodes
    # =========================================================================
    ax1 = Axis(fig[1, 1],
               xlabel = L"x",
               ylabel = L"f(x)",
               title  = "Function and Chebyshev nodes")

    N = 16
    x_cheb = [cos(pi * j / N) for j in 0:N]
    f_cheb = f.(x_cheb)

    lines!(ax1, collect(x_fine), f.(collect(x_fine)), color = NAVY, linewidth = 1.5,
           label = L"f(x) = e^x \cos(4x)")
    scatter!(ax1, x_cheb, f_cheb, color = CORAL, markersize = 8,
             label = "Chebyshev points (N=$N)")

    xlims!(ax1, -1.1, 1.1)
    axislegend(ax1, position = :lt, labelsize = 9)

    # =========================================================================
    # Panel 2: Derivative comparison for N=16
    # =========================================================================
    ax2 = Axis(fig[1, 2],
               xlabel = L"x",
               ylabel = L"f'(x)",
               title  = "Derivative comparison (N=$N)")

    f_prime_exact = f_prime.(x_cheb)
    f_prime_fft = chebfft(f_cheb)

    lines!(ax2, collect(x_fine), f_prime.(collect(x_fine)), color = NAVY, linewidth = 1.5,
           label = "Exact")
    scatter!(ax2, x_cheb, f_prime_fft, color = CORAL, markersize = 8,
             label = "FFT")

    xlims!(ax2, -1.1, 1.1)
    axislegend(ax2, position = :lt, labelsize = 9)

    # =========================================================================
    # Panel 3: Convergence with N
    # =========================================================================
    ax3 = Axis(fig[1, 3],
               xlabel = L"N",
               ylabel = "Maximum error",
               title  = "Spectral convergence",
               yscale = log10)

    N_values = 4:2:30
    errors = Float64[]

    for Nv in N_values
        xc = [cos(pi * j / Nv) for j in 0:Nv]
        fc = f.(xc)
        fp_exact = f_prime.(xc)
        fp_fft = chebfft(fc)
        err = maximum(abs.(fp_fft .- fp_exact))
        push!(errors, err)
    end

    scatterlines!(ax3, collect(N_values), errors, color = NAVY, linewidth = 1.5,
                  markersize = 7)

    # Reference line
    hlines!(ax3, [1e-14], color = CORAL, linestyle = :dash, linewidth = 0.8,
            label = "Machine precision")

    ylims!(ax3, 1e-16, 1)
    axislegend(ax3, position = :rt, labelsize = 9)

    # Save
    save(joinpath(OUTPUT_DIR, "chebfft_accuracy.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "chebfft_accuracy.png"), fig, px_per_unit = 4)
    println("Figure saved to: $(joinpath(OUTPUT_DIR, "chebfft_accuracy.pdf"))")
end

main()
