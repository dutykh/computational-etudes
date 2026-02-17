# transport_variable.jl
#
# Chapter 10: Spectral PDE Solvers with Chebyshev and Fourier Grids
#
# Variable coefficient transport equation with periodic BCs using Fourier.
#
# PDE: u_t + c(x) * u_x = 0, 0 < x < 2*pi, t > 0 (periodic)
# Speed: c(x) = 0.3 + sin^2(x - 1)
# ICs: u(x, 0) = exp(-20*(x - pi)^2)
#
# Uses FFT for spatial differentiation and RK4 for time stepping.
#
# Generates Figure 10.8: Variable speed transport waterfall plot.
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

# Output path
SCRIPT_DIR = @__DIR__
OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch10", "julia")
mkpath(OUTPUT_DIR)

# -----------------------------------------------------------------------------
# Fourier differentiation
# -----------------------------------------------------------------------------
"""
    fourier_diff(v)

Compute the derivative of a periodic function using FFT.

`v` is sampled at N equispaced points on [0, 2*pi).
"""
function fourier_diff(v::AbstractVector)
    N = length(v)
    v_hat = fft(v)

    # Wavenumbers
    if N % 2 == 0
        k = [collect(0:N÷2-1); 0; collect(1-N÷2:-1)]
    else
        k = [collect(0:(N-1)÷2); collect(-(N-1)÷2:-1)]
    end

    w = real.(ifft(1im .* k .* v_hat))
    return w
end

# -----------------------------------------------------------------------------
# Transport equation solver
# -----------------------------------------------------------------------------
"""
    transport_variable(; N=256, tmax=12.0, n_snapshots=100)

Solve variable coefficient transport equation using Fourier spectral method
with RK4 time stepping.

Returns `(x, t_save, U_save)` where `U_save` is a matrix with snapshots as rows.
"""
function transport_variable(; N=256, tmax=12.0, n_snapshots=100)
    # Grid
    h = 2pi / N
    x = h .* collect(0:N-1)

    # Variable speed: c(x) = 0.3 + sin^2(x - 1)
    c = 0.3 .+ sin.(x .- 1).^2

    # Initial condition: Gaussian pulse
    u = exp.(-20.0 .* (x .- pi).^2)

    # CFL condition for variable coefficient
    c_max = maximum(c)
    dt = 0.4 * h / c_max
    nsteps = Int(ceil(tmax / dt))

    # Store snapshots
    save_every = max(1, nsteps ÷ n_snapshots)
    U_save = [copy(u)]
    t_save = [0.0]

    # RHS function for RK4
    rhs(u) = -c .* fourier_diff(u)

    t = 0.0
    for n in 1:nsteps
        # RK4 time stepping
        k1 = rhs(u)
        k2 = rhs(u .+ 0.5 .* dt .* k1)
        k3 = rhs(u .+ 0.5 .* dt .* k2)
        k4 = rhs(u .+ dt .* k3)
        u = u .+ (dt / 6.0) .* (k1 .+ 2.0 .* k2 .+ 2.0 .* k3 .+ k4)

        t += dt

        # Save snapshot
        if n % save_every == 0
            push!(U_save, copy(u))
            push!(t_save, t)
        end
    end

    # Convert to matrix form
    U_mat = reduce(vcat, [u' for u in U_save])
    return x, t_save, U_mat
end

# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
function main()
    # Solve transport equation
    N = 256
    tmax = 8.0
    x, t_save, U_save = transport_variable(N=N, tmax=tmax, n_snapshots=120)

    # Create figure with two panels
    fig = Figure(size = (960, 400))

    # Left panel: space-time contour plot
    ax_main = Axis(fig[1, 1],
                   xlabel = L"x",
                   ylabel = L"t",
                   title  = "Variable Coefficient Transport")

    hm = heatmap!(ax_main, x, t_save, U_save, colormap = :coolwarm,
                  interpolate = true)
    Colorbar(fig[1, 1, Right()], hm, label = L"u(x, t)")

    # Right panel: waterfall view
    ax_water = Axis(fig[1, 2],
                    xlabel = L"x",
                    ylabel = L"u(x, t)" * " (offset by " * L"t" * ")",
                    title  = "Waterfall View: Pulse in Variable Speed")

    n_lines = 40
    n_total = length(t_save)
    indices = round.(Int, range(1, n_total, length=n_lines))

    for (i, idx) in enumerate(indices)
        offset = t_save[idx] * 0.05
        c_val = (i - 1) / (n_lines - 1)
        color = Makie.interpolated_getindex(Makie.to_colormap(:viridis), c_val)
        lines!(ax_water, x, U_save[idx, :] .+ offset, color = color, linewidth = 0.8)
    end

    # Add inset showing the speed profile
    ax_inset = Axis(fig[1, 2],
                    width = Relative(0.35), height = Relative(0.25),
                    halign = 0.95, valign = 0.95,
                    xlabel = L"x", ylabel = L"c(x)",
                    xlabelsize = 8, ylabelsize = 8,
                    xticklabelsize = 7, yticklabelsize = 7,
                    titlesize = 9, title = "Speed Profile",
                    backgroundcolor = :white)
    c_speed = 0.3 .+ sin.(x .- 1).^2
    lines!(ax_inset, x, c_speed, color = CORAL, linewidth = 1.5)

    hidespines!(ax_water, :t, :r)

    # Save figure
    save(joinpath(OUTPUT_DIR, "transport_variable.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "transport_variable.png"), fig, px_per_unit = 4)
    println("Figure saved to: $(joinpath(OUTPUT_DIR, "transport_variable.pdf"))")
end

main()
