# advection_fourier.jl
#
# Chapter 11: Fourier Pseudospectral Methods for Periodic PDEs
#
# Fourier pseudospectral solver for the linear advection equation:
#     u_t + c * u_x = 0 on [0, 2*pi), periodic BCs.
#
# Compares Fourier pseudospectral + RK4 with second-order upwind FD + RK4
# to demonstrate the dispersion-free nature of spectral advection.
#
# Initial condition: u(x, 0) = exp(-10*(x - pi)^2)
#
# Generates figure: advection_comparison.pdf
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
const NAVY  = colorant"#142D6E"
const SKY   = colorant"#7896D2"
const CORAL = colorant"#E74C3C"
const TEAL  = colorant"#16A085"

# Output path
SCRIPT_DIR = @__DIR__
OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch11", "julia")
mkpath(OUTPUT_DIR)

# -----------------------------------------------------------------------------
# Right-hand side functions
# -----------------------------------------------------------------------------
"""
    rhs_fourier(u, k, c)

Right-hand side for Fourier pseudospectral advection: -c * u_x.
"""
function rhs_fourier(u::Vector{Float64}, k::Vector{Float64}, c::Float64)
    u_hat = fft(u)
    ux = real.(ifft(1im .* k .* u_hat))
    return -c .* ux
end

"""
    rhs_fd_upwind2(u, dx, c)

Right-hand side for second-order upwind FD advection: -c * u_x.
"""
function rhs_fd_upwind2(u::Vector{Float64}, dx::Float64, c::Float64)
    N = length(u)
    ux = zeros(N)
    if c > 0
        # Second-order upwind (backward biased)
        for j in 1:N
            jm1 = mod(j - 2, N) + 1
            jm2 = mod(j - 3, N) + 1
            ux[j] = (3.0 * u[j] - 4.0 * u[jm1] + u[jm2]) / (2.0 * dx)
        end
    else
        for j in 1:N
            jp1 = mod(j, N) + 1
            jp2 = mod(j + 1, N) + 1
            ux[j] = (-3.0 * u[j] + 4.0 * u[jp1] - u[jp2]) / (2.0 * dx)
        end
    end
    return -c .* ux
end

"""
    rk4_step(u, dt, rhs_func)

One step of classical RK4.
"""
function rk4_step(u::Vector{Float64}, dt::Float64, rhs_func::Function)
    k1 = rhs_func(u)
    k2 = rhs_func(u .+ 0.5 .* dt .* k1)
    k3 = rhs_func(u .+ 0.5 .* dt .* k2)
    k4 = rhs_func(u .+ dt .* k3)
    return u .+ (dt / 6.0) .* (k1 .+ 2.0 .* k2 .+ 2.0 .* k3 .+ k4)
end

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
function main()
    # Parameters
    N = 128
    c = 1.0
    L = 2.0 * pi
    t_final = 2.0 * pi / c   # One full period

    # Grid
    x = L .* collect(0:N-1) ./ N
    dx = L / N

    # Wavenumbers
    k = Float64.([collect(0:N÷2-1); 0; collect(1-N÷2:-1)])

    # Initial condition (smooth localised bump)
    u0 = exp.(-10.0 .* (x .- pi).^2)

    # Time step (CFL-based)
    k_max = N / 2
    dt = 0.5 / (c * k_max)   # CFL ~ 0.5
    n_steps = Int(ceil(t_final / dt))
    dt = t_final / n_steps

    # Solve with Fourier pseudospectral
    u_fourier = copy(u0)
    for _ in 1:n_steps
        u_fourier = rk4_step(u_fourier, dt, u -> rhs_fourier(u, k, c))
    end

    # Solve with second-order upwind FD (same time step and RK4)
    u_fd = copy(u0)
    for _ in 1:n_steps
        u_fd = rk4_step(u_fd, dt, u -> rhs_fd_upwind2(u, dx, c))
    end

    # Exact solution (one full period = initial condition)
    u_exact = copy(u0)

    # -------------------------------------------------------------------------
    # Figure
    # -------------------------------------------------------------------------
    fig = Figure(size = (960, 420))

    # Left: overlay of solutions
    ax1 = Axis(fig[1, 1],
               xlabel = L"x",
               ylabel = L"u(x, t)",
               title  = "After one period (\$t = 2\\pi / c\$, \$N = $N\$)")

    lines!(ax1, x, u_exact, color = (:gray, 0.5), linewidth = 2.5,
           label = "Exact (initial)")
    lines!(ax1, x, u_fourier, color = NAVY, linewidth = 1.5,
           label = "Fourier PS + RK4")
    lines!(ax1, x, u_fd, color = CORAL, linewidth = 1.5, linestyle = :dash,
           label = "Upwind FD2 + RK4")

    xlims!(ax1, 0, L)
    axislegend(ax1, position = :rt, framevisible = true)
    hidespines!(ax1, :t, :r)

    # Right: error profiles
    err_fourier = abs.(u_fourier .- u_exact)
    err_fd      = abs.(u_fd .- u_exact)

    ax2 = Axis(fig[1, 2],
               xlabel = L"x",
               ylabel = L"|u^{\mathrm{num}} - u^{\mathrm{exact}}|",
               title  = "Pointwise Error After One Period",
               yscale = log10)

    lines!(ax2, x, err_fd, color = CORAL, linewidth = 1.5,
           label = "Upwind FD2 (max err: $(round(maximum(err_fd), sigdigits=2)))")
    lines!(ax2, x, max.(err_fourier, 1e-16), color = NAVY, linewidth = 1.5,
           label = "Fourier PS (max err: $(round(maximum(err_fourier), sigdigits=2)))")

    xlims!(ax2, 0, L)
    axislegend(ax2, position = :cb, framevisible = true, labelsize = 9)
    hidespines!(ax2, :t, :r)

    # Save
    save(joinpath(OUTPUT_DIR, "advection_comparison.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "advection_comparison.png"), fig, px_per_unit = 4)
    println("Figure saved: $(joinpath(OUTPUT_DIR, "advection_comparison.pdf"))")
end

main()
