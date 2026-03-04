# burgers_fourier.jl
#
# Chapter 11: Fourier Pseudospectral Methods for Periodic PDEs
#
# Fourier pseudospectral solver for the viscous Burgers equation:
#     u_t + u * u_x = nu * u_xx on [0, 2*pi), periodic BCs.
#
# Initial condition: u(x, 0) = sin(x)
# Features: RK4 time stepping, 2/3 dealiasing rule.
#
# Generates two figures:
#   - burgers_evolution.pdf: solution evolution showing front steepening
#   - burgers_dealiasing.pdf: comparison of energy spectra with/without dealiasing
#
# Author: Dr. Denys Dutykh
#         Mathematics Department
#         Khalifa University of Science and Technology
#         Abu Dhabi, UAE
#
# Part of the book "Computational Etudes: A Spectral Approach"
# https://github.com/dutykh/computational-etudes

using CairoMakie
using LaTeXStrings
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
# Burgers equation solver
# -----------------------------------------------------------------------------
"""
    burgers_fourier(; N=256, nu=1e-2, t_final=1.0, CFL=0.4, dealias=true)

Fourier pseudospectral solver for the viscous Burgers equation
    u_t + u u_x = nu u_xx
on [0, 2*pi) with periodic boundary conditions.

Returns `(x, t_save, U_save)` where `U_save` is a matrix with snapshots as rows.
"""
function burgers_fourier(; N::Int=256, nu::Float64=1e-2, t_final::Float64=1.0,
                          CFL::Float64=0.4, dealias::Bool=true)
    # Grid and wavenumbers
    x = 2.0 * pi .* collect(0:N-1) ./ N
    if N % 2 == 0
        k = Float64.([collect(0:N÷2-1); 0; collect(1-N÷2:-1)])
    else
        k = Float64.([collect(0:(N-1)÷2); collect(-(N-1)÷2:-1)])
    end

    # Initial condition
    u = sin.(x)

    # Time step from linearised stability
    k_max = N / 2.0
    dt_adv  = CFL / k_max
    dt_diff = nu > 0 ? 0.4 / (nu * k_max^2) : Inf
    dt = min(dt_adv, dt_diff)
    n_steps = Int(ceil(t_final / dt))
    dt = t_final / n_steps

    # Diffusive symbol
    diff_symbol = -nu .* k.^2

    # 2/3 dealiasing mask
    dealias_mask = trues(N)
    if dealias
        cutoff = N ÷ 3
        if cutoff > 0
            dealias_mask[cutoff+2:N-cutoff] .= false
        end
    end

    function rhs(u_phys)
        u_hat = fft(u_phys)
        ux = real.(ifft(1im .* k .* u_hat))

        # Nonlinear term in physical space, then dealiase
        f_nl = -u_phys .* ux
        f_nl_hat = fft(f_nl)
        f_nl_hat[.!dealias_mask] .= 0.0

        # Diffusive term in Fourier space
        f_diff_hat = diff_symbol .* u_hat

        return real.(ifft(f_nl_hat .+ f_diff_hat))
    end

    # Time integration with RK4
    n_save = min(200, n_steps)
    save_every = max(1, n_steps ÷ n_save)
    U_save = [copy(u)]
    t_save = [0.0]

    t = 0.0
    for n in 1:n_steps
        k1 = rhs(u)
        k2 = rhs(u .+ 0.5 .* dt .* k1)
        k3 = rhs(u .+ 0.5 .* dt .* k2)
        k4 = rhs(u .+ dt .* k3)
        u = u .+ (dt / 6.0) .* (k1 .+ 2.0 .* k2 .+ 2.0 .* k3 .+ k4)
        t += dt

        if n % save_every == 0
            push!(U_save, copy(u))
            push!(t_save, t)
        end
    end

    # Convert to matrix form (snapshots as rows)
    U_mat = reduce(vcat, [u' for u in U_save])
    return x, t_save, U_mat
end

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
function main()
    N = 256
    nu = 0.01
    t_final = 1.5

    # -------------------------------------------------------------------------
    # Figure 1: Solution evolution
    # -------------------------------------------------------------------------
    x, t_save, U_save = burgers_fourier(N=N, nu=nu, t_final=t_final, dealias=true)

    fig1 = Figure(size = (960, 450))

    # Left: space-time plot
    ax1 = Axis(fig1[1, 1],
               xlabel = L"x",
               ylabel = L"t",
               title  = "Burgers Equation (\$\\nu = $nu\$)")

    hm = heatmap!(ax1, x, t_save, U_save, colormap = :RdBu, interpolate = true)
    Colorbar(fig1[1, 1, Right()], hm, label = L"u(x, t)")

    # Right: snapshots
    ax2 = Axis(fig1[1, 2],
               xlabel = L"x",
               ylabel = L"u(x, t)",
               title  = "Snapshots: Front Steepening")

    snapshot_times = [0.0, 0.3, 0.6, 1.0, 1.5]
    colors_snap = [SKY, TEAL, colorant"#2ECC71", colorant"#F39C12", NAVY]

    for (i, t_target) in enumerate(snapshot_times)
        idx = argmin(abs.(t_save .- t_target))
        lines!(ax2, x, U_save[idx, :], color = colors_snap[i], linewidth = 1.5,
               label = latexstring("t = $(round(t_save[idx], digits=1))"))
    end

    xlims!(ax2, 0, 2pi)
    axislegend(ax2, position = :rt, framevisible = true, labelsize = 9)
    hidespines!(ax2, :t, :r)

    save(joinpath(OUTPUT_DIR, "burgers_evolution.pdf"), fig1)
    save(joinpath(OUTPUT_DIR, "burgers_evolution.png"), fig1, px_per_unit = 4)
    println("Figure saved: $(joinpath(OUTPUT_DIR, "burgers_evolution.pdf"))")

    # -------------------------------------------------------------------------
    # Figure 2: Dealiasing comparison
    # -------------------------------------------------------------------------
    t_compare = 1.0

    _, t1, U1 = burgers_fourier(N=N, nu=nu, t_final=t_compare, dealias=true)
    _, t2, U2 = burgers_fourier(N=N, nu=nu, t_final=t_compare, dealias=false)

    # Energy spectra at final time
    u_hat_deal   = fft(U1[end, :])
    u_hat_nodeal = fft(U2[end, :])
    spectrum_deal   = abs.(u_hat_deal[1:N÷2]).^2 ./ N^2
    spectrum_nodeal = abs.(u_hat_nodeal[1:N÷2]).^2 ./ N^2
    k_pos = collect(0:N÷2-1)

    fig2 = Figure(size = (960, 420))

    # Left: energy spectra
    ax3 = Axis(fig2[1, 1],
               xlabel = L"Wavenumber $k$",
               ylabel = L"|\hat{u}_k|^2 / N^2",
               title  = latexstring("Energy Spectrum at \$t = $t_compare\$"),
               yscale = log10)

    lines!(ax3, k_pos[2:end], spectrum_deal[2:end], color = NAVY, linewidth = 1.5,
           label = "With 2/3 dealiasing")
    lines!(ax3, k_pos[2:end], spectrum_nodeal[2:end], color = CORAL, linewidth = 1.5,
           linestyle = :dash, label = "Without dealiasing")

    axislegend(ax3, position = :rt, framevisible = true)
    hidespines!(ax3, :t, :r)

    # Right: physical space comparison
    ax4 = Axis(fig2[1, 2],
               xlabel = L"x",
               ylabel = L"u(x, t)",
               title  = latexstring("Physical Space at \$t = $t_compare\$"))

    lines!(ax4, x, U1[end, :], color = NAVY, linewidth = 1.5,
           label = "With dealiasing")
    lines!(ax4, x, U2[end, :], color = CORAL, linewidth = 1.5, linestyle = :dash,
           label = "Without dealiasing")

    xlims!(ax4, 0, 2pi)
    axislegend(ax4, position = :rt, framevisible = true)
    hidespines!(ax4, :t, :r)

    save(joinpath(OUTPUT_DIR, "burgers_dealiasing.pdf"), fig2)
    save(joinpath(OUTPUT_DIR, "burgers_dealiasing.png"), fig2, px_per_unit = 4)
    println("Figure saved: $(joinpath(OUTPUT_DIR, "burgers_dealiasing.pdf"))")
end

main()
