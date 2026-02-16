# navier_stokes_2d.jl
#
# Chapter 11: Fourier Pseudospectral Methods for Periodic PDEs
#
# Fourier pseudospectral solver for the 2D incompressible Navier-Stokes
# equations in vorticity-streamfunction form on [0, 2*pi)^2, periodic:
#
#     omega_t + J(psi, omega) = nu * Delta omega
#     Delta psi = -omega
#
# where J(psi, omega) = psi_x omega_y - psi_y omega_x is the Jacobian.
#
# Uses RK4 with 2/3 dealiasing. The Poisson equation for the streamfunction
# is solved trivially in Fourier space: psi_hat = omega_hat / (kx^2 + ky^2).
#
# Generates two figures:
#   - ns2d_vorticity.pdf: vorticity snapshots at several times
#   - ns2d_energy.pdf: energy and enstrophy decay curves
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
using Random
using Statistics

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
# 2D Navier-Stokes solver
# -----------------------------------------------------------------------------
"""
    ns2d_fourier(; N=128, nu=1e-3, t_final=20.0, dt=0.01)

Solve 2D Navier-Stokes in vorticity-streamfunction form:
    omega_t + J(psi, omega) = nu * laplacian(omega)
    laplacian(psi) = -omega
on [0, 2*pi)^2 with periodic BCs and RK4.

Returns `(x, y, t_save, omega_save, energy_save, enstrophy_save)`.
"""
function ns2d_fourier(; N::Int=128, nu::Float64=1e-3,
                        t_final::Float64=20.0, dt::Float64=0.01)
    L = 2.0 * pi

    # 1D grid and wavenumbers
    x = L .* collect(0:N-1) ./ N
    y = copy(x)
    if N % 2 == 0
        kx_1d = Float64.([collect(0:N÷2-1); 0; collect(1-N÷2:-1)])
    else
        kx_1d = Float64.([collect(0:(N-1)÷2); collect(-(N-1)÷2:-1)])
    end
    ky_1d = copy(kx_1d)

    # 2D wavenumber grids (indexing: [ix, iy])
    KX = [kx_1d[i] for i in 1:N, j in 1:N]
    KY = [ky_1d[j] for i in 1:N, j in 1:N]
    K2 = KX.^2 .+ KY.^2
    K2[1, 1] = 1.0   # Avoid division by zero for mean mode

    # 2/3 dealiasing mask (2D)
    cutoff = N ÷ 3
    mask = trues(N, N)
    if cutoff > 0
        mask_1d = trues(N)
        mask_1d[cutoff+2:N-cutoff] .= false
        for i in 1:N, j in 1:N
            mask[i, j] = mask_1d[i] && mask_1d[j]
        end
    end

    # Initial condition: random vorticity with zero mean
    Random.seed!(42)
    omega_hat = zeros(ComplexF64, N, N)
    for i in 1:N, j in 1:N
        kk = sqrt(KX[i, j]^2 + KY[i, j]^2)
        if 1.0 < kk < 10.0
            amp = 1.0 / (1.0 + kk)
            phase = 2.0 * pi * rand()
            omega_hat[i, j] = amp * exp(1im * phase)
        end
    end
    # Enforce real-valued omega
    omega = real.(ifft(omega_hat))
    omega .-= mean2d(omega)   # Ensure zero mean
    omega ./= std(vec(omega))  # Normalize to unit RMS vorticity

    # Diffusive symbol
    diff_symbol = -nu .* K2

    function rhs(omega_phys)
        """RHS: -J(psi, omega) + nu * laplacian(omega)."""
        omega_hat_loc = fft(omega_phys)

        # Poisson solve: psi_hat = omega_hat / K2
        psi_hat = omega_hat_loc ./ K2
        psi_hat[1, 1] = 0.0   # Zero mean streamfunction

        # Velocity: u = psi_y, v = -psi_x
        u_vel = real.(ifft(1im .* KY .* psi_hat))
        v_vel = real.(ifft(-1im .* KX .* psi_hat))

        # Vorticity gradients
        omega_x = real.(ifft(1im .* KX .* omega_hat_loc))
        omega_y = real.(ifft(1im .* KY .* omega_hat_loc))

        # Jacobian (advection) in physical space
        jac = u_vel .* omega_x .+ v_vel .* omega_y

        # Transform, dealiase, add diffusion
        jac_hat = fft(jac)
        jac_hat[.!mask] .= 0.0

        f_hat = -jac_hat .+ diff_symbol .* omega_hat_loc
        return real.(ifft(f_hat))
    end

    # Time stepping
    n_steps = Int(ceil(t_final / dt))
    dt = t_final / n_steps

    n_save = 200
    save_every = max(1, n_steps ÷ n_save)
    omega_save    = [copy(omega)]
    t_save        = [0.0]
    energy_save   = Float64[]
    enstrophy_save = Float64[]

    function compute_diagnostics(omega_phys)
        omega_hat_loc = fft(omega_phys)
        psi_hat = omega_hat_loc ./ K2
        psi_hat[1, 1] = 0.0
        u_vel = real.(ifft(1im .* KY .* psi_hat))
        v_vel = real.(ifft(-1im .* KX .* psi_hat))
        energy   = 0.5 * mean2d(u_vel.^2 .+ v_vel.^2)
        enstrophy = 0.5 * mean2d(omega_phys.^2)
        return energy, enstrophy
    end

    E, Z = compute_diagnostics(omega)
    push!(energy_save, E)
    push!(enstrophy_save, Z)

    t = 0.0
    for n in 1:n_steps
        k1 = rhs(omega)
        k2 = rhs(omega .+ 0.5 .* dt .* k1)
        k3 = rhs(omega .+ 0.5 .* dt .* k2)
        k4 = rhs(omega .+ dt .* k3)
        omega .= omega .+ (dt / 6.0) .* (k1 .+ 2.0 .* k2 .+ 2.0 .* k3 .+ k4)
        t += dt

        if n % save_every == 0
            push!(omega_save, copy(omega))
            push!(t_save, t)
            E, Z = compute_diagnostics(omega)
            push!(energy_save, E)
            push!(enstrophy_save, Z)
        end
    end

    return x, y, t_save, omega_save, energy_save, enstrophy_save
end

"""Helper: mean over a 2D array."""
mean2d(A) = sum(A) / length(A)

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
function main()
    N = 128
    nu = 1e-3
    t_final = 20.0

    x, y, t_save, omega_save, energy, enstrophy = ns2d_fourier(
        N=N, nu=nu, t_final=t_final, dt=0.005
    )

    # -------------------------------------------------------------------------
    # Figure 1: Vorticity snapshots
    # -------------------------------------------------------------------------
    fig1 = Figure(size = (800, 780))

    t_targets = [0.0, 5.0, 10.0, t_final]
    vmax = Float64(quantile(abs.(vec(omega_save[1])), 0.95))

    for (i, t_target) in enumerate(t_targets)
        row = (i - 1) ÷ 2 + 1
        col = (i - 1) % 2 + 1

        idx = argmin(abs.(t_save .- t_target))
        omega_snap = omega_save[idx]

        ax = Axis(fig1[row, col],
                  xlabel = L"x",
                  ylabel = L"y",
                  title  = latexstring("t = $(round(Int, t_save[idx]))"),
                  aspect = DataAspect())

        hm = heatmap!(ax, x, y, omega_snap, colormap = :RdBu,
                      colorrange = (-vmax, vmax), interpolate = true)
        Colorbar(fig1[row, col][1, 2], hm)
    end

    Label(fig1[0, :],
          "2D Navier--Stokes: Decaying Turbulence (\$N = $N\$, \$\\nu = $nu\$)",
          fontsize = 13)

    save(joinpath(OUTPUT_DIR, "ns2d_vorticity.pdf"), fig1)
    save(joinpath(OUTPUT_DIR, "ns2d_vorticity.png"), fig1, px_per_unit = 4)
    println("Figure saved: $(joinpath(OUTPUT_DIR, "ns2d_vorticity.pdf"))")

    # -------------------------------------------------------------------------
    # Figure 2: Energy and enstrophy decay
    # -------------------------------------------------------------------------
    fig2 = Figure(size = (960, 420))

    t_diag = t_save[1:length(energy)]

    ax1 = Axis(fig2[1, 1],
               xlabel = L"t",
               ylabel = L"E(t) = \frac{1}{2}\langle |\mathbf{u}|^2 \rangle",
               title  = "Kinetic Energy Decay")

    lines!(ax1, t_diag, energy, color = NAVY, linewidth = 1.5)
    hidespines!(ax1, :t, :r)

    ax2 = Axis(fig2[1, 2],
               xlabel = L"t",
               ylabel = L"Z(t) = \frac{1}{2}\langle \omega^2 \rangle",
               title  = "Enstrophy Decay")

    lines!(ax2, t_diag, enstrophy, color = CORAL, linewidth = 1.5)
    hidespines!(ax2, :t, :r)

    save(joinpath(OUTPUT_DIR, "ns2d_energy.pdf"), fig2)
    save(joinpath(OUTPUT_DIR, "ns2d_energy.png"), fig2, px_per_unit = 4)
    println("Figure saved: $(joinpath(OUTPUT_DIR, "ns2d_energy.pdf"))")
end

main()
