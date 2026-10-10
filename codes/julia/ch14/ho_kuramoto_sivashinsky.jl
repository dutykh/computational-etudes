# ho_kuramoto_sivashinsky.jl
#
# Chapter 14: Higher-Order Boundary Value Problems
#
# Computational Etude 14.9: Kuramoto-Sivashinsky equation.
#
#     u_t + u * u_x + u_xx + u_xxxx = 0
#
# on [0, 32 pi] with periodic boundary conditions. Solved using Fourier
# pseudospectral methods in space and the ETDRK4 scheme (Cox & Matthews
# 2002, with Kassam-Trefethen 2005 contour integral improvement) in time.
#
# The KS equation models laminar flame front instabilities and thin film
# flows. In Fourier space, the equation becomes:
#
#     du_hat/dt = L_k * u_hat + N_hat(u)
#
# where the linear operator L_k = k^2 - k^4 (energy production at low k,
# dissipation at high k) and the nonlinear term N_hat = -(ik/2) * FFT(u^2)
# uses the conservation form u * u_x = (u^2 / 2)_x.
#
# The ETDRK4 scheme evaluates the phi-functions via Kassam & Trefethen's
# contour integral trick on a circle in the complex plane, avoiding
# catastrophic cancellation near z = 0.
#
# The solution exhibits spatiotemporal chaos: cell merging, splitting, and
# irregular oscillations characteristic of fourth-order dissipative systems.
#
# Generates Figure 14.8: Space-time heatmap of the KS solution.
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
using LinearAlgebra
using FFTW
using Printf
using Statistics

# Publication-quality CairoMakie configuration
set_theme!(Theme(
    fontsize = 11,
    fonts = (regular = "CMU Serif", bold = "CMU Serif Bold", italic = "CMU Serif Italic"),
    Axis = (xlabelsize = 11, ylabelsize = 11, titlesize = 12,
            xticklabelsize = 10, yticklabelsize = 10,
            spinewidth = 0.8, xtickwidth = 0.8, ytickwidth = 0.8),
    Legend = (labelsize = 10,),
))

# Book color scheme
const NAVY   = colorant"#142D6E"
const SKY    = colorant"#7896D2"
const CORAL  = colorant"#E74C3C"
const TEAL   = colorant"#16A085"
const PURPLE = colorant"#8E44AD"
const ORANGE = colorant"#E67E22"

# Output paths
const SCRIPT_DIR = @__DIR__
const OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch14", "julia")

# -----------------------------------------------------------------------------
# ETDRK4 coefficients via Kassam-Trefethen contour integral trick
# -----------------------------------------------------------------------------
"""
    compute_etdrk4_coefficients(Lop, dt; M=32)

Compute ETDRK4 coefficients using the Kassam-Trefethen (2005) contour
integral method for numerical stability.

The standard ETDRK4 formulas involve ratios that suffer from catastrophic
cancellation when z = L_k * dt is near zero. The contour integral trick
evaluates these on a circle of radius 1 in the complex plane around each
L_k * dt, averaging over M roots of unity.

Returns (E, E2, Q, f1, f2, f3).
"""
function compute_etdrk4_coefficients(Lop::Vector{Float64}, dt::Float64; M::Int=32)
    Nk = length(Lop)

    # Full-step and half-step propagators
    E = exp.(dt .* Lop)
    E2 = exp.(dt .* Lop ./ 2.0)

    # Roots of unity for the contour integral
    r = exp.(2im * pi .* collect(1:M) ./ M)   # shape (M,)

    # LR[i, m] = dt * L_i + r_m, evaluated on a circle of radius 1
    # around each dt * L_i in the complex plane
    LR = dt .* Lop .+ r'   # shape (Nk, M) via broadcasting

    # Contour integral evaluations of the phi-functions:
    #
    # Q  = dt * phi_1(L*dt/2) = dt * (exp(L*dt/2) - 1) / (L*dt)
    #    ~ evaluated at half-step for the intermediate stages
    #
    # f1 = dt * (-4 - LR + exp(LR)*(4 - 3*LR + LR^2)) / LR^3
    # f2 = dt * ( 2 + LR + exp(LR)*(-2 + LR))           / LR^3
    # f3 = dt * (-4 - 3*LR - LR^2 + exp(LR)*(4 - LR))  / LR^3

    Q = dt .* real.(vec(mean((exp.(LR ./ 2.0) .- 1.0) ./ LR, dims=2)))

    f1 = dt .* real.(vec(mean(
        (-4.0 .- LR .+ exp.(LR) .* (4.0 .- 3.0 .* LR .+ LR.^2)) ./ LR.^3,
        dims=2)))

    f2 = dt .* real.(vec(mean(
        (2.0 .+ LR .+ exp.(LR) .* (-2.0 .+ LR)) ./ LR.^3,
        dims=2)))

    f3 = dt .* real.(vec(mean(
        (-4.0 .- 3.0 .* LR .- LR.^2 .+ exp.(LR) .* (4.0 .- LR)) ./ LR.^3,
        dims=2)))

    return E, E2, Q, f1, f2, f3
end

# -----------------------------------------------------------------------------
# Kuramoto-Sivashinsky solver with ETDRK4
# -----------------------------------------------------------------------------
"""
    solve_ks_etdrk4(; N=256, L=32.0*pi, t_final=150.0, dt=0.5, n_snapshots=500)

Solve the Kuramoto-Sivashinsky equation using Fourier pseudospectral
spatial discretization and ETDRK4 time integration.

The equation u_t + u*u_x + u_xx + u_xxxx = 0 is rewritten as:

    du_hat/dt = L_k * u_hat + N_hat(u)

where L_k = k^2 - k^4 and N_hat = -(ik/2) * FFT(u^2).

Returns (x, t_save, U_save).
"""
function solve_ks_etdrk4(; N::Int=256, L::Float64=32.0*pi,
                           t_final::Float64=150.0, dt::Float64=0.5,
                           n_snapshots::Int=500)
    # Spatial grid
    x = L .* collect(0:N-1) ./ N

    # Wavenumbers: k = (2*pi/L) * [0, 1, 2, ..., N/2-1, -N/2, ..., -1]
    k_base = collect(0:N-1)
    k_freq = [i <= N ÷ 2 ? i : i - N for i in k_base]
    k = (2.0 * pi / L) .* Float64.(k_freq)

    # Linear operator: L_k = k^2 - k^4
    Lop = k.^2 .- k.^4

    # 2/3 dealiasing mask for the quadratic nonlinearity u^2
    cutoff = N ÷ 3
    dealias = trues(N)
    if cutoff > 0
        dealias[cutoff+2:N-cutoff] .= false
    end

    # Time stepping setup
    n_steps = Int(ceil(t_final / dt))
    dt = t_final / n_steps  # adjust dt to hit t_final exactly
    save_every = max(1, n_steps ÷ n_snapshots)

    # Compute ETDRK4 coefficients
    E, E2, Q, f1, f2, f3 = compute_etdrk4_coefficients(Lop, dt)

    # Initial condition: u(x, 0) = cos(x/16) * (1 + sin(x/16))
    # This is the initial condition used by Kassam & Trefethen (2005)
    u = cos.(x ./ 16.0) .* (1.0 .+ sin.(x ./ 16.0))
    v = fft(u)

    # Nonlinear term in Fourier space: N_hat = -(ik/2) * FFT(u^2)
    # Uses the conservation form u*u_x = (u^2/2)_x
    function Nhat(v_hat)
        u_phys = real.(ifft(v_hat))
        nl = -0.5im .* k .* fft(u_phys.^2)
        nl[.!dealias] .= 0.0  # dealiasing
        return nl
    end

    # Storage for snapshots
    snapshots = [copy(u)]
    snap_times = [0.0]

    # ETDRK4 time integration loop
    # Scheme (Cox & Matthews 2002, Kassam & Trefethen 2005):
    #
    #   Nv = N(v_n)
    #   a  = E2 * v_n + Q * Nv
    #   Na = N(a)
    #   b  = E2 * v_n + Q * Na
    #   Nb = N(b)
    #   c  = E2 * a + Q * (2*Nb - Nv)
    #   Nc = N(c)
    #   v_{n+1} = E * v_n + f1*Nv + 2*f2*(Na + Nb) + f3*Nc

    t = 0.0
    for n in 1:n_steps
        Nv = Nhat(v)
        a = E2 .* v .+ Q .* Nv
        Na = Nhat(a)
        b = E2 .* v .+ Q .* Na
        Nb = Nhat(b)
        c = E2 .* a .+ Q .* (2.0 .* Nb .- Nv)
        Nc = Nhat(c)
        v = E .* v .+ f1 .* Nv .+ 2.0 .* f2 .* (Na .+ Nb) .+ f3 .* Nc

        t += dt

        if n % save_every == 0
            u = real.(ifft(v))
            push!(snapshots, copy(u))
            push!(snap_times, t)
        end
    end

    # Convert to matrix form: rows = time snapshots, columns = spatial points
    U_mat = reduce(vcat, [u_snap' for u_snap in snapshots])
    return x, snap_times, U_mat
end


# -----------------------------------------------------------------------------
# Main routine
# -----------------------------------------------------------------------------
function main()
    mkpath(OUTPUT_DIR)

    # =========================================================================
    # Solve the Kuramoto-Sivashinsky equation
    # =========================================================================
    N = 256                  # number of Fourier modes
    L = 32.0 * pi            # domain length
    t_final = 150.0          # final time
    dt = 0.5                 # time step (ETDRK4 allows larger steps)

    println("Kuramoto-Sivashinsky equation with ETDRK4")
    @printf("  N = %d, L = %.0f*pi, T = %.0f, dt = %.1f\n", N, L / pi, t_final, dt)
    println("  Initial condition: u(x,0) = cos(x/16) * (1 + sin(x/16))")
    println("  Computing...")

    x, t_save, U_save = solve_ks_etdrk4(
        N=N, L=L, t_final=t_final, dt=dt, n_snapshots=500)

    @printf("  %d snapshots saved, t in [%.1f, %.1f]\n",
            length(t_save), t_save[1], t_save[end])
    @printf("  max |u| = %.4f\n", maximum(abs.(U_save)))

    # =========================================================================
    # Figure: Space-time heatmap
    # =========================================================================
    fig = Figure(size = (700, 520))

    ax = Axis(fig[1, 1],
              xlabel = L"x / \pi",
              ylabel = L"t",
              title = L"Kuramoto--Sivashinsky: ETDRK4 ($L = 32\pi$, $N = 256$)")

    # Use a symmetric color range based on the 98th percentile
    # (after the initial transient has passed)
    n_transient = length(t_save) ÷ 5
    vmax = Float64(quantile(abs.(vec(U_save[n_transient:end, :])), 0.98))

    hm = heatmap!(ax, x ./ pi, t_save, U_save,
                  colormap = :RdBu, colorrange = (-vmax, vmax),
                  interpolate = true)

    Colorbar(fig[1, 2], hm, label = L"u(x, t)")

    xlims!(ax, 0, L / pi)
    ylims!(ax, 0, t_final)

    # Save
    save(joinpath(OUTPUT_DIR, "kuramoto_sivashinsky.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "kuramoto_sivashinsky.png"), fig, px_per_unit = 4)
    println("\nFigure saved to: $(joinpath(OUTPUT_DIR, "kuramoto_sivashinsky.pdf"))")
end

main()
