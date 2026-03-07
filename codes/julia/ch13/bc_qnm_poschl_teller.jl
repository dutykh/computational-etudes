# bc_qnm_poschl_teller.jl
#
# Chapter 13: Boundary Conditions in Spectral Methods
#
# Quasinormal modes (QNMs) of the Poeschl-Teller potential: the crown jewel.
#
# The QNM problem is a quadratic eigenvalue problem (QEP):
#     psi'' + (omega^2 - V(x)) psi = 0,   x in [-L, L]
#     V(x) = V0 * sech^2(x),   V0 = 2
#
# with outgoing wave boundary conditions:
#     psi'(+L) + i*omega * psi(+L) = 0   (outgoing at +L)
#     psi'(-L) - i*omega * psi(-L) = 0   (outgoing at -L)
#
# Substitution lambda = -i*omega transforms to:
#     psi'' - V(x)*psi - lambda^2 * psi = 0
#     psi'(+L) + lambda * psi(+L) = 0
#     psi'(-L) - lambda * psi(-L) = 0
#
# QEP: (M0 + lambda*M1 + lambda^2*M2) psi = 0
#   Interior rows:  M0 = D^2 - diag(V),  M1 = 0,  M2 = -I
#   Row 1 (x=+L):   M0[1,:] = D[1,:],  M1[1,1] = 1,  M2[1,:] = 0
#   Row N+1 (x=-L): M0[N+1,:] = D[N+1,:],  M1[N+1,N+1] = -1,  M2[N+1,:] = 0
#
# Companion linearization:
#   A = [M0  M1;  Z  Id],   B = [Z  -M2;  Id  Z]
#   A*w = lambda * B*w
#
# Exact fundamental QNM: omega_0 = sqrt(7)/2 - i/2
#
# Generates:
#   Figure 13.7: QNM spectrum in the complex omega-plane
#   Figure 13.8: First two QNM eigenfunctions (2 panels)
#   Figure 13.9: Convergence vs N and vs L (2 panels)
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
using Printf

# Import Chebyshev matrix function from Chapter 7
include(joinpath(@__DIR__, "..", "ch07", "cheb_matrix.jl"))

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
const OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch13", "julia")

# -----------------------------------------------------------------------------
# Exact QNM frequency
# -----------------------------------------------------------------------------
const V0 = 2.0
const OMEGA_EXACT = sqrt(7.0) / 2.0 - 0.5im  # fundamental mode

# -----------------------------------------------------------------------------
# QNM solver
# -----------------------------------------------------------------------------
"""
    compute_qnm(N, L; V0=2.0)

Compute quasinormal modes of the Poeschl-Teller potential via companion
linearization of the quadratic eigenvalue problem.

Returns `(omega_all, eigvecs_all)` where omega values are in the complex plane.
"""
function compute_qnm(N, L; V0=2.0)
    D, xi = cheb_matrix(N)
    n = N + 1

    # Scale to [-L, L]: D_phys = D / L (since x = L * xi)
    D_phys = D ./ L
    D2_phys = D_phys * D_phys
    x_phys = L .* xi  # x_phys[1]=+L, x_phys[N+1]=-L

    # Potential
    V = V0 .* sech.(x_phys).^2

    # Build QEP matrices: (M0 + lambda*M1 + lambda^2*M2) psi = 0

    # Interior: M0 = D^2 - diag(V),  M1 = 0,  M2 = -I
    M0 = copy(D2_phys)
    for i in 1:n
        M0[i, i] -= V[i]
    end
    M1 = zeros(n, n)
    M2 = -Matrix(1.0I, n, n)

    # Row 1 (x=+L): psi'(+L) + lambda * psi(+L) = 0
    M0[1, :] = D_phys[1, :]
    M1[1, :] .= 0.0
    M1[1, 1] = 1.0
    M2[1, :] .= 0.0

    # Row N+1 (x=-L): psi'(-L) - lambda * psi(-L) = 0
    M0[n, :] = D_phys[n, :]
    M1[n, :] .= 0.0
    M1[n, n] = -1.0
    M2[n, :] .= 0.0

    # Companion linearization: A*w = lambda * B*w
    # A = [M0  M1;  Z  I],   B = [Z  -M2;  I  Z]
    Z = zeros(n, n)
    Id = Matrix(1.0I, n, n)

    A = [M0 M1; Z Id]
    B = [Z -M2; Id Z]

    # Solve generalized eigenvalue problem
    lambdas = eigvals(A, B)

    # Convert to omega: omega = i * lambda (since lambda = -i*omega)
    omega_all = 1.0im .* lambdas

    # Also compute eigenvectors for mode extraction
    F = eigen(A, B)
    eigvecs_all = F.vectors

    return omega_all, eigvecs_all, x_phys
end

"""
    find_physical_modes(omega_all; omega_exact=OMEGA_EXACT, tol=2.0)

Filter for physical QNMs: modes with Im(omega) < 0 (damped)
that are close to the real axis and have moderate imaginary part.
"""
function find_physical_modes(omega_all; max_imag=-0.01, min_real=0.0, max_abs=20.0)
    mask = (imag.(omega_all) .< max_imag) .&
           (real.(omega_all) .> min_real) .&
           (abs.(omega_all) .< max_abs)
    return findall(mask)
end

# -----------------------------------------------------------------------------
# Main visualization
# -----------------------------------------------------------------------------
function main()
    mkpath(OUTPUT_DIR)

    N = 80
    L = 8.0

    println("Computing QNMs with N=$N, L=$L...")
    omega_all, eigvecs_all, x_phys = compute_qnm(N, L, V0=V0)

    # =========================================================================
    # Figure 1: QNM spectrum in complex omega-plane
    # =========================================================================
    fig1 = Figure(size = (550, 450))

    ax1 = Axis(fig1[1, 1],
               xlabel = L"\mathrm{Re}(\omega)",
               ylabel = L"\mathrm{Im}(\omega)",
               title = L"Quasinormal Mode Spectrum ($V_0 = 2$, $N = %$(N)$, $L = %$(Int(L))$)")

    # Plot all eigenvalues as light dots (all spurious except fundamental)
    scatter!(ax1, real.(omega_all), imag.(omega_all),
             color = (SKY, 0.3), markersize = 4,
             label = "Spurious")

    # Find and highlight computed fundamental mode
    phys_idx = find_physical_modes(omega_all)
    if !isempty(phys_idx)
        omega_phys = omega_all[phys_idx]
        dists = abs.(omega_phys .- OMEGA_EXACT)
        fund_local = argmin(dists)
        om_fund = omega_phys[fund_local]
        scatter!(ax1, [real(om_fund)], [imag(om_fund)],
                 color = CORAL, marker = :star5, markersize = 14,
                 strokecolor = :white, strokewidth = 0.5,
                 label = L"Fundamental $\omega_0$")
    end

    # Mark exact fundamental mode
    scatter!(ax1, [real(OMEGA_EXACT)], [imag(OMEGA_EXACT)],
             color = :black, marker = :cross, markersize = 12,
             label = @sprintf("Exact: %.4f%+.4fi", real(OMEGA_EXACT), imag(OMEGA_EXACT)))

    # Also mark the mirror mode at -Re(omega)
    scatter!(ax1, [-real(OMEGA_EXACT)], [imag(OMEGA_EXACT)],
             color = :black, marker = :cross, markersize = 12)

    xlims!(ax1, -6, 6)
    ylims!(ax1, -10, 1)
    hlines!(ax1, [0.0], color = :gray, linestyle = :dash, linewidth = 0.5)
    vlines!(ax1, [0.0], color = :gray, linestyle = :dash, linewidth = 0.5)
    axislegend(ax1, position = :rt, framevisible = true, labelsize = 9)
    hidespines!(ax1, :t, :r)

    save(joinpath(OUTPUT_DIR, "qnm_spectrum.pdf"), fig1)
    save(joinpath(OUTPUT_DIR, "qnm_spectrum.png"), fig1, px_per_unit = 4)
    println("Figure 1 saved to: $(joinpath(OUTPUT_DIR, "qnm_spectrum.pdf"))")

    # =========================================================================
    # Figure 2: First two QNM eigenfunctions
    # =========================================================================
    fig2 = Figure(size = (900, 380))

    # Find the two modes closest to exact omega_0
    dists_pos = abs.(omega_all .- OMEGA_EXACT)
    dists_neg = abs.(omega_all .+ conj(OMEGA_EXACT))  # mirror mode

    idx_mode1 = argmin(dists_pos)
    omega_mode1 = omega_all[idx_mode1]

    # Find second overtone (next mode with similar real part but more negative imag)
    # Look for modes with Im < Im(omega_0) - 0.5
    candidate_mask = (imag.(omega_all) .< imag(OMEGA_EXACT) - 0.3) .&
                     (real.(omega_all) .> 0.0) .&
                     (abs.(real.(omega_all)) .< 5.0) .&
                     (imag.(omega_all) .> -5.0)
    candidates = findall(candidate_mask)

    if !isempty(candidates)
        # Sort by imaginary part (least damped first)
        sorted_cand = sort(candidates, by = i -> -imag(omega_all[i]))
        idx_mode2 = sorted_cand[1]
    else
        # Fallback: second closest to exact
        dists_pos[idx_mode1] = Inf
        idx_mode2 = argmin(dists_pos)
    end
    omega_mode2 = omega_all[idx_mode2]

    n = N + 1

    for (panel, (idx, omega_m, title_str)) in enumerate([
        (idx_mode1, omega_mode1, "Fundamental mode"),
        (idx_mode2, omega_mode2, "First overtone")
    ])
        ax = Axis(fig2[1, panel],
                  xlabel = L"x",
                  ylabel = L"\psi(x)",
                  title = @sprintf("%s: ω = %.4f%+.4fi", title_str, real(omega_m), imag(omega_m)))

        # Extract eigenvector (first n components of the 2n companion vector)
        psi = eigvecs_all[1:n, idx]
        # Normalize
        psi ./= maximum(abs.(psi))

        lines!(ax, x_phys, real.(psi), color = NAVY, linewidth = 1.5, label = L"\mathrm{Re}(\psi)")
        lines!(ax, x_phys, imag.(psi), color = CORAL, linewidth = 1.5, linestyle = :dash, label = L"\mathrm{Im}(\psi)")
        lines!(ax, x_phys, abs.(psi), color = TEAL, linewidth = 1.0, linestyle = :dot, label = L"|\psi|")

        axislegend(ax, position = :rt, framevisible = true, labelsize = 9)
        hidespines!(ax, :t, :r)
    end

    Label(fig2[0, :], "QNM Eigenfunctions of Poeschl-Teller Potential", fontsize = 13)

    save(joinpath(OUTPUT_DIR, "qnm_eigenfunctions.pdf"), fig2)
    save(joinpath(OUTPUT_DIR, "qnm_eigenfunctions.png"), fig2, px_per_unit = 4)
    println("Figure 2 saved to: $(joinpath(OUTPUT_DIR, "qnm_eigenfunctions.pdf"))")

    # =========================================================================
    # Figure 3: Convergence (2 panels)
    # =========================================================================
    fig3 = Figure(size = (900, 380))

    # ---- Panel 1: Error vs N (fixed L=8) ----
    ax3a = Axis(fig3[1, 1],
                xlabel = L"N",
                ylabel = L"|\omega - \omega_{\mathrm{exact}}|",
                title = L"Convergence vs $N$ ($L = 8$)",
                yscale = log10)

    L_fixed = 8.0
    N_values = collect(20:5:120)
    errors_N = Float64[]

    for N_val in N_values
        omega_test, _, _ = compute_qnm(N_val, L_fixed, V0=V0)
        # Find closest to exact
        dists = abs.(omega_test .- OMEGA_EXACT)
        err = minimum(dists)
        push!(errors_N, max(err, 1e-16))
    end

    scatterlines!(ax3a, N_values, errors_N, color = CORAL, linewidth = 1.5,
                  markersize = 6, markercolor = CORAL,
                  strokecolor = :white, strokewidth = 0.5)

    hidespines!(ax3a, :t, :r)
    hlines!(ax3a, [2.2e-16], color = :gray, linestyle = :dot, linewidth = 1)

    # ---- Panel 2: Error vs L (fixed N=80) ----
    ax3b = Axis(fig3[1, 2],
                xlabel = L"L",
                ylabel = L"|\omega - \omega_{\mathrm{exact}}|",
                title = L"Convergence vs $L$ ($N = 80$)",
                yscale = log10)

    N_fixed = 80
    L_values = collect(2.0:1.0:15.0)
    errors_L = Float64[]

    for L_val in L_values
        omega_test, _, _ = compute_qnm(N_fixed, L_val, V0=V0)
        dists = abs.(omega_test .- OMEGA_EXACT)
        err = minimum(dists)
        push!(errors_L, max(err, 1e-16))
    end

    scatterlines!(ax3b, L_values, errors_L, color = NAVY, linewidth = 1.5,
                  markersize = 6, markercolor = NAVY,
                  strokecolor = :white, strokewidth = 0.5)

    hidespines!(ax3b, :t, :r)
    hlines!(ax3b, [2.2e-16], color = :gray, linestyle = :dot, linewidth = 1)

    Label(fig3[0, :], "QNM Convergence Analysis", fontsize = 13)

    save(joinpath(OUTPUT_DIR, "qnm_convergence.pdf"), fig3)
    save(joinpath(OUTPUT_DIR, "qnm_convergence.png"), fig3, px_per_unit = 4)
    println("Figure 3 saved to: $(joinpath(OUTPUT_DIR, "qnm_convergence.pdf"))")

    # =========================================================================
    # Print summary
    # =========================================================================
    println("\nQuasinormal Mode Analysis Summary")
    println("="^60)
    @printf("  Potential: V(x) = %.1f * sech^2(x)\n", V0)
    @printf("  Domain: [-%d, %d], N = %d\n", Int(L), Int(L), N)
    @printf("  Exact fundamental: omega = %.10f %+.10fi\n", real(OMEGA_EXACT), imag(OMEGA_EXACT))
    @printf("  Computed closest:  omega = %.10f %+.10fi\n", real(omega_mode1), imag(omega_mode1))
    @printf("  Error:             |delta omega| = %.2e\n", abs(omega_mode1 - OMEGA_EXACT))
    println("="^60)

    # Print convergence table
    println("\nConvergence vs N (L = $L_fixed):")
    println("-"^30)
    @printf("%6s %14s\n", "N", "Error")
    println("-"^30)
    for (N_val, err) in zip(N_values, errors_N)
        @printf("%6d %14.2e\n", N_val, err)
    end
    println("-"^30)
end

main()
