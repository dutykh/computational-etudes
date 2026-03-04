#=
disk_eigenmodes.jl - Eigenmodes of the Laplacian on the unit disk

Computational Etude 12.2: Computes eigenmodes of -Du = lu on the disk
with u = 0 on the boundary. Compares eigenvalues with Bessel function
zeros to demonstrate spectral accuracy.

Chapter 12: Spectral Methods in Polar Coordinates
Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi, UAE)
Last modified: February 2026
=#

include(joinpath(@__DIR__, "laplacian_polar.jl"))

using LinearAlgebra
using Plots
using LaTeXStrings
using SpecialFunctions
using Printf

# Output directory
fig_dir = joinpath(@__DIR__, "..", "..", "..", "textbook", "figures", "ch12", "julia")
mkpath(fig_dir)

# ---------------------------------------------------------------------------
# Bessel zero finder: compute the first n_zeros positive zeros of J_m(x)
# ---------------------------------------------------------------------------

"""
    besselj_zeros(m, n_zeros)

Find the first `n_zeros` positive zeros of the Bessel function J_m(x)
using Newton's method. Uses the derivative identity:
    J_m'(x) = (m/x)*J_m(x) - J_{m+1}(x)
"""
function besselj_zeros(m::Int, n_zeros::Int)
    zeros_found = Float64[]
    for k in 1:n_zeros
        # Initial guess: McMahon (m=0), Debye leading term for k=1 (m>=1),
        # or shift from the previous zero for k>=2.
        # The McMahon formula overestimates j_{m,1} for moderate m, causing
        # Newton's method to cross zero and converge to the trivial root at 0.
        if m == 0
            x0 = π * (k - 0.25)           # McMahon; accurate for m = 0
        elseif k == 1
            x0 = m + 1.856 * m^(1.0/3)    # Debye leading term; always underestimates
        else
            x0 = zeros_found[k-1] + π     # consecutive zeros are ~pi apart
        end
        # Newton iterations
        x = x0
        for _ in 1:50
            Jm  = besselj(m, x)
            dJm = (m / x) * Jm - besselj(m + 1, x)   # J_m'(x) = (m/x)J_m - J_{m+1}
            dx  = -Jm / dJm
            x  += dx
            abs(dx) < 1e-14 * abs(x) && break
        end
        push!(zeros_found, x)
    end
    return zeros_found
end

# ---------------------------------------------------------------------------
# Build the Laplacian and compute eigenmodes
# ---------------------------------------------------------------------------

Nr = 25   # radial points (odd)
M  = 20   # angular points (even)
N2 = (Nr - 1) ÷ 2

L, r_pos, theta = laplacian_polar(Nr, M)

# Solve the eigenvalue problem -L v = lambda v
eigenvalues_raw = eigvals(-L)
eigenvalues_all = real.(eigenvalues_raw)

# Sort by eigenvalue (get sorted eigenvalues and eigenvectors together)
eig_result = eigen(-L)
idx = sortperm(real.(eig_result.values))
eigenvalues  = real.(eig_result.values[idx])
eigenvectors = real.(eig_result.vectors[:, idx])

# ---------------------------------------------------------------------------
# Figure 12.3: Six representative eigenmodes as surface plots
# ---------------------------------------------------------------------------

# Mode indices (1-based in Julia, matching 0-based [0,1,3,5,6,8] from Python)
mode_indices = [1, 2, 4, 6, 7, 9]

fig_modes = plot(layout=(2, 3), size=(900, 600))

for (panel, mi) in enumerate(mode_indices)
    # Reshape eigenvector: the solution vector is ordered as
    # [u(r_1,th_1), u(r_1,th_2), ..., u(r_1,th_M), u(r_2,th_1), ...]
    v = real.(eigenvectors[:, mi])
    U = reshape(v, N2, M)

    # Normalise
    U = U ./ maximum(abs.(U))

    # Extend to include boundary (r=1, u=0) and wrap periodically in theta
    U_ext = zeros(N2 + 1, M + 1)
    U_ext[2:end, 1:M] = U
    U_ext[:, M+1] = U_ext[:, 1]   # periodic wrap

    r_ext     = vcat(1.0, r_pos)
    theta_ext = vcat(theta, 2π)

    # Build Cartesian coordinates
    XX = [r_ext[i] * cos(theta_ext[j]) for i in 1:length(r_ext), j in 1:length(theta_ext)]
    YY = [r_ext[i] * sin(theta_ext[j]) for i in 1:length(r_ext), j in 1:length(theta_ext)]

    ttl = @sprintf("Mode %d\nλ = %.10f", mi, eigenvalues[mi])

    surface!(fig_modes, XX, YY, U_ext, subplot=panel,
             color=:RdBu, colorbar=false, clims=(-1.1, 1.1),
             zlims=(-1.1, 1.1), xlims=(-1.1, 1.1), ylims=(-1.1, 1.1),
             camera=(300, 25), title=ttl, titlefontsize=8,
             axis=nothing, border=:none)
end

plot!(fig_modes, plot_title="Eigenmodes of the Laplacian on the unit disk",
      plot_titlefontsize=11)
savefig(fig_modes, joinpath(fig_dir, "disk_eigenmodes.pdf"))

# ---------------------------------------------------------------------------
# Figure 12.4: Eigenvalue convergence vs Bessel zeros
# ---------------------------------------------------------------------------

n_show = 25
lam_computed = eigenvalues[1:n_show]

# Collect exact Bessel eigenvalues and sort
bessel_lam = Float64[]
for m in 0:9
    z = besselj_zeros(m, 5)
    mult = m == 0 ? 1 : 2
    for zval in z
        for _ in 1:mult
            push!(bessel_lam, zval^2)
        end
    end
end
sort!(bessel_lam)
bessel_lam = bessel_lam[1:n_show]

# --- Left: computed vs exact eigenvalues ---
p1 = scatter(1:n_show, lam_computed, marker=:circle, markersize=5,
             markerstrokecolor=1, markerfillcolor=:white, label="Computed")
scatter!(p1, 1:n_show, bessel_lam, marker=:xcross, markersize=4,
         color=2, label=L"Bessel zeros $j_{m,n}^2$")
xlabel!(p1, "Mode index")
ylabel!(p1, L"Eigenvalue $\lambda$")
title!(p1, "Eigenvalue comparison", titlefontsize=10)

# --- Right: convergence of first eigenvalue with Nr ---
Nr_values = [9, 13, 17, 21, 25, 29, 33]
M_fixed   = 20
lam_exact_1 = besselj_zeros(0, 1)[1]^2   # j_{0,1}^2

errors = Float64[]
for N in Nr_values
    L_test, _, _ = laplacian_polar(N, M_fixed)
    eigs = sort(real.(eigvals(-L_test)))
    push!(errors, abs(eigs[1] - lam_exact_1))
end

p2 = plot(Nr_values, errors, yscale=:log10, marker=:circle, markersize=4,
          color=1, linewidth=1, label="")
xlabel!(p2, L"N_r" * " (radial points)")
ylabel!(p2, L"|\lambda_1 - j_{0,1}^2|")
title!(p2, "Exponential convergence", titlefontsize=10)

fig_conv = plot(p1, p2, layout=(1, 2), size=(800, 350))
savefig(fig_conv, joinpath(fig_dir, "eigenvalue_convergence.pdf"))

# Print eigenvalue comparison table
println(lpad("Mode", 5), "  ", lpad("Computed", 16), "  ",
        lpad("Exact (Bessel)", 16), "  ", lpad("Error", 12))
println("-"^58)
for i in 1:min(15, n_show)
    err = abs(lam_computed[i] - bessel_lam[i])
    @printf("  %3d  %16.10f  %16.10f  %12.2e\n", i, lam_computed[i], bessel_lam[i], err)
end

println("\nFigures saved to: ", fig_dir)
println("  disk_eigenmodes.pdf")
println("  eigenvalue_convergence.pdf")
