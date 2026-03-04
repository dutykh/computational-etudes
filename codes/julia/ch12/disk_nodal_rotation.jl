#=
disk_nodal_rotation.jl - Eigenvalue degeneracy and nodal rotations

Computational Etude 12.3: For a doubly-degenerate eigenvalue (m >= 1),
the two eigenvectors span a 2D eigenspace. Forming u_phi = cos(phi)*u1 +
sin(phi)*u2 produces continuously rotating nodal line patterns.

Chapter 12: Spectral Methods in Polar Coordinates
Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi, UAE)
Last modified: February 2026
=#

include(joinpath(@__DIR__, "laplacian_polar.jl"))

using LinearAlgebra
using Plots
using LaTeXStrings
using Printf

# ---------------------------------------------------------------------------
# Helper: interpolate polar-grid data to a regular Cartesian grid
# ---------------------------------------------------------------------------
function to_regular_grid(XX, YY, ZZ; n_fine=80)
    x_vec = collect(range(-1.0, 1.0, length=n_fine))
    y_vec = collect(range(-1.0, 1.0, length=n_fine))
    Z_grid = fill(NaN, n_fine, n_fine)
    xx_flat = vec(Float64.(XX))
    yy_flat = vec(Float64.(YY))
    zz_flat = vec(Float64.(ZZ))
    n_pts = length(xx_flat)
    for j in 1:n_fine
        yj = y_vec[j]
        for i in 1:n_fine
            xi = x_vec[i]
            xi^2 + yj^2 > 1.0 && continue
            d2_min = Inf; best_k = 1
            for k in 1:n_pts
                d2 = (xx_flat[k] - xi)^2 + (yy_flat[k] - yj)^2
                if d2 < d2_min; d2_min = d2; best_k = k; end
            end
            Z_grid[i, j] = zz_flat[best_k]
        end
    end
    x_vec, y_vec, Z_grid
end

# Output directory
fig_dir = joinpath(@__DIR__, "..", "..", "..", "textbook", "figures", "ch12", "julia")
mkpath(fig_dir)

# ---------------------------------------------------------------------------
# Build Laplacian and compute the first degenerate pair
# ---------------------------------------------------------------------------

Nr = 25
M  = 20
N2 = (Nr - 1) ÷ 2

L, r_pos, theta = laplacian_polar(Nr, M)

eig_result = eigen(-L)
idx = sortperm(real.(eig_result.values))
eigenvalues  = real.(eig_result.values[idx])
eigenvectors = real.(eig_result.vectors[:, idx])

# The first degenerate pair is modes 2 and 3 (eigenvalue j_{1,1}^2 ~ 14.682)
u1 = eigenvectors[:, 2]
u2 = eigenvectors[:, 3]

# Build the coarse grid for nodal line plotting
# Extend to include boundary at r=1 and wrap in theta
U1 = reshape(u1, N2, M)
U2 = reshape(u2, N2, M)

# Add boundary (r=1: u=0) row at top
U1_ext = zeros(N2 + 1, M + 1)
U2_ext = zeros(N2 + 1, M + 1)
U1_ext[2:end, 1:M] = U1
U2_ext[2:end, 1:M] = U2
U1_ext[:, M+1] = U1_ext[:, 1]
U2_ext[:, M+1] = U2_ext[:, 1]

r_ext     = vcat(1.0, r_pos)
theta_ext = vcat(theta, 2π)

# Build Cartesian coordinates
XX = [r_ext[i] * cos(theta_ext[j]) for i in 1:length(r_ext), j in 1:length(theta_ext)]
YY = [r_ext[i] * sin(theta_ext[j]) for i in 1:length(r_ext), j in 1:length(theta_ext)]

# ---------------------------------------------------------------------------
# Figure 12.5: Nodal line rotation for 6 values of phi
# ---------------------------------------------------------------------------

phi_values = range(0, 5π / 6, length=6)  # 0, pi/6, pi/3, pi/2, 2pi/3, 5pi/6

# Disk boundary circle points
circ_t = range(0, 2π, length=200)
circ_x = cos.(circ_t)
circ_y = sin.(circ_t)

fig = plot(layout=(2, 3), size=(800, 560))

for (idx_p, phi) in enumerate(phi_values)
    # Linear combination
    U_phi = cos(phi) .* U1_ext .+ sin(phi) .* U2_ext

    sp = idx_p

    # Draw disk boundary
    plot!(fig, circ_x, circ_y, subplot=sp, color=:black, linewidth=1.0, label="",
          aspect_ratio=:equal, xlims=(-1.1, 1.1), ylims=(-1.1, 1.1),
          axis=nothing, border=:none, framestyle=:none)

    # Filled contour for context
    xg, yg, Zg = to_regular_grid(XX, YY, U_phi)
    contourf!(fig, xg, yg, Zg,
              subplot=sp, levels=20, color=:RdBu, alpha=0.3,
              linewidth=0, colorbar=false)

    # Plot nodal lines (zero contour)
    contour!(fig, xg, yg, Zg,
             subplot=sp, levels=[0.0], color=:blue, linewidth=1.5, colorbar=false)

    ttl = idx_p > 1 ? "\$\\varphi = $(idx_p-1)\\pi/6\$" : "\$\\varphi = 0\$"
    title!(fig, subplot=sp, ttl, titlefontsize=9)
end

lam_str = @sprintf("%.4f", eigenvalues[2])
plot!(fig, plot_title="Nodal line rotation in a degenerate eigenspace (lambda ~ $lam_str)",
      plot_titlefontsize=10)
savefig(fig, joinpath(fig_dir, "nodal_rotation.pdf"))

println("Figure saved to: ", joinpath(fig_dir, "nodal_rotation.pdf"))
@printf("Degenerate eigenvalue: lambda = %.10f (modes 2 and 3)\n", eigenvalues[2])
