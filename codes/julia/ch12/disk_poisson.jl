#=
disk_poisson.jl - Poisson equation on the unit disk

Computational Etude 12.4: Solves -Du = f on the disk with u = 0 at r = 1,
using an off-centre Gaussian source. Also compares radially symmetric vs
asymmetric forcing to verify symmetry preservation.

Chapter 12: Spectral Methods in Polar Coordinates
Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi, UAE)
Last modified: February 2026
=#

include(joinpath(@__DIR__, "laplacian_polar.jl"))

using LinearAlgebra
using Plots
using LaTeXStrings
using Printf

# Output directory
fig_dir = joinpath(@__DIR__, "..", "..", "..", "textbook", "figures", "ch12", "julia")
mkpath(fig_dir)

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

# ---------------------------------------------------------------------------
# Build the Laplacian
# ---------------------------------------------------------------------------

Nr = 31   # radial points (odd)
M  = 40   # angular points (even, more than eigenmodes for smoother plots)
N2 = (Nr - 1) ÷ 2

L, r_pos, theta = laplacian_polar(Nr, M)

# Grid in (r, theta) space
RR = [r_pos[i] for i in 1:N2, j in 1:M]
TT = [theta[j] for i in 1:N2, j in 1:M]
XX = RR .* cos.(TT)
YY = RR .* sin.(TT)

# ---------------------------------------------------------------------------
# Poisson solve with off-centre Gaussian forcing
# ---------------------------------------------------------------------------

# Source: localized load at (x0, y0) = (0.4, 0.2)
x0, y0 = 0.4, 0.2
beta = 40.0
f_vals = exp.(-beta .* ((XX .- x0).^2 .+ (YY .- y0).^2))

# Flatten and solve: L u = Delta u, and we want -Delta u = f, so L u = -f
f_vec = -vec(f_vals)
u_vec = L \ f_vec
U = reshape(u_vec, N2, M)

# Extend for plotting: add boundary (r=1, u=0) and wrap in theta
U_ext = zeros(N2 + 1, M + 1)
U_ext[2:end, 1:M] = U
U_ext[:, M+1] = U_ext[:, 1]

r_ext     = vcat(1.0, r_pos)
theta_ext = vcat(theta, 2π)

RR_ext = [r_ext[i] for i in 1:length(r_ext), j in 1:length(theta_ext)]
TT_ext = [theta_ext[j] for i in 1:length(r_ext), j in 1:length(theta_ext)]
XX_ext = RR_ext .* cos.(TT_ext)
YY_ext = RR_ext .* sin.(TT_ext)

# ---------------------------------------------------------------------------
# Figure 12.6: 3D surface plot of the Poisson solution
# ---------------------------------------------------------------------------

fig_surf = surface(XX_ext, YY_ext, U_ext,
                   color=:viridis, colorbar=true,
                   xlabel=L"x", ylabel=L"y", zlabel=L"u",
                   camera=(305, 25), size=(550, 450),
                   title="Poisson equation on the disk: off-centre Gaussian load",
                   titlefontsize=10)

# Draw disk boundary circle at z=0
circ_t = range(0, 2π, length=200)
plot!(fig_surf, cos.(circ_t), sin.(circ_t), zeros(200),
      color=:black, linewidth=1.0, label="")

savefig(fig_surf, joinpath(fig_dir, "poisson_disk_surface.pdf"))

# ---------------------------------------------------------------------------
# Figure 12.7: Filled contour plot on the disk
# ---------------------------------------------------------------------------

xg, yg, Zg = to_regular_grid(XX_ext, YY_ext, U_ext)
fig_cont = contourf(xg, yg, Zg, levels=20,
                    color=:viridis, colorbar=true,
                    aspect_ratio=:equal, size=(450, 400),
                    xlims=(-1.15, 1.15), ylims=(-1.15, 1.15),
                    xlabel=L"x", ylabel=L"y",
                    title="Poisson solution: filled contour",
                    titlefontsize=10)
contour!(fig_cont, xg, yg, Zg, levels=20,
         color=:black, linewidth=0.3, alpha=0.5, colorbar=false)

# Disk boundary
plot!(fig_cont, cos.(circ_t), sin.(circ_t), color=:black, linewidth=1.2, label="")

# Mark source location
scatter!(fig_cont, [x0], [y0], marker=:star5, markersize=8,
         color=:red, markerstrokecolor=:black, markerstrokewidth=0.5, label="")
annotate!(fig_cont, x0 + 0.08, y0 + 0.08, text("source", 9, :red))

savefig(fig_cont, joinpath(fig_dir, "poisson_disk_contour.pdf"))

# ---------------------------------------------------------------------------
# Figure 12.8: Radial symmetry test
# ---------------------------------------------------------------------------

# Case 1: radially symmetric Gaussian centred at origin
f_sym = exp.(-beta .* (XX.^2 .+ YY.^2))
f_vec_sym = -vec(f_sym)
u_vec_sym = L \ f_vec_sym
U_sym = reshape(u_vec_sym, N2, M)

# Case 2: off-centre Gaussian (already computed as U)

# Left: radially symmetric case -- profiles should collapse
n_angles = 8
angle_indices = round.(Int, range(1, M, length=n_angles))

p1 = plot(xlabel=L"r", ylabel=L"u(r, \theta_k)",
          title="Centred Gaussian: radial profiles\n(all angles coincide)",
          titlefontsize=10, xflip=true)
for k in angle_indices
    plot!(p1, r_pos, U_sym[:, k], color=1, linewidth=0.8, alpha=0.7, label="")
end

# Right: off-centre case -- profiles differ
p2 = plot(xlabel=L"r", ylabel=L"u(r, \theta_k)",
          title="Off-centre Gaussian: radial profiles\n(angles differ)",
          titlefontsize=10, xflip=true)
for (cnt, k) in enumerate(angle_indices)
    lbl = cnt <= 3 ? @sprintf("\$\\theta = %.2f\$", theta[k]) : ""
    plot!(p2, r_pos, U[:, k], linewidth=0.8, alpha=0.7, label=lbl)
end

fig_sym = plot(p1, p2, layout=(1, 2), size=(800, 350))
savefig(fig_sym, joinpath(fig_dir, "radial_symmetry_test.pdf"))

println("Figures saved to: ", fig_dir)
println("  poisson_disk_surface.pdf")
println("  poisson_disk_contour.pdf")
println("  radial_symmetry_test.pdf")
@printf("\nMax deflection: %.6f at off-centre load\n", maximum(U))
@printf("Max deflection: %.6f at centred load\n", maximum(U_sym))
