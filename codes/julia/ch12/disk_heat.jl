#=
disk_heat.jl - Heat equation on the unit disk

Computational Etude 12.5: Solves the heat equation u_t = alpha*Delta*u + S(r,theta)
on the unit disk with Crank-Nicolson time stepping. Models heat dissipation
on a semiconductor wafer with a localised laser source.

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
# Parameters
# ---------------------------------------------------------------------------

Nr      = 25       # radial points (odd)
M       = 20       # angular points (even)
alpha   = 0.1      # thermal diffusivity
dt      = 0.01     # time step
t_final = 20.0
n_steps = round(Int, t_final / dt)

N2   = (Nr - 1) ÷ 2
ndof = N2 * M

# Build the Laplacian
L, r_pos, theta = laplacian_polar(Nr, M)

# Grid
RR = [r_pos[i] for i in 1:N2, j in 1:M]
TT = [theta[j] for i in 1:N2, j in 1:M]
XX = RR .* cos.(TT)
YY = RR .* sin.(TT)

# Source: localised laser spot at (0.3, 0.0)
x_src, y_src = 0.3, 0.0
S_vals = 2.0 .* exp.(-30.0 .* ((XX .- x_src).^2 .+ (YY .- y_src).^2))
S_vec  = vec(S_vals)

# ---------------------------------------------------------------------------
# Crank-Nicolson time stepping
# ---------------------------------------------------------------------------

# (I - alpha*dt/2 * L) u^{n+1} = (I + alpha*dt/2 * L) u^n + dt * S
I_mat = Matrix{Float64}(I, ndof, ndof)
A = I_mat - alpha * dt / 2 .* L     # LHS matrix
B = I_mat + alpha * dt / 2 .* L     # RHS matrix part

# LU factorise once (reuse at every step)
lu_A = lu(A)

# Initial condition: u = 0 (cold start)
u = zeros(ndof)

# Store snapshots at selected times
snapshot_times = [0.0, 0.5, 2.0, 5.0, 10.0, 20.0]
snapshots = Dict{Float64, Vector{Float64}}()
max_temp  = Float64[]
times     = Float64[]

snapshots[0.0] = copy(u)
push!(max_temp, 0.0)
push!(times, 0.0)

t = 0.0
for step in 1:n_steps
    rhs = B * u .+ dt .* S_vec
    u   = lu_A \ rhs
    t   = step * dt

    push!(max_temp, maximum(u))
    push!(times, t)

    # Check if this is a snapshot time (within tolerance)
    for ts in snapshot_times
        if abs(t - ts) < dt / 2 && !haskey(snapshots, ts)
            snapshots[ts] = copy(u)
        end
    end
end

# Ensure final snapshot is stored
if !haskey(snapshots, t_final)
    snapshots[t_final] = copy(u)
end

# ---------------------------------------------------------------------------
# Figure 12.9: Heat equation snapshots at 6 time levels
# ---------------------------------------------------------------------------

# Extend grid for plotting
r_ext     = vcat(1.0, r_pos)
theta_ext = vcat(theta, 2π)

RR_ext = [r_ext[i] for i in 1:length(r_ext), j in 1:length(theta_ext)]
TT_ext = [theta_ext[j] for i in 1:length(r_ext), j in 1:length(theta_ext)]
XX_ext = RR_ext .* cos.(TT_ext)
YY_ext = RR_ext .* sin.(TT_ext)

# Global color scale
u_max = maximum(maximum(abs.(s)) for s in values(snapshots))
if u_max < 1e-10
    u_max = 1.0
end

# Disk boundary circle points
circ_t = range(0, 2π, length=200)
circ_x = cos.(circ_t)
circ_y = sin.(circ_t)

fig_snap = plot(layout=(2, 3), size=(880, 600))

for (panel, ts) in enumerate(snapshot_times)
    U_snap = reshape(snapshots[ts], N2, M)

    # Extend for plotting
    U_snap_ext = zeros(N2 + 1, M + 1)
    U_snap_ext[2:end, 1:M] = U_snap
    U_snap_ext[:, M+1] = U_snap_ext[:, 1]

    contourf!(fig_snap, XX_ext, YY_ext, U_snap_ext,
              subplot=panel, levels=20, color=:hot,
              clims=(0, u_max), colorbar=false,
              aspect_ratio=:equal, xlims=(-1.15, 1.15), ylims=(-1.15, 1.15),
              axis=nothing, border=:none, framestyle=:none)
    plot!(fig_snap, circ_x, circ_y, subplot=panel,
          color=:black, linewidth=1.0, label="")
    title!(fig_snap, subplot=panel, @sprintf("t = %.1f", ts), titlefontsize=10)
end

plot!(fig_snap, plot_title="Heat dissipation on a disk with localised source",
      plot_titlefontsize=10)
savefig(fig_snap, joinpath(fig_dir, "heat_disk_snapshots.pdf"))

# ---------------------------------------------------------------------------
# Figure 12.10: Temperature evolution and steady-state comparison
# ---------------------------------------------------------------------------

# Left: maximum temperature vs time
p1 = plot(times, max_temp, color=1, linewidth=1.0, label="",
          xlabel=L"Time $t$",
          ylabel=L"\max_{x,y} \, u(x,y,t)",
          title="Maximum temperature evolution", titlefontsize=10)

# Right: compare steady state with direct Poisson solve
# Steady state satisfies alpha*Delta*u + S = 0, i.e., Delta*u = -S/alpha
u_poisson = L \ (-S_vec ./ alpha)

# Plot radial profiles at theta = 0
angle_idx = 1
U_final   = reshape(snapshots[t_final], N2, M)
U_poisson = reshape(u_poisson, N2, M)

p2 = plot(r_pos, U_final[:, angle_idx], color=1, linewidth=1.5,
          label="CN steady state", xflip=true)
plot!(p2, r_pos, U_poisson[:, angle_idx], color=2, linestyle=:dash,
      linewidth=1.5, label="Direct Poisson solve")
xlabel!(p2, L"r")
ylabel!(p2, L"u(r, 0)")
title!(p2, "Steady state vs Poisson solution (theta = 0)", titlefontsize=10)

fig_energy = plot(p1, p2, layout=(1, 2), size=(800, 350))
savefig(fig_energy, joinpath(fig_dir, "heat_disk_energy.pdf"))

# Print summary
err_ss = maximum(abs.(snapshots[t_final] .- u_poisson))
@printf("Steady-state vs Poisson max error: %.2e\n", err_ss)
@printf("Max temperature at t=%.1f: %.6f\n", t_final, max_temp[end])
println("\nFigures saved to: ", fig_dir)
println("  heat_disk_snapshots.pdf")
println("  heat_disk_energy.pdf")
