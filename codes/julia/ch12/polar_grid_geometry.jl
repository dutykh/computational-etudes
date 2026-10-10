#=
polar_grid_geometry.jl - Grid geometry and the cost of clustering

Computational Etude 12.1: Compares naive (r in [0,1]) and doubled
(r in [-1,1]) Chebyshev-Fourier grids on the unit disk. Visualises
grid clustering, area distribution, and CFL time-step scaling.

Chapter 12: Spectral Methods in Polar Coordinates
Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi, UAE)
Last modified: February 2026
=#

include(joinpath(@__DIR__, "..", "ch07", "cheb_matrix.jl"))

using Plots
using LaTeXStrings

# Output directory
fig_dir = joinpath(@__DIR__, "..", "..", "..", "textbook", "figures", "ch12", "julia")
mkpath(fig_dir)

# ---------------------------------------------------------------------------
# Figure 12.1: Naive vs doubled polar grids on the disk
# ---------------------------------------------------------------------------

Nr = 25   # must be odd for doubled grid
M  = 20   # angular points

# Chebyshev points on [-1, 1]
_, x_cheb = cheb_matrix(Nr)

# --- Naive grid: map x in [-1,1] to r in [0,1] via r = (x+1)/2 ---
r_naive = (x_cheb .+ 1) ./ 2   # includes r=0 and r=1
theta = range(0, 2π, length=M+1)[1:M]  # M points, endpoint=false

# Build Cartesian coordinates for naive grid
XX_n = [r * cos(t) for r in r_naive, t in theta]
YY_n = [r * sin(t) for r in r_naive, t in theta]

# --- Doubled grid: r in [-1,1], use only r > 0, odd Nr so no r=0 ---
N2 = (Nr - 1) ÷ 2
r_pos = x_cheb[2:N2+1]   # r > 0 interior points

# Build Cartesian coordinates for doubled grid
XX_d = [r * cos(t) for r in r_pos, t in theta]
YY_d = [r * sin(t) for r in r_pos, t in theta]

# Median radius (circle containing half the points)
r_naive_interior = r_naive[2:end-1]   # exclude endpoints
r_median_naive   = sort(r_naive_interior)[length(r_naive_interior) ÷ 2]
r_median_doubled = sort(r_pos)[length(r_pos) ÷ 2]

# Circle boundary points
circ_t = range(0, 2π, length=200)
circ_x = cos.(circ_t)
circ_y = sin.(circ_t)

# --- Left panel: naive grid ---
p1 = plot(circ_x, circ_y, color=:black, linewidth=1.2, label="",
          aspect_ratio=:equal, xlims=(-1.15, 1.15), ylims=(-1.15, 1.15),
          axis=nothing, border=:none, framestyle=:none)
scatter!(p1, vec(XX_n), vec(YY_n), markersize=1.5, color=:black,
         markerstrokewidth=0, label="")
# Median circle
med_x_n = r_median_naive .* cos.(circ_t)
med_y_n = r_median_naive .* sin.(circ_t)
plot!(p1, med_x_n, med_y_n, color=:red, linestyle=:dash, linewidth=1.0, label="")
area_frac_naive = round(r_median_naive^2 * 100, digits=0)
title!(p1, "Naive grid (r in [0, 1])\nMedian circle: $(Int(area_frac_naive))% of area",
       titlefontsize=9)

# --- Right panel: doubled grid ---
p2 = plot(circ_x, circ_y, color=:black, linewidth=1.2, label="",
          aspect_ratio=:equal, xlims=(-1.15, 1.15), ylims=(-1.15, 1.15),
          axis=nothing, border=:none, framestyle=:none)
scatter!(p2, vec(XX_d), vec(YY_d), markersize=1.5, color=:black,
         markerstrokewidth=0, label="")
# Median circle
med_x_d = r_median_doubled .* cos.(circ_t)
med_y_d = r_median_doubled .* sin.(circ_t)
plot!(p2, med_x_d, med_y_d, color=:red, linestyle=:dash, linewidth=1.0, label="")
area_frac_doubled = round(r_median_doubled^2 * 100, digits=0)
title!(p2, "Doubled grid (r in [-1, 1])\nMedian circle: $(Int(area_frac_doubled))% of area",
       titlefontsize=9)

fig1 = plot(p1, p2, layout=(1, 2), size=(800, 400))
savefig(fig1, joinpath(fig_dir, "polar_grids.pdf"))

# ---------------------------------------------------------------------------
# Figure 12.2: Area distribution per ring and CFL scaling
# ---------------------------------------------------------------------------

# --- Left: area per radial ring ---
r_sorted_naive = sort(r_naive)
areas_naive = Float64[]
for i in 1:length(r_sorted_naive)-1
    a = π * (r_sorted_naive[i+1]^2 - r_sorted_naive[i]^2)
    push!(areas_naive, a)
end
midpoints_naive = 0.5 .* (r_sorted_naive[1:end-1] .+ r_sorted_naive[2:end])

# Doubled grid: use sorted positive-r points, plus r=0 and r=1 as boundaries
r_sorted_doubled = sort(vcat(0.0, r_pos, 1.0))
areas_doubled = Float64[]
for i in 1:length(r_sorted_doubled)-1
    a = π * (r_sorted_doubled[i+1]^2 - r_sorted_doubled[i]^2)
    push!(areas_doubled, a)
end
midpoints_doubled = 0.5 .* (r_sorted_doubled[1:end-1] .+ r_sorted_doubled[2:end])

width = 0.012
p3 = bar(midpoints_naive .- width, areas_naive, bar_width=2*width, alpha=0.7,
         label="Naive", color=1)
bar!(p3, midpoints_doubled .+ width, areas_doubled, bar_width=2*width, alpha=0.7,
     label="Doubled", color=2)
xlabel!(p3, L"r")
ylabel!(p3, "Area of radial annulus")
title!(p3, "Area per ring", titlefontsize=10)

# --- Right: CFL time step vs Nr ---
Nr_values = 7:2:51   # odd values only
dt_naive_list   = Float64[]
dt_doubled_list = Float64[]

for N in Nr_values
    _, x = cheb_matrix(N)
    # Naive: r = (x+1)/2 in [0,1]
    r_n = (x .+ 1) ./ 2
    h_min_naive = minimum(abs.(diff(r_n)))
    push!(dt_naive_list, h_min_naive^2)   # CFL ~ h_min^2

    # Doubled: use positive interior points
    n2 = (N - 1) ÷ 2
    r_p = x[2:n2+1]
    h_min_doubled = minimum(abs.(diff(sort(r_p))))
    push!(dt_doubled_list, h_min_doubled^2)
end

Nr_ref = Float64.(collect(Nr_values))
p4 = plot(Nr_ref, dt_naive_list, xscale=:log10, yscale=:log10,
          marker=:circle, markersize=3, label="Naive", color=1, linewidth=1)
plot!(p4, Nr_ref, dt_doubled_list, marker=:square, markersize=3,
      label="Doubled", color=2, linewidth=1)
# Reference slopes
plot!(p4, Nr_ref, 5e2 .* Nr_ref .^ (-4), linestyle=:dash, color=:black,
      linewidth=0.8, alpha=0.5, label=L"O(N^{-4})")
plot!(p4, Nr_ref, 2e1 .* Nr_ref .^ (-2), linestyle=:dot, color=:black,
      linewidth=0.8, alpha=0.5, label=L"O(N^{-2})")
xlabel!(p4, L"N_r")
ylabel!(p4, L"\Delta t_{\max} \propto h_{\min}^2")
title!(p4, "CFL time-step scaling", titlefontsize=10)

fig2 = plot(p3, p4, layout=(1, 2), size=(800, 350))
savefig(fig2, joinpath(fig_dir, "polar_area_cfl.pdf"))

println("Figures saved to: ", fig_dir)
println("  polar_grids.pdf")
println("  polar_area_cfl.pdf")
