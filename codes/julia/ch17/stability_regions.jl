#!/usr/bin/env julia
#=
stability_regions.jl
Chapter 17: Stability of Spectral Methods

Computational Etude 17.2: Stability Regions and Spectral Eigenvalue Overlay.

Two figures:
  (1) stability_regions.pdf -- 2x3 grid showing stability regions of six
      classical time-stepping methods: Forward Euler, RK4, Leapfrog,
      Adams-Bashforth 2, Crank-Nicolson, Backward Euler.
  (2) spectra_overlay.pdf -- 1x2 figure overlaying spectral eigenvalues
      on the stability regions of RK4 and Forward Euler.

Generates Figures 17.2a,b: stability_regions.pdf, spectra_overlay.pdf

Author: Dr. Denys Dutykh
        Mathematics Department
        Khalifa University of Science and Technology
        Abu Dhabi, UAE

Part of the book "Computational Etudes: A Spectral Approach"
https://github.com/dutykh/computational-etudes
=#

using LinearAlgebra
using Plots
using LaTeXStrings

gr()

const NAVY       = RGB(20/255, 45/255, 110/255)
const SKY        = RGB(120/255, 150/255, 210/255)
const CORAL      = RGB(231/255, 76/255, 60/255)
const TEAL       = RGB(22/255, 160/255, 133/255)
const AMBER      = RGB(230/255, 126/255, 34/255)
const LIGHT_NAVY = RGB(197/255, 208/255, 232/255)

const SCRIPT_DIR = @__DIR__
const OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch17", "julia")
mkpath(OUTPUT_DIR)

# --- Stability functions ---
g_fe(z)  = 1.0 .+ z
g_rk4(z) = 1.0 .+ z .+ z.^2 ./2 .+ z.^3 ./6 .+ z.^4 ./24
g_be(z)  = 1.0 ./ (1.0 .- z)
g_cn(z)  = (1.0 .+ z./2) ./ (1.0 .- z./2)

"""Boundary of AB2 stability region via root locus on |g|=1."""
function ab2_boundary(n_pts=500)
    theta = range(0, 2pi, length=n_pts)
    g = exp.(1im .* theta)
    z = (g.^2 .- g) ./ (1.5 .* g .- 0.5)
    return real.(z), imag.(z)
end

"""Boundary of leapfrog stability region via root locus on |g|=1."""
function leapfrog_boundary(n_pts=500)
    theta = range(0, 2pi, length=n_pts)
    g = exp.(1im .* theta)
    z = (g.^2 .- 1.0) ./ (2.0 .* g)
    return real.(z), imag.(z)
end

"""Chebyshev differentiation matrix."""
function cheb(N)
    if N == 0; return [0.0;;], [1.0]; end
    x = [cos(pi * j / N) for j in 0:N]
    c = [j == 0 || j == N ? 2.0 : 1.0 for j in 0:N] .* [(-1.0)^j for j in 0:N]
    X = repeat(x, 1, N+1)
    dX = X - X' + I(N+1)
    D = (c * (1.0 ./ c)') ./ dX
    D .-= Diagonal(vec(sum(D, dims=2)))
    return D, x
end

"""Plot stability region for one-step methods using contour fill."""
function plot_stability_contour!(p, gfun, xlims, ylims, title_str; n_grid=400)
    xr = range(xlims[1], xlims[2], length=n_grid)
    yr = range(ylims[1], ylims[2], length=n_grid)
    Z = [complex(xi, yi) for yi in yr, xi in xr]
    G = abs.(gfun(Z))

    # Shade stable region
    contourf!(p, xr, yr, G, levels=[0, 1], c=cgrad([LIGHT_NAVY, LIGHT_NAVY]),
              alpha=0.5, colorbar=false)
    contour!(p, xr, yr, G, levels=[1], c=NAVY, lw=1.5, colorbar=false)
    hline!(p, [0], color=:gray, lw=0.5, label="")
    vline!(p, [0], color=:gray, lw=0.5, label="")
    plot!(p, xlabel=L"\mathrm{Re}(z)", ylabel=L"\mathrm{Im}(z)",
          title=title_str, aspect_ratio=:equal,
          xlims=xlims, ylims=ylims,
          grid=true, gridalpha=0.2)
end

"""Plot stability region for multistep methods using fill."""
function plot_stability_boundary!(p, bx, by, xlims, ylims, title_str)
    plot!(p, bx, by, seriestype=:shape, fillcolor=LIGHT_NAVY,
          fillalpha=0.5, linecolor=NAVY, lw=1.5, label="")
    hline!(p, [0], color=:gray, lw=0.5, label="")
    vline!(p, [0], color=:gray, lw=0.5, label="")
    plot!(p, xlabel=L"\mathrm{Re}(z)", ylabel=L"\mathrm{Im}(z)",
          title=title_str, aspect_ratio=:equal,
          xlims=xlims, ylims=ylims,
          grid=true, gridalpha=0.2)
end

function figure_stability_regions()
    # Create 6 subplots
    p1 = plot()
    plot_stability_contour!(p1, g_fe, (-3, 1.5), (-2, 2), "(a) Forward Euler")

    p2 = plot()
    plot_stability_contour!(p2, g_rk4, (-5, 2), (-4, 4), "(b) RK4")

    p3 = plot()
    bx_lf, by_lf = leapfrog_boundary(1000)
    plot_stability_boundary!(p3, bx_lf, by_lf, (-1.5, 1.5), (-1.5, 1.5), "(c) Leapfrog")

    p4 = plot()
    bx_ab, by_ab = ab2_boundary(1000)
    plot_stability_boundary!(p4, bx_ab, by_ab, (-2, 1.5), (-1.5, 1.5), "(d) Adams-Bashforth 2")

    p5 = plot()
    plot_stability_contour!(p5, g_cn, (-5, 1), (-4, 4), "(e) Crank-Nicolson")

    p6 = plot()
    plot_stability_contour!(p6, g_be, (-1.5, 5), (-4, 4), "(f) Backward Euler")

    fig = plot(p1, p2, p3, p4, p5, p6, layout=(2, 3), size=(1200, 800))

    savefig(fig, joinpath(OUTPUT_DIR, "stability_regions.pdf"))
    savefig(fig, joinpath(OUTPUT_DIR, "stability_regions.png"))
    println("Figure saved to $(joinpath(OUTPUT_DIR, "stability_regions.pdf"))")
    return fig
end

function figure_spectra_overlay()
    n_grid = 400

    # --- Left panel: RK4 + Fourier imaginary eigenvalues ---
    xr1 = range(-5, 2, length=n_grid)
    yr1 = range(-4, 4, length=n_grid)
    Z1 = [complex(xi, yi) for yi in yr1, xi in xr1]
    G1 = abs.(g_rk4(Z1))

    p1 = plot()
    contourf!(p1, xr1, yr1, G1, levels=[0, 1], c=cgrad([LIGHT_NAVY, LIGHT_NAVY]),
              alpha=0.4, colorbar=false)
    contour!(p1, xr1, yr1, G1, levels=[1], c=NAVY, lw=1.5, colorbar=false)

    # Fourier eigenvalues for advection: lambda_k = i*k
    N_f = 32
    kvals = collect(-N_f÷2+1 : N_f÷2)
    lam_f = 1im .* kvals

    # Stable dt
    dt_s = 2.8 / (N_f / 2)
    eigs_s = dt_s .* lam_f
    scatter!(p1, real.(eigs_s), imag.(eigs_s), color=TEAL, ms=3,
             markershape=:circle, label=latexstring("\\Delta t = $(round(dt_s, digits=4)) (stable)"))

    # Unstable dt
    dt_u = 3.5 / (N_f / 2)
    eigs_u = dt_u .* lam_f
    scatter!(p1, real.(eigs_u), imag.(eigs_u), color=CORAL, ms=3,
             markershape=:square, label=latexstring("\\Delta t = $(round(dt_u, digits=4)) (unstable)"))

    hline!(p1, [0], color=:gray, lw=0.5, label="")
    vline!(p1, [0], color=:gray, lw=0.5, label="")
    plot!(p1, xlabel=L"\mathrm{Re}(\Delta t \lambda)", ylabel=L"\mathrm{Im}(\Delta t \lambda)",
          title=L"(a) RK4 + Fourier advection ($N=32$)",
          aspect_ratio=:equal, grid=true, gridalpha=0.2, legendfontsize=7)

    # --- Right panel: Forward Euler + Chebyshev D^2 eigenvalues ---
    xr2 = range(-3, 1.5, length=n_grid)
    yr2 = range(-2, 2, length=n_grid)
    Z2 = [complex(xi, yi) for yi in yr2, xi in xr2]
    G2 = abs.(g_fe(Z2))

    p2 = plot()
    contourf!(p2, xr2, yr2, G2, levels=[0, 1], c=cgrad([LIGHT_NAVY, LIGHT_NAVY]),
              alpha=0.4, colorbar=false)
    contour!(p2, xr2, yr2, G2, levels=[1], c=NAVY, lw=1.5, colorbar=false)

    # Chebyshev D^2 eigenvalues (Dirichlet)
    N_c = 16
    D, xc = cheb(N_c)
    D2 = D * D
    D2_int = D2[2:N_c, 2:N_c]
    eigs_c = eigvals(D2_int)

    # Scale: choose dt so a few eigenvalues fit in the disk
    sorted_eigs = sort(abs.(real.(eigs_c)))
    dt_c = 0.8 / sorted_eigs[min(5, length(sorted_eigs))]
    scaled_eigs = dt_c .* eigs_c

    # Separate physical from outlier modes
    phys_mask = abs.(real.(eigs_c)) .< 200
    out_mask  = .!phys_mask

    scatter!(p2, real.(scaled_eigs[phys_mask]), imag.(scaled_eigs[phys_mask]),
             color=TEAL, ms=4, markershape=:circle, label="physical modes")
    if any(out_mask)
        scatter!(p2, real.(scaled_eigs[out_mask]), imag.(scaled_eigs[out_mask]),
                 color=CORAL, ms=4, markershape=:square, label="outlier modes")
    end

    hline!(p2, [0], color=:gray, lw=0.5, label="")
    vline!(p2, [0], color=:gray, lw=0.5, label="")
    plot!(p2, xlabel=L"\mathrm{Re}(\Delta t \lambda)", ylabel=L"\mathrm{Im}(\Delta t \lambda)",
          title=L"(b) Fwd Euler + Chebyshev $D^2$ ($N=16$)",
          aspect_ratio=:equal, grid=true, gridalpha=0.2, legendfontsize=7)

    fig = plot(p1, p2, layout=(1, 2), size=(1100, 500))

    savefig(fig, joinpath(OUTPUT_DIR, "spectra_overlay.pdf"))
    savefig(fig, joinpath(OUTPUT_DIR, "spectra_overlay.png"))
    println("Figure saved to $(joinpath(OUTPUT_DIR, "spectra_overlay.pdf"))")
    return fig
end

function main()
    figure_stability_regions()
    figure_spectra_overlay()
end

main()
