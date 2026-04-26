# weideman_cloot_map.jl
#
# Chapter 20: illustrative figure for the Weideman--Cloot (1990) sinh map.
#
# Three panels:
#   (a) the map y = sinh(ell t) for ell in {0.5, 1, 1.5, 2}
#   (b) sinh-mapped grids at fixed N = 64
#   (c) sech(y) in y-space and t-space (rescaled) at ell = 1
#
# Author: Dr. Denys Dutykh (Khalifa University of Science and Technology)
# Part of "Computational Etudes: A Spectral Approach"

using CairoMakie
using Colors
using FFTW
using Printf

set_theme!(Theme(
    fontsize = 11,
    fonts = (regular = "CMU Serif", bold = "CMU Serif Bold",
             italic  = "CMU Serif Italic"),
    Axis = (xlabelsize = 11, ylabelsize = 11, titlesize = 11,
            xticklabelsize = 10, yticklabelsize = 10,
            spinewidth = 0.8, xtickwidth = 0.8, ytickwidth = 0.8),
    Legend = (labelsize = 9,),
))

const NAVY   = colorant"#142D6E"
const CORAL  = colorant"#E74C3C"
const TEAL   = colorant"#16A085"
const ORANGE = colorant"#E67E22"

const SCRIPT_DIR = @__DIR__
const OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..",
                            "textbook", "figures", "ch20", "julia")
mkpath(OUTPUT_DIR)

function make_figure(; ell_values = (0.5, 1.0, 1.5, 2.0),
                     N_grid::Int = 64, ell_demo::Real = 1.0)
    fig = Figure(size = (1100, 800))
    palette_a = [NAVY, CORAL, TEAL, ORANGE]

    ax1 = Axis(fig[1, 1];
        xlabel = "computational coordinate t ∈ [-π, π]",
        ylabel = "physical coordinate y",
        title  = "(a) the map y = sinh(ℓ t)",
        limits = ((-pi - 0.1, pi + 0.1), (-30, 30)))
    t_dense = collect(range(-pi, pi, length = 401))
    for (ell, color) in zip(ell_values, palette_a)
        y_dense = sinh.(ell .* t_dense)
        y_max = sinh(ell * pi)
        lines!(ax1, t_dense, y_dense; color = color, linewidth = 1.4,
               label = "ℓ = $ell, y_max ≈ $(round(y_max, digits=1))")
    end
    hlines!(ax1, [0.0]; color = :black, linewidth = 0.4, alpha = 0.5)
    vlines!(ax1, [-pi, pi]; color = (:gray, 0.6),
            linestyle = :dot, linewidth = 0.6)
    axislegend(ax1, position = :lt, labelsize = 8, framevisible = false)

    ax2 = Axis(fig[1, 2];   # row 1, col 2
        xlabel = "physical coordinate y",
        title  = "(b) sinh-mapped grids at N = $N_grid for several ℓ",
        limits = ((-30, 30), (0, 1)),
        yticksvisible = false, yticklabelsvisible = false)
    y_axis = collect(range(-30, 30, length = 601))
    profile = 0.4 ./ cosh.(y_axis)
    lines!(ax2, y_axis, profile .+ 0.45;
           color = (:gray, 0.4), linewidth = 0.8)
    text!(ax2, 15, 0.88; text = "target sech(y)",
          fontsize = 8, color = (:gray, 0.7))
    row_y = collect(range(0.32, 0.04, length = length(ell_values)))
    j = collect(0:N_grid - 1)
    t_nodes = -pi .+ 2 .* pi .* j ./ N_grid
    for (k, (ell, color)) in enumerate(zip(ell_values, palette_a))
        y_nodes = sinh.(ell .* t_nodes)
        y_show = filter(z -> abs(z) <= 30.0, collect(y_nodes))
        scatter!(ax2, y_show, fill(row_y[k], length(y_show));
                 color = color, marker = :vline,
                 markersize = 14, strokewidth = 0,
                 label = "ℓ = $ell")
    end
    axislegend(ax2, position = :rb, labelsize = 9, framevisible = false,
               nbanks = 2)

    y_max = sinh(ell_demo * pi)
    ax3 = Axis(fig[2, 1];
        xlabel = "physical coordinate y (or scaled t)",
        ylabel = "function value",
        title  = "(c) sech(y) in two coordinate frames at ℓ = $ell_demo",
        limits = ((-y_max, y_max), (-0.1, 1.1)))
    y_phys = collect(range(-y_max, y_max, length = 801))
    sech_in_y = 1.0 ./ cosh.(y_phys)
    lines!(ax3, y_phys, sech_in_y; color = NAVY, linewidth = 1.3,
           label = "sech(y) in y-space")
    t_axis = collect(range(-pi, pi, length = 801))
    sech_in_t = 1.0 ./ cosh.(sinh.(ell_demo .* t_axis))
    t_scaled = t_axis .* (y_max / pi)
    lines!(ax3, t_scaled, sech_in_t; color = CORAL, linewidth = 1.3,
           linestyle = :dash,
           label = "same function in t-space (scaled)")
    hlines!(ax3, [0.0]; color = :black, linewidth = 0.4, alpha = 0.5)
    axislegend(ax3, position = :rt, labelsize = 9, framevisible = false)

    # ---- (d) NEW: Fourier-in-t coefficient diagnostic ----
    N_diag = 96
    j_d = collect(0:N_diag-1)
    t_d = -pi .+ 2 .* pi .* j_d ./ N_diag
    ax4 = Axis(fig[2, 2];
               xlabel = "Fourier mode k",
               ylabel = "|fhat_k| on the t-grid",
               title  = "(d) Fourier-in-t coefficients of sech(y(t)), N = $N_diag",
               yscale = log10,
               limits = (nothing, (1e-17, 1.0)))
    for (ell, color) in zip(ell_values, palette_a)
        fv = 1.0 ./ cosh.(sinh.(ell .* t_d))
        F = abs.(fft(fv)) ./ N_diag
        Fmag = F[1:div(N_diag, 2) + 1]
        ks = collect(0:length(Fmag) - 1)
        scatterlines!(ax4, ks, Fmag .+ 1e-18; color = color,
                      markercolor = :white, strokecolor = color,
                      strokewidth = 1.0, markersize = 3, linewidth = 0.9,
                      label = "ℓ = $ell")
    end
    axislegend(ax4, position = :rt, labelsize = 9, framevisible = false)

    save(joinpath(OUTPUT_DIR, "weideman_cloot_map.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "weideman_cloot_map.png"), fig, px_per_unit = 4)
    return Dict("ell_values" => collect(ell_values),
                "N_grid" => N_grid, "ell_demo" => ell_demo,
                "y_max_per_ell" => [sinh(e * pi) for e in ell_values])
end

if abspath(PROGRAM_FILE) == @__FILE__
    r = make_figure()
    println("[Section 20.x / illustrative figure]  Weideman--Cloot sinh map")
    println("  ell values shown: ", r["ell_values"])
    println("  y_max per ell: ",
            [round(y, digits=2) for y in r["y_max_per_ell"]])
    println("  figure: ", joinpath(OUTPUT_DIR, "weideman_cloot_map.pdf"))
end
