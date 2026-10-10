# rational_chebyshev_map.jl
#
# Chapter 18: Linear Spectral Eigenproblems
# Illustrative figure for the rational Chebyshev TB_n basis.
#
# Three panels:
#   (a) the map x = ell t / sqrt(1 - t^2) for ell in {1, 2, 4, 8}
#   (b) collocation grids on the real line at fixed N = 24 for the same ell
#   (c) first five basis functions TB_n(x) at ell = 2 over x in [-10, 10]
#
# Author: Dr. Denys Dutykh (Khalifa University of Science and Technology)
# Part of "Computational Etudes: A Spectral Approach"

using CairoMakie
using Colors
using Printf

include(joinpath(@__DIR__, "rational_chebyshev.jl"))
using .RationalChebyshev: rational_chebyshev_grid

set_theme!(Theme(
    fontsize = 11,
    fonts = (regular = "CMU Serif",
             bold    = "CMU Serif Bold",
             italic  = "CMU Serif Italic"),
    Axis = (xlabelsize = 11, ylabelsize = 11, titlesize = 11,
            xticklabelsize = 10, yticklabelsize = 10,
            spinewidth = 0.8, xtickwidth = 0.8, ytickwidth = 0.8),
    Legend = (labelsize = 9,),
))

const NAVY   = colorant"#142D6E"
const SKY    = colorant"#7896D2"
const CORAL  = colorant"#E74C3C"
const TEAL   = colorant"#16A085"
const ORANGE = colorant"#E67E22"
const PURPLE = colorant"#8E44AD"
const GOLD   = colorant"#D4A017"

const SCRIPT_DIR = @__DIR__
const OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..",
                            "textbook", "figures", "ch18", "julia")
mkpath(OUTPUT_DIR)

"""
    TB_n(n, x, ell)

Rational Chebyshev basis: TB_n(x) = T_n( x / sqrt(ell^2 + x^2) ),
evaluated by Chebyshev recurrence in the t-variable.
"""
function TB_n(n::Int, x::AbstractVector{<:Real}, ell::Real)
    t = x ./ sqrt.(ell^2 .+ x .^ 2)
    if n == 0
        return ones(length(t))
    elseif n == 1
        return collect(t)
    end
    Tkm2 = ones(length(t))
    Tkm1 = collect(t)
    Tk = similar(Tkm1)
    for _ in 2:n
        Tk = 2 .* t .* Tkm1 .- Tkm2
        Tkm2 = Tkm1
        Tkm1 = Tk
    end
    return Tk
end

function make_figure(; ell_values = (1.0, 2.0, 4.0, 8.0),
                     N_grid::Int = 24, ell_basis::Real = 2.0)
    fig = Figure(size = (1300, 400))
    palette_a = [NAVY, CORAL, TEAL, ORANGE]

    # ----- Panel (a): the map x(t) for several ell values -----
    ax1 = Axis(fig[1, 1];
        xlabel = "computational coordinate t",
        ylabel = "physical coordinate x",
        title  = "(a) the map x = ℓ t / √(1 - t²)",
        limits = ((-1.05, 1.05), (-15, 15)))
    t_dense = collect(range(-0.99, 0.99, length = 401))
    for (ell, color) in zip(ell_values, palette_a)
        x_dense = ell .* t_dense ./ sqrt.(1.0 .- t_dense .^ 2)
        lines!(ax1, t_dense, x_dense; color = color, linewidth = 1.4,
               label = "ℓ = $ell")
    end
    hlines!(ax1, [0.0]; color = :black, linewidth = 0.4, alpha = 0.5)
    vlines!(ax1, [-1.0, 1.0]; color = (:gray, 0.6), linestyle = :dot, linewidth = 0.6)
    axislegend(ax1, position = :lt, labelsize = 9, framevisible = false)

    # ----- Panel (b): grid clustering at N = N_grid -----
    ax2 = Axis(fig[1, 2];
        xlabel = "physical coordinate x",
        title  = "(b) collocation grids at N = $N_grid for several ℓ",
        limits = ((-15, 15), (0, 1)),
        yticksvisible = false, yticklabelsvisible = false)
    x_ref = collect(range(-15, 15, length = 601))
    profile = 0.4 .* exp.(-(x_ref ./ 4.0) .^ 2)
    lines!(ax2, x_ref, profile .+ 0.45; color = (:gray, 0.4), linewidth = 0.8)
    text!(ax2, 8.5, 0.88; text = "target exp(-(x/4)²)",
          fontsize = 8, color = (:gray, 0.7))
    row_y = collect(range(0.32, 0.04, length = length(ell_values)))
    for (k, (ell, color)) in enumerate(zip(ell_values, palette_a))
        x_int, _, _ = rational_chebyshev_grid(N_grid, ell)
        x_show = filter(z -> abs(z) <= 15.0, collect(x_int))
        scatter!(ax2, x_show, fill(row_y[k], length(x_show));
                 color = color, marker = :vline,
                 markersize = 14, strokewidth = 0,
                 label = "ℓ = $ell")
    end
    axislegend(ax2, position = :rb, labelsize = 9, framevisible = false,
               nbanks = 2)

    # ----- Panel (c): first five basis functions at ell = ell_basis -----
    ax3 = Axis(fig[1, 3];
        xlabel = "physical coordinate x",
        ylabel = "TB_n(x)",
        title  = "(c) basis functions TB_0, …, TB_4 at ℓ = $ell_basis",
        limits = ((-10, 10), (-1.2, 1.25)))
    x_basis = collect(range(-10.0, 10.0, length = 801))
    palette_c = [NAVY, CORAL, TEAL, PURPLE, GOLD]
    for (n, color) in zip(0:4, palette_c)
        lines!(ax3, x_basis, TB_n(n, x_basis, ell_basis);
               color = color, linewidth = 1.2,
               label = "TB$(n)")
    end
    hlines!(ax3, [0.0]; color = :black, linewidth = 0.4, alpha = 0.5)
    hlines!(ax3, [-1.0, 1.0]; color = (:gray, 0.6),
            linestyle = :dot, linewidth = 0.5)
    axislegend(ax3, position = :rb, labelsize = 9, framevisible = false,
               nbanks = 5)

    save(joinpath(OUTPUT_DIR, "rational_chebyshev_map.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "rational_chebyshev_map.png"), fig, px_per_unit = 4)
    return Dict("ell_values" => collect(ell_values),
                "N_grid"     => N_grid,
                "ell_basis"  => ell_basis)
end

if abspath(PROGRAM_FILE) == @__FILE__
    r = make_figure()
    println("[Etude 18.4 / illustrative figure]  rational Chebyshev map")
    println("  ell values shown: ", r["ell_values"])
    @printf("  N = %d for grid panel; ell = %g for basis panel\n",
            r["N_grid"], r["ell_basis"])
    println("  figure: ", joinpath(OUTPUT_DIR, "rational_chebyshev_map.pdf"))
end
