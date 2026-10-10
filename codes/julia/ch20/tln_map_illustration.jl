# tln_map_illustration.jl
#
# Chapter 20: Spectral Methods on Unbounded Intervals
# Illustrative figure for the rational Chebyshev TL_n basis on [0, infty).
#
# Three panels:
#   (a) the map y = ell * (1 + x) / (1 - x) for ell in {1, 2, 4, 8}
#   (b) collocation grids on [0, infty) at fixed N = 24 for the same ell
#   (c) first five basis functions TL_n(y) at ell = 2 over y in [0, 20]
#
# Author: Dr. Denys Dutykh (Khalifa University of Science and Technology)
# Part of "Computational Etudes: A Spectral Approach"

using CairoMakie
using Colors
using FFTW
using Printf

include(joinpath(@__DIR__, "..", "ch07", "cheb_matrix.jl"))

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
const PURPLE = colorant"#8E44AD"
const GOLD   = colorant"#D4A017"

const SCRIPT_DIR = @__DIR__
const OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..",
                            "textbook", "figures", "ch20", "julia")
mkpath(OUTPUT_DIR)

"TL_n(y) = T_n( (y - ell) / (y + ell) )."
function TL_n(n::Int, y::AbstractVector{<:Real}, ell::Real)
    t = (y .- ell) ./ (y .+ ell)
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

"DCT-I coefficients of sech(y) on the TL_n grid for the diagnostic panel."
function tln_coeffs_sech(N::Int, ell::Real)
    j = collect(0:N)
    x = cos.(pi .* j ./ N)
    fv = zeros(length(x))
    interior = abs.(x) .< 1 - 1e-12
    y = ell .* (1 .+ x[interior]) ./ (1 .- x[interior])
    fv[interior] .= 1 ./ cosh.(y)
    fv[x .< -1 + 1e-12] .= 1.0
    ext = vcat(fv, reverse(fv[2:N]))
    A = real.(fft(ext)) ./ N
    A[1] *= 0.5
    A[N+1] *= 0.5
    return A[1:N+1]
end


function make_figure(; ell_values = (1.0, 2.0, 4.0, 8.0),
                     N_grid::Int = 24, ell_basis::Real = 2.0)
    fig = Figure(size = (1100, 800))
    palette_a = [NAVY, CORAL, TEAL, ORANGE]

    ax1 = Axis(fig[1, 1];
        xlabel = "computational coordinate x",
        ylabel = "physical coordinate y",
        title  = "(a) the map y = ℓ (1 + x) / (1 - x)",
        limits = ((-1.05, 1.05), (0, 60)))
    x_dense = collect(range(-1.0, 0.985, length = 401))
    for (ell, color) in zip(ell_values, palette_a)
        y_dense = ell .* (1.0 .+ x_dense) ./ (1.0 .- x_dense)
        lines!(ax1, x_dense, y_dense; color = color, linewidth = 1.4,
               label = "ℓ = $ell")
    end
    hlines!(ax1, [0.0]; color = :black, linewidth = 0.4, alpha = 0.5)
    vlines!(ax1, [+1.0]; color = (:gray, 0.6), linestyle = :dot, linewidth = 0.6)
    axislegend(ax1, position = :lt, labelsize = 9, framevisible = false)

    ax2 = Axis(fig[1, 2];
        xlabel = "physical coordinate y",
        title  = "(b) collocation grids at N = $N_grid for several ℓ",
        limits = ((0, 30), (0, 1)),
        yticksvisible = false, yticklabelsvisible = false)
    y_ref = collect(range(0, 30, length = 601))
    profile = 0.4 ./ cosh.(y_ref ./ 4.0)
    lines!(ax2, y_ref, profile .+ 0.45; color = (:gray, 0.4), linewidth = 0.8)
    text!(ax2, 18, 0.88; text = "target sech(y/4)",
          fontsize = 8, color = (:gray, 0.7))
    row_y = collect(range(0.32, 0.04, length = length(ell_values)))
    _, x_int_full = cheb_matrix(N_grid)
    x_int = x_int_full[2:end-1]
    for (k, (ell, color)) in enumerate(zip(ell_values, palette_a))
        y_int = ell .* (1.0 .+ x_int) ./ (1.0 .- x_int)
        y_show = filter(z -> z <= 30.0, collect(y_int))
        scatter!(ax2, y_show, fill(row_y[k], length(y_show));
                 color = color, marker = :vline,
                 markersize = 14, strokewidth = 0,
                 label = "ℓ = $ell")
    end
    axislegend(ax2, position = :rt, labelsize = 9, framevisible = false,
               nbanks = 2)

    ax3 = Axis(fig[2, 1];
        xlabel = "physical coordinate y",
        ylabel = "TL_n(y)",
        title  = "(c) basis functions TL_0, …, TL_4 at ℓ = $ell_basis",
        limits = ((0, 20), (-1.2, 1.25)))
    y_basis = collect(range(0.001, 20.0, length = 801))
    palette_c = [NAVY, CORAL, TEAL, PURPLE, GOLD]
    for (n, color) in zip(0:4, palette_c)
        lines!(ax3, y_basis, TL_n(n, y_basis, ell_basis);
               color = color, linewidth = 1.2,
               label = "TL$(n)")
    end
    hlines!(ax3, [0.0]; color = :black, linewidth = 0.4, alpha = 0.5)
    hlines!(ax3, [-1.0, 1.0]; color = (:gray, 0.6),
            linestyle = :dot, linewidth = 0.5)
    axislegend(ax3, position = :rb, labelsize = 9, framevisible = false,
               nbanks = 5)

    # ----- (d) NEW: coefficient diagnostic for sech(y) -----
    ax4 = Axis(fig[2, 2];
        xlabel = "degree n",
        ylabel = "|a_n| of TL_n expansion",
        title  = "(d) coefficient diagnostic: TL_n expansion of sech(y)",
        yscale = log10,
        limits = (nothing, (1e-17, 10.0)))
    N_diag = 48
    diag_palette = [CORAL, NAVY, TEAL]
    for (ell, color) in zip([0.5, 2.0, 8.0], diag_palette)
        a = abs.(tln_coeffs_sech(N_diag, ell))
        scatterlines!(ax4, collect(0:N_diag), a .+ 1e-18;
                      color = color, markercolor = :white, strokecolor = color,
                      strokewidth = 1.0, markersize = 4, linewidth = 1.0,
                      label = "ℓ = $ell")
    end
    axislegend(ax4, position = :rt, labelsize = 9, framevisible = false)

    save(joinpath(OUTPUT_DIR, "tln_map_illustration.pdf"), fig)
    save(joinpath(OUTPUT_DIR, "tln_map_illustration.png"), fig, px_per_unit = 4)
    return Dict("ell_values" => collect(ell_values),
                "N_grid" => N_grid, "ell_basis" => ell_basis)
end

if abspath(PROGRAM_FILE) == @__FILE__
    r = make_figure()
    println("[Etude 20.x / illustrative figure]  rational Chebyshev TL_n map")
    println("  ell values shown: ", r["ell_values"])
    @printf("  N = %d for grid panel; ell = %g for basis panel\n",
            r["N_grid"], r["ell_basis"])
    println("  figure: ", joinpath(OUTPUT_DIR, "tln_map_illustration.pdf"))
end
