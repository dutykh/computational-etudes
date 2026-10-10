# arctan_tan_map.jl
# Chapter 19: Coordinate Transformations
# Illustration of the period-2pi arctan/tan map (eq. 19.2),
#
#     y = 2 * atan( ell * tan(x / 2) ),    x, y in [-pi, pi].
#
# Two panels: (a) the map curves y(x) for several values of ell, with the
# identity (ell = 1) drawn as a dashed grey line for reference; (b) the
# images of a uniform N = 24 x-grid under the map, shown as tick marks on
# five stacked horizontal bands -- making visible the textbook claim that
# ell < 1 clusters near y = 0 and ell > 1 clusters near y = +/- pi.
#
# Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
# Part of "Computational Etudes: A Spectral Approach"

using CairoMakie
using Printf

const EPS_END = 1.0e-12

arctan_tan_map(x, ell) = 2.0 * atan(ell * tan(clamp(x, -pi + EPS_END, pi - EPS_END) / 2))

function run()
    outdir = joinpath(@__DIR__, "..", "..", "..",
                      "textbook", "figures", "ch19", "julia")
    mkpath(outdir)

    NAVY   = colorant"#142D6E"
    CORAL  = colorant"#E74C3C"
    TEAL   = colorant"#16A085"
    ORANGE = colorant"#E67E22"
    PURPLE = colorant"#8E44AD"
    GREY   = colorant"#737373"

    ells   = [0.1, 0.3, 1.0, 3.0, 10.0]
    cols   = [CORAL, TEAL, GREY, ORANGE, PURPLE]
    styles = [:solid, :solid, :dash, :solid, :solid]

    fig = Figure(size=(1100, 440))

    # --- Panel (a): the map curves y(x) -----------------------------------
    pi_ticks = [-pi, -pi/2, 0.0, pi/2, pi]
    pi_labels = [L"-\pi", L"-\pi/2", L"0", L"\pi/2", L"\pi"]

    ax1 = Axis(fig[1, 1];
        xlabel = L"computational coordinate $x$",
        ylabel = L"physical coordinate $y$",
        title  = L"(a) the map  $y = 2\,\arctan(\ell\,\tan(x/2))$",
        limits = ((-pi, pi), (-pi, pi)),
        aspect = 1,
        xticks = (pi_ticks, pi_labels),
        yticks = (pi_ticks, pi_labels),
    )
    x_line = collect(range(-pi + EPS_END, pi - EPS_END, length=1001))
    hlines!(ax1, [0.0]; color=(:grey, 0.5), linewidth=0.4)
    vlines!(ax1, [0.0]; color=(:grey, 0.5), linewidth=0.4)
    for (k, ell) in enumerate(ells)
        y_line = arctan_tan_map.(x_line, ell)
        lines!(ax1, x_line, y_line;
               color=cols[k], linestyle=styles[k], linewidth=1.8,
               label=@sprintf("ℓ = %g", ell))
    end
    axislegend(ax1; position=:lt, framevisible=false)

    # --- Panel (b): clustering of a uniform N = 24 x-grid -----------------
    ax2 = Axis(fig[1, 2];
        xlabel = L"physical coordinate $y$",
        title  = L"(b) image of a uniform $N = 24$ $x$-grid",
        limits = ((-pi, pi), (0.4, length(ells) + 0.6)),
        xticks = (pi_ticks, pi_labels),
        yticks = ([], String[]),
    )
    N = 24
    x_uniform = [-pi + 2pi * (j + 0.5) / N for j in 0:N-1]
    bands_y = length(ells):-1:1   # stack top-down
    for (k, ell) in enumerate(ells)
        band = bands_y[k]
        col = cols[k]
        # base line
        lines!(ax2, [-pi, pi], [band, band]; color=(:grey, 0.4), linewidth=0.6)
        # mapped tick marks
        y_mapped = arctan_tan_map.(x_uniform, ell)
        for y in y_mapped
            lines!(ax2, [y, y], [band - 0.18, band + 0.18];
                   color=col, linewidth=1.5)
        end
        # left-side annotation
        text!(ax2, -pi - 0.18, band;
              text=@sprintf("ℓ = %g", ell), align=(:right, :center),
              color=NAVY, fontsize=12)
    end

    save(joinpath(outdir, "arctan_tan_map.pdf"), fig)
    save(joinpath(outdir, "arctan_tan_map.png"), fig)
    @printf("[19.0-julia] saved figure to %s\n",
            joinpath(outdir, "arctan_tan_map.pdf"))
    for ell in ells
        y_mapped = sort(arctan_tan_map.(x_uniform, ell))
        gaps = diff(y_mapped)
        @printf("  ell=%5.2f  min gap=%.4f  max gap=%.4f  ratio=%.2f\n",
                ell, minimum(gaps), maximum(gaps), maximum(gaps) / minimum(gaps))
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    run()
end
