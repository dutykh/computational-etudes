# arctan_tan_sweep.jl
# Chapter 19: Coordinate Transformations
# Computational Etude 19.7: A localised periodic pulse in the right coordinate.
#
# Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)

using CairoMakie
using FFTW
using Printf

const KAPPA = 80.0

function eval_error(N, ell, y_eval, truth)
    x  = [-pi + 2pi * k / N for k in 0:N-1]
    y  = 2.0 .* atan.(ell .* tan.(x ./ 2))
    fv = exp.(-KAPPA .* (1 .- cos.(y)))
    coeffs = fft(fv) ./ N
    ks = [0:div(N, 2)-1; -div(N, 2):-1]
    x_eval = 2.0 .* atan.(tan.(y_eval ./ 2) ./ ell)
    vals = zeros(length(y_eval))
    for m in 1:N
        vals .+= real.(coeffs[m] .* exp.(im .* ks[m] .* (x_eval .+ pi)))
    end
    return maximum(abs.(vals .- truth))
end

function run()
    outdir = joinpath(@__DIR__, "..", "..", "..",
                      "textbook", "figures", "ch19", "julia")
    mkpath(outdir)

    Ns = [12, 16, 24, 32, 48, 64, 96]
    Ls = [0.08, 0.10, 0.15, 0.20, 0.30, 0.45, 0.60, 0.80, 1.0, 1.5]
    y_eval = collect(range(-pi + 1e-9, pi - 1e-9, length=4097))
    truth  = exp.(-KAPPA .* (1 .- cos.(y_eval)))
    E = zeros(length(Ns), length(Ls))
    for (i, N) in enumerate(Ns), (j, ell) in enumerate(Ls)
        E[i, j] = eval_error(N, ell, y_eval, truth)
    end

    NAVY = colorant"#142D6E"; CORAL = colorant"#E74C3C"; TEAL = colorant"#16A085"
    ORANGE = colorant"#E67E22"; PURPLE = colorant"#8E44AD"
    fig = Figure(size=(1320, 360))

    ax1 = Axis(fig[1, 1], xlabel="y", title="(a) Mapped grids, N=32",
               limits=((-pi, pi), (-0.35, 1.1)))
    yy = collect(range(-pi, pi, length=401))
    lines!(ax1, yy, exp.(-KAPPA .* (1 .- cos.(yy))); color=NAVY, linewidth=1.2)
    Ng = 32
    for (ell, col, off) in [(0.1, CORAL, -0.08), (0.3, TEAL, -0.18), (1.0, ORANGE, -0.28)]
        x = [-pi + 2pi * k / Ng for k in 0:Ng-1]
        yc = 2.0 .* atan.(ell .* tan.(x ./ 2))
        scatter!(ax1, yc, fill(off, length(yc)); color=col, markersize=6)
    end

    ax2 = Axis(fig[1, 2], xlabel="log10 ell", ylabel="N", title="(b) Error landscape")
    heatmap!(ax2, log10.(Ls), Ns, log10.(E .+ 1e-16)')
    Colorbar(fig[1, 2][1, 2], limits=extrema(log10.(E .+ 1e-16)))

    ax3 = Axis(fig[1, 3], xlabel="ell", ylabel="max error", xscale=log10,
               yscale=log10, title="(c) Slices at fixed N")
    slice_Ns = [16, 24, 32, 48, 64]
    colours = [CORAL, ORANGE, TEAL, NAVY, PURPLE]
    for (k, N_sl) in enumerate(slice_Ns)
        i = findfirst(==(N_sl), Ns)
        scatterlines!(ax3, Ls, max.(E[i, :], 1e-16); color=colours[k], label="N=$N_sl")
    end
    axislegend(ax3; position=:rb)

    save(joinpath(outdir, "arctan_tan_sweep.pdf"), fig)
    save(joinpath(outdir, "arctan_tan_sweep.png"), fig)
    @printf("[19.7-julia] saved figure\n")
end

if abspath(PROGRAM_FILE) == @__FILE__
    run()
end
