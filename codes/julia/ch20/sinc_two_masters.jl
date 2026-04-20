# sinc_two_masters.jl
# Chapter 20: Spectral Methods on Unbounded Intervals
# Computational Etude 20.3: Sinc -- two masters, span and spacing.
#
# Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)

using CairoMakie
using Printf

target(y) = 1.0 / cosh(y)

function sinc_approx(y, N, h)
    j = (-div(N, 2)):div(N, 2)
    yj = j .* h
    fj = target.(yj)
    vals = zeros(length(y))
    for (i, yi) in enumerate(y)
        acc = 0.0
        for k in 1:length(yj)
            z = (yi - yj[k]) / h
            s = abs(z) < 1e-14 ? 1.0 : sin(pi * z) / (pi * z)
            acc += s * fj[k]
        end
        vals[i] = acc
    end
    return vals
end

function sinc_error(N, h)
    y = collect(range(-20.0, 20.0, length=4001))
    approx = sinc_approx(y, N, h)
    return maximum(abs.(approx .- target.(y)))
end

function run()
    outdir = joinpath(@__DIR__, "..", "..", "..",
                      "textbook", "figures", "ch20", "julia")
    mkpath(outdir)

    N_fix = 48
    hs = collect(range(0.15, 1.8, length=80))
    err_sweep = [sinc_error(N_fix, h) for h in hs]

    Ns = [8, 12, 16, 24, 32, 48, 64, 96, 128, 192]
    h_opts = [sqrt(pi^2 / (2 * N)) for N in Ns]
    err_opt = [sinc_error(N, h) for (N, h) in zip(Ns, h_opts)]

    NAVY = colorant"#142D6E"; CORAL = colorant"#E74C3C"; TEAL = colorant"#16A085"
    fig = Figure(size=(1360, 340))

    j_star = argmin(err_sweep)
    h_theory = sqrt(pi^2 / (2 * N_fix))
    ax1 = Axis(fig[1, 1], xlabel="h", ylabel="max error",
               yscale=log10, title="(a) Fix N=$N_fix, vary h")
    lines!(ax1, hs, err_sweep; color=CORAL, linewidth=1.2)
    vlines!(ax1, [hs[j_star]]; linestyle=:dot, color=NAVY,
            label="empirical h*")
    vlines!(ax1, [h_theory]; linestyle=:dash, color=TEAL,
            label="sqrt(pi^2/(2N))")
    axislegend(ax1; position=:rt)

    ax2 = Axis(fig[1, 2], xlabel="sqrt(N)", ylabel="max error",
               yscale=log10, title="(b) Subgeometric convergence")
    scatterlines!(ax2, sqrt.(Ns), err_opt; color=TEAL, label="sinc")
    lines!(ax2, sqrt.(Ns), exp.(-pi .* sqrt.(Ns ./ 2));
           color=NAVY, linestyle=:dot, label="e^{-pi sqrt(N/2)}")
    axislegend(ax2; position=:rt)

    ax3 = Axis(fig[1, 3], xlabel="y", title="(c) Span+spacing grow together",
               limits=((-10, 10), (0, 1.2)))
    for (k, (Nc, offset)) in enumerate([(8, 0.7), (32, 0.3)])
        hk = sqrt(pi^2 / (2 * Nc))
        j = (-div(Nc, 2)):div(Nc, 2)
        yj = j .* hk
        scatter!(ax3, yj, fill(offset, length(yj));
                 color=(k == 1 ? CORAL : TEAL), marker=:xcross, markersize=10)
    end
    y_plot = collect(range(-10.0, 10.0, length=401))
    lines!(ax3, y_plot, 0.05 .+ 0.6 .* target.(y_plot); color=NAVY, linewidth=1.0)

    save(joinpath(outdir, "sinc_two_masters.pdf"), fig)
    save(joinpath(outdir, "sinc_two_masters.png"), fig)
    @printf("[20.3-julia] saved figure\n")
end

if abspath(PROGRAM_FILE) == @__FILE__
    run()
end
