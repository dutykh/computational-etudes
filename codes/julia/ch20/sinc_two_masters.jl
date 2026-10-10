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
    fig = Figure(size=(1100, 800))

    j_star = argmin(err_sweep)
    h_theory = sqrt(pi^2 / (2 * N_fix))

    # ---- (a) V-shape: empirical error vs h ----
    ax1 = Axis(fig[1, 1], xlabel="h", ylabel="max error",
               yscale=log10, title="(a) Fix N=$N_fix, vary h")
    lines!(ax1, hs, err_sweep; color=CORAL, linewidth=1.2)
    vlines!(ax1, [hs[j_star]]; linestyle=:dot, color=NAVY,
            label="h* = $(round(hs[j_star]; digits=2)) empirical")
    vlines!(ax1, [h_theory]; linestyle=:dash, color=TEAL,
            label="sqrt(pi^2/(2N)) = $(round(h_theory; digits=2))")
    axislegend(ax1; position=:lt, labelsize=9)

    # ---- (b) NEW: two-master decomposition ----
    hs_an = collect(range(minimum(hs), maximum(hs), length=400))
    E_W  = (4 / pi) .* exp.(-pi^2 ./ (2 .* hs_an))
    E_DT = 4 .* exp.(-N_fix .* hs_an ./ 2)
    ax2 = Axis(fig[1, 2], xlabel="h", ylabel="error",
               yscale=log10,
               title="(b) Two masters: bandwidth vs grid-span",
               limits=(nothing, (1e-12, 1e1)))
    lines!(ax2, hs_an, E_W;  color=NAVY,  linestyle=:dash, linewidth=1.0,
           label="E_W (bandwidth) ~ exp(-pi^2/(2h))")
    lines!(ax2, hs_an, E_DT; color=CORAL, linestyle=:dash, linewidth=1.0,
           label="E_DT (grid span) ~ exp(-Nh/2)")
    lines!(ax2, hs_an, E_W .+ E_DT; color=TEAL, linewidth=1.4,
           label="sum: E_W + E_DT")
    lines!(ax2, hs, err_sweep; color=:gray, linewidth=0.7,
           label="empirical V")
    vlines!(ax2, [hs[j_star]]; linestyle=:dot, color=NAVY, linewidth=0.8)
    axislegend(ax2; position=:lb, labelsize=8)

    # ---- (c) Subgeometric convergence at h_opt ----
    ax3 = Axis(fig[2, 1], xlabel="sqrt(N)", ylabel="max error",
               yscale=log10, title="(c) Subgeometric convergence at h = h*")
    scatterlines!(ax3, sqrt.(Ns), err_opt; color=TEAL,
                  markercolor=:white, strokecolor=TEAL, strokewidth=1.0,
                  label="sinc at h = sqrt(pi^2/(2N))")
    lines!(ax3, sqrt.(Ns), exp.(-pi .* sqrt.(Ns ./ 2));
           color=NAVY, linestyle=:dot, label="exp(-pi sqrt(N/2))")
    axislegend(ax3; position=:rt, labelsize=9)

    # ---- (d) Grid cartoon ----
    ax4 = Axis(fig[2, 2], xlabel="y",
               title="(d) Span and spacing grow together",
               limits=((-10, 10), (0, 1.2)))
    for (k, (Nc, offset)) in enumerate([(8, 0.7), (32, 0.3)])
        hk = sqrt(pi^2 / (2 * Nc))
        j = (-div(Nc, 2)):div(Nc, 2)
        yj = j .* hk
        scatter!(ax4, yj, fill(offset, length(yj));
                 color=(k == 1 ? CORAL : TEAL), marker=:xcross, markersize=10,
                 label="N=$Nc, h=$(round(hk; digits=2))")
    end
    y_plot = collect(range(-10.0, 10.0, length=401))
    lines!(ax4, y_plot, 0.05 .+ 0.6 .* target.(y_plot); color=NAVY, linewidth=1.0)
    axislegend(ax4; position=:rt, labelsize=9)

    save(joinpath(outdir, "sinc_two_masters.pdf"), fig)
    save(joinpath(outdir, "sinc_two_masters.png"), fig)
    @printf("[20.3-julia] saved figure\n")
end

if abspath(PROGRAM_FILE) == @__FILE__
    run()
end
