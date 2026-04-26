# hermite_width_mismatch.jl
# Chapter 20: Spectral Methods on Unbounded Intervals
# Computational Etude 20.4: The cost of width mismatch for Hermite functions.
#
# Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)

using CairoMakie
using LinearAlgebra
using Printf

function hermite_psi(n, y)
    yv = collect(y)
    psi = zeros(n + 1, length(yv))
    psi[1, :] .= pi^(-0.25) .* exp.(-0.5 .* yv .^ 2)
    if n >= 1
        psi[2, :] .= sqrt(2) .* yv .* psi[1, :]
    end
    for k in 1:(n - 1)
        psi[k + 2, :] .= sqrt(2 / (k + 1)) .* yv .* psi[k + 1, :] .-
                         sqrt(k / (k + 1)) .* psi[k, :]
    end
    return psi
end

# Golub-Welsch Gauss-Hermite (physicists' normalisation)
function gauss_hermite(K)
    a = sqrt.((1:K - 1) ./ 2)
    J = diagm(1 => a, -1 => a)
    vals, vecs = eigen(J)
    x = real(vals)
    w = sqrt(pi) .* real.(vecs[1, :] .^ 2)
    return x, w
end

function hermite_expand(f, N, alpha)
    x, w = gauss_hermite(N + 32)
    y_nodes = x ./ alpha
    fv = f.(y_nodes)
    psi = hermite_psi(N, x)
    integrand = (fv .* exp.(x .^ 2)) .* permutedims(psi)
    # integrand is (K, N+1); sum with weights over first axis
    c = vec(sum(integrand .* w, dims=1)) ./ sqrt(alpha)
    return c
end

function hermite_maxerr(f, c, alpha)
    y = collect(range(-15.0, 15.0, length=4001))
    psi = hermite_psi(length(c) - 1, alpha .* y)
    approx = sqrt(alpha) .* vec(sum(c .* psi, dims=1))
    return maximum(abs.(approx .- f.(y)))
end

function run()
    outdir = joinpath(@__DIR__, "..", "..", "..",
                      "textbook", "figures", "ch20", "julia")
    mkpath(outdir)

    Ns = [8, 12, 16, 24, 32, 48, 64]
    A_list = [0.1, 0.5, 2.0, 8.0]
    err1 = zeros(length(A_list), length(Ns))
    err2 = zeros(length(A_list), length(Ns))
    for (ai, A) in enumerate(A_list)
        f = y -> exp(-A * y^2)
        alpha_m = sqrt(2A)
        for (ni, N) in enumerate(Ns)
            c1 = hermite_expand(f, N, 1.0)
            err1[ai, ni] = hermite_maxerr(f, c1, 1.0)
            c2 = hermite_expand(f, N, alpha_m)
            err2[ai, ni] = hermite_maxerr(f, c2, alpha_m)
        end
    end

    Ns_q = [1, 2, 4, 8, 16, 32]
    f_qho = y -> pi^(-0.25) * exp(-0.5 * y^2)
    err_qho = [hermite_maxerr(f_qho, hermite_expand(f_qho, N, 1.0), 1.0) for N in Ns_q]

    # NEW: coefficient decay at A = 8
    A_pick = 8.0
    f_pick = y -> exp(-A_pick * y^2)
    coeffs_unscaled = abs.(hermite_expand(f_pick, 32, 1.0))
    coeffs_matched  = abs.(hermite_expand(f_pick, 32, sqrt(2 * A_pick)))

    # NEW: optimal-alpha scan at fixed N = 16
    alphas = collect(range(0.3, 5.0, length=60))
    A_scan = [0.5, 2.0, 8.0]
    N_scan = 16
    err_scan = zeros(length(A_scan), length(alphas))
    for (ai, A) in enumerate(A_scan)
        f = y -> exp(-A * y^2)
        for (k, a) in enumerate(alphas)
            c = hermite_expand(f, N_scan, a)
            err_scan[ai, k] = hermite_maxerr(f, c, a)
        end
    end

    NAVY = colorant"#142D6E"; CORAL = colorant"#E74C3C"; TEAL = colorant"#16A085"
    PURPLE = colorant"#8E44AD"
    cols = [CORAL, NAVY, TEAL, PURPLE]
    fig = Figure(size=(1100, 1200))

    # ---- (a) Width-mismatch cartoon ----
    y_show = collect(range(-6, 6, length=401))
    psi0 = pi^(-0.25) .* exp.(-0.5 .* y_show .^ 2)
    ax1 = Axis(fig[1, 1], xlabel="y", ylabel="amplitude",
               title="(a) target exp(-Ay^2) vs basis envelope psi_0")
    lines!(ax1, y_show, psi0; color=:gray, linestyle=:dash, linewidth=1.0,
           label="psi_0 (basis envelope)")
    for (ai, A) in enumerate(A_list)
        lines!(ax1, y_show, exp.(-A .* y_show .^ 2); color=cols[ai], linewidth=1.2,
               label="A = $A")
    end
    axislegend(ax1; position=:rt, labelsize=8)

    # ---- (b) Coefficient decay at A = 8 ----
    ns = 0:length(coeffs_unscaled)-1
    ax2 = Axis(fig[1, 2], xlabel="degree n", ylabel="|a_n|",
               yscale=log10,
               title="(b) Coefficient decay at A = $A_pick")
    scatterlines!(ax2, collect(ns), coeffs_unscaled .+ 1e-20; color=CORAL,
                  markercolor=:white, strokecolor=CORAL, strokewidth=1.0,
                  markersize=4, linewidth=1.0,
                  label="alpha = 1 (slow algebraic)")
    scatter!(ax2, collect(ns), coeffs_matched .+ 1e-20; color=NAVY,
             markersize=5, label="alpha = sqrt(2A) (only a_0)")
    axislegend(ax2; position=:lb, labelsize=8)

    # ---- (c) Unscaled convergence ----
    ax3 = Axis(fig[2, 1], xlabel="N", ylabel="max error",
               yscale=log10, title="(c) alpha = 1 (unscaled)")
    for (ai, A) in enumerate(A_list)
        scatterlines!(ax3, Ns, err1[ai, :]; color=cols[ai], label="A=$A")
    end
    axislegend(ax3; position=:rt)

    # ---- (d) Matched convergence ----
    ax4 = Axis(fig[2, 2], xlabel="N", ylabel="max error",
               yscale=log10, title="(d) alpha = sqrt(2A) (matched: machine precision)")
    for (ai, A) in enumerate(A_list)
        scatterlines!(ax4, Ns, max.(err2[ai, :], 1e-18); color=cols[ai],
                      markercolor=:white, strokecolor=cols[ai], strokewidth=1.0,
                      label="A=$A")
    end
    axislegend(ax4; position=:rt)

    # ---- (e) Optimal-alpha scan ----
    ax5 = Axis(fig[3, 1], xlabel="basis scale alpha", ylabel="max error",
               yscale=log10,
               title="(e) Optimal-alpha scan at N = $N_scan (dotted: sqrt(2A))")
    cols_scan = [NAVY, TEAL, PURPLE]
    for (ai, A) in enumerate(A_scan)
        lines!(ax5, alphas, err_scan[ai, :] .+ 1e-18; color=cols_scan[ai],
               linewidth=1.0, label="A = $A")
        vlines!(ax5, [sqrt(2 * A)]; color=cols_scan[ai], linestyle=:dot,
                linewidth=0.8, alpha=0.6)
    end
    axislegend(ax5; position=:rt, labelsize=8)

    # ---- (f) QHO ground state ----
    ax6 = Axis(fig[3, 2], xlabel="N", ylabel="max error",
               yscale=log10, title="(f) QHO ground state: basis matches physics")
    scatterlines!(ax6, Ns_q, max.(err_qho, 1e-20); color=NAVY,
                  marker=:diamond, markercolor=:white, strokecolor=NAVY,
                  strokewidth=1.0, markersize=8, label="f = psi_0")
    axislegend(ax6; position=:rt)

    save(joinpath(outdir, "hermite_width_mismatch.pdf"), fig)
    save(joinpath(outdir, "hermite_width_mismatch.png"), fig)
    @printf("[20.4-julia] saved figure\n")
end

if abspath(PROGRAM_FILE) == @__FILE__
    run()
end
