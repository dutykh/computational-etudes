# J0_oscillatory.jl
# Chapter 20: Spectral Methods on Unbounded Intervals
# Computational Etude 20.11: Teach the basis the oscillation first.
#
# Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)

using CairoMakie
using LinearAlgebra
using Printf
using SpecialFunctions

target(y) = sqrt(1.0 + y) * besselj0(y)

function build_TL(y, N, ell)
    x = (y .- ell) ./ (y .+ ell)
    M = zeros(length(y), N + 1)
    M[:, 1] .= 1.0
    if N >= 1
        M[:, 2] .= x
    end
    for k in 3:N + 1
        M[:, k] .= 2 .* x .* M[:, k - 1] .- M[:, k - 2]
    end
    return M
end

function naive_fit(y_samp, f_samp, N, ell)
    M = build_TL(y_samp, N, ell)
    return M \ f_samp
end

function aug_fit(y_samp, f_samp, N, ell)
    M = build_TL(y_samp, N, ell)
    C = cos.(y_samp .- pi / 4)
    S = sin.(y_samp .- pi / 4)
    D = hcat(M .* C, M .* S)
    ab = D \ f_samp
    return ab[1:N + 1], ab[N + 2:end]
end

naive_eval(c, y, ell) = build_TL(y, length(c) - 1, ell) * c

function aug_eval(a, b, y, ell)
    M = build_TL(y, length(a) - 1, ell)
    return (M * a) .* cos.(y .- pi / 4) .+ (M * b) .* sin.(y .- pi / 4)
end

function run()
    outdir = joinpath(@__DIR__, "..", "..", "..",
                      "textbook", "figures", "ch20", "julia")
    mkpath(outdir)
    ell = 4.0
    y_fine = collect(range(0.01, 50.0, length=8001))
    truth = target.(y_fine)
    y_samp = collect(range(0.01, 80.0, length=2001))
    f_samp = target.(y_samp)

    Ns = [4, 6, 8, 10, 15, 20, 30, 40]
    err_n = zeros(length(Ns)); err_a = zeros(length(Ns))
    for (i, N) in enumerate(Ns)
        c = naive_fit(y_samp, f_samp, N, ell)
        err_n[i] = maximum(abs.(naive_eval(c, y_fine, ell) .- truth))
        a, b = aug_fit(y_samp, f_samp, N, ell)
        err_a[i] = maximum(abs.(aug_eval(a, b, y_fine, ell) .- truth))
    end

    NAVY = colorant"#142D6E"; CORAL = colorant"#E74C3C"; TEAL = colorant"#16A085"
    fig = Figure(size=(1020, 340))

    ax1 = Axis(fig[1, 1], xlabel="y",
               title="(a) sqrt(1+y) J_0(y) with amp-phase decomposition",
               limits=((0, 20), nothing))
    lines!(ax1, y_fine, truth; color=NAVY, linewidth=1.0, label="sqrt(1+y) J_0")
    a15, b15 = aug_fit(y_samp, f_samp, 15, ell)
    lines!(ax1, y_fine, abs.(build_TL(y_fine, 15, ell) * a15);
           color=CORAL, linestyle=:dash, label="|a(y)|")
    lines!(ax1, y_fine, abs.(build_TL(y_fine, 15, ell) * b15);
           color=TEAL, linestyle=:dot, label="|phi(y)|")
    axislegend(ax1; position=:rt)

    ax2 = Axis(fig[1, 2], xlabel="N", ylabel="max error",
               yscale=log10, title="(b) Asymptotic augmentation wins")
    scatterlines!(ax2, Ns, max.(err_n, 1e-18); color=CORAL, label="naive")
    scatterlines!(ax2, Ns, max.(err_a, 1e-18); color=TEAL, label="augmented")
    axislegend(ax2; position=:rt)

    save(joinpath(outdir, "J0_oscillatory.pdf"), fig)
    save(joinpath(outdir, "J0_oscillatory.png"), fig)
    @printf("[20.11-julia] saved figure\n")
end

if abspath(PROGRAM_FILE) == @__FILE__
    run()
end
