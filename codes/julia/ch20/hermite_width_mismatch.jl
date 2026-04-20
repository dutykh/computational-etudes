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

    NAVY = colorant"#142D6E"; CORAL = colorant"#E74C3C"; TEAL = colorant"#16A085"
    PURPLE = colorant"#8E44AD"
    cols = [CORAL, NAVY, TEAL, PURPLE]
    fig = Figure(size=(1340, 340))

    ax1 = Axis(fig[1, 1], xlabel="N", ylabel="max error",
               yscale=log10, title="(a) alpha = 1 (unscaled)")
    for (ai, A) in enumerate(A_list)
        scatterlines!(ax1, Ns, err1[ai, :]; color=cols[ai], label="A=$A")
    end
    axislegend(ax1; position=:rt)

    ax2 = Axis(fig[1, 2], xlabel="N", ylabel="max error",
               yscale=log10, title="(b) alpha = sqrt(2A) (matched)")
    for (ai, A) in enumerate(A_list)
        scatterlines!(ax2, Ns, max.(err2[ai, :], 1e-18); color=cols[ai], label="A=$A")
    end
    axislegend(ax2; position=:rt)

    ax3 = Axis(fig[1, 3], xlabel="N", ylabel="max error",
               yscale=log10, title="(c) QHO: psi_0 exact at N=0")
    scatterlines!(ax3, Ns_q, max.(err_qho, 1e-18); color=NAVY, marker=:diamond)

    save(joinpath(outdir, "hermite_width_mismatch.pdf"), fig)
    save(joinpath(outdir, "hermite_width_mismatch.png"), fig)
    @printf("[20.4-julia] saved figure\n")
end

if abspath(PROGRAM_FILE) == @__FILE__
    run()
end
