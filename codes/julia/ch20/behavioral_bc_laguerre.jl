# behavioral_bc_laguerre.jl
# Chapter 20: Spectral Methods on Unbounded Intervals
# Computational Etude 20.7: Behavioural boundary condition at infinity.
#
# Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)

using CairoMakie
using LinearAlgebra
using Printf

include(joinpath(@__DIR__, "..", "ch07", "cheb_matrix.jl"))

function tln_dmatrices(N, L)
    Dx, x = cheb_matrix(N)
    fp = 2L ./ (1 .- x) .^ 2
    fpp = 4L ./ (1 .- x) .^ 3
    Dy = Diagonal(1 ./ fp) * Dx
    Dy2 = Diagonal(1 ./ fp .^ 2) * (Dx * Dx) - Diagonal(fpp ./ fp .^ 3) * Dx
    y_full = L .* (1 .+ x) ./ (1 .- x)
    return y_full[2:end], Dy[2:end, 2:end], Dy2[2:end, 2:end]
end

function solve_A(N, L)
    y, Dy, Dy2 = tln_dmatrices(N, L)
    Y = Diagonal(y)
    A = Y * Dy2 + (Y + I) * Dy
    return sort(real.(eigen(-Matrix(A)).values))
end

function solve_B(N, L)
    y, Dy, Dy2 = tln_dmatrices(N, L)
    Y = Diagonal(y); M = Diagonal(0.5 .+ 0.25 .* y)
    A = Y * Dy2 + Dy - M
    return sort(real.(eigen(-Matrix(A)).values))
end

function count_good(eigs; tol=0.05)
    targets = 0:49
    taken = Set{Int}()
    good = 0
    for lam in eigs
        if lam < -0.5
            continue
        end
        order = sortperm(abs.(lam .- targets))
        for idx in order
            n = targets[idx]
            if n in taken
                continue
            end
            rel = abs(lam - n) / max(1.0, n)
            if rel < tol
                push!(taken, n); good += 1
            end
            break
        end
    end
    return good
end

function run()
    outdir = joinpath(@__DIR__, "..", "..", "..",
                      "textbook", "figures", "ch20", "julia")
    mkpath(outdir)
    L = 32.0
    Ns = [10, 20, 30, 40, 60, 80, 120]
    good_A = [count_good(solve_A(N, L)) for N in Ns]
    good_B = [count_good(solve_B(N, L)) for N in Ns]

    NAVY = colorant"#142D6E"; CORAL = colorant"#E74C3C"; TEAL = colorant"#16A085"
    fig = Figure(size=(1020, 340))

    ax1 = Axis(fig[1, 1], xlabel="N", ylabel="good eigenvalues",
               title="(a) Good-eigenvalue count")
    scatterlines!(ax1, Ns, good_A; color=CORAL, label="naive")
    scatterlines!(ax1, Ns, good_B; color=TEAL, label="behavioural")
    axislegend(ax1; position=:lt)

    eigs_A = solve_A(40, L)[1:20]
    eigs_B = solve_B(40, L)[1:20]
    ax2 = Axis(fig[1, 2], xlabel="Re(lambda)", ylabel="Im(lambda)",
               title="(b) Spectrum at N=40",
               limits=((-1, 20), (-3, 3)))
    scatter!(ax2, real.(eigs_A), imag.(eigs_A); color=CORAL, label="A naive")
    scatter!(ax2, real.(eigs_B), imag.(eigs_B); color=TEAL, label="B behavioural", marker=:rect)
    hlines!(ax2, [0]; color=:gray, linewidth=0.5)
    axislegend(ax2; position=:rt)

    save(joinpath(outdir, "behavioral_bc_laguerre.pdf"), fig)
    save(joinpath(outdir, "behavioral_bc_laguerre.png"), fig)
    @printf("[20.7-julia] saved figure\n")
end

if abspath(PROGRAM_FILE) == @__FILE__
    run()
end
