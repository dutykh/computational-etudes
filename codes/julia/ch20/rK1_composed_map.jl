# rK1_composed_map.jl
# Chapter 20: Spectral Methods on Unbounded Intervals
# Computational Etude 20.10: r K_1(r) via arcsinh(exp y) composed map.
#
# Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)

using CairoMakie
using FFTW
using Printf
using SpecialFunctions

include(joinpath(@__DIR__, "..", "ch07", "cheb_matrix.jl"))

r_of_y(y) = asinh(exp(y))
y_of_r(r) = log(sinh(r))
target_rK1(r) = r * besselk(1, r)

function dct1_coeffs(v)
    N = length(v) - 1
    V = vcat(v, reverse(v[2:N]))
    A = real.(fft(V)) ./ N
    A[1] *= 0.5; A[N + 1] *= 0.5
    return A[1:N + 1]
end

function cheb_eval(a, x, N)
    T0 = ones(length(x)); T1 = copy(x)
    val = a[1] .* T0 .+ a[2] .* T1
    for n in 2:N
        Tk = 2 .* x .* T1 .- T0
        val .+= a[n + 1] .* Tk
        T0 = T1; T1 = Tk
    end
    return val
end

function tbn_approx(N, ell)
    _, x = cheb_matrix(N)
    fv = zeros(length(x))
    ok = abs.(x) .< 1.0 - 1e-12
    y_ok = ell .* x[ok] ./ sqrt.(1 .- x[ok] .^ 2)
    r_ok = r_of_y.(y_ok)
    fv[ok] .= target_rK1.(r_ok)
    # cheb_matrix returns x descending: x[1] = 1 (y = +infty -> r K_1 -> 0),
    # x[end] = -1 (y = -infty -> r -> 0 -> r K_1 -> 1).
    fv[1] = 0.0
    fv[end] = 1.0
    return dct1_coeffs(fv)
end

function evaluate(a, r, ell)
    y = y_of_r.(r)
    x = y ./ sqrt.(ell^2 .+ y .^ 2)
    return cheb_eval(a, x, length(a) - 1)
end

function run()
    outdir = joinpath(@__DIR__, "..", "..", "..",
                      "textbook", "figures", "ch20", "julia")
    mkpath(outdir)
    ell = 4.0
    r_fine = exp10.(range(-3.0, 1.6, length=400))
    truth = target_rK1.(r_fine)
    Ns = [8, 12, 16, 20, 24, 32, 48, 64]
    errs = [maximum(abs.(evaluate(tbn_approx(N, ell), r_fine, ell) .- truth)) for N in Ns]

    NAVY = colorant"#142D6E"; CORAL = colorant"#E74C3C"; TEAL = colorant"#16A085"
    fig = Figure(size=(1020, 340))

    ax1 = Axis(fig[1, 1], xlabel="r", ylabel="f(r)", xscale=log10,
               title="(a) r K_1(r) and N=16 approximation")
    lines!(ax1, r_fine, truth; color=NAVY, linewidth=1.2, label="exact")
    lines!(ax1, r_fine, evaluate(tbn_approx(16, ell), r_fine, ell);
           color=CORAL, linestyle=:dash, label="N=16 TB_n")
    axislegend(ax1; position=:rt)

    ax2 = Axis(fig[1, 2], xlabel="N", ylabel="max error",
               yscale=log10, title="(b) Subgeometric descent, ell=$ell")
    scatterlines!(ax2, Ns, max.(errs, 1e-18); color=TEAL)

    save(joinpath(outdir, "rK1_composed_map.pdf"), fig)
    save(joinpath(outdir, "rK1_composed_map.png"), fig)
    @printf("[20.10-julia] saved figure\n")
end

if abspath(PROGRAM_FILE) == @__FILE__
    run()
end
