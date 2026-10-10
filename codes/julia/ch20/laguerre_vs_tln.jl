# laguerre_vs_tln.jl
# Chapter 20: Spectral Methods on Unbounded Intervals
# Computational Etude 20.5: Laguerre vs TL_n on [0, infty).
#
# Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)

using CairoMakie
using FFTW
using LinearAlgebra
using Printf

include(joinpath(@__DIR__, "..", "ch07", "cheb_matrix.jl"))

function laguerre_poly(n, y)
    yv = collect(y)
    ell = zeros(n + 1, length(yv))
    ell[1, :] .= 1.0
    if n >= 1
        ell[2, :] .= 1.0 .- yv
    end
    for k in 1:(n - 1)
        ell[k + 2, :] .= ((2k + 1 .- yv) .* ell[k + 1, :] .- k .* ell[k, :]) ./ (k + 1)
    end
    return ell
end

# Golub-Welsch Gauss-Laguerre
function gauss_laguerre(K)
    a = collect(2 .* (0:K - 1) .+ 1.0)
    b = collect(1.0:K - 1)
    J = diagm(0 => a, 1 => b, -1 => b)
    vals, vecs = eigen(J)
    x = real(vals)
    w = real.(vecs[1, :] .^ 2)
    return x, w
end

function laguerre_error(f, N)
    x, w = gauss_laguerre(N + 32)
    ell = laguerre_poly(N, x)
    fv = f.(x)
    # a_n = int f(y) exp(-y/2) L_n(y) dy, with Gauss-Laguerre Gauss-Laguerre integrates e^{-x} g(x)
    # so we insert e^{x/2}: a_n = sum_i w_i exp(x_i) * f(x_i) e^{-x_i/2} L_n(x_i) = sum w_i e^{x_i/2} f(x_i) L_n(x_i)
    c = vec(sum(ell .* (w .* exp.(x ./ 2) .* fv)', dims=2))
    y = collect(range(0.001, 60.0, length=4001))
    Ly = laguerre_poly(N, y)
    approx = exp.(-y ./ 2) .* vec(sum(c .* Ly, dims=1))
    return maximum(abs.(approx .- f.(y)))
end

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

function tln_error(f, N, ell)
    _, x = cheb_matrix(N)
    fv = zeros(length(x))
    ok = abs.(x) .< 1.0 - 1e-12
    y_ok = ell .* (1 .+ x[ok]) ./ (1 .- x[ok])
    fv[ok] .= f.(y_ok)
    a = dct1_coeffs(fv)
    y_fine = collect(range(0.001, 60.0, length=4001))
    x_fine = (y_fine .- ell) ./ (y_fine .+ ell)
    approx = cheb_eval(a, x_fine, N)
    return maximum(abs.(approx .- f.(y_fine)))
end

function run()
    outdir = joinpath(@__DIR__, "..", "..", "..",
                      "textbook", "figures", "ch20", "julia")
    mkpath(outdir)

    Ns = [8, 12, 16, 24, 32, 48, 64, 96]
    err_lag_e = [laguerre_error(y -> exp(-y), N) for N in Ns]
    err_lag_a = [laguerre_error(y -> 1 / (1 + y), N) for N in Ns]
    err_tln_e = [tln_error(y -> exp(-y), N, 2.0) for N in Ns]
    err_tln_a = [tln_error(y -> 1 / (1 + y), N, 5.0) for N in Ns]

    NAVY = colorant"#142D6E"; CORAL = colorant"#E74C3C"; TEAL = colorant"#16A085"
    fig = Figure(size=(1020, 340))

    ax1 = Axis(fig[1, 1], xlabel="N", ylabel="max error",
               yscale=log10, title="(a) f = e^{-y} (exponential)")
    scatterlines!(ax1, Ns, max.(err_lag_e, 1e-18); color=CORAL, label="Laguerre")
    scatterlines!(ax1, Ns, max.(err_tln_e, 1e-18); color=TEAL, label="TL_n ell=2")
    axislegend(ax1; position=:rt)

    ax2 = Axis(fig[1, 2], xlabel="N", ylabel="max error",
               yscale=log10, title="(b) f = 1/(1+y) (algebraic)")
    scatterlines!(ax2, Ns, max.(err_lag_a, 1e-18); color=CORAL, label="Laguerre")
    scatterlines!(ax2, Ns, max.(err_tln_a, 1e-18); color=TEAL, label="TL_n ell=5")
    axislegend(ax2; position=:rt)

    save(joinpath(outdir, "laguerre_vs_tln.pdf"), fig)
    save(joinpath(outdir, "laguerre_vs_tln.png"), fig)
    @printf("[20.5-julia] saved figure\n")
end

if abspath(PROGRAM_FILE) == @__FILE__
    run()
end
