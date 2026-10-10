# tbn_humiliates_hermite.jl
# Chapter 20: Spectral Methods on Unbounded Intervals
# Computational Etude 20.6: The function that humiliates Hermite.
#
# Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)

using CairoMakie
using FFTW
using LinearAlgebra
using Printf

include(joinpath(@__DIR__, "..", "ch07", "cheb_matrix.jl"))

target(y) = 1.0 / (1.0 + y^2)

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

function gauss_hermite(K)
    a = sqrt.((1:K - 1) ./ 2)
    J = diagm(1 => a, -1 => a)
    vals, vecs = eigen(J)
    return real(vals), sqrt(pi) .* real.(vecs[1, :] .^ 2)
end

function hermite_error(N)
    x, w = gauss_hermite(N + 40)
    fv = target.(x)
    psi = hermite_psi(N, x)
    integrand = (fv .* exp.(x .^ 2)) .* permutedims(psi)
    c = vec(sum(integrand .* w, dims=1))
    y = collect(range(-40.0, 40.0, length=8001))
    psi_y = hermite_psi(N, y)
    approx = vec(sum(c .* psi_y, dims=1))
    return maximum(abs.(approx .- target.(y)))
end

function sinc_error(N)
    h = sqrt(pi^2 / (2N))
    j = (-div(N, 2)):div(N, 2)
    yj = j .* h
    fj = target.(yj)
    y = collect(range(-40.0, 40.0, length=8001))
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
    return maximum(abs.(vals .- target.(y)))
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

function tbn_error(N, ell)
    _, x = cheb_matrix(N)
    fv = zeros(length(x))
    ok = abs.(x) .< 1.0 - 1e-12
    y_ok = ell .* x[ok] ./ sqrt.(1 .- x[ok] .^ 2)
    fv[ok] .= target.(y_ok)
    a = dct1_coeffs(fv)
    y_fine = collect(range(-40.0, 40.0, length=8001))
    x_fine = y_fine ./ sqrt.(ell^2 .+ y_fine .^ 2)
    approx = cheb_eval(a, x_fine, N)
    return maximum(abs.(approx .- target.(y_fine)))
end

function run()
    outdir = joinpath(@__DIR__, "..", "..", "..",
                      "textbook", "figures", "ch20", "julia")
    mkpath(outdir)

    Ns = [8, 12, 16, 24, 32, 48, 64, 96, 128]
    err_her = [hermite_error(N) for N in Ns]
    err_sinc = [sinc_error(N) for N in Ns]
    err_tbn = [tbn_error(N, 1.0) for N in Ns]

    NAVY = colorant"#142D6E"; CORAL = colorant"#E74C3C"; TEAL = colorant"#16A085"
    fig = Figure(size=(1020, 340))

    ax1 = Axis(fig[1, 1], xlabel="y", ylabel="f", title="(a) An innocuous target")
    y = collect(range(-8.0, 8.0, length=401))
    lines!(ax1, y, target.(y); color=NAVY, linewidth=1.2,
           label="f = 1/(1+y^2)")
    axislegend(ax1; position=:rt)

    ax2 = Axis(fig[1, 2], xlabel="N", ylabel="max error",
               xscale=log10, yscale=log10,
               title="(b) Algebraic / subgeometric / geometric")
    scatterlines!(ax2, Ns, err_her; color=CORAL, label="Hermite")
    scatterlines!(ax2, Ns, err_sinc; color=TEAL, label="sinc")
    scatterlines!(ax2, Ns, max.(err_tbn, 1e-18); color=NAVY, label="TB_n ell=1")
    axislegend(ax2; position=:rt)

    save(joinpath(outdir, "tbn_humiliates_hermite.pdf"), fig)
    save(joinpath(outdir, "tbn_humiliates_hermite.png"), fig)
    @printf("[20.6-julia] saved figure\n")
end

if abspath(PROGRAM_FILE) == @__FILE__
    run()
end
