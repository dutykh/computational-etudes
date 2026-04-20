# chebyshev_as_cosine.jl
# Chapter 19: Coordinate Transformations
# Computational Etude 19.2: One problem, two coordinates.
#
# Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)

using CairoMakie
using FFTW
using LinearAlgebra
using Printf

include(joinpath(@__DIR__, "..", "ch07", "cheb_matrix.jl"))

const Q = 4.0

u_ex(x) = sin(pi .* x)
f_src(x) = -(pi^2) .* sin(pi .* x) .- Q .* sin(pi .* x)

function solve_x(N)
    D, x = cheb_matrix(N)
    D2 = D * D
    A = D2[2:N, 2:N] - Q * I(N - 1)
    rhs = f_src.(x[2:N])
    u_int = A \ rhs
    u = zeros(N + 1); u[2:N] = u_int
    return x, u
end

function chebfft(v)
    N = length(v) - 1
    x = cos.(pi .* (0:N) ./ N)
    V = vcat(v, reverse(v[2:N]))
    U = real.(fft(V))
    k = 0:N
    w_hat = im .* vcat(collect(k[1:N]), 0, collect(k[2:N]) .- N) .* U
    W = real.(ifft(w_hat))
    w = zeros(N + 1)
    w[2:N] .= -W[2:N] ./ sqrt.(1 .- x[2:N] .^ 2)
    k2 = (0:N) .^ 2
    w[1] = sum(k2 .* U[1:N+1]) / N + 0.5 * N * U[N+1]
    w[N+1] = sum(((-1) .^ ((0:N) .+ 1)) .* k2 .* U[1:N+1]) / N +
             0.5 * (-1)^(N + 1) * N * U[N+1]
    return w
end

function solve_t(N)
    x = cos.(pi .* (0:N) ./ N)
    D2 = zeros(N + 1, N + 1)
    e = zeros(N + 1)
    for j in 1:(N + 1)
        fill!(e, 0.0); e[j] = 1.0
        D2[:, j] = chebfft(chebfft(e))
    end
    A = D2[2:N, 2:N] - Q * I(N - 1)
    rhs = f_src.(x[2:N])
    u_int = A \ rhs
    u = zeros(N + 1); u[2:N] = u_int
    return x, u
end

function run()
    outdir = joinpath(@__DIR__, "..", "..", "..",
                      "textbook", "figures", "ch19", "julia")
    mkpath(outdir)
    Ns = [8, 12, 16, 20, 24, 32, 40, 48]
    err_x = zeros(length(Ns)); err_t = zeros(length(Ns))
    for (i, N) in enumerate(Ns)
        x, u = solve_x(N); err_x[i] = maximum(abs.(u .- u_ex.(x)))
        x, u = solve_t(N); err_t[i] = maximum(abs.(u .- u_ex.(x)))
    end

    NAVY = colorant"#142D6E"; CORAL = colorant"#E74C3C"; TEAL = colorant"#16A085"
    fig = Figure(size=(820, 340))

    ax1 = Axis(fig[1, 1], xlabel="x", ylabel="u(x)", title="(a) Manufactured BVP")
    xN, uN = solve_x(24)
    xfine = collect(range(-1, 1, length=401))
    lines!(ax1, xfine, u_ex.(xfine); color=NAVY, linewidth=1.2, label="exact")
    scatter!(ax1, xN, uN; color=CORAL, markersize=5, label="X-form N=24")
    axislegend(ax1; position=:rt)

    ax2 = Axis(fig[1, 2], xlabel="N", ylabel="max error",
               title="(b) Same convergence, different arithmetic", yscale=log10)
    scatterlines!(ax2, Ns, err_x .+ 1e-18; color=CORAL, label="X-form (DM)")
    scatterlines!(ax2, Ns, err_t .+ 1e-18; color=TEAL, label="T-form (FFT)")
    axislegend(ax2; position=:rt)

    save(joinpath(outdir, "chebyshev_as_cosine.pdf"), fig)
    save(joinpath(outdir, "chebyshev_as_cosine.png"), fig)
    @printf("[19.2-julia] saved figure\n")
end

if abspath(PROGRAM_FILE) == @__FILE__
    run()
end
