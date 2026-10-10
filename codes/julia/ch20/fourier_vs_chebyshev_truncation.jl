# fourier_vs_chebyshev_truncation.jl
# Chapter 20: Spectral Methods on Unbounded Intervals
# Computational Etude 20.2: Fourier beats Chebyshev at the far boundary.
#
# Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)

using CairoMakie
using FFTW
using Printf

include(joinpath(@__DIR__, "..", "ch07", "cheb_matrix.jl"))

target(y) = 1.0 / cosh(y)

function dct1_coeffs(v)
    N = length(v) - 1
    V = vcat(v, reverse(v[2:N]))
    A = real.(fft(V)) ./ N
    A[1] *= 0.5; A[N + 1] *= 0.5
    return A[1:N + 1]
end

function cheb_eval(a, x, N)
    T0 = ones(length(x)); T1 = copy(x)
    val = a[1] .* T0 .+ (length(a) >= 2 ? a[2] : 0.0) .* T1
    for n in 2:N
        Tk = 2 .* x .* T1 .- T0
        val .+= a[n + 1] .* Tk
        T0 = T1; T1 = Tk
    end
    return val
end

function chebyshev_error(N, L)
    _, x = cheb_matrix(N)
    fv = target.(L .* x)
    a = dct1_coeffs(fv)
    y_fine = collect(range(-L + 1, L - 1, length=4001))
    approx = cheb_eval(a, y_fine ./ L, N)
    return maximum(abs.(approx .- target.(y_fine)))
end

function fourier_error(N, L)
    M = 2 * N
    y_nodes = [-L + 2L * j / M for j in 0:M-1]
    f_nodes = target.(y_nodes)
    coeffs = fft(f_nodes) ./ M
    ks = [0:div(M, 2) - 1; -div(M, 2):-1]
    y_fine = collect(range(-L + 1, L - 1, length=4001))
    t_fine = pi .* (y_fine .+ L) ./ L
    vals = zeros(length(y_fine))
    for m in 1:M
        vals .+= real.(coeffs[m] .* exp.(im .* ks[m] .* t_fine))
    end
    return maximum(abs.(vals .- target.(y_fine)))
end

function run()
    outdir = joinpath(@__DIR__, "..", "..", "..",
                      "textbook", "figures", "ch20", "julia")
    mkpath(outdir)

    L = 10.0
    Ns = [8, 12, 16, 24, 32, 48, 64, 96, 128]
    err_c = [chebyshev_error(N, L) for N in Ns]
    err_f = [fourier_error(N, L) for N in Ns]

    Nshow = 32
    _, x_GL = cheb_matrix(Nshow)
    y_GL = L .* x_GL
    y_F = [-L + 2L * j / (2 * Nshow) for j in 0:2 * Nshow - 1]

    NAVY = colorant"#142D6E"; CORAL = colorant"#E74C3C"; TEAL = colorant"#16A085"
    fig = Figure(size=(1020, 340))

    ax1 = Axis(fig[1, 1], xlabel="y", title="(a) Grid densities on [-L, L]",
               limits=((-L, L), (-0.25, 1.1)))
    y_plot = collect(range(-L, L, length=401))
    lines!(ax1, y_plot, target.(y_plot); color=NAVY, linewidth=1.0)
    scatter!(ax1, y_F, fill(-0.08, length(y_F)); color=CORAL, marker=:xcross)
    scatter!(ax1, y_GL, fill(-0.18, length(y_GL)); color=TEAL, markersize=6)

    ax2 = Axis(fig[1, 2], xlabel="N", ylabel="max error",
               yscale=log10, title="(b) Interior error at L=10")
    scatterlines!(ax2, Ns, err_c; color=TEAL, label="Chebyshev")
    scatterlines!(ax2, Ns, err_f; color=CORAL, label="Fourier")
    axislegend(ax2; position=:rt)

    save(joinpath(outdir, "fourier_vs_chebyshev_truncation.pdf"), fig)
    save(joinpath(outdir, "fourier_vs_chebyshev_truncation.png"), fig)
    @printf("[20.2-julia] saved figure\n")
end

if abspath(PROGRAM_FILE) == @__FILE__
    run()
end
