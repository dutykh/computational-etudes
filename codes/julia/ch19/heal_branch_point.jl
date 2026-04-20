# heal_branch_point.jl
# Chapter 19: Coordinate Transformations
# Computational Etude 19.5: Healing the square-root branch point.
#
# Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)

using CairoMakie
using FFTW
using Printf

include(joinpath(@__DIR__, "..", "ch07", "cheb_matrix.jl"))

function chebyshev_coeffs(v)
    N = length(v) - 1
    V = vcat(v, reverse(v[2:N]))
    A = real.(fft(V)) ./ N
    A[1] *= 0.5; A[N + 1] *= 0.5
    return A[1:N + 1]
end

function cheb_eval(a, xfine, N)
    T0 = ones(length(xfine)); T1 = copy(xfine)
    val = a[1] .* T0 .+ a[2] .* T1
    for n in 2:N
        Tk = 2 .* xfine .* T1 .- T0
        val .+= a[n + 1] .* Tk
        T0 = T1; T1 = Tk
    end
    return val
end

function chebyshev_error(N)
    _, x = cheb_matrix(N)
    v = sqrt.(max.(1 .- x .^ 2, 0.0))
    a = chebyshev_coeffs(v)
    xfine = collect(range(-1 + 1e-10, 1 - 1e-10, length=2001))
    val = cheb_eval(a, xfine, N)
    return maximum(abs.(val .- sqrt.(max.(1 .- xfine .^ 2, 0.0)))), abs.(a)
end

function mapped_error(N, Y)
    _, xi = cheb_matrix(N)
    y = Y .* xi
    f = 1.0 ./ cosh.(y)
    a = chebyshev_coeffs(f)
    xfine = collect(range(-1, 1, length=2001))
    yfine = Y .* xfine
    val = cheb_eval(a, xfine, N)
    return maximum(abs.(val .- 1.0 ./ cosh.(yfine))), abs.(a)
end

function run()
    outdir = joinpath(@__DIR__, "..", "..", "..",
                      "textbook", "figures", "ch19", "julia")
    mkpath(outdir)

    Ns = [8, 16, 24, 32, 48, 64, 96, 128]
    err_d = zeros(length(Ns)); err_m = zeros(length(Ns))
    for (i, N) in enumerate(Ns)
        err_d[i], _ = chebyshev_error(N)
        err_m[i], _ = mapped_error(N, 10.0)
    end
    _, a_d = chebyshev_error(64)
    _, a_m = mapped_error(64, 10.0)

    NAVY = colorant"#142D6E"; CORAL = colorant"#E74C3C"; TEAL = colorant"#16A085"
    fig = Figure(size=(1180, 340))

    ax1 = Axis(fig[1, 1], title="(a) Two views")
    Xp = collect(range(-1, 1, length=401))
    yp = collect(range(-6, 6, length=401))
    lines!(ax1, Xp, sqrt.(max.(1 .- Xp .^ 2, 0)); color=CORAL, linewidth=1.2, label="√(1-X²)")
    lines!(ax1, yp, 1.0 ./ cosh.(yp); color=TEAL, linewidth=1.2, linestyle=:dash, label="sech y")
    axislegend(ax1; position=:rt)

    ax2 = Axis(fig[1, 2], xlabel="N", ylabel="max error", xscale=log10, yscale=log10,
               title="(b) Algebraic vs geometric")
    scatterlines!(ax2, Ns, err_d; color=CORAL, label="direct")
    scatterlines!(ax2, Ns, max.(err_m, 1e-18); color=TEAL, label="tanh-mapped")
    lines!(ax2, Ns, 1.0 ./ Ns; color=NAVY, linestyle=:dot, label="1/N")
    axislegend(ax2; position=:rt)

    ax3 = Axis(fig[1, 3], xlabel="Chebyshev index", ylabel="|a_n|",
               yscale=log10, title="(c) Coefficient decay, N=64",
               limits=(nothing, (1e-17, 10.0)))
    scatter!(ax3, 0:length(a_d)-1, max.(a_d, 1e-18); color=CORAL, marker=:xcross)
    scatter!(ax3, 0:length(a_m)-1, max.(a_m, 1e-18); color=TEAL, marker=:circle)

    save(joinpath(outdir, "heal_branch_point.pdf"), fig)
    save(joinpath(outdir, "heal_branch_point.png"), fig)
    @printf("[19.5-julia] saved figure\n")
end

if abspath(PROGRAM_FILE) == @__FILE__
    run()
end
