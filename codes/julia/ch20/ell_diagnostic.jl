# L_diagnostic.jl
# Chapter 20: Spectral Methods on Unbounded Intervals
# Computational Etude 20.9: Read the coefficients before the error.
#
# Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)

using CairoMakie
using FFTW
using Printf

include(joinpath(@__DIR__, "..", "ch07", "cheb_matrix.jl"))

function dct1_coeffs(v)
    N = length(v) - 1
    V = vcat(v, reverse(v[2:N]))
    A = real.(fft(V)) ./ N
    A[1] *= 0.5; A[N + 1] *= 0.5
    return A[1:N + 1]
end

function tbn_coeffs(N, L)
    _, x = cheb_matrix(N)
    fv = zeros(length(x))
    ok = abs.(x) .< 1.0 - 1e-12
    y_ok = L .* x[ok] ./ sqrt.(1 .- x[ok] .^ 2)
    fv[ok] .= 1.0 ./ cosh.(y_ok)
    return dct1_coeffs(fv)
end

function run()
    outdir = joinpath(@__DIR__, "..", "..", "..",
                      "textbook", "figures", "ch20", "julia")
    mkpath(outdir)

    NAVY = colorant"#142D6E"; CORAL = colorant"#E74C3C"; TEAL = colorant"#16A085"
    ORANGE = colorant"#E67E22"; PURPLE = colorant"#8E44AD"
    fig = Figure(size=(1160, 400))

    N = 64
    L_list = [0.5, 1.0, 2.0, 4.0, 8.0, 16.0]
    cols = [CORAL, ORANGE, TEAL, PURPLE, NAVY, colorant"#8B4513"]
    ax1 = Axis(fig[1, 1], xlabel="n", ylabel="|a_n|",
               yscale=log10, title="(a) Coefficient decay, N=64",
               limits=(nothing, (1e-17, 10.0)))
    for (i, L) in enumerate(L_list)
        a = abs.(tbn_coeffs(N, L))
        scatterlines!(ax1, 0:N, max.(a, 1e-18); color=cols[i], label="L=$L")
    end
    axislegend(ax1; position=:lb)

    Ls = [0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 16.0, 24.0]
    ax2 = Axis(fig[1, 2], xlabel="L", ylabel="tail sum",
               xscale=log10, yscale=log10,
               title="(b) Valley of good L broadens with N")
    for (N_ref, col, label) in [(24, CORAL, "N=24"), (48, TEAL, "N=48"), (96, NAVY, "N=96")]
        errs = [sum(abs.(tbn_coeffs(N_ref, L)[(div(N_ref, 2) + 1):end])) for L in Ls]
        scatterlines!(ax2, Ls, errs; color=col, label=label)
    end
    axislegend(ax2; position=:rb)

    save(joinpath(outdir, "L_diagnostic.pdf"), fig)
    save(joinpath(outdir, "L_diagnostic.png"), fig)
    @printf("[20.9-julia] saved figure\n")
end

if abspath(PROGRAM_FILE) == @__FILE__
    run()
end
