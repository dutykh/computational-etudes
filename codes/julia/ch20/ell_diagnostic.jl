# ell_diagnostic.jl
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

function tbn_coeffs(N, ell)
    _, x = cheb_matrix(N)
    fv = zeros(length(x))
    ok = abs.(x) .< 1.0 - 1e-12
    y_ok = ell .* x[ok] ./ sqrt.(1 .- x[ok] .^ 2)
    fv[ok] .= 1.0 ./ cosh.(y_ok)
    return dct1_coeffs(fv)
end

function envelope_max(a::AbstractVector, win::Int = 3)
    n = length(a)
    out = similar(a)
    for i in 1:n
        lo = max(1, i - win); hi = min(n, i + win)
        out[i] = maximum(a[lo:hi])
    end
    return out
end

function run()
    outdir = joinpath(@__DIR__, "..", "..", "..",
                      "textbook", "figures", "ch20", "julia")
    mkpath(outdir)

    NAVY = colorant"#142D6E"; CORAL = colorant"#E74C3C"; TEAL = colorant"#16A085"
    ORANGE = colorant"#E67E22"; PURPLE = colorant"#8E44AD"
    fig = Figure(size=(1100, 800))

    N = 64
    L_full = [0.5, 1.0, 2.0, 4.0, 8.0, 16.0]
    cols_full = [CORAL, ORANGE, TEAL, PURPLE, NAVY, colorant"#8B4513"]

    # ---- (a) Three regimes ----
    ax1 = Axis(fig[1, 1], xlabel="degree n", ylabel="envelope of |a_n|",
               yscale=log10,
               title="(a) three regimes of ell (envelope view), N=64",
               limits=(nothing, (1e-17, 10.0)))
    three = [(0.5, CORAL, "small (early flatten)"),
             (2.0, NAVY, "good (clean descent)"),
             (16.0, TEAL, "large (gentle small-n slope)")]
    for (ell, color, lbl) in three
        a = abs.(tbn_coeffs(N, ell))
        env = envelope_max(a, 3)
        ns = collect(0:N)
        lines!(ax1, ns, env .+ 1e-18; color = color, linewidth = 1.2,
               label = "ell = $ell — $lbl")
        scatter!(ax1, ns[1:4:end], env[1:4:end] .+ 1e-18;
                 color = :white, strokecolor = color, strokewidth = 1.0,
                 markersize = 5)
    end
    axislegend(ax1; position = :lb, labelsize = 8)

    # ---- (b) Full sweep envelope ----
    ax2 = Axis(fig[1, 2], xlabel="degree n", ylabel="envelope of |a_n|",
               yscale=log10,
               title="(b) full ell sweep at N=64, envelope only",
               limits=(nothing, (1e-17, 10.0)))
    for (i, ell) in enumerate(L_full)
        a = abs.(tbn_coeffs(N, ell))
        env = envelope_max(a, 3)
        lines!(ax2, 0:N, env .+ 1e-18; color = cols_full[i], linewidth = 0.9,
               label = "ell = $ell")
    end
    axislegend(ax2; position = :lb, labelsize = 8, nbanks = 2)

    # ---- (c) Tail sum vs ell at three resolutions ----
    Ls = [0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 16.0, 24.0]
    ax3 = Axis(fig[2, 1], xlabel="ell", ylabel="sum_{n>N/2} |a_n| (tail size)",
               xscale=log10, yscale=log10,
               title="(c) Valley of good ell broadens with N")
    for (N_ref, col, label) in [(24, CORAL, "N=24"), (48, TEAL, "N=48"), (96, NAVY, "N=96")]
        errs = [sum(abs.(tbn_coeffs(N_ref, ell)[(div(N_ref, 2) + 1):end])) for ell in Ls]
        scatterlines!(ax3, Ls, errs; color=col, label=label)
    end
    axislegend(ax3; position = :rb)

    # ---- (d) Valley centre and width vs N ----
    Ns_extra = [24, 32, 48, 64, 96, 128]
    centres = Float64[]; widths = Float64[]
    for N_ref in Ns_extra
        errs_N = [sum(abs.(tbn_coeffs(N_ref, ell)[(div(N_ref, 2) + 1):end])) for ell in Ls]
        emin = minimum(errs_N)
        good_mask = errs_N .<= 3.0 * emin
        push!(centres, exp(sum(log.(Ls[good_mask])) / count(good_mask)))
        push!(widths,  maximum(Ls[good_mask]) / minimum(Ls[good_mask]))
    end
    ax4 = Axis(fig[2, 2], xlabel="N", ylabel="value (log-log)",
               xscale=log10, yscale=log10,
               title="(d) broad-valley centre and width vs N")
    scatterlines!(ax4, Ns_extra, centres; color = NAVY,
                  markercolor = :white, strokecolor = NAVY, strokewidth = 1.0,
                  markersize = 6, linewidth = 1.1, label = "valley centre")
    scatterlines!(ax4, Ns_extra, widths; color = TEAL,
                  marker = :rect, markercolor = :white, strokecolor = TEAL,
                  strokewidth = 1.0, markersize = 5, linewidth = 1.0,
                  label = "valley width factor")
    axislegend(ax4; position = :rt)

    save(joinpath(outdir, "ell_diagnostic.pdf"), fig)
    save(joinpath(outdir, "ell_diagnostic.png"), fig)
    @printf("[20.9-julia] saved figure\n")
end

if abspath(PROGRAM_FILE) == @__FILE__
    run()
end
