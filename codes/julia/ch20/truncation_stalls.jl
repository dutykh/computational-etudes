# truncation_stalls.jl
# Chapter 20: Spectral Methods on Unbounded Intervals
# Computational Etude 20.1: When spectral convergence stalls.
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

function cheb_trunc_err(N, L)
    _, x = cheb_matrix(N)
    y_nodes = L .* x
    fv = target.(y_nodes)
    a = dct1_coeffs(fv)
    y_fine = collect(range(-20.0, 20.0, length=4001))
    in_window = abs.(y_fine) .<= L
    approx = zeros(length(y_fine))
    approx[in_window] .= cheb_eval(a, y_fine[in_window] ./ L, N)
    return maximum(abs.(approx .- target.(y_fine)))
end

function run()
    outdir = joinpath(@__DIR__, "..", "..", "..",
                      "textbook", "figures", "ch20", "julia")
    mkpath(outdir)

    L_A = 6.0
    Ns_A = [8, 12, 16, 24, 32, 48, 64, 96, 128]
    err_A = [cheb_trunc_err(N, L_A) for N in Ns_A]

    N_B = 32
    Ls_B = [2, 3, 4, 5, 6, 8, 10, 12, 16, 20]
    err_B = [cheb_trunc_err(N_B, L) for L in Ls_B]

    pairs = [(8, 3), (12, 4), (16, 5), (24, 6), (32, 8), (48, 10), (64, 12), (96, 14)]
    err_C = [cheb_trunc_err(N, L) for (N, L) in pairs]
    Ns_C = [p[1] for p in pairs]

    NAVY = colorant"#142D6E"; CORAL = colorant"#E74C3C"; TEAL = colorant"#16A085"

    # NEW (panel d): full 2D error surface E(N, L)
    Ns_grid = [8, 12, 16, 24, 32, 48, 64, 96, 128]
    Ls_grid = [2, 3, 4, 5, 6, 8, 10, 12, 16, 20]
    err_grid = zeros(length(Ls_grid), length(Ns_grid))
    for (i, L) in enumerate(Ls_grid)
        for (j, N) in enumerate(Ns_grid)
            err_grid[i, j] = cheb_trunc_err(N, L)
        end
    end
    log_err = log10.(err_grid .+ 1e-18)

    fig = Figure(size=(1100, 800))

    # ---- (a) plateau ----
    ax1 = Axis(fig[1, 1], xlabel="N", ylabel="max error",
               yscale=log10, title="(a) Fix L=$(Int(L_A)), vary N: plateau")
    scatterlines!(ax1, Ns_A, err_A; color=CORAL,
                  label="error at L = $(Int(L_A))")
    hlines!(ax1, [exp(-L_A)]; linestyle=:dash, color=NAVY,
            label="e^{-L} ≈ $(round(exp(-L_A); sigdigits=2))")
    axislegend(ax1; position=:rb)

    # ---- (b) sweet spot ----
    ax2 = Axis(fig[1, 2], xlabel="L", ylabel="max error",
               yscale=log10, title="(b) Fix N=$N_B, vary L: sweet spot")
    scatterlines!(ax2, Ls_B, err_B; color=TEAL,
                  markercolor=:white, strokecolor=TEAL, strokewidth=1.0)

    # ---- (c) subgeometric descent ----
    ax3 = Axis(fig[2, 1], xlabel="N  (L growing with N)", ylabel="max error",
               yscale=log10,
               title="(c) Grow both: subgeometric descent")
    scatterlines!(ax3, Ns_C, err_C; color=NAVY)

    # ---- (d) NEW: full error surface E(N, L) with three sweep overlays ----
    ax4 = Axis(fig[2, 2], xlabel="N", ylabel="L",
               title="(d) error surface E(N, L) and its three slices")
    hm = heatmap!(ax4, Ns_grid, Ls_grid, log_err'; colormap=:viridis)
    Colorbar(fig[2, 2, Right()], hm; label="log10 max error", width=12)
    lines!(ax4, Ns_A, fill(L_A, length(Ns_A)); color=CORAL, linewidth=1.6,
           label="(a) L = 6, vary N")
    scatter!(ax4, Ns_A, fill(L_A, length(Ns_A)); color=:white,
             strokecolor=CORAL, strokewidth=1.0, markersize=5)
    lines!(ax4, fill(N_B, length(Ls_B)), Ls_B; color=TEAL, linewidth=1.6,
           label="(b) N = $N_B, vary L")
    scatter!(ax4, fill(N_B, length(Ls_B)), Ls_B; color=:white,
             strokecolor=TEAL, strokewidth=1.0, marker=:rect, markersize=5)
    Ls_C = [p[2] for p in pairs]
    lines!(ax4, Ns_C, Ls_C; color=:white, linewidth=2.0,
           label="(c) grow both jointly")
    scatter!(ax4, Ns_C, Ls_C; color=:white, marker=:utriangle, markersize=7)
    axislegend(ax4; position=:rt, labelsize=8, framevisible=true,
               framecolor=:black, backgroundcolor=(:black, 0.6))

    save(joinpath(outdir, "truncation_stalls.pdf"), fig)
    save(joinpath(outdir, "truncation_stalls.png"), fig)
    @printf("[20.1-julia] saved figure\n")
end

if abspath(PROGRAM_FILE) == @__FILE__
    run()
end
