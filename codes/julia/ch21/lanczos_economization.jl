# lanczos_economization.jl
# Chapter 21, Etude 21.6: Lanczos economization for an expensive determinant.
#
# Determinant D(lambda) = det(A(lambda)) of an M-by-M parameter-dependent
# matrix.  Compare:
#   (A) naive bracketing on a fine grid + 30 bisection refines per root
#   (B) Lanczos economization: K = 17 Chebyshev samples + colleague-matrix root
#
# Author: Dr. Denys Dutykh (Khalifa University of Science and Technology)

using FFTW
using LinearAlgebra
using Random
include(joinpath(@__DIR__, "tricks_common.jl"))
apply_theme!()
const OUT = output_dir(@__DIR__)


const M_SIZE = 20
const SHIFTS = Float64.(collect(1:M_SIZE)) .^ 2

function _make_T2()
    rng = MersenneTwister(0)
    A = randn(rng, M_SIZE, M_SIZE)
    F = qr(A)
    Q = Matrix(F.Q)
    diag_v = 0.05 .* randn(rng, M_SIZE)
    return Q * Diagonal(diag_v) * Q'
end
const _T2 = _make_T2()


function expensive_D(lam::Real)
    A = _T2 + Diagonal(lam .- SHIFTS)
    return det(A)
end


function roots_chebyshev_surrogate(K::Int, a::Real, b::Real)
    j = collect(0:K-1)
    t = cos.(pi .* j ./ (K - 1))
    lam = 0.5 * (a + b) .+ 0.5 * (b - a) .* t
    D_samples = [expensive_D(L) for L in lam]
    ext = vcat(D_samples, reverse(D_samples[2:K-1]))
    A = real.(fft(ext)) ./ (K - 1)
    A[1] *= 0.5; A[K] *= 0.5
    coeffs = A[1:K]

    n = K - 1
    C = zeros(n, n)
    for i in 1:n-1
        C[i, i+1] = 0.5
        C[i+1, i] = 0.5
    end
    C[1, 2] = 1.0
    last_row = -coeffs[1:n] ./ (2.0 * coeffs[n+1])
    last_row[n-1] += 0.5
    C[n, :] = last_row

    eigs = eigvals(C)
    real_in_window = sort([real(z) for z in eigs
                           if abs(imag(z)) < 1e-6 && -1.0 <= real(z) <= 1.0])
    roots = 0.5 * (a + b) .+ 0.5 * (b - a) .* real_in_window
    return roots, lam, D_samples
end


function roots_naive_brackets(K_dense::Int, a::Real, b::Real)
    lam = collect(range(a, b; length = K_dense))
    D = [expensive_D(L) for L in lam]
    roots = Float64[]
    for k in 1:length(lam)-1
        if D[k] * D[k+1] < 0
            l, r = lam[k], lam[k+1]
            fl = D[k]
            for _ in 1:30
                mid = 0.5 * (l + r)
                fm = expensive_D(mid)
                if fl * fm <= 0
                    r = mid
                else
                    l = mid; fl = fm
                end
            end
            push!(roots, 0.5 * (l + r))
        end
    end
    return roots, lam, D
end


function make_figure(; dump_path = nothing)
    a, b = 0.5, 30.0

    K_dense = 60
    roots_A, lam_A, D_A = roots_naive_brackets(K_dense, a, b)
    n_LUs_A = K_dense + 30 * length(roots_A)

    K = 17
    roots_B, lam_B, D_B = roots_chebyshev_surrogate(K, a, b)
    n_LUs_B = K

    seeds = SHIFTS[(SHIFTS .>= a) .& (SHIFTS .<= b)]
    exact_roots = Float64[]
    for seed in seeds
        x = seed
        for _ in 1:80
            d  = expensive_D(x)
            dp = (expensive_D(x + 1e-7) - expensive_D(x - 1e-7)) / 2e-7
            abs(dp) < 1e-30 && break
            x -= d / dp
        end
        push!(exact_roots, x)
    end
    sort!(exact_roots)

    err_A = [minimum(abs.(r .- exact_roots)) for r in sort(roots_A)]
    err_B = [minimum(abs.(r .- exact_roots)) for r in sort(roots_B)]

    fig = Figure(size = (1150, 440))

    ax1 = Axis(fig[1, 1]; xlabel = "lambda", ylabel = "D(lambda)",
               title = "the determinant and the two strategies' detected roots",
               limits = ((a - 1, b + 1), nothing))
    lam_dense = collect(range(a, b; length = 600))
    D_dense = [expensive_D(L) for L in lam_dense]
    lines!(ax1, lam_dense, D_dense; color = NAVY, linewidth = 1.0,
           label = "D(lambda) = det A(lambda)")
    hlines!(ax1, [0.0]; color = :gray, linewidth = 0.4, alpha = 0.5)
    scatter!(ax1, roots_A, zeros(length(roots_A));
             color = CORAL, marker = :xcross, markersize = 14,
             label = "naive: $(length(roots_A)) roots ($n_LUs_A LUs)")
    scatter!(ax1, roots_B, zeros(length(roots_B)) .+ 0.05 * maximum(abs.(D_dense));
             color = :white, strokecolor = TEAL, strokewidth = 1.4,
             markersize = 12,
             label = "Lanczos: $(length(roots_B)) roots ($n_LUs_B LUs)")
    scatter!(ax1, lam_B, D_B; color = TEAL, markersize = 6,
             label = "Lanczos sample points (K = $K)")
    axislegend(ax1, position = :rt, labelsize = 8)

    ax2 = Axis(fig[1, 2]; xlabel = "method", ylabel = "# expensive LU calls",
               title = "cost: $n_LUs_B vs $n_LUs_A LU factorisations",
               xticks = ([1, 2], ["naive\n+ 30 bisects", "Lanczos\nK=$K"]))
    barplot!(ax2, [1, 2], [n_LUs_A, n_LUs_B];
             color = [CORAL, TEAL], strokecolor = :black, strokewidth = 0.6)
    text!(ax2, [1], [n_LUs_A + 5]; text = "$n_LUs_A LUs", align = (:center, :bottom),
          fontsize = 9)
    text!(ax2, [2], [n_LUs_B + 5]; text = "$n_LUs_B LUs", align = (:center, :bottom),
          fontsize = 9)
    text!(ax2, [1], [n_LUs_A * 0.5];
          text = "max root err\n$(round(maximum(err_A); sigdigits = 2))",
          align = (:center, :center), fontsize = 8.5, color = :white)
    text!(ax2, [2], [n_LUs_B * 0.5];
          text = "max root err\n$(round(maximum(err_B); sigdigits = 2))",
          align = (:center, :center), fontsize = 8.5, color = :white)
    ylims!(ax2, 0, max(n_LUs_A, n_LUs_B) * 1.15)

    save(joinpath(OUT, "lanczos_economization.pdf"), fig)
    save(joinpath(OUT, "lanczos_economization.png"), fig, px_per_unit = 4)

    println("[Etude 21.6]  Lanczos economization for an expensive determinant")
    println("  exact roots in [$a, $b]: $exact_roots")
    println("  naive  : $(length(roots_A)) roots, $n_LUs_A LUs, max err = $(maximum(err_A))")
    println("  Lanczos: $(length(roots_B)) roots, $n_LUs_B LUs, max err = $(maximum(err_B))")
    println("  figure: ", joinpath(OUT, "lanczos_economization.pdf"))

    if dump_path !== nothing
        write_json(dump_path, Dict("exact_roots" => exact_roots,
            "roots_naive" => sort(roots_A), "roots_lanczos" => sort(roots_B),
            "n_LUs_naive" => n_LUs_A, "n_LUs_lanczos" => n_LUs_B,
            "max_err_naive" => maximum(err_A),
            "max_err_lanczos" => maximum(err_B)))
    end
end


function parse_args()
    dp = nothing; i = 1
    while i <= length(ARGS)
        if ARGS[i] == "--dump" && i < length(ARGS)
            dp = ARGS[i+1]; i += 2
        else; i += 1; end
    end
    return dp
end


if abspath(PROGRAM_FILE) == @__FILE__
    make_figure(; dump_path = parse_args())
end
