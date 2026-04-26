# rescue_naive_vs_tailored.jl
# Chapter 21, Etude 21.1: a rescue story in miniature.
#
# f(y) = sech(y) on [-L, L] by Chebyshev collocation, N = 48.
# Same N, two L values:
#   naive    L = 30 -- conservative 'large enough' choice; wastes 88% of points
#   tailored L =  8 -- just enough; ~1000x more accurate at the same cost
#
# Author: Dr. Denys Dutykh (Khalifa University of Science and Technology)
# Part of "Computational Etudes: A Spectral Approach"

using FFTW

include(joinpath(@__DIR__, "tricks_common.jl"))

apply_theme!()

const OUT = output_dir(@__DIR__)


function cgl(N::Int)
    return cos.(pi .* (0:N) ./ N)
end


function cheb_coeffs(v::AbstractVector)
    N = length(v) - 1
    ext = vcat(v, reverse(v[2:N]))   # mirror without endpoints, length 2N
    A = real.(fft(ext)) ./ N
    A[1] *= 0.5
    A[N+1] *= 0.5
    return A[1:N+1]
end


function clenshaw(a::AbstractVector, t::AbstractVector)
    N = length(a) - 1
    T0 = ones(length(t))
    if N == 0; return a[1] .* T0; end
    T1 = collect(t)
    val = a[1] .* T0 .+ a[2] .* T1
    for n in 2:N
        Tk = 2.0 .* t .* T1 .- T0
        val .+= a[n+1] .* Tk
        T0, T1 = T1, Tk
    end
    return val
end


function trunc_cheb_err(N::Int, L::Real, y_ref::AbstractVector)
    t_grid = cgl(N)
    y_grid = L .* t_grid
    a = cheb_coeffs(sech.(y_grid))
    t_ref = y_ref ./ L
    inside = abs.(t_ref) .<= 1.0
    approx = fill(NaN, length(y_ref))
    approx[inside] .= clenshaw(a, t_ref[inside])
    err = abs.(approx .- sech.(y_ref))
    err[.!inside] .= NaN
    return a, err
end


function make_figure(; dump_path = nothing)
    N = 48
    L_naive = 30.0
    L_tail  = 8.0

    y_ref = collect(range(-35.0, 35.0; length = 6001))
    a_naive, err_naive = trunc_cheb_err(N, L_naive, y_ref)
    a_tail,  err_tail  = trunc_cheb_err(N, L_tail,  y_ref)

    inside_naive = abs.(y_ref) .<= L_naive
    inside_tail  = abs.(y_ref) .<= L_tail
    e_naive = maximum(filter(!isnan, err_naive[inside_naive]))
    e_tail  = maximum(filter(!isnan, err_tail[inside_tail]))

    fig = Figure(size = (1150, 440))

    ax1 = Axis(fig[1, 1]; xlabel = "y", ylabel = "f(y)",
               title = "same f, same N = $N, two truncation lengths L",
               limits = ((-22, 22), (-0.06, 1.12)))
    lines!(ax1, y_ref, sech.(y_ref); color = NAVY, linewidth = 1.4,
           label = "f(y) = sech y")
    y_naive_grid = L_naive .* cgl(N)
    scatter!(ax1, y_naive_grid, sech.(y_naive_grid);
             color = :white, strokecolor = CORAL, strokewidth = 1.0,
             markersize = 6, label = "naive (L = $(Int(L_naive)))")
    y_tail_grid = L_tail .* cgl(N)
    scatter!(ax1, y_tail_grid, sech.(y_tail_grid);
             color = :white, strokecolor = TEAL, strokewidth = 1.0,
             markersize = 8, marker = :utriangle,
             label = "tailored (L = $(Int(L_tail)))")
    vlines!(ax1, [0.0]; color = :gray, linewidth = 0.4, alpha = 0.4)
    axislegend(ax1, position = :lt, labelsize = 9)

    ax2 = Axis(fig[1, 2]; xlabel = "y", ylabel = "|f(y) - f_N(y)|",
               yscale = log10,
               title = "pointwise error: a single tuned L buys decimals",
               limits = ((-22, 22), (1e-16, 1e0)))
    lines!(ax2, y_ref[inside_naive], err_naive[inside_naive] .+ 1e-18;
           color = CORAL, linewidth = 1.0,
           label = "naive (L=$(Int(L_naive))): max err ≈ $(round(e_naive; sigdigits=2))")
    lines!(ax2, y_ref[inside_tail], err_tail[inside_tail] .+ 1e-18;
           color = TEAL, linewidth = 1.0,
           label = "tailored (L=$(Int(L_tail))): max err ≈ $(round(e_tail; sigdigits=2))")
    hlines!(ax2, [1e-15]; color = :gray, linewidth = 0.4, alpha = 0.5)
    axislegend(ax2, position = :cb, labelsize = 9)

    save(joinpath(OUT, "rescue_naive_vs_tailored.pdf"), fig)
    save(joinpath(OUT, "rescue_naive_vs_tailored.png"), fig, px_per_unit = 4)

    println("[Etude 21.1]  rescue story")
    println("  N = $N, f(y) = sech(y)")
    println("  naive    L = $L_naive: max err = $(round(e_naive; sigdigits=4))")
    println("  tailored L = $L_tail: max err = $(round(e_tail; sigdigits=4))")
    println("  ratio (tailored / naive) = $(round(e_tail / e_naive; sigdigits=3))")
    println("  figure: ", joinpath(OUT, "rescue_naive_vs_tailored.pdf"))

    if dump_path !== nothing
        write_json(dump_path, Dict("N" => N,
            "L_naive" => L_naive, "L_tailored" => L_tail,
            "err_naive" => e_naive, "err_tailored" => e_tail))
    end
end


function parse_args()
    dp = nothing; i = 1
    while i <= length(ARGS)
        if ARGS[i] == "--dump" && i < length(ARGS)
            dp = ARGS[i+1]; i += 2
        else
            i += 1
        end
    end
    return dp
end


if abspath(PROGRAM_FILE) == @__FILE__
    make_figure(; dump_path = parse_args())
end
