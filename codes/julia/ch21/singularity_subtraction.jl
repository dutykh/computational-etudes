# singularity_subtraction.jl
# Chapter 21, Etude 21.3: singularity subtraction (1-D surrogate).
#
# f(x) = e^x + 0.1 (1+x)^(2/3) on [-1, 1].  Subtract the corner-like
# singular part; what remains is e^x, which has geometric Chebyshev
# convergence.
#
# Author: Dr. Denys Dutykh (Khalifa University of Science and Technology)

using FFTW
include(joinpath(@__DIR__, "tricks_common.jl"))
apply_theme!()
const OUT = output_dir(@__DIR__)

const CSING = 0.1
f_full(x)     = exp(x) + CSING * (1.0 + x) ^ (2.0 / 3.0)
f_singular(x) = CSING * (1.0 + x) ^ (2.0 / 3.0)
f_smooth(x)   = exp(x)


cgl(N) = cos.(pi .* (0:N) ./ N)


function cheb_coeffs(v::AbstractVector)
    N = length(v) - 1
    ext = vcat(v, reverse(v[2:N]))
    A = real.(fft(ext)) ./ N
    A[1] *= 0.5; A[N+1] *= 0.5
    return A[1:N+1]
end


function reconstruct(a::AbstractVector, t_eval::AbstractVector)
    T0 = ones(length(t_eval))
    if length(a) == 1; return a[1] .* T0; end
    T1 = collect(t_eval)
    val = a[1] .* T0 .+ a[2] .* T1
    for n in 2:length(a)-1
        Tk = 2.0 .* t_eval .* T1 .- T0
        val .+= a[n+1] .* Tk
        T0, T1 = T1, Tk
    end
    return val
end


function make_figure(; dump_path = nothing)
    N_show = 80
    x_show = cgl(N_show)
    a_naive = cheb_coeffs(f_full.(x_show))
    a_trick = cheb_coeffs(f_smooth.(x_show))

    Ns = [8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256]
    err_naive = Float64[]
    err_trick = Float64[]
    x_eval = collect(range(-1.0, 1.0; length = 5001))
    for N in Ns
        xg = cgl(N)
        a_n = cheb_coeffs(f_full.(xg))
        a_t = cheb_coeffs(f_smooth.(xg))
        approx_n = reconstruct(a_n, x_eval)
        approx_t = reconstruct(a_t, x_eval) .+ f_singular.(x_eval)
        push!(err_naive, maximum(abs.(approx_n .- f_full.(x_eval))))
        push!(err_trick, maximum(abs.(approx_t .- f_full.(x_eval))))
    end

    fig = Figure(size = (1150, 440))

    ax1 = Axis(fig[1, 1]; xlabel = "Chebyshev degree n", ylabel = "|a_n|",
               yscale = log10,
               title = "coefficient decay at N = $N_show",
               limits = (nothing, (1e-18, 5.0)))
    n_axis = collect(0:N_show)
    scatterlines!(ax1, n_axis, abs.(a_naive) .+ 1e-20; color = CORAL,
        markercolor = :white, strokecolor = CORAL, strokewidth = 1.0,
        markersize = 4, linewidth = 0.8,
        label = "naive: |a_n| for f(x)")
    scatterlines!(ax1, n_axis, abs.(a_trick) .+ 1e-20; color = TEAL,
        marker = :utriangle, markercolor = :white, strokecolor = TEAL,
        strokewidth = 1.0, markersize = 4, linewidth = 0.8,
        label = "trick: |a_n| for f - c(1+x)^(2/3)")
    axislegend(ax1, position = :rt, labelsize = 9)

    ax2 = Axis(fig[1, 2]; xlabel = "N", ylabel = "max pointwise error",
               yscale = log10, xscale = log10,
               title = "subtraction reaches machine eps at modest N",
               limits = (nothing, (1e-17, 1.0)))
    scatterlines!(ax2, Ns, err_naive; color = CORAL, markercolor = :white,
        strokecolor = CORAL, strokewidth = 1.0, markersize = 6, linewidth = 1.0,
        label = "naive: max |f - f_N|")
    scatterlines!(ax2, Ns, err_trick; color = TEAL, marker = :utriangle,
        markercolor = :white, strokecolor = TEAL, strokewidth = 1.0,
        markersize = 6, linewidth = 1.0,
        label = "trick: max |f - (smooth_N + singular)|")
    hlines!(ax2, [1e-15]; color = :gray, linewidth = 0.4, alpha = 0.5)
    axislegend(ax2, position = :lb, labelsize = 9)

    save(joinpath(OUT, "singularity_subtraction.pdf"), fig)
    save(joinpath(OUT, "singularity_subtraction.png"), fig, px_per_unit = 4)

    println("[Etude 21.3]  singularity subtraction (1-D surrogate)")
    for (k, N) in enumerate(Ns)
        println("  N = $N:  naive err = $(err_naive[k]),  trick err = $(err_trick[k])")
    end
    println("  figure: ", joinpath(OUT, "singularity_subtraction.pdf"))

    if dump_path !== nothing
        write_json(dump_path, Dict("Ns" => Ns,
            "err_naive" => err_naive, "err_trick" => err_trick))
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
