# hilbert_fourier.jl
# Chapter 21, Etude 21.8: Hilbert transform via Fourier multiplier.
#
# Test function: f(x) = exp(cos x), periodic on [-pi, pi].
# Fourier series: exp(cos x) = I_0(1) + 2 sum_k I_k(1) cos(k x).
# Hilbert transform: H{f}(y) = 2 sum_k I_k(1) sin(k y).
# Coefficients I_k(1) decay factorially; convergence is super-geometric.
#
# Author: Dr. Denys Dutykh (Khalifa University of Science and Technology)

using FFTW
using SpecialFunctions       # for besseli (= I_k)
include(joinpath(@__DIR__, "tricks_common.jl"))
apply_theme!()
const OUT = output_dir(@__DIR__)


function hilbert_via_fft(f_values::AbstractVector)
    N = length(f_values)
    F = fft(f_values)
    k = fftfreq(N, 1.0) .* N    # integer wavenumbers 0, 1, ..., -1
    multiplier = -1im .* sign.(k)
    multiplier[1] = 0.0
    return real.(ifft(multiplier .* F))
end


f_test(x) = exp(cos(x))


function hilbert_exact(y, K_max::Int = 40)
    val = zero(y)
    for k in 1:K_max
        val += 2.0 * besseli(k, 1.0) * sin(k * y)
    end
    return val
end


function make_figure(; dump_path = nothing)
    Ns = [4, 6, 8, 10, 12, 14, 16, 20, 24, 32]
    err = Float64[]
    for N in Ns
        x = collect(-pi .+ 2 * pi .* (0:N-1) ./ N)
        Hf = hilbert_via_fft(f_test.(x))
        Hf_ex = [hilbert_exact(xx) for xx in x]
        push!(err, maximum(abs.(Hf .- Hf_ex)))
    end

    N_show = 32
    x_show = collect(-pi .+ 2 * pi .* (0:N_show-1) ./ N_show)
    Hf_show = hilbert_via_fft(f_test.(x_show))
    x_dense = collect(range(-pi, pi; length = 401))
    f_dense = f_test.(x_dense)
    Hf_ex_dense = [hilbert_exact(xx) for xx in x_dense]

    fig = Figure(size = (1150, 440))

    ax1 = Axis(fig[1, 1]; xlabel = "x (or y)", ylabel = "amplitude",
               title = "f(x) = exp(cos x) and its Hilbert transform on the circle",
               limits = ((-pi - 0.1, pi + 0.1), nothing))
    lines!(ax1, x_dense, f_dense; color = NAVY, linewidth = 1.4,
           label = "f(x) = exp(cos x)")
    lines!(ax1, x_dense, Hf_ex_dense; color = TEAL, linewidth = 1.4,
           label = "H{f}(y) exact")
    scatter!(ax1, x_show, Hf_show; color = :white, strokecolor = CORAL,
             strokewidth = 1.0, markersize = 5,
             label = "H{f}_N(y_j) via FFT, N = $N_show")
    hlines!(ax1, [0.0]; color = :gray, linewidth = 0.4, alpha = 0.5)
    axislegend(ax1, position = :lb, labelsize = 9)

    ax2 = Axis(fig[1, 2]; xlabel = "N", ylabel = "max error",
               yscale = log10,
               title = "super-geometric convergence on a periodic test",
               limits = ((0, 35), (1e-17, 5.0)))
    scatterlines!(ax2, Ns, err .+ 1e-18; color = NAVY, markercolor = :white,
        strokecolor = NAVY, strokewidth = 1.0, markersize = 6, linewidth = 1.0,
        label = "max |H{f}_N - H{f}|")
    bound = 2.0 .* [besseli(N / 2, 1.0) for N in Ns]
    lines!(ax2, Ns, bound; color = CORAL, linestyle = :dash, linewidth = 0.8,
           label = "~ 2 I_(N/2)(1) (Bessel-tail bound)")
    hlines!(ax2, [1e-15]; color = :gray, linewidth = 0.4, alpha = 0.5)
    axislegend(ax2, position = :rt, labelsize = 9)

    save(joinpath(OUT, "hilbert_fourier.pdf"), fig)
    save(joinpath(OUT, "hilbert_fourier.png"), fig, px_per_unit = 4)

    println("[Etude 21.8]  Hilbert transform via Fourier multiplier")
    for (k, N) in enumerate(Ns)
        println("  N = $N: max err = $(err[k])")
    end
    println("  figure: ", joinpath(OUT, "hilbert_fourier.pdf"))

    if dump_path !== nothing
        write_json(dump_path, Dict("Ns" => Ns, "err" => err))
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
