# radiation_scattering.jl
# Chapter 21, Etude 21.5: 1-D Schroedinger scattering by a sech^2 potential
# via Boyd's radiation-basis trick.  Reproduces Boyd Table 19.1.
#
# Author: Dr. Denys Dutykh (Khalifa University of Science and Technology)

using LinearAlgebra
include(joinpath(@__DIR__, "tricks_common.jl"))
apply_theme!()
const OUT = output_dir(@__DIR__)

const V_DEPTH = 1.0
const ELL = 2.0

V(x)  = -V_DEPTH / cosh(x) ^ 2
H(x)  = 0.5 * (1.0 + tanh(x))
Hp(x) = 0.5 / cosh(x) ^ 2
Hpp(x) = -tanh(x) / cosh(x) ^ 2

reflection_exact(k) = (1.0 + cos(2 * pi * sqrt(V_DEPTH + 0.25))) /
                      (cosh(2 * pi * k) + cos(2 * pi * sqrt(V_DEPTH + 0.25)))


function collocation_grid(N::Int, ell::Real)
    i = collect(1:N)
    t = pi .* (2.0 .* i .- 1.0) ./ (2.0 * N)
    x = ell ./ tan.(t)
    return x, t
end


function basis_TB(N::Int, t::AbstractVector)
    n = collect(0:N-3)
    return [cos(n[j] * t[i]) for i in 1:N, j in 1:length(n)]
end


function basis_TB_double_prime(N::Int, t::AbstractVector, ell::Real)
    n = collect(0:N-3)
    P = zeros(N, length(n))
    for i in 1:N, j in 1:length(n)
        s, c = sin(t[i]), cos(t[i])
        nn = n[j]
        P[i, j] = -(nn / ell^2) * s^3 *
                  (nn * cos(nn * t[i]) * s + 2.0 * sin(nn * t[i]) * c)
    end
    return P
end


function solve_scattering(N::Int, ell::Real, k::Real)
    x, t = collocation_grid(N, ell)
    Phi   = zeros(N, N)
    Phipp = zeros(N, N)
    Phi[:, 1:N-2]   = basis_TB(N, t)
    Phipp[:, 1:N-2] = basis_TB_double_prime(N, t, ell)
    Phi[:, N-1] = H.(x) .* cos.(k .* x)
    Phipp[:, N-1] = (Hpp.(x) .* cos.(k .* x)
                     .- 2.0 .* k .* Hp.(x) .* sin.(k .* x)
                     .- k^2 .* H.(x) .* cos.(k .* x))
    Phi[:, N] = H.(x) .* sin.(k .* x)
    Phipp[:, N] = (Hpp.(x) .* sin.(k .* x)
                   .+ 2.0 .* k .* Hp.(x) .* cos.(k .* x)
                   .- k^2 .* H.(x) .* sin.(k .* x))

    M = Phipp .+ (k^2 .- V.(x)) .* Phi
    f = V.(x) .* cos.(k .* x)
    g = V.(x) .* sin.(k .* x)

    a = M \ f
    b = M \ g

    gamma1 = a[N-1] + 1.0
    gamma2 = a[N]
    sigma1 = b[N-1]
    sigma2 = b[N]   + 1.0

    A = [gamma1 + sigma2  sigma1 - gamma2;
         gamma2 - sigma1  gamma1 + sigma2]
    rhs = [sigma2 - gamma1, -sigma1 - gamma2]
    re_im = A \ rhs
    return complex(re_im[1], re_im[2]), sum(a[1:N-2])
end


function make_figure(; dump_path = nothing)
    N = 48
    k_axis = [0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4, 2.7, 3.0]
    R_num = Float64[]
    R_exact = Float64[]
    drifts = Float64[]
    for k in k_axis
        alpha, drift = solve_scattering(N, ELL, k)
        push!(R_num, abs(alpha)^2)
        push!(R_exact, reflection_exact(k))
        push!(drifts, abs(drift))
    end
    abs_err = abs.(R_num .- R_exact)

    fig = Figure(size = (1150, 440))

    ax1 = Axis(fig[1, 1]; xlabel = "wavenumber k", ylabel = "reflection coefficient R",
               yscale = log10,
               title = "sech^2 scattering: spectral R matches the closed form",
               limits = ((0.2, 3.1), (1e-10, 2.0)))
    lines!(ax1, k_axis, R_exact; color = NAVY, linewidth = 1.4,
           label = "R exact (Boyd Eq 19.31)")
    scatter!(ax1, k_axis, R_num; color = :white, strokecolor = CORAL,
             strokewidth = 1.0, markersize = 8,
             label = "R numerical, N = $N, ell = $(Int(ELL))")
    axislegend(ax1, position = :lb, labelsize = 10)

    ax2 = Axis(fig[1, 2]; xlabel = "wavenumber k", ylabel = "error magnitude",
               yscale = log10,
               title = "error stays well below the rapidly small R",
               limits = ((0.2, 3.1), nothing))
    scatterlines!(ax2, k_axis, abs_err .+ 1e-18; color = NAVY, marker = :rect,
        markercolor = :white, strokecolor = NAVY, strokewidth = 1.0,
        markersize = 6, linewidth = 1.0,
        label = "|R_num - R_exact|")
    scatterlines!(ax2, k_axis, drifts .+ 1e-18; color = TEAL, marker = :utriangle,
        markercolor = :white, strokecolor = TEAL, strokewidth = 1.0,
        markersize = 5, linewidth = 0.8,
        label = "|sum_n a_n| (asymptotic-constant drift)")
    axislegend(ax2, position = :lt, labelsize = 10)

    save(joinpath(OUT, "radiation_scattering.pdf"), fig)
    save(joinpath(OUT, "radiation_scattering.png"), fig, px_per_unit = 4)

    println("[Etude 21.5]  Schroedinger sech^2 scattering with radiation basis")
    println("  v = $V_DEPTH, ell = $ELL, N = $N")
    println("   k        R_numerical       R_exact         abs error")
    for (i, k) in enumerate(k_axis)
        println("  $k  $(R_num[i])   $(R_exact[i])   $(abs_err[i])")
    end
    println("  figure: ", joinpath(OUT, "radiation_scattering.pdf"))

    if dump_path !== nothing
        write_json(dump_path, Dict("N" => N, "ell" => ELL, "v" => V_DEPTH,
            "k_axis" => k_axis, "R_numerical" => R_num,
            "R_exact" => R_exact, "abs_err" => abs_err, "drifts" => drifts))
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
