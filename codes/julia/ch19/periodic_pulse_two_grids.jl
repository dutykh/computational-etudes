# periodic_pulse_two_grids.jl
# Chapter 19: Coordinate Transformations
# Computational Etude 19.1: A periodic pulse on two grids.
#
# Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
# Part of "Computational Etudes: A Spectral Approach"

using CairoMakie
using FFTW
using Printf

const KAPPA = 80.0
const L_MAP = 0.3

target(y) = exp(-KAPPA * (1.0 - cos(y)))

arctan_tan_map(x, ell) = 2.0 * atan(ell * tan(x / 2))
arctan_tan_inv(y, ell) = 2.0 * atan(tan(y / 2) / ell)

function fourier_interp(y_nodes, f_nodes, y_eval)
    N = length(y_nodes)
    coeffs = fft(f_nodes) ./ N
    ks = [0:div(N, 2)-1; -div(N, 2):-1]
    vals = zeros(length(y_eval))
    for m in 1:N
        vals .+= real.(coeffs[m] .* exp.(im .* ks[m] .* (y_eval .+ pi)))
    end
    return vals
end

function mapped_interp(x_nodes, f_nodes, y_eval, ell)
    x_eval = arctan_tan_inv.(y_eval, ell)
    return fourier_interp(x_nodes, f_nodes, x_eval)
end

function run()
    outdir = joinpath(@__DIR__, "..", "..", "..",
                      "textbook", "figures", "ch19", "julia")
    mkpath(outdir)

    Ns = [8, 12, 16, 24, 32, 48, 64, 96, 128]
    y_eval = collect(range(-pi + 1e-9, pi - 1e-9, length=4097))
    truth = target.(y_eval)
    err_U = zeros(length(Ns)); err_M = zeros(length(Ns))
    for (i, N) in enumerate(Ns)
        yU = [-pi + 2pi * k / N for k in 0:N-1]
        fU = target.(yU)
        err_U[i] = maximum(abs.(fourier_interp(yU, fU, y_eval) .- truth))

        x = [-pi + 2pi * k / N for k in 0:N-1]
        yM = arctan_tan_map.(x, L_MAP)
        fM = target.(yM)
        err_M[i] = maximum(abs.(mapped_interp(x, fM, y_eval, L_MAP) .- truth))
    end

    NAVY = colorant"#142D6E"; CORAL = colorant"#E74C3C"; TEAL = colorant"#16A085"
    fig = Figure(size=(1100, 340))

    Nshow = 32
    yU = [-pi + 2pi * k / Nshow for k in 0:Nshow-1]
    x = yU
    yM = arctan_tan_map.(x, L_MAP)
    y_line = collect(range(-pi, pi, length=401))

    ax1 = Axis(fig[1, 1], xlabel="y", title="(a) Pulse + two grids",
               limits=((-pi, pi), (-0.28, 1.1)))
    lines!(ax1, y_line, target.(y_line); color=NAVY, linewidth=1.2)
    scatter!(ax1, yU, fill(-0.08, length(yU)); color=CORAL, marker=:xcross, markersize=9)
    scatter!(ax1, yM, fill(-0.18, length(yM)); color=TEAL, marker=:circle, markersize=6)

    ax2 = Axis(fig[1, 2], xlabel="N", ylabel="max error",
               title="(b) Convergence", yscale=log10)
    scatterlines!(ax2, Ns, max.(err_U, 1e-18); color=CORAL, label="uniform")
    scatterlines!(ax2, Ns, max.(err_M, 1e-18); color=TEAL, label="arctan/tan")
    axislegend(ax2, position=:rt)

    N = 96
    yU = [-pi + 2pi * k / N for k in 0:N-1]
    cU = abs.(fftshift(fft(target.(yU)) ./ N))
    x = yU
    yM = arctan_tan_map.(x, L_MAP)
    cM = abs.(fftshift(fft(target.(yM)) ./ N))
    kfull = fftshift([0:div(N, 2)-1; -div(N, 2):-1])
    ax3 = Axis(fig[1, 3], xlabel="k", ylabel="|c_k|",
               title="(c) Coefficient decay, N=96", yscale=log10,
               limits=(nothing, (1e-16, 1)))
    scatter!(ax3, kfull[kfull .>= 0], cU[kfull .>= 0] .+ 1e-18;
             color=CORAL, marker=:xcross, markersize=7)
    scatter!(ax3, kfull[kfull .>= 0], cM[kfull .>= 0] .+ 1e-18;
             color=TEAL, marker=:circle, markersize=5)

    save(joinpath(outdir, "periodic_pulse_two_grids.pdf"), fig)
    save(joinpath(outdir, "periodic_pulse_two_grids.png"), fig)
    @printf("[19.1-julia] saved figure to %s\n",
            joinpath(outdir, "periodic_pulse_two_grids.pdf"))
end

if abspath(PROGRAM_FILE) == @__FILE__
    run()
end
