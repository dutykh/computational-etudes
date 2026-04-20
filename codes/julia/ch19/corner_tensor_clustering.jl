# corner_tensor_clustering.jl
# Chapter 19: Coordinate Transformations
# Computational Etude 19.6: A square Poisson problem with corner stress.
#
# Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)

using CairoMakie
using LinearAlgebra
using Printf
using Interpolations

include(joinpath(@__DIR__, "..", "ch07", "cheb_matrix.jl"))

function solve_unmapped(N)
    D, x = cheb_matrix(N)
    D2 = D * D
    I_N = Matrix(I, N + 1, N + 1)
    nx = N + 1
    L = kron(D2, I_N) .+ kron(I_N, D2)
    idx = reshape(1:nx*nx, nx, nx)
    interior = fill(false, nx, nx); interior[2:N, 2:N] .= true
    idx_int = idx[interior]
    A = L[idx_int, idx_int]
    b = -ones(length(idx_int))
    u_int = A \ b
    u = zeros(nx * nx); u[idx_int] = u_int
    return x, reshape(u, nx, nx)
end

function solve_mapped(N, alpha)
    D, xi = cheb_matrix(N)
    th_a = tanh(alpha)
    x_phys = tanh.(alpha .* xi) ./ th_a
    fp  = alpha ./ (cosh.(alpha .* xi) .^ 2 .* th_a)
    fpp = -2 .* alpha^2 .* tanh.(alpha .* xi) ./ (cosh.(alpha .* xi) .^ 2 .* th_a)
    D1 = Diagonal(1 ./ fp) * D
    D2 = Diagonal(1 ./ fp .^ 2) * (D * D) - Diagonal(fpp ./ fp .^ 3) * D
    I_N = Matrix(I, N + 1, N + 1); nx = N + 1
    L = kron(D2, I_N) .+ kron(I_N, D2)
    idx = reshape(1:nx*nx, nx, nx)
    interior = fill(false, nx, nx); interior[2:N, 2:N] .= true
    idx_int = idx[interior]
    A = L[idx_int, idx_int]
    u_int = A \ (-ones(length(idx_int)))
    u = zeros(nx * nx); u[idx_int] = u_int
    return x_phys, reshape(u, nx, nx)
end

function compare_err(x, U, xref, Uref)
    function sort_grid(x, U)
        if x[1] > x[end]
            return reverse(x), reverse(reverse(U; dims=1); dims=2)
        else
            return x, U
        end
    end
    xs, Us = sort_grid(x, U)
    xrs, Urs = sort_grid(xref, Uref)
    itp_m = LinearInterpolation((xs, xs), Us; extrapolation_bc=0.0)
    itp_r = LinearInterpolation((xrs, xrs), Urs; extrapolation_bc=0.0)
    xf = collect(range(-0.99, 0.99, length=121))
    emax = 0.0
    for yv in xf, xv in xf
        emax = max(emax, abs(itp_m(xv, yv) - itp_r(xv, yv)))
    end
    return emax
end

function run()
    outdir = joinpath(@__DIR__, "..", "..", "..",
                      "textbook", "figures", "ch19", "julia")
    mkpath(outdir)

    xref, Uref = solve_unmapped(96)
    Ns = [12, 16, 20, 24, 32, 40, 48, 64]
    alphas = [1.0, 2.0, 3.0]
    err_u = zeros(length(Ns))
    err_m = zeros(length(alphas), length(Ns))
    for (i, N) in enumerate(Ns)
        xu, Uu = solve_unmapped(N)
        err_u[i] = compare_err(xu, Uu, xref, Uref)
        for (j, a) in enumerate(alphas)
            xm, Um = solve_mapped(N, a)
            err_m[j, i] = compare_err(xm, Um, xref, Uref)
        end
    end

    NAVY = colorant"#142D6E"; CORAL = colorant"#E74C3C"; TEAL = colorant"#16A085"
    PURPLE = colorant"#8E44AD"
    fig = Figure(size=(1300, 340))

    ax1 = Axis(fig[1, 1], xlabel="x", ylabel="y", title="(a) -Δu=1")
    xs, Us = solve_unmapped(32)
    idx = sortperm(xs)
    contourf!(ax1, xs[idx], xs[idx], Us[idx, idx]; levels=15)

    ax2 = Axis(fig[1, 2], xlabel="physical coordinate",
               title="(b) Grids at N=24", limits=((-1.05, 1.05), (0, 1)))
    _, xi24 = cheb_matrix(24)
    x_plain = xi24; x_tanh = tanh.(2.0 .* xi24) ./ tanh(2.0)
    for xg in x_plain; vlines!(ax2, xg; color=NAVY, ymin=0.0, ymax=0.45); end
    for xg in x_tanh;  vlines!(ax2, xg; color=CORAL, ymin=0.55, ymax=1.0); end

    ax3 = Axis(fig[1, 3], xlabel="N", ylabel="error", xscale=log10, yscale=log10,
               title="(c) Corner convergence")
    scatterlines!(ax3, Ns, err_u; color=NAVY, label="unmapped")
    colours = [CORAL, TEAL, PURPLE]
    for (j, a) in enumerate(alphas)
        scatterlines!(ax3, Ns, err_m[j, :]; color=colours[j], label="tanh α=$a")
    end
    axislegend(ax3; position=:rt)

    save(joinpath(outdir, "corner_tensor_clustering.pdf"), fig)
    save(joinpath(outdir, "corner_tensor_clustering.png"), fig)
    @printf("[19.6-julia] saved figure\n")
end

if abspath(PROGRAM_FILE) == @__FILE__
    run()
end
