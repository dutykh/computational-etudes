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

function pw_error_field(x, U, xref, Uref)
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
    xf = collect(range(-0.995, 0.995, length=121))
    err = zeros(length(xf), length(xf))
    for (j, yv) in enumerate(xf), (i, xv) in enumerate(xf)
        err[i, j] = abs(itp_m(xv, yv) - itp_r(xv, yv))
    end
    return xf, err
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

    # Pointwise error fields at N_err
    N_err = 24
    xu_err, Uu_err = solve_unmapped(N_err)
    xf_err, err_um = pw_error_field(xu_err, Uu_err, xref, Uref)
    xm1, Um1 = solve_mapped(N_err, 1.0)
    _, err_t1 = pw_error_field(xm1, Um1, xref, Uref)
    xm3, Um3 = solve_mapped(N_err, 3.0)
    _, err_t3 = pw_error_field(xm3, Um3, xref, Uref)
    err_floor = max(1e-8, min(
        minimum(err_um[err_um .> 0]),
        minimum(err_t1[err_t1 .> 0]),
        minimum(err_t3[err_t3 .> 0])))
    err_ceil = max(maximum(err_um), maximum(err_t1), maximum(err_t3))

    fig = Figure(size=(1350, 800))

    # (a) Solution contour at N=32
    ax1 = Axis(fig[1, 1]; xlabel="x", ylabel="y", aspect=DataAspect(),
               title="(a) -Δu = 1, Dirichlet")
    xs, Us = solve_unmapped(32)
    idx = sortperm(xs)
    contourf!(ax1, xs[idx], xs[idx], Us[idx, idx]; levels=15)

    # (b) 1-D grids at N=24
    ax2 = Axis(fig[1, 2]; xlabel="physical coordinate",
               title="(b) Grids at N = 24",
               limits=((-1.05, 1.05), (0, 1)))
    _, xi24 = cheb_matrix(24)
    x_plain = xi24; x_tanh = tanh.(2.0 .* xi24) ./ tanh(2.0)
    for xg in x_plain
        lines!(ax2, [xg, xg], [0.0, 0.45]; color=NAVY)
    end
    for xg in x_tanh
        lines!(ax2, [xg, xg], [0.55, 1.0]; color=CORAL)
    end
    text!(ax2, -0.95, 0.22; text="standard Chebyshev", color=NAVY, fontsize=10)
    text!(ax2, -0.95, 0.78; text="tanh clustered, α = 2", color=CORAL, fontsize=10)

    # (c) Corner convergence
    ax3 = Axis(fig[1, 3]; xlabel="N", ylabel="error",
               xscale=log10, yscale=log10, title="(c) Corner convergence")
    scatterlines!(ax3, Ns, err_u; color=NAVY, label="unmapped")
    colours = [CORAL, TEAL, PURPLE]
    for (j, a) in enumerate(alphas)
        scatterlines!(ax3, Ns, err_m[j, :]; color=colours[j], label="tanh α = $a")
    end
    axislegend(ax3; position=:rt, framevisible=false)

    # (d)/(e)/(f) Pointwise-error heatmaps, common log-colour scale
    cr_lo = log10(err_floor); cr_hi = log10(err_ceil + 1e-30)

    ax4 = Axis(fig[2, 1]; xlabel="x", ylabel="y", aspect=DataAspect(),
               title="(d) |u_N - u_ref|, unmapped, N = $N_err")
    heatmap!(ax4, xf_err, xf_err, log10.(err_um .+ err_floor);
             colormap=Reverse(:hot), colorrange=(cr_lo, cr_hi))

    ax5 = Axis(fig[2, 2]; xlabel="x", ylabel="y", aspect=DataAspect(),
               title="(e) tanh α = 1, N = $N_err")
    heatmap!(ax5, xf_err, xf_err, log10.(err_t1 .+ err_floor);
             colormap=Reverse(:hot), colorrange=(cr_lo, cr_hi))

    ax6 = Axis(fig[2, 3]; xlabel="x", ylabel="y", aspect=DataAspect(),
               title="(f) tanh α = 3, N = $N_err")
    cs = heatmap!(ax6, xf_err, xf_err, log10.(err_t3 .+ err_floor);
                  colormap=Reverse(:hot), colorrange=(cr_lo, cr_hi))
    Colorbar(fig[2, 4], cs; label="log10 |u_N - u_ref|")

    save(joinpath(outdir, "corner_tensor_clustering.pdf"), fig)
    save(joinpath(outdir, "corner_tensor_clustering.png"), fig)
    @printf("[19.6-julia] saved figure\n")
end

if abspath(PROGRAM_FILE) == @__FILE__
    run()
end
