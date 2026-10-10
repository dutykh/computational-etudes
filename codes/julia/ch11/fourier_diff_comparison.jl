# fourier_diff_comparison.jl
#
# Chapter 11: Fourier Pseudospectral Methods for Periodic PDEs
#
# Compares Fourier differentiation via dense matrix multiplication (O(N^2))
# with FFT-based differentiation (O(N log N)).
#
# Test function: u(x) = exp(sin(x)) on [0, 2*pi)
# Exact derivative: u'(x) = cos(x) * exp(sin(x))
#
# Generates two figures:
#   - fourier_diff_accuracy.pdf: spectral convergence of both methods
#   - fourier_diff_timing.pdf: timing comparison (matrix vs FFT)
#
# Author: Dr. Denys Dutykh
#         Mathematics Department
#         Khalifa University of Science and Technology
#         Abu Dhabi, UAE
#
# Part of the book "Computational Etudes: A Spectral Approach"
# https://github.com/dutykh/computational-etudes

using CairoMakie
using Colors
using FFTW

# -----------------------------------------------------------------------------
# Publication-quality CairoMakie configuration
# -----------------------------------------------------------------------------
set_theme!(Theme(
    fontsize = 11,
    fonts = (regular = "CMU Serif", bold = "CMU Serif Bold", italic = "CMU Serif Italic"),
    Axis = (xlabelsize = 11, ylabelsize = 11, titlesize = 12,
            xticklabelsize = 10, yticklabelsize = 10,
            spinewidth = 0.8, xtickwidth = 0.8, ytickwidth = 0.8),
    Legend = (labelsize = 10,),
))

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
const NAVY  = colorant"#142D6E"
const SKY   = colorant"#7896D2"
const CORAL = colorant"#E74C3C"
const TEAL  = colorant"#16A085"

# Output path
SCRIPT_DIR = @__DIR__
OUTPUT_DIR = joinpath(SCRIPT_DIR, "..", "..", "..", "textbook", "figures", "ch11", "julia")
mkpath(OUTPUT_DIR)

# -----------------------------------------------------------------------------
# Fourier differentiation matrix
# -----------------------------------------------------------------------------
"""
    fourier_diff_matrix(N::Int) -> Matrix{Float64}

Build the N x N Fourier differentiation matrix for the periodic grid
x_j = 2*pi*j/N, j = 0, ..., N-1.

The (i, j) entry is the derivative of the j-th Lagrange basis function
evaluated at x_i.
"""
function fourier_diff_matrix(N::Int)
    D = zeros(N, N)
    h = 2.0 * pi / N
    for i in 0:N-1
        for j in 0:N-1
            if i != j
                D[i+1, j+1] = 0.5 * (-1.0)^(i + j) / tan((i - j) * h / 2.0)
            end
        end
    end
    return D
end

# -----------------------------------------------------------------------------
# FFT-based differentiation
# -----------------------------------------------------------------------------
"""
    fourier_diff_fft(u::AbstractVector) -> Vector{Float64}

Compute the derivative of u on a periodic grid using FFT.
u : array of length N (values at x_j = 2*pi*j/N)
"""
function fourier_diff_fft(u::AbstractVector)
    N = length(u)
    u_hat = fft(u)
    # Wavenumbers in FFT order
    if N % 2 == 0
        k = [collect(0:N÷2-1); 0; collect(1-N÷2:-1)]
    else
        k = [collect(0:(N-1)÷2); collect(-(N-1)÷2:-1)]
    end
    du_hat = 1im .* k .* u_hat
    return real.(ifft(du_hat))
end

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
function main()
    # Test function and exact derivative
    u_func(x) = exp(sin(x))
    du_exact_func(x) = cos(x) * exp(sin(x))

    # -------------------------------------------------------------------------
    # Figure 1: Accuracy comparison
    # -------------------------------------------------------------------------
    N_values = 4:2:50
    err_matrix = zeros(length(N_values))
    err_fft    = zeros(length(N_values))

    for (idx, N) in enumerate(N_values)
        x = 2.0 * pi .* collect(0:N-1) ./ N
        u = u_func.(x)
        du_exact = du_exact_func.(x)

        # Matrix method
        D = fourier_diff_matrix(N)
        du_mat = D * u
        err_matrix[idx] = maximum(abs.(du_mat .- du_exact))

        # FFT method
        du_fft_val = fourier_diff_fft(u)
        err_fft[idx] = maximum(abs.(du_fft_val .- du_exact))
    end

    fig1 = Figure(size = (600, 450))
    ax1 = Axis(fig1[1, 1],
               xlabel = L"N$ $ (number of grid points)",
               ylabel = L"\Vert u_x^{\mathrm{num}} - u_x^{\mathrm{exact}}\Vert_\infty",
               title  = L"Spectral Convergence: $u(x) = e^{\sin x}$",
               yscale = log10)

    scatterlines!(ax1, collect(N_values), err_matrix, color = NAVY, marker = :circle,
                  markersize = 7, linewidth = 1.5, label = L"Matrix $Du$")
    scatterlines!(ax1, collect(N_values), err_fft, color = CORAL, marker = :rect,
                  markersize = 7, linewidth = 1.5, linestyle = :dash,
                  label = L"FFT: $\mathcal{F}^{-1}(ik\,\hat{u})$")

    ylims!(ax1, 1e-16, nothing)
    axislegend(ax1, position = :rt, framevisible = true)
    hidespines!(ax1, :t, :r)

    save(joinpath(OUTPUT_DIR, "fourier_diff_accuracy.pdf"), fig1)
    save(joinpath(OUTPUT_DIR, "fourier_diff_accuracy.png"), fig1, px_per_unit = 4)
    println("Figure saved: $(joinpath(OUTPUT_DIR, "fourier_diff_accuracy.pdf"))")

    # -------------------------------------------------------------------------
    # Figure 2: Timing comparison
    # -------------------------------------------------------------------------
    N_timing = [2^p for p in 3:13]   # 8 to 8192
    time_matrix = Float64[]
    time_fft    = Float64[]
    n_repeats = 5

    for N in N_timing
        x = 2.0 * pi .* collect(0:N-1) ./ N
        u = u_func.(x)

        # Time matrix method (only up to N=2048)
        if N <= 2048
            D = fourier_diff_matrix(N)
            t0 = time_ns()
            for _ in 1:n_repeats
                D * u
            end
            t1 = time_ns()
            push!(time_matrix, (t1 - t0) / n_repeats / 1e9)
        else
            push!(time_matrix, NaN)
        end

        # Time FFT method
        t0 = time_ns()
        for _ in 1:n_repeats
            fourier_diff_fft(u)
        end
        t1 = time_ns()
        push!(time_fft, (t1 - t0) / n_repeats / 1e9)
    end

    fig2 = Figure(size = (600, 450))
    ax2 = Axis(fig2[1, 1],
               xlabel = L"N",
               ylabel = "Time (seconds)",
               title  = "Computational Cost: Matrix vs FFT Differentiation",
               xscale = log10, yscale = log10)

    # Matrix timing (only where valid)
    valid = .!isnan.(time_matrix)
    scatterlines!(ax2, N_timing[valid], time_matrix[valid], color = NAVY, marker = :circle,
                  markersize = 7, linewidth = 1.5, label = L"Matrix $O(N^2)$")
    scatterlines!(ax2, N_timing, time_fft, color = CORAL, marker = :rect,
                  markersize = 7, linewidth = 1.5, label = L"FFT $O(N \log N)$")

    # Reference slopes
    N_ref = Float64.([64, 2048])
    N_valid = N_timing[valid]
    scale_n2 = time_matrix[valid][end] .* (N_ref ./ N_valid[end]).^2
    lines!(ax2, N_ref, scale_n2, color = SKY, linewidth = 1.0, linestyle = :dot,
           label = L"\sim N^2")

    N_ref2 = Float64.([64, 8192])
    scale_nlogn = time_fft[end] .* (N_ref2 .* log.(N_ref2)) ./
                  (N_timing[end] * log(N_timing[end]))
    lines!(ax2, N_ref2, scale_nlogn, color = colorant"#E67E22", linewidth = 1.0,
           linestyle = :dot, label = L"\sim N \log N")

    axislegend(ax2, position = :lt, framevisible = true)
    hidespines!(ax2, :t, :r)

    save(joinpath(OUTPUT_DIR, "fourier_diff_timing.pdf"), fig2)
    save(joinpath(OUTPUT_DIR, "fourier_diff_timing.png"), fig2, px_per_unit = 4)
    println("Figure saved: $(joinpath(OUTPUT_DIR, "fourier_diff_timing.pdf"))")
end

main()
