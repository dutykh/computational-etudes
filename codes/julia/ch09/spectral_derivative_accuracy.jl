# spectral_derivative_accuracy.jl
#
# Chapter 9: Physical and Fourier Space on Grids
#
# Demonstrates spectral accuracy of FFT-based differentiation.
# Tests f(x) = exp(sin(x)) whose exact derivative is cos(x)*exp(sin(x)).
# This is Etude 3 from the chapter.
#
# Author: Dr. Denys Dutykh
#         Mathematics Department
#         Khalifa University of Science and Technology
#         Abu Dhabi, UAE
#
# Part of the book "Computational Etudes: A Spectral Approach"
# https://github.com/dutykh/computational-etudes

using FFTW
using Printf

function spectral_derivative(v)
    N = length(v)
    v_hat = fft(v)
    k = [0:N÷2-1; 0; -N÷2+1:-1]  # Wavenumbers with zeroed Nyquist
    w_hat = 1im .* k .* v_hat
    return real.(ifft(w_hat))
end

function main()
    # Grid sizes to test
    N_values = [4, 8, 12, 16, 20, 24, 28, 32]

    @printf("%4s   %12s\n", "N", "Max error")
    println(repeat("-", 20))

    for N in N_values
        x = 2π .* collect(0:N-1) ./ N

        # Test function and exact derivative
        v = exp.(sin.(x))
        w_exact = cos.(x) .* exp.(sin.(x))

        # Spectral derivative
        w = spectral_derivative(v)

        # Maximum error
        err = maximum(abs.(w .- w_exact))
        @printf("%4d   %12.2e\n", N, err)
    end
end

# Run if executed as script
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
