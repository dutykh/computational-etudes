% spectral_derivative_accuracy.m
% Chapter 9: Physical and Fourier Space on Grids
%
% Demonstrates spectral accuracy of FFT-based differentiation.
% Tests f(x) = exp(sin(x)) whose exact derivative is cos(x)*exp(sin(x)).
% This is Etude 3 from the chapter.
%
% Author: Dr. Denys Dutykh
% Date: February 2026

function spectral_derivative_accuracy()
    % Grid sizes to test
    N_values = [4, 8, 12, 16, 20, 24, 28, 32];

    fprintf('%4s   %12s\n', 'N', 'Max error');
    fprintf('%s\n', repmat('-', 1, 20));

    for i = 1:length(N_values)
        N = N_values(i);
        x = 2*pi*(0:N-1)/N;

        % Test function and exact derivative
        v = exp(sin(x));
        w_exact = cos(x) .* exp(sin(x));

        % Spectral derivative
        w = spectral_derivative(v);

        % Maximum error
        err = max(abs(w - w_exact));
        fprintf('%4d   %12.2e\n', N, err);
    end
end

function w = spectral_derivative(v)
    N = length(v);
    v_hat = fft(v);
    k = [0:N/2-1, 0, -N/2+1:-1];  % Wavenumbers with zeroed Nyquist
    w_hat = 1i * k .* v_hat;
    w = real(ifft(w_hat));
end
