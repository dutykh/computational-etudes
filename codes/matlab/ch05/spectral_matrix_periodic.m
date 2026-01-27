function [D, x] = spectral_matrix_periodic(N)
%SPECTRAL_MATRIX_PERIODIC Construct the periodic spectral differentiation matrix.
%
%   [D, x] = spectral_matrix_periodic(N) constructs the N x N spectral
%   differentiation matrix for periodic functions on [0, 2*pi) with N
%   equispaced nodes.
%
%   Input:
%       N   - Number of grid points (should be even for best results)
%
%   Output:
%       D   - N x N differentiation matrix
%       x   - N x 1 vector of grid points
%
%   The matrix entries are:
%       D(j,k) = 0.5 * (-1)^(j-k) * cot((j-k)*pi/N)  for j ~= k
%       D(j,j) = 0
%
%   Properties:
%       - Skew-symmetric: D' = -D
%       - Toeplitz (and circulant due to periodicity)
%       - Dense: all off-diagonal entries are nonzero
%       - Exact for trigonometric polynomials of degree <= N/2
%
%   Example:
%       N = 16;
%       [D, x] = spectral_matrix_periodic(N);
%       u = sin(x);           % Test function
%       du_exact = cos(x);    % Exact derivative
%       du_approx = D * u;    % Spectral approximation
%       max_error = max(abs(du_approx - du_exact))  % Should be ~eps
%
%   Reference:
%       Trefethen, L.N. (2000). "Spectral Methods in MATLAB", SIAM.
%
%   Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
%   Part of "Computational Études: A Spectral Approach"
%   https://github.com/dutykh/computational-etudes

    h = 2 * pi / N;
    x = h * (0:N-1)';
    D = zeros(N, N);

    for i = 1:N
        for j = 1:N
            if i ~= j
                diff = i - j;
                D(i, j) = 0.5 * ((-1)^diff) / tan(diff * pi / N);
            end
            % D(i,i) = 0 (already initialized)
        end
    end
end
