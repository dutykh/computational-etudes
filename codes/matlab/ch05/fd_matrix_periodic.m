function [D, x] = fd_matrix_periodic(N, order)
%FD_MATRIX_PERIODIC Construct a periodic finite difference differentiation matrix.
%
%   [D, x] = fd_matrix_periodic(N, order) constructs the N x N finite
%   difference differentiation matrix for periodic functions on [0, 2*pi)
%   with N equispaced nodes.
%
%   Input:
%       N     - Number of grid points
%       order - Order of accuracy (2, 4, 6, ...)
%
%   Output:
%       D   - N x N differentiation matrix (sparse circulant)
%       x   - N x 1 vector of grid points
%
%   The stencil width is (order + 1), centered around each point.
%   Periodic boundary conditions are enforced via wrap-around indexing.
%
%   Example:
%       N = 32;
%       [D2, x] = fd_matrix_periodic(N, 2);   % 2nd order
%       [D4, x] = fd_matrix_periodic(N, 4);   % 4th order
%       [D6, x] = fd_matrix_periodic(N, 6);   % 6th order
%
%   Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
%   Part of "Computational Études: A Spectral Approach"

    h = 2 * pi / N;
    x = h * (0:N-1)';

    % Determine stencil size
    stencil_half = order / 2;
    stencil_nodes = h * (-stencil_half:stencil_half)';

    % Compute FD weights for first derivative at center
    weights = fdweights(0, stencil_nodes, 1);

    % Build circulant matrix
    D = zeros(N, N);
    for i = 1:N
        for k = 1:length(weights)
            j = mod(i - 1 + k - 1 - stencil_half, N) + 1;  % Periodic wrap-around
            D(i, j) = D(i, j) + weights(k);
        end
    end
end
