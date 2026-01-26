function [D, x] = cheb_matrix(N)
%% cheb_matrix - Construct Chebyshev differentiation matrix
%
% Constructs the Chebyshev differentiation matrix and grid points using
% the Chebyshev-Gauss-Lobatto points: x_j = cos(j*pi/N), j = 0, 1, ..., N
%
% Uses the "negative sum trick" for diagonal entries to ensure numerical
% stability (derivative of a constant is exactly zero).
%
% Inputs:
%   N - Number of intervals (matrix size is (N+1) x (N+1))
%
% Outputs:
%   D - Chebyshev differentiation matrix, size (N+1) x (N+1)
%   x - Chebyshev-Gauss-Lobatto points on [-1, 1], ordered from x_0=1 to x_N=-1
%
% Notes:
%   The matrix entries (before applying negative sum trick) are:
%
%   Corner entries:
%       D(1,1) = (2N^2 + 1) / 6
%       D(N+1,N+1) = -(2N^2 + 1) / 6
%
%   Diagonal entries (j = 2, ..., N in MATLAB indexing):
%       D(j,j) = -x(j) / (2(1 - x(j)^2))
%
%   Off-diagonal entries:
%       D(i,j) = (c(i) / c(j)) * (-1)^(i+j) / (x(i) - x(j))
%
%   where c(1) = c(N+1) = 2 and c(j) = 1 for j = 2, ..., N.
%
%   For numerical stability, we compute diagonals via:
%       D(j,j) = -sum_{k != j} D(j,k)
%
%   This ensures that D * ones(N+1,1) = zeros(N+1,1) to machine precision.
%
% Example:
%   [D, x] = cheb_matrix(4);
%   max(abs(D * ones(5, 1)))  % Should be ~eps
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
% Part of "Computational Etudes: A Spectral Approach"
% https://github.com/dutykh/computational-etudes

    if N == 0
        D = 0;
        x = 1;
        return;
    end

    % Chebyshev-Gauss-Lobatto points (column vector)
    x = cos(pi * (0:N)' / N);

    % Coefficient vector c: c(1) = c(N+1) = 2, others = 1
    c = ones(N + 1, 1);
    c(1) = 2;
    c(N + 1) = 2;

    % Construct the differentiation matrix using vectorized operations
    % X(i,j) = x(j), dX(i,j) = x(i) - x(j)
    X = repmat(x', N + 1, 1);
    dX = X' - X;

    % Off-diagonal entries: D(i,j) = c(i)/c(j) * (-1)^(i+j) / (x(i) - x(j))
    C = c * (1 ./ c');  % C(i,j) = c(i) / c(j)

    % Sign matrix: (-1)^(i+j) using 0-based indices
    i_idx = (0:N)';
    sign_vec = (-1) .^ i_idx;
    sign_mat = sign_vec * sign_vec';

    % Compute off-diagonal entries
    D = C .* sign_mat ./ dX;

    % Fix the diagonal (currently NaN from 0/0)
    % Apply negative sum trick: D(j,j) = -sum_{k != j} D(j,k)
    D(1:N+2:end) = 0;  % Set diagonal to zero first
    D = D - diag(sum(D, 2));  % Subtract row sums from diagonal

end
