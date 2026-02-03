function [D, x] = cheb_matrix(N)
%CHEB_MATRIX Chebyshev differentiation matrix
%
%   [D, x] = cheb_matrix(N) constructs the Chebyshev differentiation matrix
%   D_N of size (N+1) x (N+1) and the Chebyshev-Gauss-Lobatto points.
%
%   Input:
%       N : number of intervals (matrix is (N+1) x (N+1))
%
%   Output:
%       D : Chebyshev differentiation matrix
%       x : Chebyshev points x_j = cos(j*pi/N), j=0,...,N
%
%   Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
%   Part of "Computational Etudes: A Spectral Approach"
%   https://github.com/dutykh/computational-etudes

if N == 0
    D = 0;
    x = 1;
    return
end

% Chebyshev-Gauss-Lobatto points
x = cos(pi * (0:N)' / N);

% Weights: c_0 = c_N = 2, others = 1
c = ones(N+1, 1);
c(1) = 2;
c(N+1) = 2;

% Build the matrix
X = repmat(x, 1, N+1);
dX = X - X';

% Create c_i / c_j matrix
C = (c * (1./c'));

% Sign matrix: (-1)^{i+j}
ii = (0:N)';
Sign = (-1).^(ii + ii');

% Off-diagonal entries
D = C .* Sign ./ dX;

% Replace NaN on diagonal (from 0/0) with 0
D(1:N+2:end) = 0;

% Fix diagonal using negative sum trick
D = D - diag(sum(D, 2));

end
