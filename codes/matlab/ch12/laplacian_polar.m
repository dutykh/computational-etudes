function [L, r_pos, theta] = laplacian_polar(Nr, M)
% LAPLACIAN_POLAR  Discrete polar Laplacian on the unit disk
%
%   [L, r_pos, theta] = laplacian_polar(Nr, M)
%
%   Uses Chebyshev discretisation in r (doubling trick, r in [-1,1])
%   and Fourier discretisation in theta. Imposes u = 0 at r = 1.
%
%   Inputs:
%     Nr    - number of Chebyshev radial points (must be odd, >= 3)
%     M     - number of Fourier angular points (must be even, >= 4)
%
%   Outputs:
%     L     - discrete Laplacian matrix, size (N2*M) x (N2*M)
%     r_pos - positive radial grid points (N2 x 1, decreasing)
%     theta - angular grid points (M x 1)
%
% Chapter 12: Spectral Methods in Polar Coordinates
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology,
%         Abu Dhabi, UAE)
% Last modified: February 2026

    if mod(Nr, 2) == 0
        error('Nr must be odd to avoid a grid point at r = 0.');
    end
    if mod(M, 2) ~= 0
        error('M must be even for the Fourier angular grid.');
    end

    N2 = (Nr - 1) / 2;

    % --- Radial direction: Chebyshev on [-1, 1] ---
    [D, r] = cheb_matrix(Nr);
    D2full = D^2;

    % Block decomposition with reversed negative-r columns
    D1 = D2full(2:N2+1, 2:N2+1);
    D2b = D2full(2:N2+1, Nr:-1:N2+2);
    E1 = D(2:N2+1, 2:N2+1);
    E2 = D(2:N2+1, Nr:-1:N2+2);

    r_pos = r(2:N2+1);
    R = diag(1 ./ r_pos);

    % --- Angular direction: Fourier on [0, 2*pi) ---
    dt = 2 * pi / M;
    theta = dt * (0:M-1)';
    M2 = M / 2;

    % Fourier second derivative Toeplitz matrix
    col = zeros(M, 1);
    col(1) = -pi^2 / (3*dt^2) - 1/6;
    for j = 1:M-1
        col(j+1) = 0.5 * (-1)^(j+1) / sin(dt*j/2)^2;
    end
    D2t = toeplitz(col);

    % --- Kronecker product assembly ---
    Z = zeros(M2);
    I_half = eye(M2);
    S = [Z, I_half; I_half, Z];

    L = kron(D1 + R*E1, eye(M)) + kron(D2b + R*E2, S) + kron(R^2, D2t);
end
