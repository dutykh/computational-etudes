function w = chebfft(v)
%CHEBFFT Chebyshev differentiation via FFT
%
%   w = chebfft(v) computes the first derivative of a function sampled at
%   Chebyshev points using the FFT.
%
%   The algorithm transforms from the x-variable to the theta-variable
%   (where x = cos(theta)), performs differentiation via FFT in theta-space,
%   and transforms back.
%
%   Input:
%       v : column vector of length N+1
%           Function values at Chebyshev points x_j = cos(j*pi/N), j=0,...,N
%
%   Output:
%       w : column vector of length N+1
%           Approximate derivative values at the same Chebyshev points
%
%   Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
%   Part of "Computational Etudes: A Spectral Approach"
%   https://github.com/dutykh/computational-etudes

v = v(:);  % Ensure column vector
N = length(v) - 1;

if N == 0
    w = 0;
    return
end

x = cos((0:N)' * pi / N);

% Extend to even periodic function in theta
% V has length 2N: V = [v_0, v_1, ..., v_N, v_{N-1}, ..., v_1]
V = [v; flipud(v(2:N))];

% FFT to get Fourier coefficients
U = real(fft(V));

% Wavenumber vector for differentiation
% k = [0, 1, 2, ..., N-1, 0, -(N-1), -(N-2), ..., -1]
ii = (0:N-1)';
k = [ii; 0; (1-N:-1)'];

% Differentiate in theta: multiply by ik
W = real(ifft(1i * k .* fft(V)));

% Transform from theta derivative to x derivative
% q'(x) = Q'(theta) / (dx/dtheta) = -Q'(theta) / sin(theta)
%       = -Q'(theta) / sqrt(1 - x^2)
w = zeros(N+1, 1);
w(2:N) = -W(2:N) ./ sqrt(1 - x(2:N).^2);

% Boundary values via summation formulas (l'Hopital's rule at x = +/- 1)
w(1) = sum(ii.^2 .* U(ii+1)) / N + 0.5 * N * U(N+1);
w(N+1) = sum((-1).^(ii+1) .* ii.^2 .* U(ii+1)) / N + 0.5 * (-1)^(N+1) * N * U(N+1);

end
