function [D1_x, D2_x, x_int] = rational_chebyshev(N, L)
%% rational_chebyshev - Chapter 18 shared utility.
%
% Rational Chebyshev first and second derivative matrices on (-inf, +inf)
% via the algebraic map
%
%     x = L * t / sqrt(1 - t^2),
%
% with t on the Chebyshev-Gauss-Lobatto grid in [-1, 1]. Endpoints
% t = +- 1 map to x = +- inf; the (N - 1) interior nodes provide the
% real-line collocation grid.
%
% Chain-rule factors:
%     dt/dx         = (1 - t^2)^(3/2) / L
%     d^2 t / d x^2 = -3 t (1 - t^2)^2 / L^2.
%
% Inputs:
%   N - number of intervals on the t-grid (matrix size N-1)
%   L - map parameter (default 4)
%
% Outputs:
%   D1_x  - (N-1) x (N-1) first  derivative matrix in x
%   D2_x  - (N-1) x (N-1) second derivative matrix in x
%   x_int - (N-1) column vector of interior real-line collocation points
%
% When invoked with no output arguments, runs the Gaussian self-test.
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
% Part of "Computational Etudes: A Spectral Approach"
% https://github.com/dutykh/computational-etudes

    if nargin < 2; L = 4.0; end
    if nargin < 1 || isempty(N); self_test(); return; end

    script_dir = fileparts(mfilename('fullpath'));
    addpath(fullfile(script_dir, '..', 'ch07'));

    [D_t, t] = cheb_matrix(N);
    t_int = t(2:N);
    x_int = L .* t_int ./ sqrt(1 - t_int.^2);

    one_minus_t2 = 1 - t_int.^2;
    dt_dx   = (one_minus_t2 .^ 1.5) / L;
    d2t_dx2 = -3 .* t_int .* (one_minus_t2 .^ 2) / (L ^ 2);

    D_t_int  = D_t(2:N, 2:N);
    D_t2_int = D_t * D_t;
    D_t2_int = D_t2_int(2:N, 2:N);

    D1_x = diag(dt_dx)      * D_t_int;
    D2_x = diag(dt_dx .^ 2) * D_t2_int + diag(d2t_dx2) * D_t_int;
end

function self_test()
    for N = [16, 32, 64]
        [~, D2, x] = rational_chebyshev(N, 4.0);
        u = exp(-x.^2 / 2);
        u_xx_num = D2 * u;
        u_xx_ex  = (x.^2 - 1) .* u;
        err = max(abs(u_xx_num - u_xx_ex));
        fprintf('N = %3d  L = 4   max error in u_xx (Gaussian) = %.3e\n', N, err);
    end
end
