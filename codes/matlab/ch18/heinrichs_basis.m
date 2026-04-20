function varargout = heinrichs_basis(action, N)
%% heinrichs_basis - Chapter 18 shared utility (dispatch-style).
%
% Multiple constructors dispatched by the `action` string:
%
%   [A, xi]    = heinrichs_basis('dirichlet_heinrichs', N)
%   [A, M, xi] = heinrichs_basis('clamped_heinrichs',   N)
%   [A, xi]    = heinrichs_basis('dirichlet_naive',     N)
%   [A, B]     = heinrichs_basis('clamped_naive',       N)
%
% When called with no arguments, runs the condition-number self-test.
%
% See Python docstring for derivation details.
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
% Part of "Computational Etudes: A Spectral Approach"

    script_dir = fileparts(mfilename('fullpath'));
    addpath(fullfile(script_dir, '..', 'ch07'));

    if nargin < 1
        self_test();
        return;
    end

    switch action
        case 'dirichlet_heinrichs'
            [varargout{1}, varargout{2}] = dirichlet_heinrichs(N);
        case 'clamped_heinrichs'
            [varargout{1}, varargout{2}, varargout{3}] = clamped_heinrichs(N);
        case 'dirichlet_naive'
            [varargout{1}, varargout{2}] = dirichlet_naive(N);
        case 'clamped_naive'
            [varargout{1}, varargout{2}] = clamped_naive(N);
        otherwise
            error('heinrichs_basis:unknown_action', 'Unknown action: %s', action);
    end
end


function [A, xi] = dirichlet_heinrichs(N)
    [D, x] = cheb_matrix(N);
    D2 = D * D;
    Di  = D(2:N, 2:N);
    D2i = D2(2:N, 2:N);
    xi  = x(2:N);
    S = diag(1 - xi.^2);
    X = diag(xi);
    L = -S * D2i + 4.0 * X * Di + 2.0 * eye(N - 1);
    A = S \ L;
end


function [A, M, xi] = clamped_heinrichs(N)
    [D, x] = cheb_matrix(N);
    D2 = D * D;
    D3 = D2 * D;
    D4 = D3 * D;
    Di  = D(2:N, 2:N);
    D2i = D2(2:N, 2:N);
    D3i = D3(2:N, 2:N);
    D4i = D4(2:N, 2:N);
    xi  = x(2:N);
    I_int = eye(N - 1);
    S  = diag(1 - xi.^2);
    S2 = diag((1 - xi.^2).^2);
    X  = diag(xi);
    X2 = diag(xi.^2);
    A = S2 * D4i ...
        - 16.0 * X * S * D3i ...
        + (48.0 * X2 - 24.0 * I_int) * D2i ...
        + 96.0 * X * Di ...
        + 24.0 * I_int;
    M = S2;
end


function [A, xi] = dirichlet_naive(N)
    [D, x] = cheb_matrix(N);
    D2 = D * D;
    A  = -D2(2:N, 2:N);
    xi = x(2:N);
end


function [A, B] = clamped_naive(N)
    [D, ~] = cheb_matrix(N);
    D4 = (D * D) * (D * D);
    A = D4;
    B = eye(N + 1);
    A(1,   :) = 0; A(1,   1)   = 1; B(1,   :) = 0;
    A(2,   :) = D(1, :);             B(2,   :) = 0;
    A(N,   :) = D(N + 1, :);         B(N,   :) = 0;
    A(N+1, :) = 0; A(N+1, N+1) = 1; B(N+1, :) = 0;
end


function self_test()
    fprintf('Condition-number comparison for second-order operator:\n');
    fprintf('%6s %16s %20s\n', 'N', 'kappa(naive)', 'kappa(Heinrichs)');
    for N = [16, 32, 64, 96]
        [A_naive, ~] = dirichlet_naive(N);
        [A_hein,  ~] = dirichlet_heinrichs(N);
        fprintf('%6d %16.3e %20.3e\n', N, cond(A_naive), cond(A_hein));
    end
end
