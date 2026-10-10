% quad_exactness_table.m
% Chapter 15: Quadrature in Spectral Methods
%
% Computational Etude 15.3: Exactness of Newton--Cotes Rules.
%
% Demonstrates the stark difference in exactness between the monomial basis
% {1, x, x^2, ...} and the Chebyshev basis {T_0, T_1, T_2, ...} for
% Newton--Cotes quadrature with n = 30 equispaced nodes.
%
% For even k from 26 to 38, this script computes and prints:
%
%     |E_n(x^k)| = |integral_{-1}^{1} x^k dx - sum_j w_j x_j^k|
%     |E_n(T_k)| = |integral_{-1}^{1} T_k(x) dx - sum_j w_j T_k(x_j)|
%
% The exact integrals are:
%     integral x^k dx  = 2/(k+1)  for even k,  0 for odd k,
%     integral T_k dx   = 2/(1-k^2) for even k, 0 for odd k.
%
% This script prints a formatted table to stdout (no figure output).
%
% Author: Dr. Denys Dutykh
%         Mathematics Department
%         Khalifa University of Science and Technology
%         Abu Dhabi, UAE
%
% Part of the book "Computational Etudes: A Spectral Approach"
% https://github.com/dutykh/computational-etudes

clear; close all; clc;

% -------------------------------------------------------------------------
% Compute Newton--Cotes weights for n = 30
% -------------------------------------------------------------------------
n = 30;
[x, w] = newton_cotes_weights(n);

% Header
fprintf('\nNewton--Cotes exactness table (n = %d)\n', n);
fprintf('%s\n', repmat('=', 1, 58));
fprintf('%4s   %16s   %16s   %10s\n', 'k', '|E_n(x^k)|', '|E_n(T_k)|', 'Ratio');
fprintf('%s\n', repmat('-', 1, 58));

for k = 26:2:38
    % Monomial x^k
    % Exact integral of x^k over [-1, 1] for even k
    I_mono_exact = 2.0 / (k + 1);
    I_mono_quad = w' * (x.^k);
    err_mono = abs(I_mono_quad - I_mono_exact);

    % Chebyshev polynomial T_k
    % Exact integral of T_k(x) over [-1, 1] for even k: 2/(1-k^2)
    I_cheb_exact = 2.0 / (1.0 - k^2);
    T_vals = cos(k * acos(x));
    I_cheb_quad = w' * T_vals;
    err_cheb = abs(I_cheb_quad - I_cheb_exact);

    % Ratio
    if err_mono > 0 && err_cheb > 0
        ratio = err_cheb / err_mono;
    else
        ratio = Inf;
    end

    fprintf('%4d   %16.6e   %16.6e   %10.2e\n', k, err_mono, err_cheb, ratio);
end

fprintf('%s\n', repmat('-', 1, 58));
fprintf('\nFor k <= n = %d, Newton--Cotes is exact for both bases.\n', n);
fprintf('For k > n, errors in T_k are typically much larger than in x^k,\n');
fprintf('because T_k concentrates all its degree-k content in one mode.\n\n');

% =========================================================================
% Local functions
% =========================================================================

function [x, w] = newton_cotes_weights(n)
    x = linspace(-1, 1, n+1)';
    V = zeros(n+1);
    for k = 0:n
        V(k+1,:) = x.^k;
    end
    rhs = arrayfun(@(k) (1 - (-1)^(k+1))/(k+1), (0:n)');
    w = V \ rhs;
end
