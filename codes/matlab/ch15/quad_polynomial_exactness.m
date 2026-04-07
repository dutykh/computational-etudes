% quad_polynomial_exactness.m
% Chapter 15: Quadrature in Spectral Methods
%
% Computational Etude 15.1: Polynomial Exactness Test.
%
% For a fixed number of nodes n + 1 = 33, computes the absolute monomial
% quadrature error
%
%     |E_n(x^k)| = | sum_j w_j x_j^k - integral_{-1}^{1} x^k dx |
%
% for k = 0, 1, ..., 2n = 64, using three classical interpolatory rules:
%
%     1. Newton--Cotes (equispaced nodes, Vandermonde-based weights).
%     2. Gauss--Legendre (Golub--Welsch algorithm).
%     3. Clenshaw--Curtis (Chebyshev nodes, FFT-based weights).
%
% The exact integral is 2/(k+1) for even k and 0 for odd k.
%
% Theoretical degrees of precision:
%
%     Newton--Cotes:    n   (here, k <= 32)
%     Clenshaw--Curtis: n   (here, k <= 32)
%     Gauss--Legendre:  2n + 1   (here, k <= 65, the entire test range)
%
% The figure shows that NC and CC sit at machine precision for k <= n and
% lift off for larger k, while Gauss stays at the floor across the whole
% range. Odd k are zero by symmetry and are clamped to the floor.
%
% Generates Figure 15.1b: polynomial_exactness.pdf
%
% Author: Dr. Denys Dutykh
%         Mathematics Department
%         Khalifa University of Science and Technology
%         Abu Dhabi, UAE
%
% Part of the book "Computational Etudes: A Spectral Approach"
% https://github.com/dutykh/computational-etudes

clear; close all; clc;

% Publication-quality figure settings
set(0, 'DefaultAxesFontSize', 11);
set(0, 'DefaultTextFontSize', 11);
set(0, 'DefaultLineLineWidth', 1.5);
set(0, 'DefaultAxesLineWidth', 0.8);

% Book color scheme
NAVY   = [20, 45, 110] / 255;
SKY    = [120, 150, 210] / 255;
CORAL  = [231, 76, 60] / 255;
TEAL   = [22, 160, 133] / 255;
PURPLE = [142, 68, 173] / 255;
ORANGE = [230, 126, 34] / 255;

% Output path
output_dir = fullfile(fileparts(mfilename('fullpath')), ...
    '..', '..', '..', 'textbook', 'figures', 'ch15', 'matlab');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

% -------------------------------------------------------------------------
% Compute monomial-exactness errors for the three quadrature rules
% -------------------------------------------------------------------------
n = 32;                  % n + 1 = 33 quadrature points
k_vals = (0:2*n)';       % 0, 1, ..., 64

% Build the three rules once
[x_nc, w_nc] = newton_cotes_weights(n);
[x_gl, w_gl] = gauss_legendre(n + 1);
[x_cc, w_cc] = clenshaw_curtis_weights(n);

err_nc = zeros(size(k_vals));
err_gl = zeros(size(k_vals));
err_cc = zeros(size(k_vals));

for i = 1:length(k_vals)
    k = k_vals(i);
    if mod(k, 2) == 0
        I_exact = 2 / (k + 1);
    else
        I_exact = 0;
    end
    err_nc(i) = abs(w_nc' * (x_nc .^ k) - I_exact);
    err_gl(i) = abs(w_gl' * (x_gl .^ k) - I_exact);
    err_cc(i) = abs(w_cc' * (x_cc .^ k) - I_exact);
end

% Clamp errors at the machine-epsilon floor for log plotting
eps_floor = 1e-17;
err_nc = max(err_nc, eps_floor);
err_gl = max(err_gl, eps_floor);
err_cc = max(err_cc, eps_floor);

% Quick sanity prints
fprintf('n + 1 = %d quadrature points; testing k = 0 .. %d\n', n + 1, 2 * n);
fprintf('max NC error for k <= n: %.2e\n', max(err_nc(k_vals <= n)));
fprintf('max CC error for k <= n: %.2e\n', max(err_cc(k_vals <= n)));
fprintf('max GL error for entire range: %.2e\n', max(err_gl));
fprintf('max NC error for k > n:  %.2e\n', max(err_nc(k_vals > n)));
fprintf('max CC error for k > n:  %.2e\n', max(err_cc(k_vals > n)));

% -------------------------------------------------------------------------
% Create figure
% -------------------------------------------------------------------------
fig = figure('Position', [100, 100, 800, 480]);

semilogy(k_vals, err_nc, 'o-', 'Color', CORAL, 'MarkerSize', 4, ...
    'LineWidth', 1.1, 'MarkerFaceColor', CORAL, 'DisplayName', 'Newton-Cotes');
hold on;
semilogy(k_vals, err_gl, 's-', 'Color', NAVY, 'MarkerSize', 4, ...
    'LineWidth', 1.1, 'MarkerFaceColor', NAVY, 'DisplayName', 'Gauss-Legendre');
semilogy(k_vals, err_cc, 'd-', 'Color', TEAL, 'MarkerSize', 4, ...
    'LineWidth', 1.1, 'MarkerFaceColor', TEAL, 'DisplayName', 'Clenshaw-Curtis');

% Vertical guides at k = n and k = 2n+1
xline(n, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8, ...
    'Label', 'NC/CC exactness (k = n)', 'LabelHorizontalAlignment', 'left', ...
    'LabelOrientation', 'aligned', 'FontSize', 9, ...
    'HandleVisibility', 'off');
xline(2 * n + 1, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8, ...
    'Label', 'Gauss exactness (k = 2n+1)', 'LabelHorizontalAlignment', 'left', ...
    'LabelOrientation', 'aligned', 'FontSize', 9, ...
    'HandleVisibility', 'off');
hold off;

xlabel('Monomial degree $k$', 'Interpreter', 'latex');
ylabel('Absolute error $|E_n(x^k)|$', 'Interpreter', 'latex');
title('Polynomial exactness of three quadrature rules ($n + 1 = 33$ points)', ...
    'Interpreter', 'latex', 'FontSize', 12);
legend('Location', 'southeast', 'EdgeColor', [0.5 0.5 0.5]);
grid on;
set(gca, 'GridAlpha', 0.3);
xlim([-1, 2 * n + 2]);
ylim([1e-18, 1e1]);
xticks(0:8:2*n);

% Save output
exportgraphics(gcf, fullfile(output_dir, 'polynomial_exactness.pdf'), 'ContentType', 'vector');
exportgraphics(gcf, fullfile(output_dir, 'polynomial_exactness.png'), 'Resolution', 300);
fprintf('Figure saved to %s\n', fullfile(output_dir, 'polynomial_exactness.pdf'));

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

function [x, w] = gauss_legendre(n)
    % n-point Gauss-Legendre nodes and weights via Golub-Welsch.
    if n == 1
        x = 0; w = 2; return;
    end
    j = (1:n-1)';
    beta = j ./ sqrt(4*j.^2 - 1);
    J = diag(beta, 1) + diag(beta, -1);
    [V, D] = eig(J);
    [x, idx] = sort(diag(D));
    w = 2 * V(1, idx)'.^2;
end

function [x, w] = clenshaw_curtis_weights(n)
    % Clenshaw-Curtis nodes and weights via DCT-I / FFT.
    if n == 0
        x = 0; w = 2; return;
    end
    if n == 1
        x = [1; -1]; w = [1; 1]; return;
    end
    x = cos(pi * (0:n)' / n);
    c = zeros(n+1, 1);
    c(1:2:end) = 2 ./ (1 - (0:2:n)'.^2);
    v = [c; c(n:-1:2)];  % length 2n, mirror
    f = real(fft(v));
    w = f(1:n+1) / n;
    w(1) = w(1) / 2;
    w(end) = w(end) / 2;
end
