% quad_newton_cotes_runge.m
% Chapter 15: Quadrature in Spectral Methods
%
% Computational Etude 15.2: The Runge Phenomenon in Quadrature.
%
% Integrates the Runge function f(x) = 1/(1 + 25x^2) over [-1, 1] using
% three quadrature rules of increasing order n:
%
%     1. Newton--Cotes (equispaced nodes, Vandermonde-based weights).
%     2. Gauss--Legendre (Golub--Welsch algorithm).
%     3. Clenshaw--Curtis (Chebyshev nodes, FFT-based weights).
%
% The exact integral is (2/5) arctan(5).
%
% For Newton--Cotes, the weights develop large alternating values as the
% Vandermonde matrix becomes ill-conditioned: the quadrature analogue of
% the Runge phenomenon. Gauss--Legendre and Clenshaw--Curtis converge.
%
% Reference: L. N. Trefethen, "Is Gauss Quadrature Better than
%            Clenshaw--Curtis?", SIAM Review 50(1):67--87, 2008.
%
% Generates Figure 15.2: newton_cotes_runge.pdf
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
% The Runge function and its exact integral
% -------------------------------------------------------------------------
f = @(x) 1.0 ./ (1.0 + 25.0 * x.^2);
I_exact = (2.0 / 5.0) * atan(5.0);

% Range of n values
n_values = 4:2:50;
nn = length(n_values);

err_nc = zeros(nn, 1);
err_gl = zeros(nn, 1);
err_cc = zeros(nn, 1);

for i = 1:nn
    n = n_values(i);

    % Newton--Cotes
    [x_nc, w_nc] = newton_cotes_weights(n);
    I_nc = w_nc' * f(x_nc);
    err_nc(i) = abs(I_nc - I_exact);

    % Gauss--Legendre (n+1 points)
    [x_gl, w_gl] = gauss_legendre(n + 1);
    I_gl = w_gl' * f(x_gl);
    err_gl(i) = abs(I_gl - I_exact);

    % Clenshaw--Curtis
    [x_cc, w_cc] = clenshaw_curtis_weights(n);
    I_cc = w_cc' * f(x_cc);
    err_cc(i) = abs(I_cc - I_exact);
end

% Clamp zeros for log plot
err_nc = max(err_nc, 1e-16);
err_gl = max(err_gl, 1e-16);
err_cc = max(err_cc, 1e-16);

% -------------------------------------------------------------------------
% Create figure
% -------------------------------------------------------------------------
fig = figure('Position', [100, 100, 700, 500]);

semilogy(n_values, err_nc, 'o-', 'Color', CORAL, 'MarkerSize', 4, ...
    'LineWidth', 1.2, 'MarkerFaceColor', CORAL, 'DisplayName', 'Newton-Cotes');
hold on;
semilogy(n_values, err_gl, 's-', 'Color', NAVY, 'MarkerSize', 4, ...
    'LineWidth', 1.2, 'MarkerFaceColor', NAVY, 'DisplayName', 'Gauss-Legendre');
semilogy(n_values, err_cc, 'd-', 'Color', TEAL, 'MarkerSize', 4, ...
    'LineWidth', 1.2, 'MarkerFaceColor', TEAL, 'DisplayName', 'Clenshaw-Curtis');
hold off;

xlabel('Number of nodes $n$', 'Interpreter', 'latex');
ylabel('Absolute error $|I - I_n|$', 'Interpreter', 'latex');
title('Quadrature of $f(x) = 1/(1 + 25x^2)$ on $[-1, 1]$', ...
    'Interpreter', 'latex', 'FontSize', 12);
legend('Location', 'best', 'EdgeColor', [0.5 0.5 0.5]);
grid on;
set(gca, 'GridAlpha', 0.3);
xlim([n_values(1) - 1, n_values(end) + 1]);

% Save output
exportgraphics(gcf, fullfile(output_dir, 'newton_cotes_runge.pdf'), 'ContentType', 'vector');
exportgraphics(gcf, fullfile(output_dir, 'newton_cotes_runge.png'), 'Resolution', 300);
fprintf('Figure saved to %s\n', fullfile(output_dir, 'newton_cotes_runge.pdf'));

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
