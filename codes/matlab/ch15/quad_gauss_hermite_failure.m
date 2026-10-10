% quad_gauss_hermite_failure.m
% Chapter 15: Quadrature in Spectral Methods
%
% Computational Etude 15.9: Gauss--Hermite vs Truncated Trapezoidal Rule.
%
% Compares two quadrature strategies for the oscillatory integral
%
%     I = int_{-inf}^{inf} e^{-x^2} cos(x^3) dx
%
% Gauss--Hermite quadrature: The weight function e^{-x^2} is built into
% the rule, so I_n = sum_k w_k cos(x_k^3). Theory predicts root-exponential
% convergence: errors decay like O(exp(-C sqrt(n))).
%
% Truncated trapezoidal rule: The full integrand g(x) = e^{-x^2} cos(x^3)
% is integrated on [-L, L] using the composite trapezoidal rule. The
% trapezoidal rule converges exponentially in n, dramatically outperforming
% Gauss--Hermite.
%
% This demonstrates that Gauss rules matched to the weight function are not
% always the best choice; a simpler rule on a truncated domain can be far
% more effective when the integrand decays rapidly.
%
% Generates Figure 15.9: gauss_hermite_failure.pdf
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
% Reference value via high-accuracy adaptive quadrature
% -------------------------------------------------------------------------
integrand = @(x) exp(-x.^2) .* cos(x.^3);
I_exact = integral(integrand, -Inf, Inf, 'AbsTol', 1e-15, 'RelTol', 1e-15);
fprintf('Reference integral value: %.16e\n', I_exact);

% -------------------------------------------------------------------------
% Gauss--Hermite quadrature: I_n = sum w_k cos(x_k^3)
% The weight function e^{-x^2} is absorbed into the weights.
% -------------------------------------------------------------------------
n_gh = 10:10:1000;
err_gh = zeros(size(n_gh));

for i = 1:length(n_gh)
    n = n_gh(i);
    [x_h, w_h] = gauss_hermite(n);
    I_gh = w_h' * cos(x_h.^3);
    err_gh(i) = abs(I_gh - I_exact);
end

% -------------------------------------------------------------------------
% Truncated trapezoidal rule on [-L, L]
% g(x) = e^{-x^2} cos(x^3), with L = 6 (truncation error < 10^{-15})
% -------------------------------------------------------------------------
L = 6.0;
n_trap = 10:2:200;
err_trap = zeros(size(n_trap));

for i = 1:length(n_trap)
    n = n_trap(i);
    h = 2.0 * L / n;
    x = linspace(-L, L, n + 1)';
    g = exp(-x.^2) .* cos(x.^3);
    % Composite trapezoidal rule
    I_trap = h * (0.5 * g(1) + sum(g(2:end-1)) + 0.5 * g(end));
    err_trap(i) = abs(I_trap - I_exact);
end

% -------------------------------------------------------------------------
% Create figure
% -------------------------------------------------------------------------
fl = 1e-16;

fig = figure('Position', [100, 100, 700, 500]);

% Gauss--Hermite errors
semilogy(n_gh, max(err_gh, fl), 'o', 'Color', NAVY, ...
    'MarkerSize', 3, 'MarkerEdgeColor', NAVY, 'MarkerFaceColor', NAVY, ...
    'DisplayName', 'Gauss-Hermite');
hold on;

% Trapezoidal errors
semilogy(n_trap, max(err_trap, fl), 's', 'Color', CORAL, ...
    'MarkerSize', 3.5, 'MarkerFaceColor', 'w', ...
    'MarkerEdgeColor', CORAL, 'LineWidth', 0.7, ...
    'DisplayName', 'Trapezoidal on $[-6, 6]$');

% Machine epsilon reference
yline(eps, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);
text(950, 4e-16, '$\epsilon_{\mathrm{mach}}$', 'Interpreter', 'latex', ...
    'FontSize', 9, 'Color', [0.5 0.5 0.5], 'HorizontalAlignment', 'right');

hold off;

xlabel('Number of quadrature points $n$', 'Interpreter', 'latex');
ylabel('Absolute error $|I - I_n|$', 'Interpreter', 'latex');
title('Quadrature for $\int_{-\infty}^{\infty} e^{-x^2} \cos(x^3)\, dx$', ...
    'Interpreter', 'latex');
legend('Interpreter', 'latex', 'Location', 'northeast', ...
    'EdgeColor', [0.5 0.5 0.5], 'FontSize', 10);
grid on; set(gca, 'GridAlpha', 0.3);
xlim([0, 1050]);
ylim([1e-16, 1e1]);

% Save output
exportgraphics(gcf, fullfile(output_dir, 'gauss_hermite_failure.pdf'), 'ContentType', 'vector');
exportgraphics(gcf, fullfile(output_dir, 'gauss_hermite_failure.png'), 'Resolution', 300);
fprintf('Figure saved to %s\n', fullfile(output_dir, 'gauss_hermite_failure.pdf'));

% =========================================================================
% Local functions
% =========================================================================

function [x, w] = gauss_hermite(n)
    % Gauss-Hermite nodes and weights via Golub-Welsch.
    j = (1:n-1)';
    beta = sqrt(j / 2);
    J = diag(beta, 1) + diag(beta, -1);
    [V, D] = eig(J);
    [x, idx] = sort(diag(D));
    w = sqrt(pi) * V(1, idx)'.^2;
end
