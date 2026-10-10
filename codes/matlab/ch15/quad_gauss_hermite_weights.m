% quad_gauss_hermite_weights.m
% Chapter 15: Quadrature in Spectral Methods
%
% Computational Etude 15.8: Gauss--Hermite Weight Distribution.
%
% Examines the pathological behaviour of Gauss--Hermite quadrature weights
% for large n, which limits the practical utility of this classical rule
% for high-accuracy integration.
%
% Left panel: For n = 100, the Gauss--Hermite weights w_k are plotted
% on a semilogarithmic scale against the corresponding nodes x_k. The
% weights span a vast dynamic range: the central weights are of O(1),
% while the outermost weights fall far below machine epsilon.
%
% Right panel: The fraction of weights below machine epsilon is plotted
% as a function of n for n = 10, 20, ..., 500. This fraction increases
% monotonically and approaches 1, showing that Gauss--Hermite quadrature
% becomes numerically useless in standard double-precision for large n.
%
% Generates Figure 15.8: gauss_hermite_weights.pdf
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

% Machine epsilon
eps_mach = eps;  % ~2.22e-16

% -------------------------------------------------------------------------
% Left panel: weights for n = 100
% -------------------------------------------------------------------------
n_demo = 100;
[x_h, w_h] = gauss_hermite(n_demo);

% -------------------------------------------------------------------------
% Right panel: fraction of weights below eps vs n
% -------------------------------------------------------------------------
n_values = 10:10:500;
fraction_below = zeros(size(n_values));

for i = 1:length(n_values)
    n = n_values(i);
    [~, w] = gauss_hermite(n);
    fraction_below(i) = sum(abs(w) < eps_mach) / n;
end

% -------------------------------------------------------------------------
% Create figure
% -------------------------------------------------------------------------
fig = figure('Position', [100, 100, 1100, 450]);

% --- Left panel: weight distribution ---
subplot(1, 2, 1);

semilogy(x_h, abs(w_h), 'o', 'Color', NAVY, 'MarkerSize', 2.5, ...
    'MarkerEdgeColor', NAVY, 'MarkerFaceColor', NAVY, ...
    'DisplayName', '$|w_k|$');
hold on;

% Machine epsilon reference line
yline(eps_mach, '--', 'Color', CORAL, 'LineWidth', 1.0, ...
    'DisplayName', 'Machine $\epsilon \approx 2.2 \times 10^{-16}$');

hold off;

xlabel('Node $x_k$', 'Interpreter', 'latex');
ylabel('Weight $|w_k|$', 'Interpreter', 'latex');
title('(a) Gauss-Hermite weights ($n = 100$)', 'Interpreter', 'latex');
legend('Interpreter', 'latex', 'Location', 'south', ...
    'EdgeColor', [0.5 0.5 0.5], 'FontSize', 9);
grid on; set(gca, 'GridAlpha', 0.3);
ylim([1e-300, 1e2]);

% --- Right panel: fraction below epsilon ---
subplot(1, 2, 2);

plot(n_values, fraction_below, 'o-', 'Color', NAVY, 'MarkerSize', 3, ...
    'LineWidth', 1.2, 'MarkerEdgeColor', NAVY, 'MarkerFaceColor', NAVY);
hold on;
yline(1.0, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);
hold off;

xlabel('Number of points $n$', 'Interpreter', 'latex');
ylabel('Fraction of $|w_k| < \epsilon$', 'Interpreter', 'latex');
title('(b) Fraction of underflowed weights', 'Interpreter', 'latex');
xlim([0, 510]);
ylim([-0.02, 1.05]);
grid on; set(gca, 'GridAlpha', 0.3);

% Save output
exportgraphics(gcf, fullfile(output_dir, 'gauss_hermite_weights.pdf'), 'ContentType', 'vector');
exportgraphics(gcf, fullfile(output_dir, 'gauss_hermite_weights.png'), 'Resolution', 300);
fprintf('Figure saved to %s\n', fullfile(output_dir, 'gauss_hermite_weights.pdf'));

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
