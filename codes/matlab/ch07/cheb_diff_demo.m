%% cheb_diff_demo.m - Chebyshev spectral differentiation demonstration
%
% Demonstrates Chebyshev spectral differentiation using the "Witch of Agnesi"
% test function u(x) = 1/(1 + 4x^2). Shows the function, its derivative,
% and the spectral accuracy for different values of N.
%
% This script generates Figure 6.4 for Chapter 6.
%
% Author: Dr. Denys Dutykh
%         Mathematics Department
%         Khalifa University of Science and Technology
%         Abu Dhabi, UAE
%
% Part of the book "Computational Etudes: A Spectral Approach"
% https://github.com/dutykh/computational-etudes
%
% Last modified: January 2026

clear;
close all;
clc;

%% Publication-quality figure settings
set(groot, 'DefaultAxesFontSize', 10);
set(groot, 'DefaultAxesFontName', 'CMU Serif');
set(groot, 'DefaultTextFontSize', 10);
set(groot, 'DefaultTextFontName', 'CMU Serif');
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');
set(groot, 'DefaultAxesLineWidth', 0.8);
set(groot, 'DefaultLineLineWidth', 1.5);

%% Color scheme
NAVY = [20, 45, 110] / 255;
SKY = [120, 150, 210] / 255;
CORAL = [231, 76, 60] / 255;
TEAL = [22, 160, 133] / 255;
PURPLE = [142, 68, 173] / 255;
ORANGE = [230, 126, 34] / 255;

%% Output directory
script_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch06', 'matlab');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% Test functions: Witch of Agnesi
u_func = @(x) 1 ./ (1 + 4 * x.^2);
u_prime_exact = @(x) -8 * x ./ (1 + 4 * x.^2).^2;

%% Fine grid for exact solution
N_FINE = 500;
x_fine = linspace(-1, 1, N_FINE)';

%% Create 2x2 figure
fig = figure('Position', [100, 100, 1000, 800]);

N_values = [10, 20];

for col = 1:2
    N = N_values(col);

    % Construct differentiation matrix
    [D, x] = cheb_matrix(N);

    % Evaluate function at Chebyshev points
    v = u_func(x);

    % Compute numerical derivative
    w = D * v;

    % Exact derivative at Chebyshev points
    w_exact = u_prime_exact(x);

    % Error
    error_val = max(abs(w - w_exact));

    % Top row: Function and numerical points
    subplot(2, 2, col);
    plot(x_fine, u_func(x_fine), '-', 'Color', NAVY, 'LineWidth', 1.5, ...
        'DisplayName', '$u(x) = 1/(1+4x^2)$');
    hold on;
    plot(x, v, 'o', 'Color', CORAL, 'MarkerSize', 6, ...
        'MarkerFaceColor', CORAL, 'MarkerEdgeColor', 'white', ...
        'DisplayName', 'Chebyshev points');
    hold off;

    xlabel('$x$');
    ylabel('$u(x)$');
    title(sprintf('Function with $N = %d$ points', N));
    legend('Location', 'northeast', 'FontSize', 9);
    xlim([-1.05, 1.05]);
    box off;
    grid on;
    set(gca, 'GridAlpha', 0.2);

    % Bottom row: Derivative comparison
    subplot(2, 2, col + 2);
    plot(x_fine, u_prime_exact(x_fine), '-', 'Color', NAVY, 'LineWidth', 1.5, ...
        'DisplayName', 'Exact $u''(x)$');
    hold on;
    plot(x, w, 'o', 'Color', TEAL, 'MarkerSize', 6, ...
        'MarkerFaceColor', TEAL, 'MarkerEdgeColor', 'white', ...
        'DisplayName', 'Spectral');
    hold off;

    xlabel('$x$');
    ylabel('$u''(x)$');
    title(sprintf('Derivative with $N = %d$ (max error: %.2e)', N, error_val));
    legend('Location', 'northeast', 'FontSize', 9);
    xlim([-1.05, 1.05]);
    box off;
    grid on;
    set(gca, 'GridAlpha', 0.2);
end

% Main title
sgtitle('Chebyshev Differentiation of $u(x) = 1/(1+4x^2)$ (Witch of Agnesi)', ...
    'FontSize', 13);

%% Save figure
output_file = fullfile(output_dir, 'cheb_diff_demo');
exportgraphics(fig, [output_file, '.pdf'], 'ContentType', 'vector');
exportgraphics(fig, [output_file, '.png'], 'Resolution', 300);
fprintf('Figure saved to: %s.pdf\n', output_file);

%% Print convergence table
fprintf('\nChebyshev differentiation accuracy for u(x) = 1/(1+4x^2):\n');
fprintf('%s\n', repmat('-', 1, 50));
fprintf('%6s %14s %14s\n', 'N', 'Max Error', 'Convergence');
fprintf('%s\n', repmat('-', 1, 50));

prev_error = [];
for N = [4, 8, 10, 12, 16, 20, 24, 32, 40, 48, 64]
    [D, x] = cheb_matrix(N);
    v = u_func(x);
    w = D * v;
    w_exact = u_prime_exact(x);
    error_val = max(abs(w - w_exact));

    if ~isempty(prev_error) && error_val > 0
        rate = log10(prev_error / error_val) / log10(2);
        fprintf('%6d %14.6e %14.2f\n', N, error_val, rate);
    else
        fprintf('%6d %14.6e %14s\n', N, error_val, '---');
    end

    prev_error = error_val;
end

fprintf('%s\n', repmat('-', 1, 50));
fprintf('Note: The Witch of Agnesi is analytic in [-1,1] with poles at x = +/-i/2,\n');
fprintf('so exponential convergence is expected until machine precision is reached.\n');
