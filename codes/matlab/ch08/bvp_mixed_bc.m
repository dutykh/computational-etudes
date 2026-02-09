%% bvp_mixed_bc.m - Mixed boundary conditions with Chebyshev spectral collocation
%
% Solves a BVP with mixed boundary conditions:
%
%     u_xx = -cos(pi*x/2),  x in (-1, 1)
%     u(-1) = 0              (Dirichlet: fixed temperature)
%     u'(1) = 0              (Neumann: insulated boundary)
%
% Exact solution: u(x) = (4/pi^2)*cos(pi*x/2) + (2/pi)*(x + 1)
%
% This demonstrates the "boundary bordering" technique: instead of extracting
% the interior system (matrix stripping), we keep the full (N+1)x(N+1) system
% and replace boundary rows with the appropriate conditions.
%
% Author: Dr. Denys Dutykh
%         Mathematics Department
%         Khalifa University of Science and Technology
%         Abu Dhabi, UAE
%
% Part of the book "Computational Etudes: A Spectral Approach"
% https://github.com/dutykh/computational-etudes
%
% Last modified: February 2026

clear;
close all;
clc;

% Add path to Chapter 7 functions
addpath('../ch07');

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
output_dir = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch08', 'matlab');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% Problem functions
rhs_function = @(x) -cos(pi * x / 2);
exact_solution = @(x) (4 / pi^2) * cos(pi * x / 2) + (2 / pi) * (x + 1);

%% Solve with N = 16
N = 16;
[x, u_num] = solve_mixed_bc(N, rhs_function);
u_exact_grid = exact_solution(x);
error_N16 = max(abs(u_num - u_exact_grid));

%% Create figure
fig = figure('Position', [100, 100, 1100, 450]);

% Fine grid for exact solution
x_fine = linspace(-1, 1, 500)';
u_exact_fine = exact_solution(x_fine);

% Panel 1: Solution
subplot(1, 2, 1);
plot(x_fine, u_exact_fine, '-', 'Color', NAVY, 'LineWidth', 1.5, 'DisplayName', 'Exact');
hold on;
plot(x, u_num, 'o', 'Color', TEAL, 'MarkerSize', 6, ...
    'MarkerFaceColor', TEAL, 'MarkerEdgeColor', 'white', 'DisplayName', 'Spectral');

% Annotate boundary conditions
text(-0.7, 0.2, '$u(-1) = 0$', 'FontSize', 10, 'Color', PURPLE);
u_right = exact_solution(1.0);
plot([0.9, 1.05], [u_right, u_right], '-', 'Color', CORAL, 'LineWidth', 2);
text(0.45, 1.1, "$u'(1) = 0$", 'FontSize', 10, 'Color', CORAL);
hold off;

xlabel('$x$');
ylabel('$u(x)$');
title(sprintf('Solution ($N = %d$, max error: %.2e)', N, error_N16));
legend('Location', 'northwest', 'FontSize', 9);
xlim([-1.05, 1.05]);
box off;
grid on;
set(gca, 'GridAlpha', 0.2);

% Panel 2: Convergence study
subplot(1, 2, 2);
N_values = [4, 6, 8, 10, 12, 16, 20, 24, 28, 32, 36, 40, 48];
errors = zeros(size(N_values));

for idx = 1:length(N_values)
    N_val = N_values(idx);
    [x_n, u_n] = solve_mixed_bc(N_val, rhs_function);
    u_exact_n = exact_solution(x_n);
    errors(idx) = max(max(abs(u_n - u_exact_n)), 1e-16);
end

semilogy(N_values, errors, 'o-', 'Color', CORAL, 'LineWidth', 1.5, ...
    'MarkerSize', 6, 'MarkerFaceColor', CORAL, 'MarkerEdgeColor', 'white');

xlabel('$N$');
ylabel('Max Error');
title('Convergence');
xlim([0, 52]);
ylim([1e-16, 1e1]);
grid on;
box off;

% Machine precision line
hold on;
yline(2.2e-16, ':', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1);
text(50, 4e-16, 'Machine precision', 'FontSize', 9, 'Color', [0.5, 0.5, 0.5], ...
     'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
hold off;

% Main title
sgtitle("Mixed BCs: $u'' = -\cos(\pi x/2)$, $u(-1)=0$, $u'(1)=0$", 'FontSize', 13);

%% Save figure
output_file = fullfile(output_dir, 'mixed_bc');
exportgraphics(fig, [output_file, '.pdf'], 'ContentType', 'vector');
exportgraphics(fig, [output_file, '.png'], 'Resolution', 300);
fprintf('Figure saved to: %s.pdf\n', output_file);

%% Print convergence table
fprintf('\nConvergence Table: Mixed BC Problem\n');
fprintf('%s\n', repmat('-', 1, 40));
fprintf('%6s %14s\n', 'N', 'Max Error');
fprintf('%s\n', repmat('-', 1, 40));
for idx = 1:length(N_values)
    fprintf('%6d %14.2e\n', N_values(idx), errors(idx));
end
fprintf('%s\n', repmat('-', 1, 40));

%% Local function: solve BVP with mixed BCs using boundary bordering
function [x, u] = solve_mixed_bc(N, rhs_func)
    % Get Chebyshev matrix and grid
    [D, x] = cheb_matrix(N);
    D2 = D * D;  % Second derivative matrix

    % Build full system: L u = rhs
    L = D2;
    rhs = rhs_func(x);

    % Neumann condition at x(1) = 1: u'(1) = 0
    L(1, :) = D(1, :);
    rhs(1) = 0;

    % Dirichlet condition at x(N+1) = -1: u(-1) = 0
    L(N + 1, :) = 0;
    L(N + 1, N + 1) = 1;
    rhs(N + 1) = 0;

    % Solve
    u = L \ rhs;
end
