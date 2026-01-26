%% bvp_linear.m - 1D Poisson equation with Chebyshev spectral collocation
%
% Solves the 1D Poisson equation:
%
%     u_xx = sin(pi*x) + 2*cos(2*pi*x),  x in (-1, 1),  u(+/-1) = 0
%
% The exact solution is determined by integrating twice and applying
% boundary conditions.
%
% This script generates Figure 7.3 for Chapter 7.
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

% Add path to Chapter 6 functions
addpath('../ch06_chebyshev_differentiation');

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
output_dir = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch07', 'matlab');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% Problem functions
rhs_function = @(x) sin(pi * x) + 2 * cos(2 * pi * x);
exact_solution = @(x) -sin(pi * x) / pi^2 + (1 - cos(2 * pi * x)) / (2 * pi^2);

%% Solve Poisson equation
N = 16;
[x, u_num] = solve_poisson_dirichlet(N, rhs_function);
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
hold off;

xlabel('$x$');
ylabel('$u(x)$');
title(sprintf('Solution ($N = %d$, max error: %.2e)', N, error_N16));
legend('Location', 'northeast', 'FontSize', 9);
xlim([-1.05, 1.05]);
box off;
grid on;
set(gca, 'GridAlpha', 0.2);

% Panel 2: Convergence study
subplot(1, 2, 2);
N_values = [4, 6, 8, 10, 12, 16, 20, 24, 32];
errors = zeros(size(N_values));

for idx = 1:length(N_values)
    N_val = N_values(idx);
    [x_n, u_n] = solve_poisson_dirichlet(N_val, rhs_function);
    u_exact_n = exact_solution(x_n);
    errors(idx) = max(abs(u_n - u_exact_n));
end

semilogy(N_values, errors, 'o-', 'Color', CORAL, 'LineWidth', 1.5, ...
    'MarkerSize', 6, 'MarkerFaceColor', CORAL, 'MarkerEdgeColor', 'white');

xlabel('$N$');
ylabel('Max Error');
title('Convergence');
xlim([0, 35]);
ylim([1e-15, 1]);
grid on;
box off;

% Machine epsilon line
hold on;
yline(1e-14, '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1);
text(30, 3e-14, '$\approx$ machine $\epsilon$', 'FontSize', 9, 'Color', [0.5, 0.5, 0.5]);
hold off;

% Main title
sgtitle('1D Poisson Equation: $u_{xx} = \sin(\pi x) + 2\cos(2\pi x)$, $u(\pm 1) = 0$', ...
    'FontSize', 13);

%% Save figure
output_file = fullfile(output_dir, 'poisson_1d');
exportgraphics(fig, [output_file, '.pdf'], 'ContentType', 'vector');
exportgraphics(fig, [output_file, '.png'], 'Resolution', 300);
fprintf('Figure saved to: %s.pdf\n', output_file);

%% Print convergence table
fprintf('\nConvergence Table: 1D Poisson Equation\n');
fprintf('%s\n', repmat('-', 1, 40));
fprintf('%6s %14s\n', 'N', 'Max Error');
fprintf('%s\n', repmat('-', 1, 40));
for idx = 1:length(N_values)
    fprintf('%6d %14.2e\n', N_values(idx), errors(idx));
end
fprintf('%s\n', repmat('-', 1, 40));

%% Local function: solve Poisson equation with Dirichlet BCs
function [x, u] = solve_poisson_dirichlet(N, rhs_func)
    % Get Chebyshev matrix and grid
    [D, x] = cheb_matrix(N);
    D2 = D * D;  % Second derivative matrix

    % Extract interior points (remove first and last rows/columns)
    % x(1) = 1 (right boundary), x(N+1) = -1 (left boundary)
    D2_int = D2(2:N, 2:N);
    f_int = rhs_func(x(2:N));

    % Solve the linear system
    u_int = D2_int \ f_int;

    % Assemble full solution with boundary conditions
    u = zeros(N + 1, 1);
    u(2:N) = u_int;
    u(1) = 0;      % u(1) = 0
    u(N + 1) = 0;  % u(-1) = 0
end
