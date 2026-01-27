%% bvp_variable_coeff.m - Variable coefficient BVP with Chebyshev collocation
%
% Solves an Airy-type equation with variable coefficients:
%
%     u_xx - (1 + x^2)*u = 1,  x in (-1, 1),  u(+/-1) = 0
%
% This demonstrates that variable coefficient problems require no additional
% complexity with spectral methods - the variable coefficient becomes a
% diagonal matrix multiplication.
%
% This script generates Figure 7.4 for Chapter 7.
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

%% Exact solution for constant coefficient case
exact_constant_coeff = @(x) cosh(x) / cosh(1) - 1;

%% Solve both problems
N = 32;

% Variable coefficient: u_xx - (1 + x^2)*u = 1
[x_var, u_var] = solve_variable_coeff_bvp(N);

% Constant coefficient: u_xx - u = 1
[x_const, u_const] = solve_constant_coeff_bvp(N);
u_const_exact = exact_constant_coeff(x_const);

% Fine grid for exact solution
x_fine = linspace(-1, 1, 500)';
u_const_exact_fine = exact_constant_coeff(x_fine);

%% Create figure
fig = figure('Position', [100, 100, 1100, 450]);

% Panel 1: Compare solutions
subplot(1, 2, 1);
plot(x_var, u_var, 'o-', 'Color', TEAL, 'LineWidth', 1.5, 'MarkerSize', 4, ...
    'MarkerFaceColor', TEAL, 'MarkerEdgeColor', 'white', ...
    'DisplayName', 'Variable: $u_{xx} - (1+x^2)u = 1$');
hold on;
plot(x_const, u_const, 's-', 'Color', CORAL, 'LineWidth', 1.5, 'MarkerSize', 4, ...
    'MarkerFaceColor', CORAL, 'MarkerEdgeColor', 'white', ...
    'DisplayName', 'Constant: $u_{xx} - u = 1$');
hold off;

xlabel('$x$');
ylabel('$u(x)$');
title(sprintf('Comparison of Solutions ($N = %d$)', N));
legend('Location', 'south', 'FontSize', 9);
xlim([-1.05, 1.05]);
box off;
grid on;
set(gca, 'GridAlpha', 0.2);

% Add annotation
text(0.4, -0.4, {'Variable coefficient', 'reduces amplitude'}, ...
    'FontSize', 9, 'Color', [0.4, 0.4, 0.4]);

% Panel 2: Verify constant coefficient case
subplot(1, 2, 2);
plot(x_fine, u_const_exact_fine, '-', 'Color', NAVY, 'LineWidth', 1.5, ...
    'DisplayName', 'Exact');
hold on;
plot(x_const, u_const, 'o', 'Color', CORAL, 'MarkerSize', 6, ...
    'MarkerFaceColor', CORAL, 'MarkerEdgeColor', 'white', ...
    'DisplayName', 'Spectral');
hold off;

error_const = max(abs(u_const - u_const_exact));

xlabel('$x$');
ylabel('$u(x)$');
title(sprintf('Constant Coeff. Verification (error: %.2e)', error_const));
legend('Location', 'south', 'FontSize', 9);
xlim([-1.05, 1.05]);
box off;
grid on;
set(gca, 'GridAlpha', 0.2);

% Main title
sgtitle('Variable Coefficient BVP: $u_{xx} - (1+x^2)u = 1$, $u(\pm 1) = 0$', ...
    'FontSize', 13);

%% Save figure
output_file = fullfile(output_dir, 'variable_coeff');
exportgraphics(fig, [output_file, '.pdf'], 'ContentType', 'vector');
exportgraphics(fig, [output_file, '.png'], 'Resolution', 300);
fprintf('Figure saved to: %s.pdf\n', output_file);

%% Print solution properties
fprintf('\nSolution comparison:\n');
fprintf('%s\n', repmat('-', 1, 50));
idx_mid = floor((N + 1) / 2) + 1;  % Index for x = 0 approximately
fprintf('Variable coefficient u(0) = %.6f\n', u_var(idx_mid));
fprintf('Constant coefficient u(0) = %.6f\n', u_const(idx_mid));
fprintf('Constant coeff exact u(0) = %.6f\n', exact_constant_coeff(0));
fprintf('Numerical error (const coeff): %.2e\n', error_const);
fprintf('%s\n', repmat('-', 1, 50));

%% Local function: solve variable coefficient BVP
function [x, u] = solve_variable_coeff_bvp(N)
    % Get Chebyshev matrix and grid
    [D, x] = cheb_matrix(N);
    D2 = D * D;  % Second derivative matrix

    % Variable coefficient: a(x) = 1 + x^2
    a = 1 + x.^2;

    % Build the operator: L = D^2 - diag(a)
    L = D2 - diag(a);

    % Extract interior system
    L_int = L(2:N, 2:N);
    rhs_int = ones(N - 1, 1);  % f(x) = 1

    % Solve
    u_int = L_int \ rhs_int;

    % Assemble with boundary conditions
    u = zeros(N + 1, 1);
    u(2:N) = u_int;
end

%% Local function: solve constant coefficient BVP
function [x, u] = solve_constant_coeff_bvp(N)
    % Get Chebyshev matrix and grid
    [D, x] = cheb_matrix(N);
    D2 = D * D;  % Second derivative matrix

    % Constant coefficient: a = 1
    L = D2 - eye(N + 1);

    % Extract interior system
    L_int = L(2:N, 2:N);
    rhs_int = ones(N - 1, 1);

    % Solve
    u_int = L_int \ rhs_int;

    % Assemble with boundary conditions
    u = zeros(N + 1, 1);
    u(2:N) = u_int;
end
