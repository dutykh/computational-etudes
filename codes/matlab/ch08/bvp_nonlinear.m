%% bvp_nonlinear.m - Bratu equation (nonlinear BVP) using Newton iteration
%
% Solves the Bratu equation:
%
%     u_xx + lambda * exp(u) = 0,  x in (-1, 1),  u(+/-1) = 0
%
% with lambda = 0.5 (subcritical case, unique solution exists).
%
% IMPORTANT: The critical value for [-1, 1] domain is lambda_c ~ 0.878.
% For lambda > lambda_c, no solution exists! We use lambda = 0.5 for stability.
%
% This demonstrates that nonlinear BVPs are easily handled with spectral
% methods using Newton iteration.
%
% This script generates Figure 7.5 and Table 7.1 for Chapter 7.
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

%% Create figure
fig = figure('Position', [100, 100, 1100, 450]);

% Panel 1: Solution for various N
subplot(1, 2, 1);

N_values = [8, 16, 32];
colors = {CORAL, TEAL, PURPLE};

for idx = 1:length(N_values)
    N = N_values(idx);
    color = colors{idx};
    [x, u, residuals, n_iter] = solve_bratu_newton(N);
    plot(x, u, 'o-', 'Color', color, 'LineWidth', 1.2, 'MarkerSize', 4, ...
        'MarkerFaceColor', color, 'MarkerEdgeColor', 'white', ...
        'DisplayName', sprintf('$N = %d$ (%d iter)', N, n_iter));
    hold on;
end
hold off;

xlabel('$x$');
ylabel('$u(x)$');
title('Bratu Equation: $u_{xx} + e^u = 0$');
legend('Location', 'south', 'FontSize', 9, 'NumColumns', 3);
xlim([-1.05, 1.05]);
box off;
grid on;
set(gca, 'GridAlpha', 0.2);

% Panel 2: Convergence history
subplot(1, 2, 2);

% Compare fixed-point vs Newton for N=16
N = 16;
[~, ~, res_fp, n_fp] = solve_bratu_fixedpoint(N, 0.5, 1e-14, 50);
[~, ~, res_newton, n_newton] = solve_bratu_newton(N, 0.5, 1e-14, 20);

semilogy(1:length(res_fp), res_fp, 'o-', 'Color', CORAL, 'LineWidth', 1.5, ...
    'MarkerSize', 5, 'MarkerFaceColor', CORAL, 'DisplayName', 'Fixed-point');
hold on;
semilogy(1:length(res_newton), res_newton, 's-', 'Color', TEAL, 'LineWidth', 1.5, ...
    'MarkerSize', 5, 'MarkerFaceColor', TEAL, 'DisplayName', 'Newton');

% Tolerance line
yline(1e-12, '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1);
text(2, 3e-12, 'tol', 'FontSize', 9, 'Color', [0.5, 0.5, 0.5]);
hold off;

xlabel('Iteration');
ylabel('Residual norm');
title(sprintf('Convergence History ($N = %d$)', N));
legend('Location', 'northeast', 'FontSize', 9);
grid on;
box off;
ylim([1e-16, 10]);

% Main title
sgtitle('Nonlinear BVP: Bratu Equation with $\lambda = 0.5$', 'FontSize', 13);

%% Save figure
output_file = fullfile(output_dir, 'bratu');
exportgraphics(fig, [output_file, '.pdf'], 'ContentType', 'vector');
exportgraphics(fig, [output_file, '.png'], 'Resolution', 300);
fprintf('Figure saved to: %s.pdf\n', output_file);

%% Generate convergence table
fprintf('\n%s\n', repmat('=', 1, 60));
fprintf('Table 7.1: Bratu Equation Convergence (Newton)\n');
fprintf('%s\n', repmat('=', 1, 60));
fprintf('%6s %12s %16s\n', 'N', 'Iterations', 'u(0)');
fprintf('%s\n', repmat('-', 1, 60));

for N = [4, 8, 16, 32, 64]
    [x, u, ~, n_iter] = solve_bratu_newton(N);
    % Find index closest to x=0
    [~, idx_0] = min(abs(x));
    fprintf('%6d %12d %16.10f\n', N, n_iter, u(idx_0));
end

fprintf('%s\n', repmat('=', 1, 60));
fprintf('\nNote: Newton iteration converges quadratically (in 5-8 iterations)\n');
fprintf('while fixed-point converges linearly (in ~18 iterations).\n');

%% Local function: solve Bratu equation with Newton iteration
function [x, u, residuals, n_iter] = solve_bratu_newton(N, lam, tol, max_iter)
    % Default parameters
    if nargin < 2, lam = 0.5; end
    if nargin < 3, tol = 1e-12; end
    if nargin < 4, max_iter = 30; end

    % Get Chebyshev matrix and grid
    [D, x] = cheb_matrix(N);
    D2 = D * D;  % Second derivative matrix

    % Extract interior system
    D2_int = D2(2:N, 2:N);

    % Initial guess: u = 0
    u = zeros(N + 1, 1);

    residuals = [];

    for it = 1:max_iter
        % Residual: F(u) = D^2*u + lam*exp(u) (interior points)
        exp_u = exp(u(2:N));
        F = D2_int * u(2:N) + lam * exp_u;

        res_norm = max(abs(F));
        residuals(end + 1) = res_norm;

        if res_norm < tol
            n_iter = it;
            return;
        end

        % Jacobian: J = D^2 + lam*diag(exp(u))
        J = D2_int + lam * diag(exp_u);

        % Newton step
        delta_u = J \ (-F);

        % Backtracking line search
        alpha = 1.0;
        u_trial = u;
        for bt = 1:10
            u_trial(2:N) = u(2:N) + alpha * delta_u;
            exp_u_trial = exp(u_trial(2:N));
            F_trial = D2_int * u_trial(2:N) + lam * exp_u_trial;
            res_trial = max(abs(F_trial));

            if res_trial < res_norm
                break;
            end
            alpha = alpha * 0.5;
        end

        % Update interior points
        u(2:N) = u(2:N) + alpha * delta_u;
    end

    n_iter = max_iter;
end

%% Local function: solve Bratu equation with fixed-point iteration
function [x, u, residuals, n_iter] = solve_bratu_fixedpoint(N, lam, tol, max_iter)
    % Default parameters
    if nargin < 2, lam = 0.5; end
    if nargin < 3, tol = 1e-10; end
    if nargin < 4, max_iter = 100; end

    % Get Chebyshev matrix and grid
    [D, x] = cheb_matrix(N);
    D2 = D * D;
    D2_int = D2(2:N, 2:N);

    % Initial guess
    u = zeros(N + 1, 1);

    residuals = [];

    for it = 1:max_iter
        % RHS
        rhs = -lam * exp(u(2:N));

        % Solve
        u_new_int = D2_int \ rhs;

        % Compute change
        change = max(abs(u_new_int - u(2:N)));
        residuals(end + 1) = change;

        u(2:N) = u_new_int;

        if change < tol
            n_iter = it;
            return;
        end
    end

    n_iter = max_iter;
end
