%% bc_helmholtz_robin.m - Helmholtz equation with Robin boundary conditions
%
% Chapter 13: Boundary Conditions in Spectral Methods
%
% Helmholtz equation with Robin boundary conditions (tau method / boundary
% bordering):
%
%   -u_xx + k^2 u = f(x),  x in (-1, 1)
%   u'(-1) - alpha * u(-1) = 0   (Robin at left)
%   u(1) = 0                      (Dirichlet at right)
%
% Parameters: k = 2, f(x) = exp(-x^2)
%
% Method II (boundary bordering): the full (N+1) x (N+1) system is built
% and boundary rows are overwritten:
%   Row 1   (x = 1):   identity row for Dirichlet u(1) = 0
%   Row N+1 (x = -1):  D(N+1,:) - alpha * I(N+1,:) for Robin condition
%
% Generates:
%   Figure 1: helmholtz_robin.pdf  - solutions for alpha = 0,1,5,20 + convergence
%   Figure 2: robin_conditioning.pdf - condition number vs alpha
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
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'ch07'));

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
NAVY   = [20,  45, 110] / 255;
SKY    = [120, 150, 210] / 255;
CORAL  = [231,  76,  60] / 255;
TEAL   = [22,  160, 133] / 255;
PURPLE = [142,  68, 173] / 255;
ORANGE = [230, 126,  34] / 255;

%% Output directory
script_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch13', 'matlab');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% Parameters
k = 2.0;

%% =====================================================================
%  Figure 1: Solutions for various alpha + convergence
%  =====================================================================
fig1 = figure('Units', 'inches', 'Position', [1, 1, 11, 4.5]);

% ---- Left panel: solutions for different alpha ----
subplot(1, 2, 1);

alpha_values = [0, 1, 5, 20];
colors_alpha = {NAVY, TEAL, CORAL, PURPLE};
N_plot = 32;

hold on;
for j = 1:length(alpha_values)
    alpha = alpha_values(j);
    [x, u] = solve_helmholtz_robin(N_plot, k, alpha);
    plot(x, u, '-o', 'Color', colors_alpha{j}, 'LineWidth', 1.5, ...
         'MarkerSize', 3, 'MarkerEdgeColor', 'w', 'MarkerFaceColor', colors_alpha{j}, ...
         'DisplayName', sprintf('$\\alpha = %d$', alpha));
end
hold off;

xlabel('$x$');
ylabel('$u(x)$');
title(sprintf('Solutions ($N = %d$, $k = %d$)', N_plot, k));
legend('Location', 'best', 'FontSize', 9);
xlim([-1.05, 1.05]);
grid on;
set(gca, 'GridAlpha', 0.3);
box on;

% ---- Right panel: convergence for alpha = 5 ----
subplot(1, 2, 2);

alpha_conv = 5.0;

% Reference solution with high N
N_ref = 100;
[x_ref, u_ref] = solve_helmholtz_robin(N_ref, k, alpha_conv);

N_values = 4:2:40;
errors = zeros(size(N_values));

for j = 1:length(N_values)
    N_val = N_values(j);
    [x_n, u_n] = solve_helmholtz_robin(N_val, k, alpha_conv);
    u_ref_n = bary_interp(x_ref, u_ref, x_n);
    errors(j) = max(max(abs(u_n - u_ref_n)), 1e-16);
end

semilogy(N_values, errors, 'o-', 'Color', CORAL, 'LineWidth', 1.5, ...
         'MarkerSize', 5, 'MarkerFaceColor', CORAL, 'MarkerEdgeColor', 'w');
hold on;
yline(2.2e-16, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
text(40, 4e-16, 'Machine precision', 'FontSize', 9, 'Color', [0.5 0.5 0.5], ...
     'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
hold off;

xlabel('$N$');
ylabel('Max Error');
title(sprintf('Convergence ($\\alpha = %d$)', alpha_conv));
xlim([0, 44]);
ylim([1e-16, 1e1]);
grid on;
set(gca, 'GridAlpha', 0.3);
box on;

sgtitle('Helmholtz with Robin BC: $-u_{xx} + k^2 u = e^{-x^2}$', 'FontSize', 12);

exportgraphics(fig1, fullfile(output_dir, 'helmholtz_robin.pdf'), ...
               'ContentType', 'vector');
fprintf('Figure 1 saved to: %s\n', fullfile(output_dir, 'helmholtz_robin.pdf'));
close(fig1);

%% =====================================================================
%  Figure 2: Condition number vs alpha
%  =====================================================================
fig2 = figure('Units', 'inches', 'Position', [1, 1, 6, 4.5]);

alpha_range = logspace(-2, 2, 80);
N_cond_values = [16, 32, 64];
colors_cond = {TEAL, CORAL, PURPLE};

hold on;
set(gca, 'XScale', 'log', 'YScale', 'log');
for j = 1:length(N_cond_values)
    N_c = N_cond_values(j);
    cond_numbers = zeros(size(alpha_range));

    for m = 1:length(alpha_range)
        alpha_val = alpha_range(m);
        [D, x] = cheb_matrix(N_c);
        D2 = D * D;
        I_mat = eye(N_c + 1);

        % Build operator: L = -D^2 + k^2 I
        L = -D2 + k^2 * I_mat;

        % Row 1 (x(1) = 1): Dirichlet u(1) = 0
        L(1, :) = 0.0;
        L(1, 1) = 1.0;

        % Row N+1 (x(N+1) = -1): Robin u'(-1) - alpha * u(-1) = 0
        L(N_c + 1, :) = D(N_c + 1, :) - alpha_val * I_mat(N_c + 1, :);

        cond_numbers(m) = cond(L);
    end

    loglog(alpha_range, cond_numbers, '-', 'Color', colors_cond{j}, ...
           'LineWidth', 1.5, 'DisplayName', sprintf('$N = %d$', N_c));
end
hold off;

xlabel('$\alpha$');
ylabel('$\mathrm{cond}(L)$');
title('Condition Number vs Robin Parameter $\alpha$', 'FontSize', 11);
legend('Location', 'best', 'FontSize', 10);
grid on;
set(gca, 'GridAlpha', 0.3);
box on;

exportgraphics(fig2, fullfile(output_dir, 'robin_conditioning.pdf'), ...
               'ContentType', 'vector');
fprintf('Figure 2 saved to: %s\n', fullfile(output_dir, 'robin_conditioning.pdf'));
close(fig2);

%% Print convergence table
fprintf('\nConvergence Table: Helmholtz with Robin BC (alpha = %g)\n', alpha_conv);
fprintf('%s\n', repmat('-', 1, 30));
fprintf('%6s  %14s\n', 'N', 'Max Error');
fprintf('%s\n', repmat('-', 1, 30));
for j = 1:length(N_values)
    fprintf('%6d  %14.2e\n', N_values(j), errors(j));
end
fprintf('%s\n', repmat('-', 1, 30));


%% =====================================================================
%  Local functions
%  =====================================================================

function [x, u] = solve_helmholtz_robin(N, k, alpha)
% SOLVE_HELMHOLTZ_ROBIN  Solve -u_xx + k^2*u = f(x) with Robin/Dirichlet BCs.
%
%   Row 1   (x = 1, Dirichlet):  L(1,:)   = e_1^T,       rhs(1)   = 0
%   Row N+1 (x = -1, Robin):     L(N+1,:) = D(N+1,:) - alpha*I(N+1,:), rhs(N+1) = 0

    [D, x] = cheb_matrix(N);
    D2 = D * D;
    I_mat = eye(N + 1);

    % Build operator: L = -D^2 + k^2 I
    L   = -D2 + k^2 * I_mat;
    rhs = exp(-x.^2);   % f(x) = exp(-x^2)

    % Row 1 (x(1) = 1): Dirichlet u(1) = 0
    L(1, :) = 0.0;
    L(1, 1) = 1.0;
    rhs(1)  = 0.0;

    % Row N+1 (x(N+1) = -1): Robin u'(-1) - alpha * u(-1) = 0
    L(N + 1, :) = D(N + 1, :) - alpha * I_mat(N + 1, :);
    rhs(N + 1)  = 0.0;

    % Solve
    u = L \ rhs;
end

function p = bary_interp(x_nodes, f_nodes, x_eval)
% BARY_INTERP  Barycentric interpolation at Chebyshev nodes.

    N = length(x_nodes) - 1;

    % Barycentric weights for Chebyshev-Gauss-Lobatto points
    w = ones(N + 1, 1);
    w(1)   = 0.5;
    w(N+1) = 0.5;
    w = w .* (-1).^(0:N)';

    p = zeros(size(x_eval));
    for j = 1:length(x_eval)
        xj = x_eval(j);
        diff = xj - x_nodes;

        [min_diff, idx] = min(abs(diff));
        if min_diff < 1e-15
            p(j) = f_nodes(idx);
        else
            terms = w ./ diff;
            p(j) = sum(terms .* f_nodes) / sum(terms);
        end
    end
end
