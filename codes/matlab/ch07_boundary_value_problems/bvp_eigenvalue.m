%% bvp_eigenvalue.m - Eigenvalue problem with Chebyshev collocation
%
% Solves the eigenvalue problem:
%
%     u_xx = lambda * u,  x in (-1, 1),  u(+/-1) = 0
%
% Exact eigenvalues: lambda_n = -(n*pi/2)^2 for n = 1, 2, 3, ...
% Exact eigenfunctions: u_n(x) = sin(n*pi*(x+1)/2)
%
% This demonstrates the resolution limits of spectral methods - higher modes
% require more points per wavelength (ppw) for accuracy.
%
% This script generates Figure 7.6 for Chapter 7.
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

colors = {CORAL, TEAL, PURPLE, ORANGE, SKY, NAVY};

%% Output directory
script_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch07', 'matlab');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% Exact solutions
exact_eigenvalue = @(n) -(n * pi / 2)^2;
exact_eigenfunction = @(x, n) sin(n * pi * (x + 1) / 2);

%% Solve eigenvalue problem
N = 36;
[eigenvalues, eigenvectors, x] = solve_eigenvalue_problem(N);

%% Create figure with 6 panels (2x3)
fig = figure('Position', [100, 100, 1200, 700]);

% Modes to display (1-indexed)
modes = [5, 10, 15, 20, 25, 30];

% Fine grid for exact solutions
x_fine = linspace(-1, 1, 500)';

for idx = 1:length(modes)
    subplot(2, 3, idx);

    mode = modes(idx);
    color = colors{idx};

    % Numerical eigenvalue and eigenfunction
    lam_num = eigenvalues(mode);
    u_num = eigenvectors(:, mode);

    % Exact values
    lam_exact = exact_eigenvalue(mode);
    u_exact_fine = exact_eigenfunction(x_fine, mode);
    u_exact_grid = exact_eigenfunction(x(2:N), mode);

    % Normalize numerical eigenfunction to match sign of exact
    if sign(u_num(1)) ~= sign(u_exact_grid(1))
        u_num = -u_num;
    end

    % Normalize amplitude
    u_num = u_num / max(abs(u_num));
    u_exact_fine = u_exact_fine / max(abs(u_exact_fine));

    % Eigenvalue error
    lam_error = abs(lam_num - lam_exact) / abs(lam_exact);

    % Points per wavelength
    wavelength = 4 / mode;  % wavelength = 2L/n = 4/n for L=2
    avg_spacing = 2 / N;
    ppw = wavelength / avg_spacing;

    % Plot
    plot(x_fine, u_exact_fine, '-', 'Color', NAVY, 'LineWidth', 1.5, ...
        'DisplayName', 'Exact');
    hold on;
    plot(x(2:N), u_num, 'o', 'Color', color, 'MarkerSize', 4, ...
        'MarkerFaceColor', color, 'MarkerEdgeColor', 'white', ...
        'DisplayName', 'Spectral');
    hold off;

    xlabel('$x$');
    ylabel('$u(x)$');
    title(sprintf('Mode %d: ppw = %.1f, $\\lambda$ error = %.1e', mode, ppw, lam_error));
    xlim([-1.05, 1.05]);
    ylim([-1.3, 1.3]);
    box off;

    if idx == 1
        legend('Location', 'northeast', 'FontSize', 8);
    end

    % Add reference line
    yline(0, '-', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5);
end

% Main title
sgtitle(sprintf('Eigenvalue Problem: $u_{xx} = \\lambda u$, $u(\\pm 1) = 0$ (N = %d)', N), ...
    'FontSize', 13);

%% Save figure
output_file = fullfile(output_dir, 'eigenvalue_problem');
exportgraphics(fig, [output_file, '.pdf'], 'ContentType', 'vector');
exportgraphics(fig, [output_file, '.png'], 'Resolution', 300);
fprintf('Figure saved to: %s.pdf\n', output_file);

%% Print eigenvalue table
fprintf('\nEigenvalue accuracy for N = %d:\n', N);
fprintf('%s\n', repmat('-', 1, 70));
fprintf('%6s %14s %14s %12s %8s\n', 'Mode', 'lambda (exact)', 'lambda (num)', 'Rel. Error', 'ppw');
fprintf('%s\n', repmat('-', 1, 70));

for mode = 1:min(N - 1, 35)
    lam_exact = exact_eigenvalue(mode);
    lam_num = eigenvalues(mode);
    rel_error = abs(lam_num - lam_exact) / abs(lam_exact);
    ppw = (4 / mode) / (2 / N);

    fprintf('%6d %14.4f %14.4f %12.2e %8.2f\n', mode, lam_exact, lam_num, rel_error, ppw);
end

fprintf('%s\n', repmat('-', 1, 70));
fprintf('\nNote: Accuracy degrades when ppw < pi (approx 3.14).\n');
fprintf('The rule of thumb is: need >= pi points per wavelength for spectral accuracy.\n');

%% Local function: solve eigenvalue problem
function [eigenvalues, eigenvectors, x] = solve_eigenvalue_problem(N)
    % Get Chebyshev matrix and grid
    [D, x] = cheb_matrix(N);
    D2 = D * D;  % Second derivative matrix

    % Extract interior system
    D2_int = D2(2:N, 2:N);

    % Solve eigenvalue problem
    [V, Lambda] = eig(D2_int);
    eigenvalues = diag(Lambda);

    % Sort by eigenvalue (largest = least negative first)
    [eigenvalues, idx] = sort(eigenvalues, 'descend');
    eigenvalues = real(eigenvalues);
    eigenvectors = real(V(:, idx));
end
