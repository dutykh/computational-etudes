%% bc_vibrating_string.m - Vibrating string with Dirichlet-Neumann BCs
%
% Chapter 13: Boundary Conditions in Spectral Methods
%
% Eigenvalue problem for a vibrating string with mixed BCs:
%
%   -u_xx = lambda * u,   x in (0, 1)
%   u(0) = 0              (Dirichlet, fixed end)
%   u'(1) = 0             (Neumann, free end)
%
% Mapping to Chebyshev domain [-1, 1]: x = (xi + 1)/2, D_phys = 2 * D_xi
%
% In Chebyshev ordering:
%   xi(1) = +1  -> physical x = 1 (Neumann end)
%   xi(N+1) = -1 -> physical x = 0 (Dirichlet end)
%
% The eigenvalue problem A * u = lambda * u is set up as:
%   A = -D_phys^2 (interior)
%   Row 1: D_phys(1,:)       (Neumann at physical x = 1)
%   Row N+1: identity row    (Dirichlet at physical x = 0)
%
% Exact eigenvalues: lambda_n = ((2n-1)*pi/2)^2,  n = 1, 2, ...
% Exact eigenmodes:  u_n(x) = sin((2n-1)*pi*x/2)
%
% Generates Figure: vibrating_string.pdf (2 panels)
%   Left:  first 6 eigenmodes (computed vs exact)
%   Right: eigenvalue errors vs mode number
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
N = 36;

%% Setup grid and operators
[D_xi, xi] = cheb_matrix(N);

% Physical derivative: D_phys = 2 * D_xi (mapping [0,1] -> [-1,1])
D_phys = 2 * D_xi;
D2_phys = D_phys * D_phys;

% Physical coordinates
x_phys = (xi + 1) / 2;

%% Build eigenvalue problem matrix
A = -D2_phys;

% Row 1 (xi = +1, physical x = 1): Neumann u'(1) = 0
A(1, :) = D_phys(1, :);

% Row N+1 (xi = -1, physical x = 0): Dirichlet u(0) = 0
A(N + 1, :)     = 0;
A(N + 1, N + 1) = 1;

%% Construct corresponding B matrix for generalized problem A*u = lambda*B*u
B = eye(N + 1);
B(1, :)   = 0;    % BC row: not an eigenvalue equation
B(N+1, :) = 0;    % BC row: not an eigenvalue equation

%% Solve generalized eigenvalue problem
[V, D_eig] = eig(A, B);
eigenvalues = diag(D_eig);

% Filter out infinite and NaN eigenvalues (from zero rows in B)
valid = isfinite(eigenvalues) & (real(eigenvalues) > 0);
eigenvalues = real(eigenvalues(valid));
V = V(:, valid);

% Sort by eigenvalue
[eigenvalues, sort_idx] = sort(eigenvalues);
V = V(:, sort_idx);

%% Exact eigenvalues and eigenfunctions
n_modes = min(20, length(eigenvalues));
n_exact = (1:n_modes)';
lambda_exact = ((2 * n_exact - 1) * pi / 2).^2;

% Fine grid for exact eigenfunctions
x_fine = linspace(0, 1, 200)';

%% Print eigenvalue comparison
fprintf('Eigenvalue Comparison: Vibrating String with Free End\n');
fprintf('%s\n', repmat('-', 1, 60));
fprintf('%5s  %16s  %16s  %12s\n', 'Mode', 'Computed', 'Exact', 'Error');
fprintf('%s\n', repmat('-', 1, 60));
for k = 1:min(12, n_modes)
    err = abs(eigenvalues(k) - lambda_exact(k));
    fprintf('%5d  %16.10f  %16.10f  %12.2e\n', k, eigenvalues(k), lambda_exact(k), err);
end
fprintf('%s\n', repmat('-', 1, 60));

%% Create figure
fig = figure('Units', 'inches', 'Position', [1, 1, 11, 4.5]);

% ---- Left panel: first 6 eigenmodes ----
subplot(1, 2, 1);

colors_mode = {NAVY, TEAL, CORAL, PURPLE, ORANGE, SKY};
n_show = 6;

hold on;
for k = 1:n_show
    % Computed eigenmode
    psi = real(V(:, k));
    % Normalize: match sign of exact mode
    u_exact_k = sin((2*k - 1) * pi * x_phys / 2);
    if sum(psi .* u_exact_k) < 0
        psi = -psi;
    end
    psi = psi / max(abs(psi));

    % Exact eigenmode on fine grid
    u_exact_fine = sin((2*k - 1) * pi * x_fine / 2);

    % Plot exact (thin line) and computed (circles)
    plot(x_fine, u_exact_fine, '-', 'Color', colors_mode{k}, 'LineWidth', 0.8, ...
         'HandleVisibility', 'off');
    plot(x_phys, psi, 'o', 'Color', colors_mode{k}, 'MarkerSize', 4, ...
         'MarkerFaceColor', colors_mode{k}, 'MarkerEdgeColor', 'w', ...
         'LineWidth', 0.3, 'DisplayName', sprintf('Mode %d', k));
end
hold off;

xlabel('$x$');
ylabel('$u_n(x)$');
title('Eigenmodes (lines: exact, circles: spectral)', 'FontSize', 10);
legend('Location', 'best', 'FontSize', 8, 'NumColumns', 2);
xlim([0, 1]);
ylim([-1.15, 1.15]);
grid on;
set(gca, 'GridAlpha', 0.3);
box on;

% ---- Right panel: eigenvalue errors ----
subplot(1, 2, 2);

errors = abs(eigenvalues(1:n_modes) - lambda_exact);
errors = max(errors, 1e-16);

semilogy(1:n_modes, errors, 'o-', 'Color', NAVY, 'LineWidth', 1.5, ...
         'MarkerSize', 5, 'MarkerFaceColor', NAVY, 'MarkerEdgeColor', 'w');
hold on;
yline(2.2e-16, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
text(n_modes, 4e-16, 'Machine precision', 'FontSize', 9, 'Color', [0.5 0.5 0.5], ...
     'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
hold off;

xlabel('Mode number $n$');
ylabel('$|\lambda_n - \lambda_n^{\mathrm{exact}}|$');
title(sprintf('Eigenvalue Errors ($N = %d$)', N), 'FontSize', 11);
xlim([0, n_modes + 1]);
grid on;
set(gca, 'GridAlpha', 0.3);
box on;

sgtitle('Vibrating String: $-u_{xx} = \lambda u$, $u(0) = 0$, $u''(1) = 0$', ...
        'FontSize', 12);

%% Save figure
exportgraphics(fig, fullfile(output_dir, 'vibrating_string.pdf'), ...
               'ContentType', 'vector');
fprintf('\nFigure saved to: %s\n', fullfile(output_dir, 'vibrating_string.pdf'));

close(fig);
