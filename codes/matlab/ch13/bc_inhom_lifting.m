%% bc_inhom_lifting.m - Inhomogeneous Poisson equation with lifting
%
% Chapter 13: Boundary Conditions in Spectral Methods
%
% Solves the inhomogeneous Poisson equation via the lifting technique:
%
%   u_xx = exp(4x),  x in (-1, 1),  u(-1) = 0,  u(1) = 1
%
% The lifting function w(x) = (x+1)/2 satisfies the boundary conditions:
%   w(-1) = 0,  w(1) = 1
%
% Since w'' = 0, the homogeneous variable v = u - w satisfies:
%   v_xx = exp(4x),  v(-1) = 0,  v(1) = 0
%
% After solving for v with homogeneous Dirichlet BCs, we reconstruct u = v + w.
%
% Generates Figure: inhom_lifting.pdf
%   Left panel:  solution (circles) vs reference (line) for N = 16
%   Right panel: convergence (max error vs N) for N = 4, 6, ..., 40
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

%% Reference solution (high resolution)
N_ref = 100;
[x_ref, u_ref] = solve_poisson_lifting(N_ref);

% Fine grid for smooth plotting
x_fine = linspace(-1, 1, 500)';
u_fine = bary_interp(x_ref, u_ref, x_fine);

%% Panel 1 data: N = 16 solution
N_plot = 16;
[x16, u16] = solve_poisson_lifting(N_plot);
u_ref_on16 = bary_interp(x_ref, u_ref, x16);
error16 = max(abs(u16 - u_ref_on16));

%% Panel 2 data: convergence study
N_values = 4:2:40;
errors = zeros(size(N_values));

for k = 1:length(N_values)
    N_val = N_values(k);
    [x_n, u_n] = solve_poisson_lifting(N_val);
    u_ref_n = bary_interp(x_ref, u_ref, x_n);
    errors(k) = max(max(abs(u_n - u_ref_n)), 1e-16);
end

%% Create figure
fig = figure('Units', 'inches', 'Position', [1, 1, 11, 4.5]);

% ---- Left panel: solution for N = 16 ----
subplot(1, 2, 1);
plot(x_fine, u_fine, '-', 'Color', NAVY, 'LineWidth', 1.5, ...
     'DisplayName', 'Reference');
hold on;
plot(x16, u16, 'o', 'Color', TEAL, 'MarkerSize', 6, ...
     'MarkerFaceColor', TEAL, 'MarkerEdgeColor', 'w', ...
     'LineWidth', 0.5, 'DisplayName', 'Spectral');
hold off;
xlabel('$x$');
ylabel('$u(x)$');
title(sprintf('Solution ($N = %d$, max error: %.2e)', N_plot, error16));
legend('Location', 'northwest', 'FontSize', 9);
xlim([-1.05, 1.05]);
grid on;
set(gca, 'GridAlpha', 0.3);
box on;

% ---- Right panel: convergence ----
subplot(1, 2, 2);
semilogy(N_values, errors, 'o-', 'Color', CORAL, 'LineWidth', 1.5, ...
         'MarkerSize', 5, 'MarkerFaceColor', CORAL, ...
         'MarkerEdgeColor', 'w');
hold on;
yline(2.2e-16, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
text(40, 4e-16, 'Machine precision', 'FontSize', 9, 'Color', [0.5 0.5 0.5], ...
     'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
hold off;
xlabel('$N$');
ylabel('Max Error');
title('Convergence');
xlim([0, 44]);
ylim([1e-16, 1e1]);
grid on;
set(gca, 'GridAlpha', 0.3);
box on;

sgtitle('Inhomogeneous Poisson: $u_{xx} = e^{4x}$, $u(-1)=0$, $u(1)=1$ (Lifting)', ...
        'FontSize', 12);

%% Save figure
exportgraphics(fig, fullfile(output_dir, 'inhom_lifting.pdf'), ...
               'ContentType', 'vector');

fprintf('Figure saved to: %s\n', fullfile(output_dir, 'inhom_lifting.pdf'));

%% Print convergence table
fprintf('\nConvergence Table: Inhomogeneous Poisson with Lifting\n');
fprintf('%s\n', repmat('-', 1, 30));
fprintf('%6s  %14s\n', 'N', 'Max Error');
fprintf('%s\n', repmat('-', 1, 30));
for k = 1:length(N_values)
    fprintf('%6d  %14.2e\n', N_values(k), errors(k));
end
fprintf('%s\n', repmat('-', 1, 30));

close(fig);


%% =====================================================================
%  Local functions
%  =====================================================================

function [x, u] = solve_poisson_lifting(N)
% SOLVE_POISSON_LIFTING  Solve u_xx = exp(4x) with u(-1)=0, u(1)=1 via lifting.
%
%   Lifting: w(x) = (x+1)/2, v = u - w satisfies v_xx = exp(4x), v(+-1) = 0.
%   Since w'' = 0, the RHS is unchanged.

    [D, x] = cheb_matrix(N);
    D2 = D * D;

    % RHS: f(x) = exp(4x)
    f = exp(4.0 * x);

    % Interior system (strip boundary rows): rows 2..N, cols 2..N
    D2_int = D2(2:N, 2:N);
    f_int  = f(2:N);

    % Solve v_xx = f with homogeneous Dirichlet BCs
    v_int = D2_int \ f_int;

    % Assemble full v
    v    = zeros(N + 1, 1);
    v(2:N) = v_int;

    % Reconstruct u = v + w(x)
    w = (x + 1.0) / 2.0;
    u = v + w;
end

function p = bary_interp(x_nodes, f_nodes, x_eval)
% BARY_INTERP  Barycentric interpolation at Chebyshev nodes.
%
%   Uses the barycentric formula with Chebyshev weights.

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

        % Check if xj coincides with a node
        [min_diff, idx] = min(abs(diff));
        if min_diff < 1e-15
            p(j) = f_nodes(idx);
        else
            terms = w ./ diff;
            p(j) = sum(terms .* f_nodes) / sum(terms);
        end
    end
end
