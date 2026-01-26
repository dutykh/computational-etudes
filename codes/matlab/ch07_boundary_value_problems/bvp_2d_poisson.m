%% bvp_2d_poisson.m - 2D Poisson equation with tensor product Chebyshev methods
%
% Solves the 2D Poisson equation:
%
%     u_xx + u_yy = f(x,y),  (x,y) in (-1,1)^2,  u = 0 on boundary
%
% Test problem: f(x,y) = -2*pi^2 * sin(pi*x) * sin(pi*y)
% Exact solution: u(x,y) = sin(pi*x) * sin(pi*y)
%
% This demonstrates tensor product grids and Kronecker product operators.
%
% This script generates Figures 7.7, 7.8, and 7.9 for Chapter 7.
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
rhs_function = @(X, Y) -2 * pi^2 * sin(pi * X) .* sin(pi * Y);
exact_solution = @(X, Y) sin(pi * X) .* sin(pi * Y);

%% Solve problem
N = 16;
[X, Y, U, L] = solve_2d_poisson(N, rhs_function);
U_exact = exact_solution(X, Y);
error_val = max(max(abs(U - U_exact)));

%% Figure 7.7: Tensor product grid
fig1 = figure('Position', [100, 100, 600, 600]);

% Plot grid points
plot(X(:), Y(:), 'o', 'Color', TEAL, 'MarkerSize', 5, ...
    'MarkerFaceColor', TEAL, 'MarkerEdgeColor', 'white');
hold on;

% Add grid lines
for i = 1:N+1
    plot(X(i, :), Y(i, :), '-', 'Color', SKY, 'LineWidth', 0.5);
    plot(X(:, i), Y(:, i), '-', 'Color', SKY, 'LineWidth', 0.5);
end
hold off;

xlabel('$x$');
ylabel('$y$');
title(sprintf('Chebyshev Tensor Product Grid ($N = %d$)', N));
xlim([-1.1, 1.1]);
ylim([-1.1, 1.1]);
axis equal;
box off;

% Add annotation about clustering
text(0.3, 0.5, {'Grid clusters', 'near boundaries'}, ...
    'FontSize', 9, 'Color', [0.4, 0.4, 0.4]);

output_file = fullfile(output_dir, 'tensor_grid');
exportgraphics(fig1, [output_file, '.pdf'], 'ContentType', 'vector');
exportgraphics(fig1, [output_file, '.png'], 'Resolution', 300);
fprintf('Figure saved to: %s.pdf\n', output_file);

%% Figure 7.8: 2D Poisson solution
fig2 = figure('Position', [100, 100, 1100, 500]);

% 3D surface plot
subplot(1, 2, 1);

% Create finer grid for plotting
x_fine = linspace(-1, 1, 100);
[X_fine, Y_fine] = meshgrid(x_fine, x_fine);
U_exact_fine = exact_solution(X_fine, Y_fine);

surf(X_fine, Y_fine, U_exact_fine, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
hold on;
scatter3(X(:), Y(:), U(:), 20, CORAL, 'filled', 'MarkerEdgeColor', 'white');
hold off;

colormap(viridis(256));
xlabel('$x$');
ylabel('$y$');
zlabel('$u(x,y)$');
title(sprintf('Solution (max error: %.2e)', error_val));
view(25, 45);

% Contour plot
subplot(1, 2, 2);
levels = linspace(-1, 1, 21);
contourf(X_fine, Y_fine, U_exact_fine, levels);
hold on;
contour(X_fine, Y_fine, U_exact_fine, levels, 'LineColor', 'white', ...
    'LineWidth', 0.5);

% Mark numerical points
plot(X(:), Y(:), 'o', 'Color', CORAL, 'MarkerSize', 3, ...
    'MarkerFaceColor', CORAL, 'MarkerEdgeColor', 'none');
hold off;

colormap(viridis(256));
colorbar('Label', '$u(x,y)$');
xlabel('$x$');
ylabel('$y$');
title('Contour Plot with Grid Points');
axis equal;

sgtitle('2D Poisson: $\nabla^2 u = -2\pi^2 \sin(\pi x)\sin(\pi y)$', 'FontSize', 13);

output_file = fullfile(output_dir, 'poisson_2d');
exportgraphics(fig2, [output_file, '.pdf'], 'ContentType', 'vector');
exportgraphics(fig2, [output_file, '.png'], 'Resolution', 300);
fprintf('Figure saved to: %s.pdf\n', output_file);

%% Figure 7.9: Laplacian sparsity pattern
fig3 = figure('Position', [100, 100, 600, 600]);

n_int = N - 1;
spy(L, 1);
colormap([NAVY; NAVY]);

title(sprintf('2D Laplacian Structure ($N = %d$, size $%d \\times %d$)', N, n_int^2, n_int^2));
xlabel('Column index');
ylabel('Row index');

% Add annotation
text(0.5, -0.08, '$L = I \otimes D^2 + D^2 \otimes I$', ...
    'Units', 'normalized', 'FontSize', 11, 'HorizontalAlignment', 'center');

output_file = fullfile(output_dir, 'laplacian_sparsity');
exportgraphics(fig3, [output_file, '.pdf'], 'ContentType', 'vector');
exportgraphics(fig3, [output_file, '.png'], 'Resolution', 300);
fprintf('Figure saved to: %s.pdf\n', output_file);

%% Print convergence table
fprintf('\nConvergence Table: 2D Poisson Equation\n');
fprintf('%s\n', repmat('-', 1, 50));
fprintf('%6s %12s %14s\n', 'N', 'Grid Size', 'Max Error');
fprintf('%s\n', repmat('-', 1, 50));

for N_val = [4, 8, 12, 16, 20, 24, 32]
    [X_n, Y_n, U_n, ~] = solve_2d_poisson(N_val, rhs_function);
    U_exact_n = exact_solution(X_n, Y_n);
    err = max(max(abs(U_n - U_exact_n)));
    grid_size = (N_val + 1)^2;
    fprintf('%6d %12d %14.2e\n', N_val, grid_size, err);
end

fprintf('%s\n', repmat('-', 1, 50));

%% Local function: solve 2D Poisson equation
function [X, Y, U, L] = solve_2d_poisson(N, rhs_func)
    % Get Chebyshev matrix and grid
    [D, x] = cheb_matrix(N);
    D2 = D * D;  % Second derivative matrix

    % Interior points only
    n_int = N - 1;
    D2_int = D2(2:N, 2:N);
    x_int = x(2:N);

    % Build 2D Laplacian using Kronecker products
    % L = I kron D^2 + D^2 kron I
    I = eye(n_int);
    L = kron(I, D2_int) + kron(D2_int, I);

    % Create meshgrid for interior points
    [X_int, Y_int] = meshgrid(x_int, x_int);

    % Right-hand side
    F = rhs_func(X_int, Y_int);
    f_vec = F(:);  % Vectorize column-major (MATLAB default)

    % Solve the linear system
    u_vec = L \ f_vec;

    % Reshape solution
    U_int = reshape(u_vec, [n_int, n_int]);

    % Embed in full grid with boundary conditions
    U = zeros(N + 1, N + 1);
    U(2:N, 2:N) = U_int;

    % Create full meshgrid
    [X, Y] = meshgrid(x, x);
end

%% Viridis colormap (for compatibility)
function cmap = viridis(n)
    % Approximate viridis colormap
    if nargin < 1
        n = 256;
    end

    % Key colors from viridis
    c = [0.267004, 0.004874, 0.329415;
         0.282327, 0.140926, 0.457517;
         0.253935, 0.265254, 0.529983;
         0.206756, 0.371758, 0.553117;
         0.163625, 0.471133, 0.558148;
         0.127568, 0.566949, 0.550556;
         0.134692, 0.658636, 0.517649;
         0.266941, 0.748751, 0.440573;
         0.477504, 0.821444, 0.318195;
         0.741388, 0.873449, 0.149561;
         0.993248, 0.906157, 0.143936];

    t = linspace(0, 1, size(c, 1));
    t_interp = linspace(0, 1, n);
    cmap = interp1(t, c, t_interp);
end
