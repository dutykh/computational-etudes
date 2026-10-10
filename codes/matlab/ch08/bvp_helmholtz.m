%% bvp_helmholtz.m - Helmholtz equation with Chebyshev collocation
%
% Solves the Helmholtz equation:
%
%     u_xx + u_yy + k^2 * u = f(x,y),  (x,y) in (-1,1)^2,  u = 0 on boundary
%
% with localized Gaussian forcing:
%     f(x,y) = exp(-20*[(x-0.3)^2 + (y+0.4)^2])
%
% Uses k = 7 to demonstrate near-resonance behavior with the (2,4) eigenmode.
% (Eigenvalue for mode (m,n) is k^2 = (pi/2)^2*(m^2 + n^2), so (2,4) gives k ~ 7.02)
%
% This script generates Figure 7.10 for Chapter 7.
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

%% Forcing function
forcing_function = @(X, Y) exp(-20 * ((X - 0.3).^2 + (Y + 0.4).^2));

%% Theoretical resonance wavenumbers
fprintf('Resonance wavenumbers for first few modes:\n');
fprintf('%s\n', repmat('-', 1, 50));
fprintf('%15s %12s\n', 'Mode (m,n)', 'k_mn');
fprintf('%s\n', repmat('-', 1, 50));
for m = 1:5
    for n = 1:5
        k_mn = (pi / 2) * sqrt(m^2 + n^2);
        if k_mn < 10
            fprintf('%15s %12.4f\n', sprintf('(%d,%d)', m, n), k_mn);
        end
    end
end
fprintf('%s\n', repmat('-', 1, 50));

%% Solve Helmholtz equation at the showcase k = 7
N = 32;
k = 7.0;  % Near resonance with (2,4) mode

[X, Y, U] = solve_helmholtz(N, k, forcing_function);

%% NEW (panel c): resonance amplitude sweep |u(force)| vs k
[~, idx_x] = min(abs(X(1, :) - 0.3));
[~, idx_y] = min(abs(Y(:, 1) - (-0.4)));
k_axis = linspace(3.0, 12.0, 60);
amp_at_force = zeros(size(k_axis));
fprintf('\nSweeping k...\n');
for ki = 1:numel(k_axis)
    [~, ~, U_test] = solve_helmholtz(N, k_axis(ki), forcing_function);
    amp_at_force(ki) = abs(U_test(idx_y, idx_x));
end
k_24 = (pi / 2) * sqrt(4 + 16);  % sqrt(20) approx 7.0248

%% NEW (panel d): the (2,4) Dirichlet eigenmode
x_fine = linspace(-1, 1, 200);
y_fine = linspace(-1, 1, 200);
[X_fine, Y_fine] = meshgrid(x_fine, y_fine);
eigenmode_24 = sin(2 * pi * (X_fine + 1) / 2) .* sin(4 * pi * (Y_fine + 1) / 2);

%% Create 2x2 figure
fig = figure('Position', [100, 100, 1100, 900]);

% Panel (a): 3D surface
subplot(2, 2, 1);
surf(X, Y, U, 'EdgeColor', 'none', 'FaceAlpha', 0.9);
colormap(redblue(256));
xlabel('$x$'); ylabel('$y$'); zlabel('$u(x,y)$');
title(sprintf('(a) Solution surface ($k = %.1f$)', k));
view(25, 45);

% Panel (b): Contour plot of solution
subplot(2, 2, 2);
x_sorted = sort(X(1, :));
y_sorted = sort(Y(:, 1));
[~, sort_x_idx] = sort(X(1, :));
[~, sort_y_idx] = sort(Y(:, 1));
U_sorted = U(sort_y_idx, sort_x_idx);
U_fine = interp2(x_sorted, y_sorted, U_sorted, X_fine, Y_fine, 'spline');
vmax = max(abs(U_fine(:)));
levels = linspace(-vmax, vmax, 31);
contourf(X_fine, Y_fine, U_fine, levels);
hold on;
contour(X_fine, Y_fine, U_fine, levels(1:2:end), 'LineColor', 'k', ...
    'LineWidth', 0.3);
plot(0.3, -0.4, 'p', 'Color', ORANGE, 'MarkerSize', 15, ...
    'MarkerFaceColor', ORANGE, 'MarkerEdgeColor', 'white', ...
    'DisplayName', 'forcing centre');
hold off;
colormap(redblue(256));
cb = colorbar; cb.Label.String = '$u(x,y)$'; cb.Label.Interpreter = 'latex';
xlabel('$x$'); ylabel('$y$');
title('(b) Contour of solution');
axis equal;
legend('Location', 'northeast', 'FontSize', 9);

% Panel (c): NEW resonance sweep
subplot(2, 2, 3);
semilogy(k_axis, amp_at_force + 1e-18, '-', 'Color', NAVY, ...
         'LineWidth', 1.2, 'DisplayName', '$|u(\mathrm{forcing\ centre})|$');
hold on;
xline(k_24, '--', 'Color', CORAL, 'LineWidth', 0.8, ...
      'DisplayName', sprintf('$(2, 4)$ resonance, $k_{24} = %.2f$', k_24));
xline(k, ':', 'Color', TEAL, 'LineWidth', 0.8, ...
      'DisplayName', sprintf('showcase $k = %.1f$', k));
% Other resonance lines for context
modes_extra = {[1 1] [1 2] [2 2] [1 3] [2 3] [3 3] [1 4] [3 4] [2 5]};
for kk_idx = 1:numel(modes_extra)
    mn = modes_extra{kk_idx};
    kk_val = (pi / 2) * sqrt(mn(1)^2 + mn(2)^2);
    if kk_val >= 3.0 && kk_val <= 12.0
        xline(kk_val, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.4, ...
              'Alpha', 0.4, 'HandleVisibility', 'off');
    end
end
hold off;
xlabel('wavenumber $k$'); ylabel('$|u|$ at forcing centre');
title('(c) resonance amplification of $|u|$ vs $k$');
grid on; box on;
legend('Location', 'northwest', 'FontSize', 9);

% Panel (d): NEW (2,4) eigenmode
subplot(2, 2, 4);
em_max = max(abs(eigenmode_24(:)));
levels_e = linspace(-em_max, em_max, 31);
contourf(X_fine, Y_fine, eigenmode_24, levels_e);
hold on;
contour(X_fine, Y_fine, eigenmode_24, levels_e(1:4:end), 'LineColor', 'k', ...
    'LineWidth', 0.3);
hold off;
colormap(redblue(256));
colorbar;
xlabel('$x$'); ylabel('$y$');
title('(d) $(2,4)$ Dirichlet eigenmode shape');
axis equal;

% Main title
sgtitle('Helmholtz Equation: $\nabla^2 u + k^2 u = f(x,y)$, $(2,4)$ near-resonance', ...
    'FontSize', 12);

%% Save figure
output_file = fullfile(output_dir, 'helmholtz');
exportgraphics(fig, [output_file, '.pdf'], 'ContentType', 'vector');
exportgraphics(fig, [output_file, '.png'], 'Resolution', 300);
fprintf('\nFigure saved to: %s.pdf\n', output_file);

%% Compare solutions at different k values
fprintf('\nSolution amplitude at forcing location for different k:\n');
fprintf('%s\n', repmat('-', 1, 50));
fprintf('%8s %16s %15s\n', 'k', '|u| at center', 'Near mode');
fprintf('%s\n', repmat('-', 1, 50));

% Find forcing point in grid
[~, idx_x] = min(abs(X(1, :) - 0.3));
[~, idx_y] = min(abs(Y(:, 1) - (-0.4)));

for k_test = [3.0, 5.0, 7.0, 7.02, 8.0, 10.0]
    [~, ~, U_test] = solve_helmholtz(N, k_test, forcing_function);
    u_at_force = abs(U_test(idx_y, idx_x));

    % Find nearest eigenvalue
    nearest_mode = '';
    min_dist = inf;
    for m = 1:9
        for n = 1:9
            k_mn = (pi / 2) * sqrt(m^2 + n^2);
            if abs(k_mn - k_test) < min_dist
                min_dist = abs(k_mn - k_test);
                nearest_mode = sprintf('(%d,%d)', m, n);
            end
        end
    end

    fprintf('%8.2f %16.6f %15s\n', k_test, u_at_force, nearest_mode);
end

fprintf('%s\n', repmat('-', 1, 50));
fprintf('\nNote: Amplitude increases dramatically near resonance values.\n');

%% Local function: solve Helmholtz equation
function [X, Y, U] = solve_helmholtz(N, k, forcing_func)
    % Get Chebyshev matrix and grid
    [D, x] = cheb_matrix(N);
    D2 = D * D;  % Second derivative matrix

    % Interior points only
    n_int = N - 1;
    D2_int = D2(2:N, 2:N);
    x_int = x(2:N);

    % Build 2D Laplacian + k^2*I using Kronecker products
    % L = I kron D^2 + D^2 kron I + k^2*(I kron I)
    I = eye(n_int);
    L = kron(I, D2_int) + kron(D2_int, I) + k^2 * eye(n_int^2);

    % Create meshgrid for interior points
    [X_int, Y_int] = meshgrid(x_int, x_int);

    % Right-hand side
    F = forcing_func(X_int, Y_int);
    f_vec = F(:);

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

%% Red-blue diverging colormap
function cmap = redblue(n)
    % Red-blue diverging colormap
    if nargin < 1
        n = 256;
    end

    % Create diverging colormap from blue to white to red
    half = floor(n / 2);

    % Blue to white
    r1 = linspace(0.2, 1, half);
    g1 = linspace(0.2, 1, half);
    b1 = linspace(0.7, 1, half);

    % White to red
    r2 = linspace(1, 0.7, n - half);
    g2 = linspace(1, 0.2, n - half);
    b2 = linspace(1, 0.2, n - half);

    cmap = [r1(:), g1(:), b1(:); r2(:), g2(:), b2(:)];
end
