%% poisson2d_cheb.m
%
% Figure 10.7: 2D Poisson equation on Chebyshev tensor grid
%
% PDE: -Lap(u) = f, -1 < x, y < 1
% BCs: u = 0 on boundary (Dirichlet)
%
% Uses manufactured solution:
%   u(x, y) = (1 - x^2)(1 - y^2) cos(pi*x/2) cos(pi*y/2)
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
% Part of "Computational Etudes: A Spectral Approach"
% https://github.com/dutykh/computational-etudes

clear; close all; clc;

%% Configuration
N = 32;

% Colors
NAVY = [0.078 0.176 0.431];

% Output path
script_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch10', 'matlab');
output_file = fullfile(output_dir, 'poisson2d_solution.pdf');

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% Exact solution (manufactured)
exact_solution = @(x, y) (1 - x.^2) .* (1 - y.^2) .* cos(pi*x/2) .* cos(pi*y/2);

%% Setup grid and operators
[D, x] = cheb_matrix(N);
D2 = D * D;
[xx, yy] = meshgrid(x, x);

% Exact solution
U_exact = exact_solution(xx, yy);

%% Interior block
D2_int = D2(2:N, 2:N);
I_int = eye(N-1);

% Kronecker sum Laplacian (represents Δ)
L = kron(D2_int, I_int) + kron(I_int, D2_int);

% For -Δu = f, we have -L * u = f, so L * u = -f
% Or equivalently, use A = -L so that A * u = f
A = -L;

%% RHS: apply discrete Laplacian to exact solution
Lap_U = D2 * U_exact + U_exact * D2';
f_int = -Lap_U(2:N, 2:N);

%% Solve
u_vec = A \ f_int(:);

%% Reconstruct full solution
U = zeros(N+1, N+1);
U(2:N, 2:N) = reshape(u_vec, N-1, N-1);

%% Compute error
err = abs(U - U_exact);
max_error = max(err(:));
fprintf('N = %d, Max error = %.2e\n', N, max_error);

%% Create figure with 3 panels
fig = figure('Units', 'inches', 'Position', [1, 1, 14, 4.5]);

% Common color limits
vmin = min(U_exact(:));
vmax = max(U_exact(:));

% Panel 1: Exact solution
subplot(1, 3, 1);
contourf(xx, yy, U_exact, 30, 'LineStyle', 'none');
hold on;
contour(xx, yy, U_exact, 10, 'k', 'LineWidth', 0.3);
xlabel('x', 'FontSize', 11);
ylabel('y', 'FontSize', 11);
title('Exact Solution', 'FontSize', 12);
axis equal;
xlim([-1, 1]); ylim([-1, 1]);
caxis([vmin, vmax]);
colorbar;
colormap(gca, 'jet');

% Panel 2: Numerical solution
subplot(1, 3, 2);
contourf(xx, yy, U, 30, 'LineStyle', 'none');
hold on;
contour(xx, yy, U, 10, 'k', 'LineWidth', 0.3);
xlabel('x', 'FontSize', 11);
ylabel('y', 'FontSize', 11);
title('Numerical Solution', 'FontSize', 12);
axis equal;
xlim([-1, 1]); ylim([-1, 1]);
caxis([vmin, vmax]);
colorbar;
colormap(gca, 'jet');

% Panel 3: Error (log scale)
subplot(1, 3, 3);
err_plot = log10(max(err, 1e-16));
contourf(xx, yy, err_plot, 20, 'LineStyle', 'none');
xlabel('x', 'FontSize', 11);
ylabel('y', 'FontSize', 11);
title(sprintf('log_{10} Error (max = %.2e)', max_error), 'FontSize', 12);
axis equal;
xlim([-1, 1]); ylim([-1, 1]);
cb = colorbar;
ylabel(cb, 'log_{10}|u - u_{exact}|', 'FontSize', 10);
colormap(gca, 'parula');

sgtitle('2D Poisson Equation: Spectral Solution', 'FontSize', 14);

%% Save figure
exportgraphics(fig, output_file, 'ContentType', 'vector');
exportgraphics(fig, strrep(output_file, '.pdf', '.png'), 'Resolution', 300);

fprintf('Figure saved to: %s\n', output_file);

%% Convergence study
fprintf('\nConvergence Study:\n');
fprintf('----------------------------------------\n');
fprintf('%6s  %14s  %10s\n', 'N', 'Max Error', 'Ratio');
fprintf('----------------------------------------\n');

prev_error = NaN;
N_values = [8, 12, 16, 20, 24, 28, 32, 40, 48];

for i = 1:length(N_values)
    Ni = N_values(i);
    [Di, xi] = cheb_matrix(Ni);
    D2i = Di * Di;
    [xxi, yyi] = meshgrid(xi, xi);

    U_exact_i = exact_solution(xxi, yyi);
    D2_int_i = D2i(2:Ni, 2:Ni);
    I_int_i = eye(Ni-1);
    Li = kron(D2_int_i, I_int_i) + kron(I_int_i, D2_int_i);
    Ai = -Li;  % For -Δu = f

    Lap_Ui = D2i * U_exact_i + U_exact_i * D2i';
    f_int_i = -Lap_Ui(2:Ni, 2:Ni);

    u_vec_i = Ai \ f_int_i(:);
    Ui = zeros(Ni+1, Ni+1);
    Ui(2:Ni, 2:Ni) = reshape(u_vec_i, Ni-1, Ni-1);

    erri = max(abs(Ui(:) - U_exact_i(:)));

    if ~isnan(prev_error) && erri > 1e-15
        ratio = prev_error / erri;
    else
        ratio = NaN;
    end

    fprintf('%6d  %14.6e  %10.2f\n', Ni, erri, ratio);
    prev_error = erri;
end
fprintf('----------------------------------------\n');

close(fig);
