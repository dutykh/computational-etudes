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

%% RHS from analytical Laplacian (not discrete, to test convergence)
f_full = -exact_laplacian(xx, yy);
f_int = f_full(2:N, 2:N);

%% Solve
u_vec = A \ f_int(:);

%% Reconstruct full solution
U = zeros(N+1, N+1);
U(2:N, 2:N) = reshape(u_vec, N-1, N-1);

%% Compute error
err = abs(U - U_exact);
max_error = max(err(:));
fprintf('N = %d, Max error = %.2e\n', N, max_error);

%% Convergence study upfront so panel (d) can plot it
fprintf('\nConvergence Study:\n');
fprintf('----------------------------------------\n');
fprintf('%6s  %14s  %10s\n', 'N', 'Max Error', 'Ratio');
fprintf('----------------------------------------\n');

N_values = [8, 12, 16, 20, 24, 28, 32, 40, 48];
errs_conv = zeros(size(N_values));
prev_error = NaN;

for i = 1:length(N_values)
    Ni = N_values(i);
    [Di, xi] = cheb_matrix(Ni);
    D2i = Di * Di;
    [xxi, yyi] = meshgrid(xi, xi);

    U_exact_i = exact_solution(xxi, yyi);
    D2_int_i = D2i(2:Ni, 2:Ni);
    I_int_i = eye(Ni-1);
    Li = kron(D2_int_i, I_int_i) + kron(I_int_i, D2_int_i);
    Ai = -Li;

    f_full_i = -exact_laplacian(xxi, yyi);
    f_int_i = f_full_i(2:Ni, 2:Ni);

    u_vec_i = Ai \ f_int_i(:);
    Ui = zeros(Ni+1, Ni+1);
    Ui(2:Ni, 2:Ni) = reshape(u_vec_i, Ni-1, Ni-1);

    erri = max(abs(Ui(:) - U_exact_i(:)));
    errs_conv(i) = erri;

    if ~isnan(prev_error) && erri > 1e-15
        ratio = prev_error / erri;
    else
        ratio = NaN;
    end

    fprintf('%6d  %14.6e  %10.2f\n', Ni, erri, ratio);
    prev_error = erri;
end
fprintf('----------------------------------------\n');

%% Create figure with 4 panels (2x2)
fig = figure('Units', 'inches', 'Position', [1, 1, 11, 9]);

vmin = min(U_exact(:));
vmax = max(U_exact(:));

% Panel (a): Exact solution
subplot(2, 2, 1);
contourf(xx, yy, U_exact, 30, 'LineStyle', 'none');
hold on;
contour(xx, yy, U_exact, 10, 'k', 'LineWidth', 0.3);
xlabel('x', 'FontSize', 11);
ylabel('y', 'FontSize', 11);
title('(a) Exact solution u_{exact}(x, y)', 'FontSize', 11);
axis equal;
xlim([-1, 1]); ylim([-1, 1]);
caxis([vmin, vmax]);
colorbar;
colormap(gca, 'jet');

% Panel (b): Numerical solution
subplot(2, 2, 2);
contourf(xx, yy, U, 30, 'LineStyle', 'none');
hold on;
contour(xx, yy, U, 10, 'k', 'LineWidth', 0.3);
xlabel('x', 'FontSize', 11);
ylabel('y', 'FontSize', 11);
title(sprintf('(b) Numerical solution u_N(x, y) at N = %d', N), 'FontSize', 11);
axis equal;
xlim([-1, 1]); ylim([-1, 1]);
caxis([vmin, vmax]);
colorbar;
colormap(gca, 'jet');

% Panel (c): Pointwise log10 error
subplot(2, 2, 3);
err_plot = log10(max(err, 1e-16));
contourf(xx, yy, err_plot, 20, 'LineStyle', 'none');
xlabel('x', 'FontSize', 11);
ylabel('y', 'FontSize', 11);
title(sprintf('(c) log_{10}|u_N - u_{exact}| (max = %.2e)', max_error), 'FontSize', 11);
axis equal;
xlim([-1, 1]); ylim([-1, 1]);
cb = colorbar;
ylabel(cb, 'log_{10}|u_N - u_{exact}|', 'FontSize', 9);
colormap(gca, 'parula');

% Panel (d): NEW spectral convergence
subplot(2, 2, 4);
semilogy(N_values, errs_conv + 1e-18, '-o', 'Color', NAVY, ...
         'LineWidth', 1.2, 'MarkerSize', 6, 'MarkerFaceColor', 'w', ...
         'DisplayName', 'max |u_N - u_{exact}|');
hold on;
yline(1e-15, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.6, ...
      'DisplayName', 'machine epsilon');
hold off;
xlabel('N (Chebyshev degree, each direction)');
ylabel('max-norm error');
title('(d) spectral convergence: error vs N', 'FontSize', 11);
grid on; box on;
legend('Location', 'northeast', 'FontSize', 9);

sgtitle('2D Poisson Equation: Spectral Solution', 'FontSize', 13);

%% Save figure
exportgraphics(fig, output_file, 'ContentType', 'vector');
exportgraphics(fig, strrep(output_file, '.pdf', '.png'), 'Resolution', 300);

fprintf('Figure saved to: %s\n', output_file);

close(fig);

%% Local functions

function lap_u = exact_laplacian(xx, yy)
% Analytical Laplacian of the manufactured solution.
% u = g(x)*h(y), g(x) = (1-x^2)*cos(pi*x/2), h(y) = (1-y^2)*cos(pi*y/2).
% Lap(u) = g''(x)*h(y) + g(x)*h''(y).
    p = pi;
    cx = cos(p*xx/2); cy = cos(p*yy/2);
    sx = sin(p*xx/2); sy = sin(p*yy/2);
    gx  = (1 - xx.^2) .* cx;
    hy  = (1 - yy.^2) .* cy;
    gpp = -(2 + p^2/4 * (1 - xx.^2)) .* cx + 2*p*xx .* sx;
    hpp = -(2 + p^2/4 * (1 - yy.^2)) .* cy + 2*p*yy .* sy;
    lap_u = gpp .* hy + gx .* hpp;
end
