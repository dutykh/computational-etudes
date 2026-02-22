%% ho_pseudospectra.m - Pseudospectra of the Orr-Sommerfeld operator
%
% Chapter 14: Higher-Order Boundary Value Problems
%
% Computational Etude 14.7: Pseudospectra of the Orr-Sommerfeld operator
% for plane Poiseuille flow at R = 5772, alpha = 1, N = 100.
%
% The epsilon-pseudospectrum of the matrix pencil (A, B) is:
%
%   Lambda_eps = { z in C : sigma_min(A - z*B) < eps }
%
% For normal operators the pseudospectra are simply eps-balls around the
% eigenvalues. The Orr-Sommerfeld operator is highly non-normal, so its
% pseudospectra extend far from the eigenvalues, revealing extreme
% sensitivity to perturbations.
%
% At R = 5772, all eigenvalues lie barely in the stable half-plane.
% The contour plot of log10(sigma_min) shows that eps-perturbations of
% size as small as 10^{-8} can push eigenvalues into the unstable
% half-plane, explaining subcritical transition to turbulence.
%
% Reference: Trefethen, "Spectral Methods in MATLAB", Program 40.
%            Trefethen & Embree, "Spectra and Pseudospectra" (2005).
%
% Generates Figure: pseudospectra.pdf - Filled contour plot
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
output_dir = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch14', 'matlab');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% Parameters
R     = 5772;      % Critical Reynolds number
alpha = 1.0;       % Streamwise wavenumber
N     = 100;       % Number of Chebyshev intervals

fprintf('Pseudospectra of the Orr-Sommerfeld Operator\n');
fprintf('%s\n', repmat('=', 1, 60));
fprintf('  Reynolds number:      R = %d\n', R);
fprintf('  Wavenumber:           alpha = %g\n', alpha);
fprintf('  Chebyshev points:     N = %d\n', N);
fprintf('  Matrix size:          %d x %d\n', N - 1, N - 1);
fprintf('\n');

%% Build Orr-Sommerfeld matrices
fprintf('  Building Orr-Sommerfeld matrices ...\n');
[A, B] = build_os_matrices(N, R, alpha);
fprintf('  Matrix shapes: A = [%d x %d], B = [%d x %d]\n', ...
        size(A, 1), size(A, 2), size(B, 1), size(B, 2));

%% Compute eigenvalues
fprintf('  Computing eigenvalues ...\n');
eigenvalues = eig(A, B);

% Filter out infinite and NaN eigenvalues
finite_mask = isfinite(eigenvalues) & (abs(eigenvalues) < 1e8);
eigenvalues = eigenvalues(finite_mask);

% Sort by imaginary part (most unstable first)
[~, sort_idx] = sort(-imag(eigenvalues));
eigenvalues = eigenvalues(sort_idx);

fprintf('  Number of finite eigenvalues: %d\n', length(eigenvalues));
if ~isempty(eigenvalues)
    lam_crit = eigenvalues(1);
    fprintf('  Most unstable eigenvalue:\n');
    fprintf('    lambda = %.10f + %.10fi\n', real(lam_crit), imag(lam_crit));
end

%% Form C = B\A for standard eigenvalue pseudospectra
% Compute sigma_min(C - z*I) on a grid
fprintf('\n  Forming C = B\\A for pseudospectra computation ...\n');
C = B \ A;

%% Define grid for pseudospectra
nx = 120;
ny = 120;
x_re = linspace(-0.8, 0.2, nx);
x_im = linspace(-1.0, 0.05, ny);

fprintf('  Computing pseudospectra on %d x %d grid ...\n', nx, ny);
fprintf('  Grid: Re in [%.1f, %.1f], Im in [%.2f, %.2f]\n', ...
        x_re(1), x_re(end), x_im(1), x_im(end));
fprintf('  This requires %d SVD computations of %dx%d matrices.\n', ...
        nx * ny, size(C, 1), size(C, 2));
fprintf('\n');

%% Compute sigma_min on the grid
M_size = size(C, 1);
I_mat = eye(M_size);
sigma_grid = zeros(ny, nx);

tic;
for j = 1:ny
    for i = 1:nx
        z = x_re(i) + 1i * x_im(j);
        sv = svd(C - z * I_mat);
        sigma_grid(j, i) = log10(max(sv(end), 1e-20));
    end
    % Progress report every 10 rows
    if mod(j, 10) == 0
        pct = 100.0 * j / ny;
        fprintf('    ... %.0f%% complete\n', pct);
    end
end
t_elapsed = toc;
fprintf('\n  Pseudospectra computed in %.1f seconds.\n\n', t_elapsed);

%% Create figure: pseudospectra contour plot
fig = figure('Units', 'inches', 'Position', [1, 1, 8, 6.5]);

% Contour levels: log10(eps) = -8, -7, ..., -1
levels = -8:-1;

% Filled contours
contourf(x_re, x_im, sigma_grid, levels, 'LineColor', 'k', ...
         'LineWidth', 0.5);
hold on;

% Custom colormap (YlOrRd_r equivalent)
colormap(ylorrd_r_colormap(256));

% Contour line labels
[C_cont, h_cont] = contour(x_re, x_im, sigma_grid, levels, 'k', ...
                            'LineWidth', 0.5);
clabel(C_cont, h_cont, 'FontSize', 7, 'Interpreter', 'latex');

% Overlay eigenvalues (only those within the plot region)
in_region = (real(eigenvalues) >= x_re(1)) & (real(eigenvalues) <= x_re(end)) & ...
            (imag(eigenvalues) >= x_im(1)) & (imag(eigenvalues) <= x_im(end));
eigs_plot = eigenvalues(in_region);

plot(real(eigs_plot), imag(eigs_plot), 'o', 'Color', 'white', ...
     'MarkerSize', 4, 'MarkerFaceColor', 'white', ...
     'MarkerEdgeColor', NAVY, 'LineWidth', 0.8, ...
     'DisplayName', 'Eigenvalues');

% Stability boundary: vertical dashed line at Re = 0
xline(0, '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.8);

hold off;

% Colorbar
cb = colorbar;
cb.Label.String = '$\log_{10}\,\sigma_{\min}(C - z\,I)$';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 11;
cb.Ticks = levels;

% Labels and title
xlabel('$\mathrm{Re}(\lambda)$');
ylabel('$\mathrm{Im}(\lambda)$');
title(sprintf('$\\varepsilon$-Pseudospectra of Orr--Sommerfeld: $R = %d$, $\\alpha = %g$, $N = %d$', ...
      R, alpha, N), 'FontSize', 12);

legend('Location', 'northwest', 'FontSize', 9);
xlim([x_re(1), x_re(end)]);
ylim([x_im(1), x_im(end)]);
box on;

%% Save figure
exportgraphics(fig, fullfile(output_dir, 'pseudospectra.pdf'), ...
               'ContentType', 'vector');
fprintf('Figure saved to: %s\n', fullfile(output_dir, 'pseudospectra.pdf'));

close(fig);

fprintf('\n%s\n', repmat('=', 1, 60));


%% =====================================================================
%  Local functions
%  =====================================================================

function [A, B] = build_os_matrices(N, R, alpha)
% BUILD_OS_MATRICES  Build the generalized eigenvalue matrices A, B
% for the Orr-Sommerfeld equation for plane Poiseuille flow.
%
%   Boundary conditions: v(+-1) = v'(+-1) = 0 (no-slip).
%   Dirichlet conditions are enforced by restriction to interior points.
%   Neumann conditions are imposed by replacing the first and last rows.

    [D_cheb, x] = cheb_matrix(N);

    % Second and fourth derivative matrices
    D2_full = D_cheb * D_cheb;
    D4_full = D2_full * D2_full;

    % Interior grid (remove boundary points x_0=1, x_N=-1)
    x_int = x(2:N);

    % Base flow: U(y) = 1 - y^2, U''(y) = -2
    U_diag = diag(1.0 - x_int.^2);
    U2 = -2.0 * eye(N - 1);

    % Interior blocks of D^2 and D^4
    D2_int = D2_full(2:N, 2:N);
    D4_int = D4_full(2:N, 2:N);

    % Shift operator: S = D^2 - alpha^2 * I
    I_int = eye(N - 1);
    S = D2_int - alpha^2 * I_int;

    % S^2 = D4 - 2*alpha^2*D2 + alpha^4*I
    S2 = D4_int - 2.0 * alpha^2 * D2_int + alpha^4 * I_int;

    % A = S^2 / (i*alpha*R) - U*S + U''
    % B = -S
    A = S2 / (1i * alpha * R) - U_diag * S + U2;
    B = -S;

    % Impose Neumann boundary conditions via row replacement:
    % v'(+1) = D(1, 2:N) . v_int = 0 => row 1
    % v'(-1) = D(N+1, 2:N) . v_int = 0 => row N-1
    A(1, :)     = D_cheb(1, 2:N);
    B(1, :)     = 0.0;

    A(end, :)   = D_cheb(N + 1, 2:N);
    B(end, :)   = 0.0;
end


function cmap = ylorrd_r_colormap(m)
% YLORRD_R_COLORMAP  Yellow-Orange-Red reversed colormap.
    if nargin < 1, m = 256; end
    nodes = [0.5020 0.0000 0.1490;   % dark red
             0.7412 0.0000 0.1490;   % red
             0.8902 0.1020 0.1098;   % red-orange
             0.9882 0.3059 0.1647;   % orange-red
             0.9922 0.5529 0.2353;   % orange
             0.9961 0.6980 0.2980;   % yellow-orange
             0.9961 0.8510 0.4627;   % light yellow-orange
             1.0000 0.9294 0.6275;   % pale yellow
             1.0000 1.0000 0.8000;   % very pale yellow
             1.0000 1.0000 0.8980];  % near white
    x_nodes = linspace(0, 1, size(nodes, 1));
    xi = linspace(0, 1, m);
    cmap = interp1(x_nodes, nodes, xi);
end
