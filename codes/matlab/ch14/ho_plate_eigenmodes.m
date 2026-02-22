%% ho_plate_eigenmodes.m - Eigenmodes of a clamped square plate
%
% Chapter 14: Higher-Order Boundary Value Problems
%
% Computational Etude 14.4: Eigenmodes of a clamped square plate.
%
% Solves the biharmonic eigenvalue problem:
%
%   Delta^2 u = lambda u,   (x, y) in (-1, 1)^2
%   u = 0,  du/dn = 0       on the boundary
%
% The clamped plate problem models the vibrations of a thin elastic plate
% with all four edges rigidly fixed.  The biharmonic operator Delta^2 is
% the composition of two Laplacians.
%
% Discretization follows Trefethen's Program 39 (Spectral Methods in MATLAB):
%   1. Build 1D clamped biharmonic operator D4 via a polynomial trick
%      that automatically enforces u = du/dn = 0 on the boundary.
%   2. Assemble the 2D operator via Kronecker products:
%      L = I x D4 + D4 x I + 2 * D2int x D2int
%   3. Solve the eigenvalue problem L v = lambda v.
%   4. Display the first 25 eigenmodes as nodal-line contour plots.
%
% Parameters: N = 17 (default for figure generation).
%
% Generates Figure: plate_eigenmodes.pdf
%   5 x 5 grid of nodal-line contour plots with normalized eigenvalues.
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
N = 17;
n_modes = 25;

fprintf('Eigenmodes of a Clamped Square Plate\n');
fprintf('%s\n', repmat('=', 1, 55));
fprintf('  Chebyshev grid:     N = %d  (%d points per direction)\n', N, N + 1);
fprintf('  Interior points:    %d x %d = %d\n', N - 1, N - 1, (N - 1)^2);
fprintf('  Eigenvalue problem: %d x %d matrix\n', (N - 1)^2, (N - 1)^2);
fprintf('  Modes to compute:   %d\n\n', n_modes);

%% Build 1D clamped biharmonic operator
[D4, D2int, x] = build_clamped_biharmonic(N);

M = N - 1;   % number of interior points per direction
I_mat = eye(M);

%% Assemble 2D biharmonic via Kronecker products
L = kron(I_mat, D4) + kron(D4, I_mat) + 2.0 * kron(D2int, D2int);

%% Solve eigenvalue problem
[V, D_eig] = eig(L);
eigenvalues = real(diag(D_eig));

% Sort by ascending eigenvalue (smallest first)
[eigenvalues, idx] = sort(eigenvalues);
V = V(:, idx);

% Keep only positive eigenvalues
pos_mask = eigenvalues > 0;
eigenvalues = eigenvalues(pos_mask);
V = V(:, pos_mask);

% Keep first n_modes
eigenvalues = eigenvalues(1:n_modes);
V = V(:, 1:n_modes);

% Normalize: Lam_k = sqrt(lambda_k / lambda_1)
Lam = sqrt(eigenvalues / eigenvalues(1));

%% Reshape eigenvectors to 2D grids with zero-padded boundaries
modes = cell(n_modes, 1);
for k = 1:n_modes
    v = real(V(:, k));
    VV = reshape(v, [M, M]);

    % Pad with zeros on all boundaries
    U = zeros(N + 1, N + 1);
    U(2:N, 2:N) = VV;
    modes{k} = U;
end

%% Print eigenvalue table
fprintf('  Normalized eigenvalues  sqrt(lambda_k / lambda_1):\n');
fprintf('  %s\n', repmat('-', 1, 50));
for k = 1:n_modes
    fprintf('    Mode %2d:  %.6f\n', k, Lam(k));
end
fprintf('  %s\n\n', repmat('-', 1, 50));

%% Create figure: 5 x 5 grid of nodal-line contour plots
fig = figure('Units', 'inches', 'Position', [1, 1, 10, 10]);

[XX, YY] = meshgrid(x, x);
nrows = 5;
ncols = 5;

for k = 1:n_modes
    subplot(nrows, ncols, k);

    U = modes{k};

    % Plot nodal lines (contour at level 0)
    contour(XX, YY, U, [0, 0], 'Color', NAVY, 'LineWidth', 0.8);
    hold on;

    % Draw domain boundary
    plot([-1, 1, 1, -1, -1], [-1, -1, 1, 1, -1], 'k-', 'LineWidth', 1.0);
    hold off;

    % Label with normalized eigenvalue
    title(sprintf('$%.4f$', Lam(k)), 'FontSize', 9);

    % Clean up axes
    xlim([-1.05, 1.05]);
    ylim([-1.05, 1.05]);
    axis equal;
    set(gca, 'XTick', [], 'YTick', []);
end

sgtitle(['Eigenmodes of a Clamped Square Plate: ', ...
         '$\Delta^2 u = \lambda\, u$, ', ...
         '$u = \partial_n u = 0$ on $\partial\Omega$'], ...
        'FontSize', 13);

%% Save figure
exportgraphics(fig, fullfile(output_dir, 'plate_eigenmodes.pdf'), ...
               'ContentType', 'vector');
fprintf('  Figure saved to: %s\n', fullfile(output_dir, 'plate_eigenmodes.pdf'));
fprintf('%s\n', repmat('=', 1, 55));

close(fig);


%% =====================================================================
%  Local functions
%  =====================================================================

function [D4, D2int, x] = build_clamped_biharmonic(N)
% BUILD_CLAMPED_BIHARMONIC  1D clamped biharmonic operator via polynomial trick.
%
%   The clamped boundary conditions u(+-1) = u'(+-1) = 0 are enforced
%   implicitly: represent u(x) = (1 - x^2) * p(x), so that u and u'
%   vanish at x = +-1.
%
%   The fourth derivative of u = (1 - x^2) p is:
%     u'''' = (1 - x^2) p'''' - 8x p''' - 12 p''
%
%   so the effective operator on p at interior nodes is:
%     D4 = [diag(1 - x^2) D^4 - 8 diag(x) D^3 - 12 D^2] S
%
%   where S = diag(1/(1 - x_j^2)) at interior points (zeros at boundary).
%
%   Outputs:
%     D4    - Clamped biharmonic operator on interior points, (N-1) x (N-1)
%     D2int - Standard D^2 restricted to interior points, (N-1) x (N-1)
%     x     - Full Chebyshev grid (N+1 points)

    [D, x] = cheb_matrix(N);

    % Powers of D
    D2 = D * D;
    D3 = D2 * D;
    D4full = D3 * D;

    % Diagonal weight matrices
    x2 = x.^2;
    s = zeros(N + 1, 1);
    s(2:N) = 1.0 ./ (1.0 - x2(2:N));
    S = diag(s);

    % Build D4 on full grid, then restrict to interior
    D4mat = (diag(1.0 - x2) * D4full - 8.0 * diag(x) * D3 - 12.0 * D2) * S;

    % Restrict to interior points (indices 2 to N)
    D4 = D4mat(2:N, 2:N);
    D2int = D2(2:N, 2:N);
end
