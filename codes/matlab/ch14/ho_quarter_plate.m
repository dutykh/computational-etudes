%% ho_quarter_plate.m - Eigenmodes of a clamped plate on the quarter-square
%
% Chapter 14: Higher-Order Boundary Value Problems
%
% Computational Etude 14.6: Quarter-square symmetry reduction on [0,1]^2.
%
% Eigenmodes of the biharmonic operator Delta^2 u = lambda u on the
% quarter plate [0,1]^2 using symmetry reduction. Neumann (symmetry)
% conditions at x=0, y=0 select even-even modes, while clamped
% conditions are imposed at x=1, y=1.
%
% The Chebyshev grid on [-1,1] is mapped to [0,1] via x_phys = (x+1)/2,
% scaling all derivative operators by a factor of 2.
%
% The 2D biharmonic operator is assembled via Kronecker products:
%
%   Delta^2 = D4 (x) I + 2 D2 (x) D2 + I (x) D4
%
% Boundary conditions (fourth-order requires two per edge):
%   At x=0, y=0 (symmetry): u' = 0  and  u''' = 0
%   At x=1, y=1 (clamped):  u = 0   and  u' = 0
%
% Generates Figure: quarter_plate.pdf - 4x3 grid of mode contour plots
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
N = 13;
n_modes = 12;

%% Setup grid and operators on [0,1]
[D, x_cheb] = cheb_matrix(N);

% Map [-1,1] -> [0,1]: x_phys = (x+1)/2, D_phys = 2*D
D_phys = 2.0 * D;
D2 = D_phys * D_phys;
D3 = D2 * D_phys;
D4 = D2 * D2;

x_phys = (x_cheb + 1.0) / 2.0;  % x_phys(1)=1, x_phys(N+1)=0

% Index ordering:
%   i = 1   -> x_phys = 1  (clamped boundary)
%   i = N+1 -> x_phys = 0  (symmetry boundary)

n = N + 1;   % number of grid points per direction
I1d = eye(n);

%% Assemble 2D biharmonic operator via Kronecker products
L = kron(I1d, D4) + 2.0 * kron(D2, D2) + kron(D4, I1d);

% Mass matrix for the generalized eigenvalue problem
M = eye(n * n);

%% Apply boundary conditions by boundary bordering
% Grid ordering: u(i + (j-1)*n) corresponds to point (x_phys(i), x_phys(j))
% In MATLAB: i goes 1..n, j goes 1..n
bc_rows = false(n * n, 1);

% Pass 1: Dirichlet u = 0 at clamped boundaries (x=1 and y=1)
for j = 1:n
    for i = 1:n
        k = i + (j - 1) * n;

        if i == 1
            % x_phys = 1 (clamped): u(1, y) = 0
            L(k, :) = 0.0;
            L(k, k) = 1.0;
            M(k, :) = 0.0;
            bc_rows(k) = true;

        elseif j == 1 && ~bc_rows(k)
            % y_phys = 1 (clamped): u(x, 1) = 0
            L(k, :) = 0.0;
            L(k, k) = 1.0;
            M(k, :) = 0.0;
            bc_rows(k) = true;
        end
    end
end

% Pass 2: Normal derivative u_n = 0 at clamped boundaries
for j = 1:n
    for i = 1:n
        k = i + (j - 1) * n;
        if bc_rows(k), continue; end

        if i == 2
            % du/dx = 0 at x_phys = 1 (use row for i=2, next to clamped edge)
            L(k, :) = 0.0;
            for ii = 1:n
                L(k, ii + (j - 1) * n) = D_phys(1, ii);
            end
            M(k, :) = 0.0;
            bc_rows(k) = true;

        elseif j == 2 && ~bc_rows(k)
            % du/dy = 0 at y_phys = 1 (use row for j=2)
            L(k, :) = 0.0;
            for jj = 1:n
                L(k, i + (jj - 1) * n) = D_phys(1, jj);
            end
            M(k, :) = 0.0;
            bc_rows(k) = true;
        end
    end
end

% Pass 3: Neumann u_n = 0 at symmetry boundaries (x=0 and y=0)
for j = 1:n
    for i = 1:n
        k = i + (j - 1) * n;
        if bc_rows(k), continue; end

        if i == n
            % du/dx = 0 at x_phys = 0 (symmetry)
            L(k, :) = 0.0;
            for ii = 1:n
                L(k, ii + (j - 1) * n) = D_phys(n, ii);
            end
            M(k, :) = 0.0;
            bc_rows(k) = true;

        elseif j == n && ~bc_rows(k)
            % du/dy = 0 at y_phys = 0 (symmetry)
            L(k, :) = 0.0;
            for jj = 1:n
                L(k, i + (jj - 1) * n) = D_phys(n, jj);
            end
            M(k, :) = 0.0;
            bc_rows(k) = true;
        end
    end
end

% Pass 4: Third derivative = 0 at symmetry boundaries
% For even functions: u_xxx(0) = 0 and u_yyy(0) = 0
for j = 1:n
    for i = 1:n
        k = i + (j - 1) * n;
        if bc_rows(k), continue; end

        if i == n - 1
            % d^3u/dx^3 = 0 at x_phys = 0 (symmetry)
            L(k, :) = 0.0;
            for ii = 1:n
                L(k, ii + (j - 1) * n) = D3(n, ii);
            end
            M(k, :) = 0.0;
            bc_rows(k) = true;

        elseif j == n - 1 && ~bc_rows(k)
            % d^3u/dy^3 = 0 at y_phys = 0 (symmetry)
            L(k, :) = 0.0;
            for jj = 1:n
                L(k, i + (jj - 1) * n) = D3(n, jj);
            end
            M(k, :) = 0.0;
            bc_rows(k) = true;
        end
    end
end

%% Solve the generalized eigenvalue problem L v = lambda M v
[V, D_eig] = eig(L, M);
eigenvalues = diag(D_eig);

% Filter: keep only finite, real, positive eigenvalues
finite_mask = isfinite(eigenvalues);
eigenvalues = eigenvalues(finite_mask);
V = V(:, finite_mask);

% Take real parts (self-adjoint problem)
eigenvalues_real = real(eigenvalues);

% Discard spurious near-zero and negative eigenvalues
pos_mask = eigenvalues_real > 1.0;
eigenvalues_real = eigenvalues_real(pos_mask);
V = V(:, pos_mask);

% Sort by ascending eigenvalue
[eigenvalues_real, idx] = sort(eigenvalues_real);
V = V(:, idx);

% Keep the requested number of modes
n_avail = min(n_modes, length(eigenvalues_real));
eigenvalues_real = eigenvalues_real(1:n_avail);
V = V(:, 1:n_avail);

% Normalized eigenvalues: Lam_k = sqrt(lambda_k / lambda_1)
Lam = sqrt(eigenvalues_real / eigenvalues_real(1));

% Reshape eigenvectors to 2D grids
modes = cell(n_avail, 1);
for k = 1:n_avail
    v = real(V(:, k));
    U = reshape(v, [n, n]);  % column-major: U(i,j) at (x_phys(i), x_phys(j))
    modes{k} = U;
end

%% Print eigenvalue table
fprintf('\nQuarter-Plate Symmetry Reduction: Even-Even Eigenmodes\n');
fprintf('%s\n', repmat('=', 1, 60));
fprintf('  Domain:             [0, 1]^2  (quarter of [-1, 1]^2)\n');
fprintf('  Chebyshev grid:     N = %d  (%d points per direction)\n', N, n);
fprintf('  System size:        %d x %d\n', n^2, n^2);
fprintf('  Modes computed:     %d\n', n_avail);
fprintf('  BCs at x=0, y=0:   Neumann (symmetry)\n');
fprintf('  BCs at x=1, y=1:   Clamped (u = du/dn = 0)\n');
fprintf('\n  Even-even eigenvalues of the clamped square plate:\n');
fprintf('  %s\n', repmat('-', 1, 55));
fprintf('  %4s  %14s  %18s\n', 'Mode', 'lambda_k', 'sqrt(lam_k/lam_1)');
fprintf('  %s\n', repmat('-', 1, 55));
for k = 1:n_avail
    fprintf('  %4d  %14.6f  %18.6f\n', k, eigenvalues_real(k), Lam(k));
end
fprintf('  %s\n', repmat('-', 1, 55));

%% Create figure: 4x3 grid of contour plots
[XX, YY] = meshgrid(x_phys, x_phys);

fig = figure('Units', 'inches', 'Position', [1, 1, 8, 10]);

nrows = 4;
ncols = 3;

for k = 1:n_avail
    row = ceil(k / ncols);
    col = k - (row - 1) * ncols;
    subplot(nrows, ncols, k);

    U = modes{k}';  % transpose so rows=y, cols=x for contour

    % Filled contours for the mode shape
    max_val = max(abs(U(:)));
    levels = linspace(-max_val, max_val, 21);
    contourf(XX, YY, U, levels, 'LineColor', 'none');
    hold on;

    % Custom RdBu_r colormap
    colormap(gca, rdbu_colormap(256));
    caxis([-max_val, max_val]);

    % Nodal lines (contour at level 0)
    contour(XX, YY, U, [0, 0], 'Color', NAVY, 'LineWidth', 0.9);

    % Draw domain boundary
    plot([0, 1, 1, 0, 0], [0, 0, 1, 1, 0], 'k-', 'LineWidth', 1.0);

    % Mark symmetry edges with dashed lines
    plot([0, 0], [0, 1], '--', 'Color', TEAL, 'LineWidth', 1.0);
    plot([0, 1], [0, 0], '--', 'Color', TEAL, 'LineWidth', 1.0);

    hold off;

    % Label with normalized eigenvalue
    title(sprintf('$%.4f$', Lam(k)), 'FontSize', 9);

    % Clean up axes
    xlim([-0.05, 1.05]);
    ylim([-0.05, 1.05]);
    axis equal;
    set(gca, 'XTick', [0, 0.5, 1], 'YTick', [0, 0.5, 1]);
    set(gca, 'FontSize', 7);
end

sgtitle({['Even-Even Eigenmodes of a Clamped Plate (Quarter Domain)'], ...
         ['$\Delta^2 u = \lambda\, u$ on $[0,1]^2$, ', ...
          '$\partial_n u = 0$ (symmetry), ', ...
          '$u = \partial_n u = 0$ (clamped)']}, ...
        'FontSize', 11, 'Interpreter', 'latex');

%% Save figure
exportgraphics(fig, fullfile(output_dir, 'quarter_plate.pdf'), ...
               'ContentType', 'vector');
fprintf('\nFigure saved to: %s\n', fullfile(output_dir, 'quarter_plate.pdf'));

close(fig);

fprintf('%s\n', repmat('=', 1, 60));


%% =====================================================================
%  Local functions
%  =====================================================================

function cmap = rdbu_colormap(m)
% RDBU_COLORMAP  Red-Blue diverging colormap similar to matplotlib's RdBu_r.
    if nargin < 1, m = 256; end
    nodes = [0.0196 0.1882 0.3804;   % dark blue
             0.1294 0.4000 0.6745;   % blue
             0.2627 0.5765 0.7647;   % light blue
             0.5725 0.7725 0.8706;   % very light blue
             0.8196 0.8980 0.9412;   % pale blue
             0.9686 0.9686 0.9686;   % near white
             0.9922 0.8588 0.7804;   % pale red
             0.9569 0.6471 0.5098;   % light red
             0.8392 0.3765 0.3020;   % red
             0.6980 0.0941 0.1686;   % dark red
             0.4039 0.0000 0.1216];  % very dark red
    x_nodes = linspace(0, 1, size(nodes, 1));
    xi = linspace(0, 1, m);
    cmap = interp1(x_nodes, nodes, xi);
end
