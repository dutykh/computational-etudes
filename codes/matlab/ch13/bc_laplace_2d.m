%% bc_laplace_2d.m - 2D Laplace equation with piecewise boundary data
%
% Chapter 13: Boundary Conditions in Spectral Methods
%
% Laplace equation with piecewise boundary data on [-1, 1]^2:
%
%   u_xx + u_yy = 0,  (x, y) in (-1, 1)^2
%
% Boundary conditions:
%   u(x, 1)  = sin^4(pi*x)  for x < 0  (top, left half)
%   u(1, y)  = (1/5)*sin(3*pi*y)        (right boundary)
%   u = 0    everywhere else on the boundary
%
% Solved via boundary bordering (Method II):
%   1. Build full 2D Laplacian L = I kron D^2 + D^2 kron I
%   2. Replace boundary rows with identity rows and set boundary data
%   3. Solve the global linear system
%
% Generates Figure: laplace_2d_mixed.pdf - 3D surface plot
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
N = 24;

%% Setup grid and operators
[D, x] = cheb_matrix(N);
D2 = D * D;
n = N + 1;   % Points per direction

% 2D tensor product grid
[XX, YY] = meshgrid(x, x);
xx = XX(:);
yy = YY(:);

% Build 2D Laplacian: L = I kron D2 + D2 kron I
I_mat = eye(n);
L = kron(I_mat, D2) + kron(D2, I_mat);

% Right-hand side (zero for Laplace equation)
rhs = zeros(n * n, 1);

%% Apply boundary conditions
for k = 1:n * n
    is_boundary = (abs(abs(xx(k)) - 1.0) < 1e-14) || ...
                  (abs(abs(yy(k)) - 1.0) < 1e-14);

    if is_boundary
        % Replace row with identity row
        L(k, :) = 0.0;
        L(k, k) = 1.0;

        xk = xx(k);
        yk = yy(k);

        % Top boundary: y = 1
        if abs(yk - 1.0) < 1e-14
            if xk < 0
                rhs(k) = sin(pi * xk)^4;
            else
                rhs(k) = 0.0;
            end
        % Right boundary: x = 1
        elseif abs(xk - 1.0) < 1e-14
            rhs(k) = sin(3 * pi * yk) / 5.0;
        % All other boundaries: u = 0
        else
            rhs(k) = 0.0;
        end
    end
end

%% Solve the linear system
u_vec = L \ rhs;
U = reshape(u_vec, [n, n]);

%% Create figure: 3D surface plot
fig = figure('Units', 'inches', 'Position', [1, 1, 8, 6]);

% Custom colormap
cmap = book_colormap(NAVY, TEAL, SKY, 256);

surf(XX, YY, U, 'EdgeColor', [0.5 0.5 0.5], 'LineWidth', 0.3, ...
     'FaceAlpha', 0.9);
colormap(cmap);
xlabel('$x$');
ylabel('$y$');
zlabel('$u(x, y)$');
title(sprintf('Laplace Equation: $\\nabla^2 u = 0$ with Piecewise BCs ($N = %d$)', N), ...
      'FontSize', 12);
view(25, -55);
cb = colorbar;
cb.Label.String = '$u(x, y)$';
cb.Label.Interpreter = 'latex';

%% Save figure
exportgraphics(fig, fullfile(output_dir, 'laplace_2d_mixed.pdf'), ...
               'ContentType', 'vector');
fprintf('Figure saved to: %s\n', fullfile(output_dir, 'laplace_2d_mixed.pdf'));

%% Print summary
fprintf('\n2D Laplace Equation with Piecewise Boundary Data\n');
fprintf('%s\n', repmat('=', 1, 55));
fprintf('  Grid size:          N = %d  (%d unknowns)\n', N, n^2);
fprintf('  Domain:             [-1, 1]^2\n');
fprintf('  Solution range:     [%.6f, %.6f]\n', min(U(:)), max(U(:)));
fprintf('  BCs applied:\n');
fprintf('    Top (y=1, x<0):   sin^4(pi*x)\n');
fprintf('    Right (x=1):      (1/5)*sin(3*pi*y)\n');
fprintf('    Elsewhere:        u = 0\n');
fprintf('%s\n', repmat('=', 1, 55));

% Verify boundary conditions
u_top   = U(1, :)';  % y = x(1) = 1
bc_top_exact = zeros(n, 1);
mask = x < 0;
bc_top_exact(mask) = sin(pi * x(mask)).^4;
top_err = max(abs(u_top - bc_top_exact));

u_right = U(:, 1);   % x = x(1) = 1
bc_right_exact = sin(3 * pi * x) / 5.0;
right_err = max(abs(u_right - bc_right_exact));

u_bottom = U(n, :)';    % y = x(N+1) = -1
bottom_err = max(abs(u_bottom));

u_left = U(:, n);    % x = x(N+1) = -1
left_err = max(abs(u_left));

fprintf('\n  BC verification:\n');
fprintf('    Top boundary error:    %.2e\n', top_err);
fprintf('    Right boundary error:  %.2e\n', right_err);
fprintf('    Bottom boundary error: %.2e\n', bottom_err);
fprintf('    Left boundary error:   %.2e\n', left_err);

close(fig);


%% =====================================================================
%  Helper function
%  =====================================================================

function cmap = book_colormap(c1, c2, c3, m)
% BOOK_COLORMAP  Custom colormap using book colors transitioning to white.
    if nargin < 4, m = 256; end
    nodes = [c1; c2; c3; 1 1 1];
    x_nodes = linspace(0, 1, size(nodes, 1));
    xi = linspace(0, 1, m);
    cmap = interp1(x_nodes, nodes, xi);
end
