%% heat2d_cheb.m
%
% Figure 10.6: 2D heat equation on Chebyshev tensor grid
%
% PDE: u_t = u_xx + u_yy, -1 < x, y < 1, t > 0
% BCs: u = 0 on boundary (Dirichlet)
% ICs: u(x, y, 0) = exp(-25*((x+0.3)^2 + (y-0.1)^2))
%
% Uses backward Euler with Kronecker sum formulation.
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
% Part of "Computational Etudes: A Spectral Approach"
% https://github.com/dutykh/computational-etudes

clear; close all; clc;

%% Configuration
N = 32;           % Chebyshev intervals
tmax = 0.6;       % Final time
dt = 0.002;       % Time step

% Colors
NAVY = [0.078 0.176 0.431];

% Output path
script_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch10', 'matlab');
output_file = fullfile(output_dir, 'heat2d_snapshots.pdf');

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% Setup grid and operators
[D, x] = cheb_matrix(N);
D2 = D * D;
[xx, yy] = meshgrid(x, x);

% Interior block
D2_int = D2(2:N, 2:N);
I_int = eye(N-1);

% Kronecker sum Laplacian for interior
L = kron(D2_int, I_int) + kron(I_int, D2_int);
A_mat = eye((N-1)^2) - dt * L;

% Pre-factorize A (constant across all time steps)
dA = decomposition(A_mat, 'lu');

%% Initial condition
U = exp(-25 * ((xx + 0.3).^2 + (yy - 0.1).^2));
U(1, :) = 0; U(N+1, :) = 0;
U(:, 1) = 0; U(:, N+1) = 0;

%% Snapshot times
t_snap = [0.0, 0.05, 0.15, 0.5];
snapshots = cell(4, 1);
snapshots{1} = U;
current_snap = 2;

%% Time stepping
nsteps = ceil(tmax / dt);
t = 0;

for n = 1:nsteps
    % Extract interior and solve
    u_int = U(2:N, 2:N);
    u_vec = u_int(:);
    u_vec_new = dA \ u_vec;

    % Reconstruct
    U_new = zeros(N+1, N+1);
    U_new(2:N, 2:N) = reshape(u_vec_new, N-1, N-1);
    U = U_new;

    t = t + dt;

    % Save snapshot
    if current_snap <= length(t_snap) && t >= t_snap(current_snap)
        snapshots{current_snap} = U;
        current_snap = current_snap + 1;
    end
end

%% Create figure with 4 panels
fig = figure('Units', 'inches', 'Position', [1, 1, 10, 9]);

% Plot each snapshot with individual color scaling
for i = 1:4
    subplot(2, 2, i);

    % Individual color scaling for each panel
    vmin = 0;
    vmax = max(snapshots{i}(:));

    contourf(xx, yy, snapshots{i}, 30, 'LineStyle', 'none');
    hold on;
    contour(xx, yy, snapshots{i}, 8, 'k', 'LineWidth', 0.3);

    % Mark initial center on first panel
    if i == 1
        plot(-0.3, 0.1, 'w+', 'MarkerSize', 10, 'LineWidth', 2);
    end

    xlabel('x', 'FontSize', 11);
    ylabel('y', 'FontSize', 11);
    title(sprintf('t = %.3f', t_snap(i)), 'FontSize', 12);
    axis equal;
    xlim([-1, 1]);
    ylim([-1, 1]);
    caxis([vmin, vmax]);
    colorbar;
    colormap('hot');
end

sgtitle('2D Heat Equation: Diffusion on a Square', 'FontSize', 14);

%% Save figure
exportgraphics(fig, output_file, 'ContentType', 'vector');
exportgraphics(fig, strrep(output_file, '.pdf', '.png'), 'Resolution', 300);

fprintf('Figure saved to: %s\n', output_file);

close(fig);
