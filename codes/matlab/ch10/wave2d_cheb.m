%% wave2d_cheb.m
%
% Figure 10.4: 2D wave equation on Chebyshev tensor grid
%
% PDE: u_tt = c^2 (u_xx + u_yy), -1 < x, y < 1, t > 0
% BCs: u = 0 on boundary (Dirichlet)
% ICs: u(x, y, 0) = exp(-30*((x-0.2)^2 + (y+0.3)^2)), u_t = 0
%
% Uses leapfrog time stepping with tensor product Laplacian.
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
% Part of "Computational Etudes: A Spectral Approach"
% https://github.com/dutykh/computational-etudes

clear; close all; clc;

%% Configuration
N = 40;           % Chebyshev intervals in each direction
c = 1.0;          % Wave speed
tmax = 1.5;       % Final time
alpha = 3.0;      % CFL parameter

% Colors
NAVY = [0.078 0.176 0.431];

% Output path
script_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch10', 'matlab');
output_file = fullfile(output_dir, 'wave2d_snapshots.pdf');

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% Setup grid
x = cos(pi * (0:N)' / N);
[xx, yy] = meshgrid(x, x);

%% Initial condition: offset Gaussian
U = exp(-30 * ((xx - 0.2).^2 + (yy + 0.3).^2));
U_prev = U;

%% Time stepping setup
dt = alpha / N^2;
nsteps = ceil(tmax / dt);

% Snapshot times
t_snap = [0.0, 0.3, 0.7, 1.2];
snapshots = cell(4, 1);
snapshots{1} = U;
current_snap = 2;

% Taylor start
Lap_U = laplacian_chebfft(U, N);
U_prev = U - 0.5 * (c * dt)^2 * Lap_U;

t = 0;
for n = 1:nsteps
    % 2D Laplacian via chebfft along rows and columns
    Lap_U = laplacian_chebfft(U, N);

    % Leapfrog step
    U_new = 2*U - U_prev + (c * dt)^2 * Lap_U;

    % Boundary conditions
    U_new(1, :) = 0;
    U_new(N+1, :) = 0;
    U_new(:, 1) = 0;
    U_new(:, N+1) = 0;

    U_prev = U;
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

% Color limits
vmin = -0.5;
vmax = 1.0;

for i = 1:4
    subplot(2, 2, i);

    contourf(xx, yy, snapshots{i}, 30, 'LineStyle', 'none');
    hold on;
    contour(xx, yy, snapshots{i}, 10, 'k', 'LineWidth', 0.3);

    % Mark initial center on first panel
    if i == 1
        plot(0.2, -0.3, 'k+', 'MarkerSize', 10, 'LineWidth', 2);
    end

    xlabel('x', 'FontSize', 11);
    ylabel('y', 'FontSize', 11);
    title(sprintf('t = %.2f', t_snap(i)), 'FontSize', 12);
    axis equal;
    xlim([-1, 1]);
    ylim([-1, 1]);
    caxis([vmin, vmax]);
    colorbar;
    colormap('jet');
end

sgtitle('2D Wave Equation: Vibrating Membrane', 'FontSize', 14);

%% Save figure
exportgraphics(fig, output_file, 'ContentType', 'vector');
exportgraphics(fig, strrep(output_file, '.pdf', '.png'), 'Resolution', 300);

fprintf('Figure saved to: %s\n', output_file);

close(fig);

%% -----------------------------------------------------------------------
function Lap = laplacian_chebfft(U, N)
%LAPLACIAN_CHEBFFT  Compute 2D Laplacian using chebfft along rows/columns.
    Lap = zeros(N+1, N+1);
    for i = 1:N+1
        Lap(i, :) = Lap(i, :) + chebfft(chebfft(U(i, :)'))';  % d^2/dx^2
    end
    for j = 1:N+1
        Lap(:, j) = Lap(:, j) + chebfft(chebfft(U(:, j)));    % d^2/dy^2
    end
end
