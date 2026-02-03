%% wave1d_cheb.m
%
% Figure 10.3: 1D wave equation on Chebyshev grid
%
% PDE: u_tt = c^2 u_xx, -1 < x < 1, t > 0
% BCs: u(-1, t) = u(1, t) = 0 (Dirichlet)
% ICs: u(x, 0) = sin(pi/2 * (1 + x)), u_t(x, 0) = 0
%
% Uses leapfrog time stepping with Taylor start.
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
% Part of "Computational Etudes: A Spectral Approach"
% https://github.com/dutykh/computational-etudes

clear; close all; clc;

%% Configuration
N = 80;           % Chebyshev intervals
c = 1.0;          % Wave speed
tmax = 6.0;       % Final time
alpha = 4.0;      % CFL parameter

% Colors
NAVY = [0.078 0.176 0.431];

% Output path
script_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch10', 'matlab');
output_file = fullfile(output_dir, 'wave1d_waterfall.pdf');

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% Setup grid and operators
x = cos(pi * (0:N)' / N);
[D, ~] = cheb_matrix(N);
D2 = D * D;

%% Initial condition
u = sin(0.5 * pi * (1 + x));

%% Time stepping
dt = alpha / N^2;
nsteps = ceil(tmax / dt);

% Storage for waterfall plot
n_save = 80;
save_every = max(1, floor(nsteps / n_save));
U_save = zeros(N+1, n_save+1);
t_save = zeros(1, n_save+1);
U_save(:, 1) = u;
t_save(1) = 0;
save_idx = 2;

% Taylor start: u_prev = u + 0.5 * dt^2 * u_tt with u_t = 0
u_tt0 = c^2 * (D2 * u);
u_tt0(1) = 0;
u_tt0(N+1) = 0;
u_prev = u - 0.5 * dt^2 * u_tt0;

t = 0;
for n = 1:nsteps
    % Compute Laplacian
    Lap_u = c^2 * (D2 * u);
    Lap_u(1) = 0;
    Lap_u(N+1) = 0;

    % Leapfrog step
    u_new = 2*u - u_prev + dt^2 * Lap_u;

    % Boundary conditions
    u_new(1) = 0;
    u_new(N+1) = 0;

    u_prev = u;
    u = u_new;
    t = t + dt;

    % Save snapshot
    if mod(n, save_every) == 0 && save_idx <= n_save + 1
        U_save(:, save_idx) = u;
        t_save(save_idx) = t;
        save_idx = save_idx + 1;
    end
end

%% Create figure
fig = figure('Units', 'inches', 'Position', [1, 1, 10, 7]);

% 3D surface plot
[X, T] = meshgrid(x, t_save);
surf(X, T, U_save', 'EdgeColor', 'none', 'FaceAlpha', 0.9);
colormap('jet');
hold on;

% Add wireframe
mesh(X(1:4:end, 1:4:end), T(1:4:end, 1:4:end), U_save(1:4:end, 1:4:end)', ...
     'EdgeColor', 'k', 'FaceColor', 'none', 'LineWidth', 0.2);

xlabel('x', 'FontSize', 12);
ylabel('t', 'FontSize', 12);
zlabel('u(x, t)', 'FontSize', 12);
title('1D Wave Equation: Vibrating String', 'FontSize', 12);

view([-60, 25]);
colorbar;

%% Save figure
exportgraphics(fig, output_file, 'ContentType', 'vector');
exportgraphics(fig, strrep(output_file, '.pdf', '.png'), 'Resolution', 300);

fprintf('Figure saved to: %s\n', output_file);

close(fig);
