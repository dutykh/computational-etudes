%% wave_equation_waterfall.m
%
% Waterfall (3D surface) plot showing the complete space-time evolution of the
% wave equation solution. This visualization demonstrates the standing wave
% oscillations of a plucked string over time.
%
% Wave equation (Dirichlet BCs on [0, L]):
%     u_tt = c^2 u_xx
%
% Analytical solution:
%     u(x,t) = sum (a_n cos(omega_n t) + b_n sin(omega_n t)) sin(n*pi*x/L)
%
% where omega_n = c*n*pi/L
%
% Plucked string initial condition (triangle):
%     f(x) = (2h/L)*x           for 0 <= x <= L/2
%     f(x) = 2h*(1 - x/L)       for L/2 <= x <= L
%     g(x) = 0                  (initially at rest)
%
% Author: Dr. Denys Dutykh
%         Mathematics Department
%         Khalifa University of Science and Technology
%         Abu Dhabi, UAE
%
% Part of the book "Computational Études: A Spectral Approach"
% https://github.com/dutykh/computational-etudes

clear; clc; close all;

%% Configuration
N_MODES = 50;           % Number of Fourier modes in truncation
NX = 150;               % Spatial grid points
NT = 100;               % Time grid points
L = pi;                 % Domain length [0, L]
C = 1.0;                % Wave speed
H = 1.0;                % Height of plucked string at center

% Period of oscillation: T = 2L/c
T_PERIOD = 2 * L / C;
T_MAX = T_PERIOD;       % Show one full period

% Output path (relative to this script's location)
script_dir = fileparts(mfilename('fullpath'));
output_file = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch02', 'matlab', 'wave_waterfall.pdf');

%% Compute Fourier coefficients for the plucked string
a_n = zeros(1, N_MODES + 1);
b_n = zeros(1, N_MODES + 1);

for n = 1:N_MODES
    if mod(n, 2) == 1  % odd n only
        if mod(n, 4) == 1
            sgn = 1;
        else
            sgn = -1;
        end
        a_n(n + 1) = sgn * 8.0 * H / (n^2 * pi^2);
    end
end

%% Create 2D grid for space-time
x = linspace(0, L, NX);
t = linspace(0, T_MAX, NT);
[X, T] = meshgrid(x, t);

%% Evaluate solution on grid
U = zeros(size(X));
for i = 1:NT
    ti = t(i);
    u = zeros(1, NX);
    for n = 1:N_MODES
        omega_n = C * n * pi / L;
        spatial = sin(n * pi * x / L);
        temporal = a_n(n + 1) * cos(omega_n * ti) + b_n(n + 1) * sin(omega_n * ti);
        u = u + temporal * spatial;
    end
    U(i, :) = u;
end

%% Create figure
fig = figure('Units', 'inches', 'Position', [1, 1, 7, 5]);

% Custom colormap: navy -> white -> sky blue
n_colors = 256;
navy_rgb = [20, 45, 110] / 255;
white_rgb = [1, 1, 1];
sky_rgb = [120, 150, 210] / 255;

cmap = zeros(n_colors, 3);
for i = 1:n_colors/2
    t_val = (i - 1) / (n_colors/2 - 1);
    cmap(i, :) = navy_rgb * (1 - t_val) + white_rgb * t_val;
end
for i = 1:n_colors/2
    t_val = (i - 1) / (n_colors/2 - 1);
    cmap(n_colors/2 + i, :) = white_rgb * (1 - t_val) + sky_rgb * t_val;
end
colormap(cmap);

% Plot the surface
surf(X, T, U, 'EdgeColor', [0.5, 0.5, 0.5], 'EdgeAlpha', 0.3, ...
     'FaceAlpha', 0.9, 'FaceLighting', 'gouraud');
hold on;

% Labels
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 11);
ylabel('$t$', 'Interpreter', 'latex', 'FontSize', 11);
zlabel('$u(x,t)$', 'Interpreter', 'latex', 'FontSize', 11);

% Custom x-ticks with pi notation
set(gca, 'XTick', [0, L/4, L/2, 3*L/4, L]);
set(gca, 'XTickLabel', {'$0$', '$\pi/4$', '$\pi/2$', '$3\pi/4$', '$\pi$'}, ...
    'TickLabelInterpreter', 'latex');

% Custom t-ticks as fractions of period
set(gca, 'YTick', [0, T_PERIOD/4, T_PERIOD/2, 3*T_PERIOD/4, T_PERIOD]);
set(gca, 'YTickLabel', {'$0$', '$T/4$', '$T/2$', '$3T/4$', '$T$'}, ...
    'TickLabelInterpreter', 'latex');

% Set viewing angle
view(-60, 25);

% Add lighting for better 3D perception
camlight('headlight');
lighting gouraud;

%% Save figure
% Ensure output directory exists
[output_dir, ~, ~] = fileparts(output_file);
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Export to PDF
exportgraphics(fig, output_file, 'ContentType', 'vector', 'Resolution', 300);
fprintf('Figure saved to: %s\n', output_file);

close(fig);
