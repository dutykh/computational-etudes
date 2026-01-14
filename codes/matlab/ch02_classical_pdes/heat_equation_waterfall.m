%% heat_equation_waterfall.m
%
% Waterfall (3D surface) plot showing the complete space-time evolution of the
% heat equation solution. This visualization demonstrates how the initial
% triangle wave smooths out over time as higher frequency modes decay faster.
%
% Heat equation (periodic BCs on [0, 2*pi]):
%     u_t = u_xx
%
% Analytical solution:
%     u(x,t) = a_0 + sum (a_n cos(nx) + b_n sin(nx)) exp(-n^2 t)
%
% Triangle wave initial condition:
%     f(x) = pi - |x - pi|  on [0, 2*pi]
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
T_MAX = 0.3;            % Maximum time

% Output path (relative to this script's location)
script_dir = fileparts(mfilename('fullpath'));
output_file = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch02', 'matlab', 'heat_waterfall.pdf');

%% Compute Fourier coefficients for the triangle wave
a0 = pi / 2;
a_n = zeros(1, N_MODES + 1);
b_n = zeros(1, N_MODES + 1);

for n = 1:N_MODES
    if mod(n, 2) == 1  % odd n only
        a_n(n + 1) = 4.0 / (pi * n^2);
    end
end

%% Create 2D grid for space-time
x = linspace(0, 2*pi, NX);
t = linspace(0, T_MAX, NT);
[X, T] = meshgrid(x, t);

%% Evaluate solution on grid
U = zeros(size(X));
for i = 1:NT
    ti = t(i);
    u = a0 * ones(1, NX);
    for n = 1:N_MODES
        decay = exp(-n^2 * ti);
        u = u + (a_n(n + 1) * cos(n * x) + b_n(n + 1) * sin(n * x)) * decay;
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
set(gca, 'XTick', [0, pi/2, pi, 3*pi/2, 2*pi]);
set(gca, 'XTickLabel', {'$0$', '$\pi/2$', '$\pi$', '$3\pi/2$', '$2\pi$'}, ...
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
