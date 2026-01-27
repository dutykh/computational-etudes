%% heat_equation_evolution.m
%
% Visualizes the time evolution of the heat equation solution using the
% truncated Fourier series derived in Chapter 2. The initial condition is
% a triangle wave, which demonstrates how higher frequency modes decay
% faster than lower ones.
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
% Fourier coefficients (computed analytically):
%     a_0 = pi/2
%     a_n = (4/pi) / n^2  for odd n, 0 for even n
%     b_n = 0             (symmetric function)
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
N_POINTS = 500;         % Spatial grid points for plotting
TIMES = [0.0, 0.01, 0.05, 0.1, 0.5];  % Time snapshots

% Output path (relative to this script's location)
script_dir = fileparts(mfilename('fullpath'));
output_file = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch02', 'matlab', 'heat_evolution.pdf');

%% Compute Fourier coefficients for the triangle wave
% The Fourier series is:
%     f(x) = pi/2 + (4/pi) * sum_{k=1}^inf cos((2k-1)x) / (2k-1)^2

a0 = pi / 2;
a_n = zeros(1, N_MODES + 1);  % a_n(n) = coefficient of cos(nx), 1-indexed
b_n = zeros(1, N_MODES + 1);  % b_n(n) = coefficient of sin(nx)

for n = 1:N_MODES
    if mod(n, 2) == 1  % odd n only
        a_n(n + 1) = 4.0 / (pi * n^2);
    end
end

%% Spatial grid
x = linspace(0, 2*pi, N_POINTS);

%% Create figure
fig = figure('Units', 'inches', 'Position', [1, 1, 6, 4]);
hold on;

% Color gradient from dark to light for time progression
colors = [linspace(0.08, 0.7, length(TIMES))', ...
          linspace(0.18, 0.8, length(TIMES))', ...
          linspace(0.43, 0.95, length(TIMES))'];

%% Plot solution at each time
legends = cell(1, length(TIMES));
for i = 1:length(TIMES)
    t = TIMES(i);

    % Evaluate the heat equation solution
    % u(x,t) = a_0 + sum_{n=1}^N (a_n cos(nx) + b_n sin(nx)) exp(-n^2 t)
    u = a0 * ones(size(x));
    for n = 1:N_MODES
        decay = exp(-n^2 * t);
        u = u + (a_n(n + 1) * cos(n * x) + b_n(n + 1) * sin(n * x)) * decay;
    end

    plot(x, u, 'Color', colors(i, :), 'LineWidth', 1.5);
    legends{i} = sprintf('$t = %.2g$', t);
end

%% Formatting
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 11);
ylabel('$u(x,t)$', 'Interpreter', 'latex', 'FontSize', 11);
xlim([0, 2*pi]);
ylim([-0.2, 3.5]);
set(gca, 'XTick', [0, pi/2, pi, 3*pi/2, 2*pi]);
set(gca, 'XTickLabel', {'$0$', '$\pi/2$', '$\pi$', '$3\pi/2$', '$2\pi$'}, ...
    'TickLabelInterpreter', 'latex');
leg = legend(legends, 'Location', 'north', 'Interpreter', 'latex', 'NumColumns', 3);
set(leg, 'Box', 'off');
box off;

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
