%% wave_equation_evolution.m
%
% Visualizes the time evolution of the wave equation solution using the
% truncated Fourier sine series derived in Chapter 2. The initial condition is
% a plucked string (triangular displacement), which demonstrates standing wave
% oscillations.
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
N_POINTS = 500;         % Spatial grid points for plotting
L = pi;                 % Domain length [0, L]
C = 1.0;                % Wave speed
H = 1.0;                % Height of plucked string at center

% Period of oscillation: T = 2L/c
T_PERIOD = 2 * L / C;

% Time snapshots: show half a period of oscillation
TIMES = [0.0, T_PERIOD/8, T_PERIOD/4, 3*T_PERIOD/8, T_PERIOD/2];

% Output path (relative to this script's location)
script_dir = fileparts(mfilename('fullpath'));
output_file = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch02', 'matlab', 'wave_evolution.pdf');

%% Compute Fourier coefficients for the plucked string
% The Fourier sine series is:
%     f(x) = sum a_n sin(n*pi*x/L)
%
% where a_n = (8h/(n^2*pi^2)) * sin(n*pi/2)
%           = 8h/(n^2*pi^2) for n = 1, 5, 9, ...  (n mod 4 = 1)
%           = -8h/(n^2*pi^2) for n = 3, 7, 11, ... (n mod 4 = 3)
%           = 0 for even n

a_n = zeros(1, N_MODES + 1);  % a_n(n) = coefficient of sin(n*pi*x/L), 1-indexed
b_n = zeros(1, N_MODES + 1);  % b_n(n) = 0 since initially at rest

for n = 1:N_MODES
    if mod(n, 2) == 1  % odd n only
        % sin(n*pi/2) = 1 for n = 1, 5, 9, ... and -1 for n = 3, 7, 11, ...
        if mod(n, 4) == 1
            sgn = 1;
        else
            sgn = -1;
        end
        a_n(n + 1) = sgn * 8.0 * H / (n^2 * pi^2);
    end
end

%% Spatial grid
x = linspace(0, L, N_POINTS);

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

    % Evaluate the wave equation solution
    % u(x,t) = sum_{n=1}^N (a_n cos(omega_n t) + b_n sin(omega_n t)) sin(n*pi*x/L)
    u = zeros(size(x));
    for n = 1:N_MODES
        omega_n = C * n * pi / L;
        spatial = sin(n * pi * x / L);
        temporal = a_n(n + 1) * cos(omega_n * t) + b_n(n + 1) * sin(omega_n * t);
        u = u + temporal * spatial;
    end

    plot(x, u, 'Color', colors(i, :), 'LineWidth', 1.5);

    % Create time label as fraction of period
    if t == 0
        legends{i} = '$t = 0$';
    elseif abs(t - T_PERIOD/8) < 0.001
        legends{i} = '$t = T/8$';
    elseif abs(t - T_PERIOD/4) < 0.001
        legends{i} = '$t = T/4$';
    elseif abs(t - 3*T_PERIOD/8) < 0.001
        legends{i} = '$t = 3T/8$';
    elseif abs(t - T_PERIOD/2) < 0.001
        legends{i} = '$t = T/2$';
    else
        legends{i} = sprintf('$t = %.3f$', t);
    end
end

%% Add reference line at y=0
plot([0, L], [0, 0], 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5, 'LineStyle', '-');

%% Formatting
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 11);
ylabel('$u(x,t)$', 'Interpreter', 'latex', 'FontSize', 11);
xlim([0, L]);
ylim([-1.2, 1.2]);
set(gca, 'XTick', [0, L/4, L/2, 3*L/4, L]);
set(gca, 'XTickLabel', {'$0$', '$\pi/4$', '$\pi/2$', '$3\pi/4$', '$\pi$'}, ...
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
