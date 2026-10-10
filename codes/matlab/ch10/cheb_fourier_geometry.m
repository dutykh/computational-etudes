%% cheb_fourier_geometry.m
%
% Figure 10.1: Chebyshev--Fourier geometry visualization
%
% Left panel: Unit circle with angle theta, projection to x-axis showing
% Chebyshev point clustering near boundaries.
%
% Right panel: Wrapped cosine interpretation of Chebyshev polynomials.
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
% Part of "Computational Etudes: A Spectral Approach"
% https://github.com/dutykh/computational-etudes

clear; close all; clc;

%% Configuration
N = 16;  % Number of Chebyshev points

% Colors
NAVY = [0.078 0.176 0.431];
CORAL = [0.906 0.298 0.235];
TEAL = [0.086 0.627 0.522];
GRAY = [0.5 0.5 0.5];

% Output path
script_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch10', 'matlab');
output_file = fullfile(output_dir, 'cheb_fourier_geometry.pdf');

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% Chebyshev points
theta = pi * (0:N)' / N;
x_cheb = cos(theta);

%% Create figure
fig = figure('Units', 'inches', 'Position', [1, 1, 12, 5]);

%% Panel 1: Circle projection
subplot(1, 2, 1);
hold on;

% Draw unit semicircle
theta_circle = linspace(0, pi, 200);
plot(cos(theta_circle), sin(theta_circle), '-', 'Color', GRAY, 'LineWidth', 1.5);

% Draw diameter
plot([-1, 1], [0, 0], '-', 'Color', GRAY, 'LineWidth', 1);

% Draw projection lines and points
for j = 1:length(theta)
    % Point on circle
    xc = cos(theta(j));
    yc = sin(theta(j));

    % Plot projection line
    plot([xc, xc], [yc, 0], ':', 'Color', TEAL, 'LineWidth', 0.8);

    % Plot point on circle
    plot(xc, yc, 'o', 'Color', CORAL, 'MarkerSize', 8, 'MarkerFaceColor', CORAL);

    % Plot point on x-axis (Chebyshev point)
    plot(xc, 0, 'o', 'Color', NAVY, 'MarkerSize', 8, 'MarkerFaceColor', NAVY);
end

% Draw angle arc for one point
j_show = 5;
theta_arc = linspace(0, theta(j_show), 30);
r_arc = 0.3;
plot(r_arc * cos(theta_arc), r_arc * sin(theta_arc), '-', 'Color', CORAL, 'LineWidth', 1.5);
text(r_arc * cos(theta(j_show)/2) + 0.1, r_arc * sin(theta(j_show)/2) + 0.05, ...
     '\theta_j', 'FontSize', 12, 'Color', CORAL);

% Labels
text(0, 1.15, 'x = cos(\theta)', 'FontSize', 11, 'HorizontalAlignment', 'center');
xlabel('x', 'FontSize', 11);

axis equal;
xlim([-1.3, 1.3]);
ylim([-0.2, 1.3]);
title('Chebyshev Points from Circle Projection', 'FontSize', 12);
box off;

%% Panel 2: Chebyshev polynomial as wrapped cosine
subplot(1, 2, 2);
hold on;

% T_n(x) = cos(n * arccos(x)) for n = 4
n = 4;

% In theta space
theta_fine = linspace(0, pi, 200);
y_theta = cos(n * theta_fine);
plot(theta_fine, y_theta, '-', 'Color', NAVY, 'LineWidth', 2);

% Mark the Chebyshev points
for j = 1:length(theta)
    y_point = cos(n * theta(j));
    plot(theta(j), y_point, 'o', 'Color', CORAL, 'MarkerSize', 8, 'MarkerFaceColor', CORAL);
end

% Grid lines
for k = 0:N
    xline(k * pi / N, ':', 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5);
end
yline(0, '-', 'Color', GRAY, 'LineWidth', 0.5);

xlabel('\theta', 'FontSize', 11);
ylabel(sprintf('T_%d(cos(\\theta)) = cos(%d\\theta)', n, n), 'FontSize', 11);
title(sprintf('Chebyshev Polynomial T_%d as Cosine Wave', n), 'FontSize', 12);

xlim([0, pi]);
ylim([-1.2, 1.2]);
box on;
grid on;

% Add text annotation
text(pi/2, -1.5, 'T_n(x) = cos(n \cdot arccos(x))', 'FontSize', 10, ...
     'HorizontalAlignment', 'center', 'Color', NAVY);

%% Save figure
exportgraphics(fig, output_file, 'ContentType', 'vector');
exportgraphics(fig, strrep(output_file, '.pdf', '.png'), 'Resolution', 300);

fprintf('Figure saved to: %s\n', output_file);

close(fig);
