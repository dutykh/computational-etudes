%% chebyshev_points_circle.m
%
% Visualizes the geometric construction of Chebyshev points from the unit circle.
%
% Chebyshev-Gauss-Lobatto points are obtained by:
% 1. Placing N+1 equally spaced points on the upper half of the unit circle
% 2. Projecting these points vertically onto the x-axis
%
% This reveals why Chebyshev points cluster near ±1.
%
% Author: Dr. Denys Dutykh
%         Mathematics Department
%         Khalifa University of Science and Technology
%         Abu Dhabi, UAE
%
% Part of the book "Computational Études: A Spectral Approach"

clear; close all; clc;

%% Configuration
N = 8;  % Number of intervals

% Book color scheme
NAVY = [20, 45, 110] / 255;
SKY = [120, 150, 210] / 255;
CORAL = [231, 76, 60] / 255;
TEAL = [22, 160, 133] / 255;
LIGHT_GRAY = [0.8, 0.8, 0.8];

% Output path
script_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch04', 'matlab');
output_file = fullfile(output_dir, 'chebyshev_points_circle.pdf');

%% Create figure
fig = figure('Units', 'inches', 'Position', [1, 1, 8, 5]);
hold on;

% Draw unit circle (upper half)
theta_circle = linspace(0, pi, 200);
x_circle = cos(theta_circle);
y_circle = sin(theta_circle);
plot(x_circle, y_circle, 'Color', NAVY, 'LineWidth', 1.5);

% Draw x-axis
plot([-1.3, 1.3], [0, 0], 'Color', [0.4, 0.4, 0.4], 'LineWidth', 0.8);

% Chebyshev points: angles equally spaced on [0, π]
j = 0:N;
theta_j = j * pi / N;

% Points on the circle
x_circle_pts = cos(theta_j);
y_circle_pts = sin(theta_j);

% Points on x-axis (Chebyshev nodes)
x_cheb = cos(theta_j);
y_cheb = zeros(size(x_cheb));

% Draw projection lines
for i = 1:N+1
    plot([x_circle_pts(i), x_cheb(i)], [y_circle_pts(i), y_cheb(i)], ...
         '--', 'Color', LIGHT_GRAY, 'LineWidth', 0.8);
end

% Draw points on circle
scatter(x_circle_pts, y_circle_pts, 60, SKY, 'filled', ...
        'MarkerEdgeColor', NAVY, 'LineWidth', 1);

% Draw Chebyshev points on axis
scatter(x_cheb, y_cheb, 80, CORAL, 's', 'filled', ...
        'MarkerEdgeColor', [0.7, 0.2, 0.15], 'LineWidth', 1.5);

% Draw angle arcs
for i = [2, 3, 4]
    if theta_j(i) > 0.1
        theta_arc = linspace(0, theta_j(i), 50);
        plot(0.15 * cos(theta_arc), 0.15 * sin(theta_arc), 'Color', TEAL, 'LineWidth', 1);
    end
end

% Annotate angle
text(0.25, 0.15, '$\theta_j = \frac{j\pi}{N}$', 'Interpreter', 'latex', ...
     'FontSize', 10, 'Color', TEAL);

% Label some points
indices_to_label = [1, 3, N/2+1, N-1, N+1];
for idx = indices_to_label
    j_val = idx - 1;
    if y_circle_pts(idx) > 0.1
        offset_y = 0.12;
    else
        offset_y = -0.18;
    end
    text(x_cheb(idx), offset_y, sprintf('$j=%d$', j_val), ...
         'Interpreter', 'latex', 'FontSize', 8, 'HorizontalAlignment', 'center', ...
         'Color', [0.3, 0.3, 0.3]);
end

% Add formula box
text(0.98, 0.95, {'$x_j = \cos\left(\frac{j\pi}{N}\right)$', '$j = 0, 1, \ldots, N$'}, ...
     'Units', 'normalized', 'Interpreter', 'latex', 'FontSize', 11, ...
     'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
     'BackgroundColor', 'white', 'EdgeColor', NAVY, 'LineWidth', 0.5);

% Annotations
text(-0.5, 0.8, 'Equal spacing on circle', 'FontSize', 9, 'Color', SKY, ...
     'HorizontalAlignment', 'center', 'FontAngle', 'italic');
text(0.7, -0.25, 'Clustering near boundaries', 'FontSize', 9, 'Color', CORAL, ...
     'HorizontalAlignment', 'center', 'FontAngle', 'italic');

% Styling
axis equal;
xlim([-1.4, 1.4]);
ylim([-0.4, 1.2]);
xlabel('$x$', 'Interpreter', 'latex');
title(sprintf('Geometric Construction of Chebyshev Points ($N = %d$)', N), 'FontSize', 11);

set(gca, 'FontSize', 10, 'TickLabelInterpreter', 'latex');
set(gca, 'XTick', [-1, -0.5, 0, 0.5, 1]);
set(gca, 'YTick', [0, 0.5, 1]);

% Legend
h1 = scatter(nan, nan, 60, SKY, 'filled', 'MarkerEdgeColor', NAVY);
h2 = scatter(nan, nan, 80, CORAL, 's', 'filled');
legend([h1, h2], {'Points on circle', 'Chebyshev points'}, ...
       'Location', 'northwest', 'Interpreter', 'latex');

%% Save figure
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

exportgraphics(fig, output_file, 'ContentType', 'vector');
fprintf('Figure saved to: %s\n', output_file);

png_file = strrep(output_file, '.pdf', '.png');
exportgraphics(fig, png_file, 'Resolution', 300);

%% Print node positions
fprintf('\nChebyshev-Gauss-Lobatto points for N = %d:\n', N);
fprintf('%s\n', repmat('-', 1, 40));
fprintf('%4s %12s %12s\n', 'j', 'theta_j', 'x_j');
fprintf('%s\n', repmat('-', 1, 40));
for i = 1:N+1
    j_val = i - 1;
    fprintf('%4d %12s %12.6f\n', j_val, sprintf('%d*pi/%d', j_val, N), x_cheb(i));
end
fprintf('%s\n', repmat('-', 1, 40));

% Show clustering effect
fprintf('\nSpacing between consecutive points:\n');
spacings = diff(sort(x_cheb));
for i = 1:length(spacings)
    fprintf('Interval %d: dx = %.6f\n', i, spacings(i));
end
fprintf('Ratio (boundary/center): %.2f\n', spacings(1)/spacings(round(N/2)));
