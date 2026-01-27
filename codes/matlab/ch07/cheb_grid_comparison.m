%% cheb_grid_comparison.m - Equispaced vs Chebyshev grid comparison
%
% Visualizes the comparison between equispaced and Chebyshev-Gauss-Lobatto grids.
% Shows why Chebyshev points cluster near the boundaries and how this relates
% to the projection from a circle.
%
% This script generates Figure 6.1 for Chapter 6.
%
% Author: Dr. Denys Dutykh
%         Mathematics Department
%         Khalifa University of Science and Technology
%         Abu Dhabi, UAE
%
% Part of the book "Computational Etudes: A Spectral Approach"
% https://github.com/dutykh/computational-etudes
%
% Last modified: January 2026

clear;
close all;
clc;

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
NAVY = [20, 45, 110] / 255;
SKY = [120, 150, 210] / 255;
CORAL = [231, 76, 60] / 255;
TEAL = [22, 160, 133] / 255;
PURPLE = [142, 68, 173] / 255;
ORANGE = [230, 126, 34] / 255;

%% Output directory
script_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch06', 'matlab');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% Configuration
N = 16;  % Number of intervals

%% Grid points
x_equi = linspace(-1, 1, N + 1)';  % Equispaced
x_cheb = cos(pi * (0:N)' / N);     % Chebyshev-GL

%% Create figure with four panels
fig = figure('Position', [100, 100, 1000, 700]);

% Panel 1: Equispaced grid (top)
subplot(3, 2, [1, 2]);
plot(x_equi, zeros(size(x_equi)), 'o', 'Color', CORAL, 'MarkerSize', 8, ...
    'MarkerFaceColor', CORAL, 'MarkerEdgeColor', 'white');
hold on;
yline(0, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5);
hold off;

xlim([-1.1, 1.1]);
ylim([-0.2, 0.2]);
xlabel('$x$');
title(sprintf('Equispaced Grid ($N = %d$ intervals)', N));
set(gca, 'YTick', []);
box off;

% Panel 2: Chebyshev grid (middle)
subplot(3, 2, [3, 4]);
plot(x_cheb, zeros(size(x_cheb)), 'o', 'Color', TEAL, 'MarkerSize', 8, ...
    'MarkerFaceColor', TEAL, 'MarkerEdgeColor', 'white');
hold on;
yline(0, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5);
hold off;

xlim([-1.1, 1.1]);
ylim([-0.2, 0.2]);
xlabel('$x$');
title(sprintf('Chebyshev-Gauss-Lobatto Grid ($N = %d$ intervals)', N));
set(gca, 'YTick', []);
box off;

% Panel 3: Circle projection (bottom-left)
subplot(3, 2, 5);

% Draw circle
theta = linspace(0, pi, 200);
plot(cos(theta), sin(theta), '-', 'Color', SKY, 'LineWidth', 2);
hold on;

% Draw projection lines and points for Chebyshev
theta_cheb = pi * (0:N)' / N;
for i = 1:length(theta_cheb)
    t = theta_cheb(i);
    x = cos(t);
    y = sin(t);
    % Point on circle
    plot(x, y, 'o', 'Color', TEAL, 'MarkerSize', 6, 'MarkerFaceColor', TEAL);
    % Projection line
    plot([x, x], [y, 0], '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5);
    % Point on x-axis
    plot(x, 0, 's', 'Color', NAVY, 'MarkerSize', 5, 'MarkerFaceColor', NAVY);
end

xline(0, 'Color', [0.3, 0.3, 0.3], 'LineWidth', 0.5);
yline(0, 'Color', [0.3, 0.3, 0.3], 'LineWidth', 0.8);
hold off;

xlim([-1.3, 1.3]);
ylim([-0.2, 1.3]);
axis equal;
xlabel('$x = \cos\theta$');
ylabel('$\sin\theta$');
title('Circle Projection: $x_j = \cos(j\pi/N)$');
box off;

% Panel 4: Boundary spacing comparison (bottom-right)
subplot(3, 2, 6);

% Compute spacing for both grids
spacing_equi = diff(sort(x_equi));
spacing_cheb = diff(sort(x_cheb));

% Bar plot of spacings near boundary
n_show = 8;  % Show first 8 spacings
x_pos = 1:n_show;
width = 0.35;

bar(x_pos - width/2, spacing_equi(1:n_show), width, 'FaceColor', CORAL, ...
    'EdgeColor', 'white', 'FaceAlpha', 0.8, 'DisplayName', 'Equispaced');
hold on;
bar(x_pos + width/2, spacing_cheb(1:n_show), width, 'FaceColor', TEAL, ...
    'EdgeColor', 'white', 'FaceAlpha', 0.8, 'DisplayName', 'Chebyshev');
hold off;

xlabel('Interval index (from left boundary)');
ylabel('Spacing $\Delta x$');
title('Grid Spacing Near Boundary');
legend('Location', 'northeast', 'FontSize', 9);
box off;

% Main title
sgtitle('Equispaced vs. Chebyshev-Gauss-Lobatto Grids', 'FontSize', 13);

%% Save figure
output_file = fullfile(output_dir, 'grid_comparison');
exportgraphics(fig, [output_file, '.pdf'], 'ContentType', 'vector');
exportgraphics(fig, [output_file, '.png'], 'Resolution', 300);
fprintf('Figure saved to: %s.pdf\n', output_file);

%% Print spacing information
fprintf('\nGrid spacing analysis (N = %d):\n', N);
fprintf('%s\n', repmat('-', 1, 50));
fprintf('%-12s %12s %12s %12s\n', 'Spacing', 'Equispaced', 'Chebyshev', 'Ratio');
fprintf('%s\n', repmat('-', 1, 50));
fprintf('%-12s %12.6f %12.6f %12.2f\n', 'Minimum', min(spacing_equi), min(spacing_cheb), ...
    min(spacing_equi)/min(spacing_cheb));
fprintf('%-12s %12.6f %12.6f %12.2f\n', 'Maximum', max(spacing_equi), max(spacing_cheb), ...
    max(spacing_equi)/max(spacing_cheb));
fprintf('%-12s %12.6f %12.6f %12.2f\n', 'At boundary', spacing_equi(1), spacing_cheb(1), ...
    spacing_equi(1)/spacing_cheb(1));
fprintf('%s\n', repmat('-', 1, 50));
fprintf('Theoretical boundary spacing: pi^2/(2*N^2) = %.6f\n', pi^2/(2*N^2));
