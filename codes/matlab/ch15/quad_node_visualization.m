% quad_node_visualization.m
% Chapter 15: Quadrature in Spectral Methods
%
% Computational Etude 15.1: Quadrature Node Distributions.
%
% Visualises three classical sets of quadrature nodes on [-1, 1] for n = 32
% points, displayed on parallel horizontal lines:
%
%     1. Newton--Cotes (equispaced) nodes: x_j = -1 + 2j/n, j = 0, ..., n.
%     2. Gauss--Legendre nodes: zeros of P_{n+1}(x).
%     3. Clenshaw--Curtis (Chebyshev) nodes: x_j = cos(j pi / n), j = 0, ..., n.
%
% The figure highlights the fundamental difference in node clustering near
% the endpoints: equispaced nodes remain uniformly distributed, whereas both
% Gauss--Legendre and Clenshaw--Curtis nodes cluster near x = +-1.
%
% Reference: L. N. Trefethen, "Is Gauss Quadrature Better than
%            Clenshaw--Curtis?", SIAM Review 50(1):67--87, 2008, Figure 1.
%
% Generates Figure 15.1: node_visualization.pdf
%
% Author: Dr. Denys Dutykh
%         Mathematics Department
%         Khalifa University of Science and Technology
%         Abu Dhabi, UAE
%
% Part of the book "Computational Etudes: A Spectral Approach"
% https://github.com/dutykh/computational-etudes

clear; close all; clc;

% Publication-quality figure settings
set(0, 'DefaultAxesFontSize', 11);
set(0, 'DefaultTextFontSize', 11);
set(0, 'DefaultLineLineWidth', 1.5);
set(0, 'DefaultAxesLineWidth', 0.8);

% Book color scheme
NAVY   = [20, 45, 110] / 255;
SKY    = [120, 150, 210] / 255;
CORAL  = [231, 76, 60] / 255;
TEAL   = [22, 160, 133] / 255;
PURPLE = [142, 68, 173] / 255;
ORANGE = [230, 126, 34] / 255;

% Output path
output_dir = fullfile(fileparts(mfilename('fullpath')), ...
    '..', '..', '..', 'textbook', 'figures', 'ch15', 'matlab');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

% -------------------------------------------------------------------------
% Number of intervals (n+1 points)
% -------------------------------------------------------------------------
n = 32;

% Newton--Cotes (equispaced) nodes
nc_nodes = linspace(-1, 1, n + 1);

% Gauss--Legendre nodes (n+1 points)
[gl_nodes, ~] = gauss_legendre(n + 1);

% Clenshaw--Curtis nodes
cc_nodes = cos(pi * (0:n)' / n);

% Vertical positions for the three rows
y_nc = 3.0;
y_gl = 2.0;
y_cc = 1.0;

% -------------------------------------------------------------------------
% Create figure
% -------------------------------------------------------------------------
fig = figure('Position', [100, 100, 800, 300]);

hold on;

% Plot horizontal base lines
for y = [y_nc, y_gl, y_cc]
    plot([-1, 1], [y, y], 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5);
end

% Plot nodes as filled circles
ms = 5;

plot(nc_nodes, y_nc * ones(size(nc_nodes)), 'o', ...
    'Color', CORAL, 'MarkerSize', ms, 'MarkerFaceColor', CORAL, ...
    'MarkerEdgeColor', CORAL, 'LineWidth', 0.5);

plot(gl_nodes, y_gl * ones(size(gl_nodes)), 'o', ...
    'Color', NAVY, 'MarkerSize', ms, 'MarkerFaceColor', NAVY, ...
    'MarkerEdgeColor', NAVY, 'LineWidth', 0.5);

plot(cc_nodes, y_cc * ones(size(cc_nodes)), 'o', ...
    'Color', TEAL, 'MarkerSize', ms, 'MarkerFaceColor', TEAL, ...
    'MarkerEdgeColor', TEAL, 'LineWidth', 0.5);

% Labels on the left
label_x = -1.25;
text(label_x, y_nc, 'Newton-Cotes', 'FontSize', 11, ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
    'Color', CORAL, 'FontWeight', 'bold');
text(label_x, y_gl, 'Gauss-Legendre', 'FontSize', 11, ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
    'Color', NAVY, 'FontWeight', 'bold');
text(label_x, y_cc, 'Clenshaw-Curtis', 'FontSize', 11, ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
    'Color', TEAL, 'FontWeight', 'bold');

% Mark the interval endpoints
for y = [y_nc, y_gl, y_cc]
    plot([-1, 1], [y, y], '|', 'Color', 'k', 'MarkerSize', 8, ...
        'LineWidth', 1.0);
end

hold off;

% Axis formatting
xlim([-1.8, 1.15]);
ylim([0.4, 3.6]);
xlabel('$x$', 'Interpreter', 'latex');
set(gca, 'XTick', [-1, -0.5, 0, 0.5, 1]);
set(gca, 'YTick', []);

% Remove y-axis and top/right spines
box off;
ax = gca;
ax.YAxis.Visible = 'off';

% Title
title('Quadrature nodes on $[-1, 1]$ with $n = 32$', ...
    'Interpreter', 'latex', 'FontSize', 12);

% Save output
exportgraphics(gcf, fullfile(output_dir, 'node_visualization.pdf'), 'ContentType', 'vector');
exportgraphics(gcf, fullfile(output_dir, 'node_visualization.png'), 'Resolution', 300);
fprintf('Figure saved to %s\n', fullfile(output_dir, 'node_visualization.pdf'));

% =========================================================================
% Local functions
% =========================================================================

function [x, w] = gauss_legendre(n)
    % n-point Gauss-Legendre nodes and weights via Golub-Welsch.
    if n == 1
        x = 0; w = 2; return;
    end
    j = (1:n-1)';
    beta = j ./ sqrt(4*j.^2 - 1);
    J = diag(beta, 1) + diag(beta, -1);
    [V, D] = eig(J);
    [x, idx] = sort(diag(D));
    w = 2 * V(1, idx)'.^2;
end
