%% cheb_cardinal.m - Chebyshev cardinal functions visualization
%
% Visualizes Chebyshev cardinal functions (Lagrange interpolation basis) and
% shows how the differentiation matrix entries are the derivatives of these
% cardinal functions evaluated at the grid points.
%
% This script generates Figure 6.3 for Chapter 6.
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
output_dir = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch07', 'matlab');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% Configuration
N = 10;  % Number of intervals
N_FINE = 500;  % Fine grid for plotting

%% Chebyshev grid points
[D, x_nodes] = cheb_matrix(N);
x_fine = linspace(-1, 1, N_FINE)';

%% Create figure
fig = figure('Position', [100, 100, 1100, 500]);

% Panel 1: Several cardinal functions
subplot(1, 2, 1);

% Select cardinal functions to display
indices = [1, 4, 6, N+1];  % 0-indexed: 0, 3, 5, N
colors = {CORAL, TEAL, PURPLE, ORANGE};
linestyles = {'-', '--', '-.', ':'};

for k = 1:length(indices)
    j = indices(k);  % 1-indexed
    L_j = chebyshev_cardinal(x_fine, x_nodes, j);
    plot(x_fine, L_j, linestyles{k}, 'Color', colors{k}, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('$L_{%d}(x)$', j-1));
    hold on;
    % Mark the cardinal point where L_j = 1
    plot(x_nodes(j), 1.0, 'o', 'Color', colors{k}, 'MarkerSize', 7, ...
        'MarkerFaceColor', colors{k}, 'MarkerEdgeColor', 'white', ...
        'HandleVisibility', 'off');
end

% Mark all nodes on x-axis
plot(x_nodes, zeros(N+1, 1), 'o', 'Color', NAVY, 'MarkerSize', 5, ...
    'MarkerFaceColor', NAVY, 'MarkerEdgeColor', 'white', ...
    'HandleVisibility', 'off');

% Reference lines
yline(0, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5);
yline(1, '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5);
hold off;

xlabel('$x$');
ylabel('$L_j(x)$');
title(sprintf('Chebyshev Cardinal Functions ($N = %d$)', N));
legend('Location', 'northeast', 'FontSize', 9, 'NumColumns', 2);
xlim([-1.05, 1.05]);
ylim([-0.3, 1.2]);
box off;
grid on;
set(gca, 'GridAlpha', 0.2);

% Add annotation
text(0.5, 0.85, '$L_j(x_k) = \delta_{jk}$', 'Units', 'normalized', ...
    'FontSize', 10, 'Color', [0.4, 0.4, 0.4]);

% Panel 2: Cardinal function with derivative slopes at nodes
subplot(1, 2, 2);

% Plot cardinal function L_5 (index 6 in MATLAB)
j = 6;  % L_5
L_j = chebyshev_cardinal(x_fine, x_nodes, j);
plot(x_fine, L_j, '-', 'Color', TEAL, 'LineWidth', 2, ...
    'DisplayName', sprintf('$L_{%d}(x)$', j-1));
hold on;

% Mark the cardinal point
plot(x_nodes(j), 1.0, 'o', 'Color', TEAL, 'MarkerSize', 10, ...
    'MarkerFaceColor', TEAL, 'MarkerEdgeColor', 'white', ...
    'HandleVisibility', 'off');

% Draw tangent lines at grid points (slopes from D matrix)
for k = 1:2:N+1  % Every other point for clarity
    slope = D(k, j);  % This is L_j'(x_k)
    x_k = x_nodes(k);
    y_k = chebyshev_cardinal(x_k, x_nodes, j);

    % Draw short tangent line
    dx = 0.15;
    x_tan = [x_k - dx, x_k + dx];
    y_tan = y_k + slope * (x_tan - x_k);

    plot(x_tan, y_tan, '-', 'Color', CORAL, 'LineWidth', 1.5, ...
        'HandleVisibility', 'off');
    plot(x_k, y_k, 's', 'Color', NAVY, 'MarkerSize', 4, ...
        'MarkerFaceColor', NAVY, 'HandleVisibility', 'off');
end

% Reference line
yline(0, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5);
hold off;

xlabel('$x$');
ylabel('$L_j(x)$');
title(sprintf('Cardinal Function $L_{%d}$ with Derivative Slopes', j-1));
xlim([-1.05, 1.05]);
ylim([-0.5, 1.2]);
box off;
grid on;
set(gca, 'GridAlpha', 0.2);

% Add annotation
text(0.02, 0.15, sprintf('Slopes = column %d of $D_N$:\n$D_{kj} = L_j''(x_k)$', j-1), ...
    'Units', 'normalized', 'FontSize', 9, 'Color', [0.4, 0.4, 0.4], ...
    'BackgroundColor', 'white', 'EdgeColor', [0.8, 0.8, 0.8]);

%% Save figure
output_file = fullfile(output_dir, 'cheb_cardinal');
exportgraphics(fig, [output_file, '.pdf'], 'ContentType', 'vector');
exportgraphics(fig, [output_file, '.png'], 'Resolution', 300);
fprintf('Figure saved to: %s.pdf\n', output_file);

%% Verify cardinal property
fprintf('\nVerification of cardinal property (N = %d):\n', N);
fprintf('%s\n', repmat('-', 1, 50));
fprintf('%4s %12s %20s\n', 'j', 'L_j(x_j)', 'max|L_j(x_k)|, k!=j');
fprintf('%s\n', repmat('-', 1, 50));
for j = 1:N+1
    L_j_at_nodes = chebyshev_cardinal(x_nodes, x_nodes, j);
    val_at_own = L_j_at_nodes(j);
    others = L_j_at_nodes;
    others(j) = [];
    max_at_others = max(abs(others));
    fprintf('%4d %12.6f %20.2e\n', j-1, val_at_own, max_at_others);
end
fprintf('%s\n', repmat('-', 1, 50));

%% Local function: Chebyshev cardinal function using barycentric interpolation
function L_j = chebyshev_cardinal(x_eval, x_nodes, j)
    N = length(x_nodes) - 1;

    % Barycentric weights for Chebyshev points
    % w_k = (-1)^k * (1/2 if k=0,N else 1)
    w = ones(N + 1, 1);
    w(1) = 0.5;
    w(N + 1) = 0.5;
    w = w .* ((-1).^(0:N)');

    L_j = zeros(size(x_eval));

    for i = 1:length(x_eval)
        x = x_eval(i);
        diffs = x - x_nodes;

        if any(abs(diffs) < 1e-14)
            % x is at a node
            [~, k] = min(abs(diffs));
            L_j(i) = double(k == j);
        else
            % Barycentric formula
            terms = w ./ diffs;
            L_j(i) = (w(j) / diffs(j)) / sum(terms);
        end
    end
end
