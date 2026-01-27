%% cheb_matrix_structure.m - Chebyshev differentiation matrix structure
%
% Visualizes the structure of the Chebyshev differentiation matrix: heatmap
% showing the matrix entries and row profile showing the decay pattern.
%
% This script generates Figure 6.2 for Chapter 6.
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
N = 16;

%% Construct differentiation matrix
[D, x] = cheb_matrix(N);

%% Create figure with two panels
fig = figure('Position', [100, 100, 1100, 500]);

% Panel 1: Heatmap of D matrix
subplot(1, 2, 1);
vmax = max(abs(D(:)));

imagesc(D);
colormap(redblue(256));
caxis([-vmax, vmax]);
cb = colorbar;
cb.Label.String = 'Matrix entry value';

title(sprintf('Chebyshev Differentiation Matrix $D_N$ ($N = %d$)', N));
xlabel('Column index $j$');
ylabel('Row index $i$');
axis equal tight;

% Set tick labels
tick_pos = [1, N/4+1, N/2+1, 3*N/4+1, N+1];
tick_labels = {'0', num2str(N/4), num2str(N/2), num2str(3*N/4), num2str(N)};
set(gca, 'XTick', tick_pos, 'XTickLabel', tick_labels);
set(gca, 'YTick', tick_pos, 'YTickLabel', tick_labels);

% Add annotations for corner entries
text(2, 4, sprintf('$D_{00} = %.1f$', D(1,1)), 'FontSize', 9, 'Color', 'white', ...
    'BackgroundColor', NAVY, 'Margin', 2);
text(N-3, N-2, sprintf('$D_{N,N} = %.1f$', D(N+1,N+1)), 'FontSize', 9, 'Color', 'white', ...
    'BackgroundColor', NAVY, 'Margin', 2);

% Panel 2: Row profiles
subplot(1, 2, 2);

j_indices = 0:N;

% First row (boundary)
stem(j_indices, D(1, :), 'filled', 'Color', CORAL, 'MarkerFaceColor', CORAL, ...
    'MarkerSize', 4, 'DisplayName', 'Row $i=0$ (boundary)');
hold on;

% Middle row
mid_row = N/2 + 1;
stem(j_indices + 0.15, D(mid_row, :), 'filled', 'Color', TEAL, 'MarkerFaceColor', TEAL, ...
    'MarkerSize', 4, 'DisplayName', sprintf('Row $i=%d$ (interior)', mid_row-1));

yline(0, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5);
hold off;

xlabel('Column index $j$');
ylabel('Matrix entry $D_{ij}$');
title('Row Profiles of $D_N$');
legend('Location', 'northeast', 'FontSize', 9);
box off;
grid on;
set(gca, 'GridAlpha', 0.2);
xlim([-0.5, N + 0.5]);

%% Save figure
output_file = fullfile(output_dir, 'cheb_matrix_structure');
exportgraphics(fig, [output_file, '.pdf'], 'ContentType', 'vector');
exportgraphics(fig, [output_file, '.png'], 'Resolution', 300);
fprintf('Figure saved to: %s.pdf\n', output_file);

%% Print matrix properties
fprintf('\nChebyshev Matrix Properties (N = %d):\n', N);
fprintf('%s\n', repmat('-', 1, 50));
fprintf('  Matrix size: (%d, %d)\n', N+1, N+1);
fprintf('  Corner entry D[0,0]: %.6f\n', D(1,1));
fprintf('  Corner entry D[N,N]: %.6f\n', D(N+1,N+1));
fprintf('  Theoretical D[0,0] = (2N^2+1)/6: %.6f\n', (2*N^2+1)/6);
fprintf('  Max entry: %.6f\n', max(D(:)));
fprintf('  Min entry: %.6f\n', min(D(:)));
fprintf('  Negative sum trick error: %.2e\n', max(abs(D * ones(N+1, 1))));
fprintf('%s\n', repmat('-', 1, 50));

%% Red-blue diverging colormap
function cmap = redblue(n)
    if nargin < 1
        n = 256;
    end

    % Create diverging colormap from blue to white to red
    half = floor(n / 2);

    % Blue to white
    r1 = linspace(0, 1, half);
    g1 = linspace(0, 1, half);
    b1 = linspace(1, 1, half);

    % White to red
    r2 = linspace(1, 1, n - half);
    g2 = linspace(1, 0, n - half);
    b2 = linspace(1, 0, n - half);

    cmap = [r1(:), g1(:), b1(:); r2(:), g2(:), b2(:)];
end
