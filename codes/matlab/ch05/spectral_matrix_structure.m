%% spectral_matrix_structure.m
%
% Visualizes the structure of the periodic spectral differentiation matrix.
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
% Part of "Computational Études: A Spectral Approach"
% https://github.com/dutykh/computational-etudes

clear; close all; clc;

%% Configuration
N = 16;  % Number of grid points

% Output path
script_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch05', 'matlab');
output_file = fullfile(output_dir, 'spectral_matrix_structure.pdf');

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% Construct matrix
[D, x] = spectral_matrix_periodic(N);

%% Create figure
fig = figure('Units', 'inches', 'Position', [1, 1, 10, 4.5]);

% Left panel: Heatmap
subplot(1, 2, 1);
vmax = max(abs(D(:)));
imagesc(D, [-vmax, vmax]);
colormap(redblue(256));
colorbar;
title(sprintf('Spectral Matrix D (N = %d)', N), 'FontSize', 12);
xlabel('Column index j', 'FontSize', 11);
ylabel('Row index i', 'FontSize', 11);
axis square;
set(gca, 'XTick', [1, N/2, N], 'XTickLabel', {'0', num2str(N/2-1), num2str(N-1)});
set(gca, 'YTick', [1, N/2, N], 'YTickLabel', {'0', num2str(N/2-1), num2str(N-1)});

% Add annotation
text(0.05, 0.95, 'D^T = -D', 'Units', 'normalized', 'FontSize', 11, ...
     'FontWeight', 'bold', 'Color', [0.078 0.176 0.431], ...
     'BackgroundColor', 'white');

% Right panel: First row profile
subplot(1, 2, 2);
row = D(1, :);
j_indices = 0:N-1;

% Shift indices for centered display
j_shifted = j_indices;
j_shifted(j_indices > N/2) = j_indices(j_indices > N/2) - N;
[j_sorted, sort_idx] = sort(j_shifted);
row_sorted = row(sort_idx);

stem(j_sorted, row_sorted, 'filled', 'Color', [0.078 0.176 0.431], 'MarkerSize', 5);
hold on;

% Theoretical curve
j_theory = linspace(-N/2 + 0.1, N/2 - 0.1, 500);
j_theory = j_theory(abs(j_theory) > 0.1);
D_theory = 0.5 * ((-1).^round(j_theory)) ./ tan(j_theory * pi / N);
plot(j_theory, D_theory, '--', 'Color', [0.91 0.30 0.24], 'LineWidth', 1.5);

xlabel('Column offset k = j - i', 'FontSize', 11);
ylabel('Matrix entry D_{0,k}', 'FontSize', 11);
title('First Row Profile (Toeplitz Structure)', 'FontSize', 12);
xlim([-N/2 - 0.5, N/2 + 0.5]);
grid on;
legend({'Entries', 'cot formula'}, 'Location', 'northeast');

%% Save figure
exportgraphics(fig, output_file, 'ContentType', 'vector');
exportgraphics(fig, strrep(output_file, '.pdf', '.png'), 'Resolution', 300);

fprintf('Figure saved to: %s\n', output_file);

%% Print properties
fprintf('\nSpectral Matrix Properties (N = %d):\n', N);
fprintf('--------------------------------------------------\n');
fprintf('  Skew-symmetry error ||D + D^T||_max: %.2e\n', max(abs(D + D'), [], 'all'));
fprintf('  Diagonal sum: %.2e (should be 0)\n', sum(diag(D)));
fprintf('  Max entry: %.4f\n', max(D(:)));
fprintf('  Min entry: %.4f\n', min(D(:)));
fprintf('--------------------------------------------------\n');

close(fig);

%% Helper function: red-blue colormap
function cmap = redblue(m)
    if nargin < 1, m = 256; end
    % Create a red-white-blue diverging colormap
    r = [linspace(0, 1, m/2), ones(1, m/2)];
    g = [linspace(0, 1, m/2), linspace(1, 0, m/2)];
    b = [ones(1, m/2), linspace(1, 0, m/2)];
    cmap = [r', g', b'];
end
