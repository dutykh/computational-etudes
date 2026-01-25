%% fd_matrix_bandwidth.m
%
% Visualizes the sparsity patterns of finite difference and spectral
% differentiation matrices, showing how the bandwidth expands as the
% order of accuracy increases.
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
% Part of "Computational Études: A Spectral Approach"
% https://github.com/dutykh/computational-etudes

clear; close all; clc;

%% Configuration
N = 20;  % Number of grid points

% Output path
script_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch05', 'matlab');
output_file = fullfile(output_dir, 'fd_matrix_bandwidth.pdf');

% Create output directory if needed
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% Construct matrices
[D_fd2, ~] = fd_matrix_periodic(N, 2);   % 2nd order: tridiagonal
[D_fd4, ~] = fd_matrix_periodic(N, 4);   % 4th order: pentadiagonal
[D_spec, ~] = spectral_matrix_periodic(N);  % Spectral: dense

%% Create figure
fig = figure('Units', 'inches', 'Position', [1, 1, 11, 3.5]);

matrices = {D_fd2, D_fd4, D_spec};
titles = {'2nd Order FD', '4th Order FD', 'Spectral'};
bandwidths = {'Bandwidth = 3', 'Bandwidth = 5', sprintf('Dense (N \\times N)')};

for k = 1:3
    subplot(1, 3, k);

    % Create binary sparsity pattern
    threshold = 1e-14;
    sparsity = abs(matrices{k}) > threshold;

    % Plot using imagesc
    imagesc(sparsity);
    colormap([1 1 1; 0.078 0.176 0.431]);  % White and Navy

    title(titles{k}, 'FontSize', 12, 'FontWeight', 'bold');
    % Combine column label and bandwidth info
    xlabel({['Column index ' char(106)], bandwidths{k}}, 'FontSize', 10);
    if k == 1
        ylabel('Row index i', 'FontSize', 10);
    end

    axis square;
    set(gca, 'XTick', [1, N/2, N], 'XTickLabel', {'0', num2str(N/2-1), num2str(N-1)});
    set(gca, 'YTick', [1, N/2, N], 'YTickLabel', {'0', num2str(N/2-1), num2str(N-1)});
end

sgtitle('Finite Difference Matrix Sparsity: From Local to Global', 'FontSize', 13);

%% Save figure
exportgraphics(fig, output_file, 'ContentType', 'vector');
exportgraphics(fig, strrep(output_file, '.pdf', '.png'), 'Resolution', 300);

fprintf('Figure saved to: %s\n', output_file);

%% Print statistics
fprintf('\nMatrix Statistics (N = %d):\n', N);
fprintf('--------------------------------------------------\n');
for k = 1:3
    nnz_count = sum(abs(matrices{k}(:)) > 1e-14);
    density = nnz_count / numel(matrices{k}) * 100;
    fprintf('%15s: %4d nonzeros (%5.1f%% dense)\n', titles{k}, nnz_count, density);
end
fprintf('--------------------------------------------------\n');

close(fig);
