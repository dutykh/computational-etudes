%% aliasing_demo.m
%
% Visualization of the aliasing phenomenon (Theorem 2: Poisson summation).
%
% When a function is sampled at N equispaced points, frequencies differing
% by multiples of N become indistinguishable. High frequencies "fold" onto
% low frequencies.
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
% Part of "Computational Etudes: A Spectral Approach"
% https://github.com/dutykh/computational-etudes

clear; close all; clc;

%% Configuration
% Colors
NAVY = [0.078 0.176 0.431];
CORAL = [0.906 0.298 0.235];
TEAL = [0.086 0.627 0.522];
PURPLE = [0.608 0.349 0.714];

% Output path
script_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch06', 'matlab');
output_file = fullfile(output_dir, 'aliasing_visualization.pdf');

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% Create figure with 3 panels
fig = figure('Units', 'inches', 'Position', [1, 1, 12, 4]);

%% Panel 1: Coarse sampling misses oscillations
subplot(1, 3, 1);

% Function with low and high frequency content
f = @(x) sin(x) + 0.5 * sin(10*x);

x_fine = linspace(0, 2*pi, 500);
N_coarse = 12;
x_coarse = linspace(0, 2*pi, N_coarse+1);
x_coarse = x_coarse(1:end-1);

% Aliased interpolant: sin(10x) -> sin(10x - 12x) = sin(-2x)
f_aliased = @(x) sin(x) - 0.5 * sin(2*x);

plot(x_fine, f(x_fine), '-', 'Color', NAVY, 'LineWidth', 1.5, ...
     'DisplayName', 'True function');
hold on;
plot(x_coarse, f(x_coarse), 'o', 'Color', CORAL, 'MarkerSize', 8, ...
     'MarkerFaceColor', CORAL, 'DisplayName', sprintf('Sampled (N=%d)', N_coarse));
plot(x_fine, f_aliased(x_fine), '--', 'Color', TEAL, 'LineWidth', 1.5, ...
     'DisplayName', 'Aliased interpolant');

xlabel('x');
ylabel('f(x)');
title('Coarse Sampling Misses High Frequencies');
legend('Location', 'northeast', 'FontSize', 9);
xlim([0, 2*pi]);
ylim([-1.8, 1.8]);
box on;

%% Panel 2: Frequency folding schematic
subplot(1, 3, 2);

N = 12;
nyquist = N / 2;

% Frequency axis
plot([-18, 18], [0, 0], 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
hold on;

% Nyquist boundaries
xline(-nyquist, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
xline(nyquist, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);

% Shade resolved region
fill([-nyquist, nyquist, nyquist, -nyquist], [-0.8, -0.8, 0.8, 0.8], ...
     TEAL, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
text(0, 0.6, {'Resolved', 'frequencies'}, 'HorizontalAlignment', 'center', ...
     'FontSize', 10, 'Color', TEAL);

% k=10 folds to k=-2
k_original = 10;
k_aliased = k_original - N;  % = -2
annotation('arrow', [0.55, 0.47], [0.65, 0.65], 'Color', CORAL, 'LineWidth', 2);
plot(k_original, 0.3, 'o', 'Color', CORAL, 'MarkerSize', 8, 'MarkerFaceColor', CORAL);
text(k_original, 0.45, sprintf('k=%d', k_original), 'HorizontalAlignment', 'center', ...
     'FontSize', 10, 'Color', CORAL);
text((k_original + k_aliased)/2, 0.5, 'folds to', 'HorizontalAlignment', 'center', ...
     'FontSize', 9, 'Color', CORAL);

% k=-10 folds to k=2
k_original2 = -10;
k_aliased2 = k_original2 + N;  % = 2
annotation('arrow', [0.45, 0.53], [0.35, 0.35], 'Color', PURPLE, 'LineWidth', 2);
plot(k_original2, -0.3, 'o', 'Color', PURPLE, 'MarkerSize', 8, 'MarkerFaceColor', PURPLE);
text(k_original2, -0.45, sprintf('k=%d', k_original2), 'HorizontalAlignment', 'center', ...
     'FontSize', 10, 'Color', PURPLE);

% Labels for boundaries
text(-nyquist, -0.7, {'-N/2', sprintf('=%d', -nyquist)}, 'HorizontalAlignment', 'center', ...
     'FontSize', 10, 'Color', [0.5 0.5 0.5]);
text(nyquist, -0.7, {'N/2', sprintf('=%d', nyquist)}, 'HorizontalAlignment', 'center', ...
     'FontSize', 10, 'Color', [0.5 0.5 0.5]);

xlabel('Wavenumber k');
title(sprintf('Frequency Folding (N=%d)', N));
xlim([-18, 18]);
ylim([-1, 1]);
set(gca, 'YTick', []);
box on;

%% Panel 3: True vs aliased spectrum
subplot(1, 3, 3);

% Function with significant high-frequency content
g = @(x) 1 ./ (1.5 + cos(x));

% Fine spectrum (truth)
N_fine = 128;
x_fine = linspace(0, 2*pi, N_fine+1);
x_fine = x_fine(1:end-1);
g_hat_fine = abs(fft(g(x_fine))) / N_fine;

% Coarse spectrum (aliased)
N_coarse = 16;
x_coarse = linspace(0, 2*pi, N_coarse+1);
x_coarse = x_coarse(1:end-1);
g_hat_coarse = abs(fft(g(x_coarse))) / N_coarse;

k_fine = 0:N_fine/2;
k_coarse = 0:N_coarse/2;

semilogy(k_fine, g_hat_fine(1:N_fine/2+1), 'o-', 'Color', NAVY, ...
         'MarkerSize', 3, 'LineWidth', 1, 'DisplayName', 'True spectrum (N=128)');
hold on;
semilogy(k_coarse, g_hat_coarse(1:N_coarse/2+1), 's-', 'Color', CORAL, ...
         'MarkerSize', 6, 'LineWidth', 1.5, 'DisplayName', sprintf('Aliased spectrum (N=%d)', N_coarse));

xline(N_coarse/2, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'HandleVisibility', 'off');
text(N_coarse/2 + 1, 1e-1, {'Nyquist', sprintf('(N/2=%d)', N_coarse/2)}, ...
     'FontSize', 9, 'Color', [0.5 0.5 0.5]);

xlabel('Wavenumber k');
ylabel('|\^g_k|');
title('True vs. Aliased Spectrum');
legend('Location', 'northeast', 'FontSize', 9);
xlim([-1, 50]);
ylim([1e-10, 1e1]);
box on;

%% Save figure
exportgraphics(fig, output_file, 'ContentType', 'vector');
exportgraphics(fig, strrep(output_file, '.pdf', '.png'), 'Resolution', 300);

fprintf('Figure saved to: %s\n', output_file);

close(fig);
