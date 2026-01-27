%% fourier_decay.m
%
% Demonstration of the decay hierarchy for Fourier coefficients.
% Illustrates Theorem 1: smoothness determines spectral decay.
%
% Three test functions on [0, 2pi]:
%   1. |sin(x)|^3         - Finite regularity, O(k^{-4}) decay
%   2. 1/(1+sin^2(x/2))   - Analytic in strip, geometric decay
%   3. exp(sin(x))        - Entire function, super-geometric decay
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
% Part of "Computational Etudes: A Spectral Approach"
% https://github.com/dutykh/computational-etudes

clear; close all; clc;

%% Configuration
N = 128;  % Number of sample points

% Colors
NAVY = [0.078 0.176 0.431];
CORAL = [0.906 0.298 0.235];
TEAL = [0.086 0.627 0.522];

% Output path
script_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch06', 'matlab');
output_file = fullfile(output_dir, 'decay_hierarchy.pdf');

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% Test functions
f1 = @(x) abs(sin(x)).^3;               % Finite regularity
f2 = @(x) 1 ./ (1 + sin(x/2).^2);       % Analytic in strip
f3 = @(x) exp(sin(x));                   % Entire

%% Compute Fourier coefficients
x = linspace(0, 2*pi, N+1);
x = x(1:end-1);  % Remove endpoint for periodicity

f1_vals = f1(x);
f2_vals = f2(x);
f3_vals = f3(x);

f1_hat = abs(fft(f1_vals)) / N;
f2_hat = abs(fft(f2_vals)) / N;
f3_hat = abs(fft(f3_vals)) / N;

% Keep only positive frequencies
k = 0:N/2;
f1_hat = f1_hat(1:N/2+1);
f2_hat = f2_hat(1:N/2+1);
f3_hat = f3_hat(1:N/2+1);

%% Create figure
fig = figure('Units', 'inches', 'Position', [1, 1, 8, 5.5]);

% Plot coefficients (skip k=0)
semilogy(k(2:end), f1_hat(2:end), 'o-', 'Color', TEAL, 'LineWidth', 1.5, ...
         'MarkerSize', 4, 'DisplayName', '|sin(x)|^3 (finite regularity)');
hold on;
semilogy(k(2:end), f2_hat(2:end), 's-', 'Color', CORAL, 'LineWidth', 1.5, ...
         'MarkerSize', 4, 'DisplayName', '1/(1+sin^2(x/2)) (analytic in strip)');
semilogy(k(2:end), f3_hat(2:end), 'd-', 'Color', NAVY, 'LineWidth', 1.5, ...
         'MarkerSize', 4, 'DisplayName', 'exp(sin(x)) (entire)');

% Reference lines
k_ref = 5:60;

% O(k^{-4}) reference
C1 = f1_hat(11) * 10^4;
TEAL_LIGHT = TEAL * 0.5 + [1 1 1] * 0.5;
semilogy(k_ref, C1 ./ k_ref.^4, '--', 'Color', TEAL_LIGHT, 'LineWidth', 1, ...
         'HandleVisibility', 'off');
text(55, C1 / 55^4 * 2, 'O(k^{-4})', 'FontSize', 10, 'Color', TEAL);

% Geometric decay reference (fit the rate)
k_fit = 10:40;
log_f2 = log(f2_hat(k_fit+1) + 1e-20);
coeffs = polyfit(k_fit, log_f2, 1);
decay_rate = -coeffs(1);
CORAL_LIGHT = CORAL * 0.5 + [1 1 1] * 0.5;
semilogy(k_ref, exp(coeffs(2)) * exp(-decay_rate * k_ref), '--', ...
         'Color', CORAL_LIGHT, 'LineWidth', 1, 'HandleVisibility', 'off');
text(55, exp(coeffs(2)) * exp(-decay_rate * 55) * 3, ...
     sprintf('O(e^{-%.2f k})', decay_rate), 'FontSize', 10, 'Color', CORAL);

% Machine precision line
yline(1e-16, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'HandleVisibility', 'off');
text(62, 2e-16, 'Machine precision', 'FontSize', 9, 'Color', [0.5 0.5 0.5], ...
     'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

% Labels
xlabel('Wavenumber k', 'FontSize', 11);
ylabel('Fourier coefficient magnitude |\^f_k|', 'FontSize', 11);
title('The Decay Hierarchy: Smoothness Determines Spectral Decay', 'FontSize', 12);

xlim([0, 65]);
ylim([1e-17, 1e1]);
legend('Location', 'northeast', 'FontSize', 10);
grid on;
box on;

%% Save figure
exportgraphics(fig, output_file, 'ContentType', 'vector');
exportgraphics(fig, strrep(output_file, '.pdf', '.png'), 'Resolution', 300);

fprintf('Figure saved to: %s\n', output_file);

%% Print coefficient values
fprintf('\nFourier Coefficient Magnitudes:\n');
fprintf('============================================================\n');
fprintf('%4s  %14s  %14s  %14s\n', 'k', '|sin(x)|^3', '1/(1+sin^2)', 'exp(sin)');
fprintf('------------------------------------------------------------\n');
for kk = [1, 2, 5, 10, 20, 30, 40, 50]
    if kk <= length(f1_hat)
        fprintf('%4d  %14.6e  %14.6e  %14.6e\n', kk, f1_hat(kk+1), f2_hat(kk+1), f3_hat(kk+1));
    end
end
fprintf('============================================================\n');

close(fig);
