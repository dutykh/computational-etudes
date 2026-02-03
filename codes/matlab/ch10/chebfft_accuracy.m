%% chebfft_accuracy.m
%
% Figure 10.2: Demonstration of chebfft spectral accuracy
%
% Tests the Chebyshev FFT differentiation on f(x) = exp(x) * cos(4x)
% and shows spectral convergence to machine precision.
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

% Output path
script_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch10', 'matlab');
output_file = fullfile(output_dir, 'chebfft_accuracy.pdf');

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% Test function
f = @(x) exp(x) .* cos(4*x);
f_prime = @(x) exp(x) .* (cos(4*x) - 4*sin(4*x));

%% Panel 1: Function with Chebyshev nodes (N = 16)
N_demo = 16;
x_demo = cos(pi * (0:N_demo)' / N_demo);
y_demo = f(x_demo);

%% Panel 2: Derivative comparison
x_fine = linspace(-1, 1, 200)';
y_exact = f_prime(x_fine);
y_spectral = chebfft(y_demo);  % Derivative at Chebyshev points

%% Panel 3: Convergence study
N_values = 4:2:40;
errors = zeros(size(N_values));

for i = 1:length(N_values)
    N = N_values(i);
    x = cos(pi * (0:N)' / N);
    v = f(x);
    v_prime_exact = f_prime(x);
    v_prime_fft = chebfft(v);
    errors(i) = max(abs(v_prime_fft - v_prime_exact));
end

%% Create figure
fig = figure('Units', 'inches', 'Position', [1, 1, 14, 4.5]);

% Panel 1: Function and nodes
subplot(1, 3, 1);
plot(x_fine, f(x_fine), '-', 'Color', NAVY, 'LineWidth', 1.5);
hold on;
plot(x_demo, y_demo, 'o', 'Color', CORAL, 'MarkerSize', 8, 'MarkerFaceColor', CORAL);
xlabel('x', 'FontSize', 11);
ylabel('f(x)', 'FontSize', 11);
title('f(x) = e^x cos(4x) with Chebyshev Nodes', 'FontSize', 12);
legend({'f(x)', sprintf('N = %d nodes', N_demo)}, 'Location', 'northwest');
xlim([-1, 1]);
grid on;
box on;

% Panel 2: Derivative comparison
subplot(1, 3, 2);
plot(x_fine, y_exact, '-', 'Color', NAVY, 'LineWidth', 1.5);
hold on;
plot(x_demo, y_spectral, 'o', 'Color', CORAL, 'MarkerSize', 8, 'MarkerFaceColor', CORAL);
xlabel('x', 'FontSize', 11);
ylabel("f'(x)", 'FontSize', 11);
title(sprintf("Exact vs Spectral Derivative (N = %d)", N_demo), 'FontSize', 12);
legend({'Exact', 'chebfft'}, 'Location', 'northwest');
xlim([-1, 1]);
grid on;
box on;

% Panel 3: Convergence
subplot(1, 3, 3);
semilogy(N_values, errors, 'o-', 'Color', NAVY, 'LineWidth', 1.5, 'MarkerSize', 6, ...
         'MarkerFaceColor', NAVY);
hold on;

% Machine precision line
yline(eps, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
text(38, eps * 3, 'Machine precision', 'FontSize', 9, 'Color', [0.5 0.5 0.5], ...
     'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

xlabel('N (number of intervals)', 'FontSize', 11);
ylabel("Max error ||f' - Df||_\infty", 'FontSize', 11);
title('Spectral Convergence of chebfft', 'FontSize', 12);
xlim([0, 42]);
ylim([1e-16, 1e0]);
grid on;
box on;

%% Save figure
exportgraphics(fig, output_file, 'ContentType', 'vector');
exportgraphics(fig, strrep(output_file, '.pdf', '.png'), 'Resolution', 300);

fprintf('Figure saved to: %s\n', output_file);

%% Print convergence table
fprintf('\nConvergence Study:\n');
fprintf('---------------------------\n');
fprintf('%6s  %14s\n', 'N', 'Max Error');
fprintf('---------------------------\n');
for i = 1:length(N_values)
    fprintf('%6d  %14.6e\n', N_values(i), errors(i));
end
fprintf('---------------------------\n');

close(fig);
