%% convergence_rates.m
%
% Demonstration of spectral differentiation convergence rates for functions
% of varying smoothness, illustrating Theorems 3 and 4.
%
% Three test functions on [0, 2pi]:
%   1. |sin(x)|^3         - Algebraic convergence O(N^{-3})
%   2. 1/(1+sin^2(x/2))   - Geometric convergence O(c^{-N})
%   3. exp(sin(x))        - Super-geometric convergence
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
% Part of "Computational Etudes: A Spectral Approach"
% https://github.com/dutykh/computational-etudes

clear; close all; clc;

%% Configuration
N_values = [6, 8, 10, 12, 14, 16, 20, 24, 28, 32, 40, 48, 56, 64];

% Colors
NAVY = [0.078 0.176 0.431];
CORAL = [0.906 0.298 0.235];
TEAL = [0.086 0.627 0.522];

% Output path
script_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch06_smoothness', 'matlab');
output_file = fullfile(output_dir, 'convergence_rates.pdf');

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% Test functions and derivatives
f1 = @(x) abs(sin(x)).^3;
f1_deriv = @(x) 3 * sin(x) .* abs(sin(x)) .* cos(x);

f2 = @(x) 1 ./ (1 + sin(x/2).^2);
f2_deriv = @(x) -sin(x/2) .* cos(x/2) ./ (1 + sin(x/2).^2).^2;

f3 = @(x) exp(sin(x));
f3_deriv = @(x) cos(x) .* exp(sin(x));

%% Compute errors
errors1 = zeros(1, length(N_values));
errors2 = zeros(1, length(N_values));
errors3 = zeros(1, length(N_values));

for i = 1:length(N_values)
    N = N_values(i);
    [D, x] = spectral_diff_periodic(N);

    % Function 1: |sin(x)|^3
    u1 = f1(x);
    du1_exact = f1_deriv(x);
    du1_spectral = D * u1;
    errors1(i) = max(abs(du1_spectral - du1_exact));

    % Function 2: 1/(1 + sin^2(x/2))
    u2 = f2(x);
    du2_exact = f2_deriv(x);
    du2_spectral = D * u2;
    errors2(i) = max(abs(du2_spectral - du2_exact));

    % Function 3: exp(sin(x))
    u3 = f3(x);
    du3_exact = f3_deriv(x);
    du3_spectral = D * u3;
    errors3(i) = max(abs(du3_spectral - du3_exact));
end

%% Create figure
fig = figure('Units', 'inches', 'Position', [1, 1, 8, 5.5]);

semilogy(N_values, errors1, 'o-', 'Color', TEAL, 'LineWidth', 1.5, ...
         'MarkerSize', 6, 'DisplayName', '|sin(x)|^3 (finite regularity)');
hold on;
semilogy(N_values, errors2, 's-', 'Color', CORAL, 'LineWidth', 1.5, ...
         'MarkerSize', 6, 'DisplayName', '1/(1+sin^2(x/2)) (analytic in strip)');
semilogy(N_values, errors3, 'd-', 'Color', NAVY, 'LineWidth', 1.5, ...
         'MarkerSize', 6, 'DisplayName', 'exp(sin(x)) (entire)');

% Reference lines
N_ref = linspace(8, 64, 100);

% O(N^{-3}) reference
C1 = errors1(4) * N_values(4)^3;
TEAL_LIGHT = TEAL * 0.5 + [1 1 1] * 0.5;
semilogy(N_ref, C1 ./ N_ref.^3, '--', 'Color', TEAL_LIGHT, 'LineWidth', 1, ...
         'HandleVisibility', 'off');
text(66, C1 / 66^3 * 1.5, 'O(N^{-3})', 'FontSize', 10, 'Color', TEAL, ...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');

% Geometric decay reference (fit the rate)
idx_fit = (N_values >= 10) & (N_values <= 40) & (errors2 > 1e-15);
if sum(idx_fit) >= 2
    coeffs = polyfit(N_values(idx_fit), log(errors2(idx_fit)), 1);
    decay_rate = -coeffs(1);
    C2 = exp(coeffs(2));
    CORAL_LIGHT = CORAL * 0.5 + [1 1 1] * 0.5;
    semilogy(N_ref, C2 * exp(-decay_rate * N_ref), '--', 'Color', CORAL_LIGHT, ...
             'LineWidth', 1, 'HandleVisibility', 'off');
    text(50, C2 * exp(-decay_rate * 50) * 3, sprintf('O(e^{-%.2f N})', decay_rate), ...
         'FontSize', 10, 'Color', CORAL);
end

% Machine precision line
yline(2.2e-16, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'HandleVisibility', 'off');
text(66, 4e-16, 'Machine precision', 'FontSize', 9, 'Color', [0.5 0.5 0.5], ...
     'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

% Super-geometric annotation
text(35, 1e-8, {'Super-geometric', 'convergence'}, 'FontSize', 10, 'Color', NAVY, ...
     'HorizontalAlignment', 'center');
annotation('arrow', [0.52, 0.42], [0.35, 0.25], 'Color', NAVY, 'LineWidth', 1.5);

% Labels
xlabel('Number of grid points N', 'FontSize', 11);
ylabel('Differentiation error ||f'' - Df||_\infty', 'FontSize', 11);
title('Spectral Differentiation Convergence: Smoothness Matters', 'FontSize', 12);

xlim([0, 70]);
ylim([1e-16, 1e1]);
legend('Location', 'northeast', 'FontSize', 10);
grid on;
box on;

%% Save figure
exportgraphics(fig, output_file, 'ContentType', 'vector');
exportgraphics(fig, strrep(output_file, '.pdf', '.png'), 'Resolution', 300);

fprintf('Figure saved to: %s\n', output_file);

%% Print results table
fprintf('\nConvergence Study Results:\n');
fprintf('=================================================================\n');
fprintf('%4s  %14s  %14s  %14s\n', 'N', '|sin(x)|^3', '1/(1+sin^2)', 'exp(sin)');
fprintf('-----------------------------------------------------------------\n');
for i = 1:length(N_values)
    fprintf('%4d  %14.6e  %14.6e  %14.6e\n', N_values(i), errors1(i), errors2(i), errors3(i));
end
fprintf('=================================================================\n');

close(fig);


%% Helper function: periodic spectral differentiation matrix
function [D, x] = spectral_diff_periodic(N)
    h = 2 * pi / N;
    x = h * (0:N-1)';
    D = zeros(N, N);

    for i = 1:N
        for j = 1:N
            if i ~= j
                D(i, j) = 0.5 * ((-1)^(i-j)) / tan((i-j) * pi / N);
            end
        end
    end
end
