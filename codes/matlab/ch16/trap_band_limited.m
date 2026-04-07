% trap_band_limited.m
% Chapter 16: Integration of Periodic Functions
%
% Computational Etude 16.2: Band-Limited Exactness and One-Mode Aliasing.
%
% (a) The N-point periodic trapezoidal rule is exact for any
%     trigonometric polynomial of degree m < N.
% (b) For a single mode cos(k*theta) with N fixed, the rule still
%     gives the correct value of zero unless k is an integer multiple
%     of N -- even when k > N/2 (the aliased Nyquist regime).
%
% Generates Figure 16.2: band_limited.pdf  (two-panel figure)
%
% Author: Dr. Denys Dutykh
%         Mathematics Department
%         Khalifa University of Science and Technology
%         Abu Dhabi, UAE
%
% Part of the book "Computational Etudes: A Spectral Approach"
% https://github.com/dutykh/computational-etudes

clear; close all; clc;

set(0, 'DefaultAxesFontSize', 11);
set(0, 'DefaultTextFontSize', 11);
set(0, 'DefaultLineLineWidth', 1.5);
set(0, 'DefaultAxesLineWidth', 0.8);

NAVY  = [20, 45, 110] / 255;
CORAL = [231, 76, 60] / 255;

output_dir = fullfile(fileparts(mfilename('fullpath')), ...
    '..', '..', '..', 'textbook', 'figures', 'ch16', 'matlab');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

% -------------------------------------------------------------------------
% Part (a): random trigonometric polynomial of degree m = 10
% -------------------------------------------------------------------------
rng(42);
m = 10;
a = randn(m + 1, 1);
b = randn(m + 1, 1);

% trig_poly(theta) = a(1) + sum_{k=1}^m [a(k+1)*cos(k*theta) + b(k+1)*sin(k*theta)]
trig_poly = @(t) a(1) + sum(a(2:end) .* cos((1:m)' * t(:)') + ...
                            b(2:end) .* sin((1:m)' * t(:)'), 1);

I_exact = 2 * pi * a(1);
N_values_a = 4:30;
errors_a = zeros(size(N_values_a));
for i = 1:length(N_values_a)
    N = N_values_a(i);
    theta = 2*pi*(0:N-1)/N;
    errors_a(i) = abs((2*pi/N) * sum(trig_poly(theta)) - I_exact);
end
errors_a = max(errors_a, 1e-17);

fprintf('Part (a): random trig polynomial of degree m = %d\n', m);
fprintf('Exact integral: 2*pi*a(1) = %.10f\n', I_exact);
fprintf('Error at N = m = %d: %.4e\n', m, errors_a(m - 3));
fprintf('Error at N = m+1 = %d: %.4e\n\n', m + 1, errors_a(m - 2));

% -------------------------------------------------------------------------
% Part (b): single mode cos(k*theta), N = 16 fixed, sweep k
% -------------------------------------------------------------------------
N_fixed = 16;
k_values = 0:32;
errors_b = zeros(size(k_values));
for i = 1:length(k_values)
    k = k_values(i);
    theta = 2*pi*(0:N_fixed-1)/N_fixed;
    if k == 0
        I_exact_k = 2*pi;
    else
        I_exact_k = 0;
    end
    errors_b(i) = abs((2*pi/N_fixed) * sum(cos(k * theta)) - I_exact_k);
end
errors_b = max(errors_b, 1e-17);

% -------------------------------------------------------------------------
% Two-panel plot
% -------------------------------------------------------------------------
fig = figure('Position', [100, 100, 1100, 450]);

subplot(1, 2, 1);
semilogy(N_values_a, errors_a, 'o-', 'Color', CORAL, 'MarkerSize', 5, ...
    'LineWidth', 1.1, 'MarkerFaceColor', CORAL);
hold on;
xline(m, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.9);
hold off;
xlabel('Number of nodes $N$', 'Interpreter', 'latex');
ylabel('Absolute error $|I_N - I|$', 'Interpreter', 'latex');
title('(a) trigonometric polynomial of degree $m = 10$', 'Interpreter', 'latex');
grid on; set(gca, 'GridAlpha', 0.3);
xlim([3, 31]); ylim([1e-18, 1e2]);

subplot(1, 2, 2);
% Bar plot with separate colours for the failure modes
mask_fail = errors_b > 1e-10;
bar(k_values(~mask_fail), errors_b(~mask_fail), 'FaceColor', NAVY, 'EdgeColor', 'none');
hold on;
bar(k_values(mask_fail), errors_b(mask_fail), 'FaceColor', CORAL, 'EdgeColor', 'none');
xline(N_fixed, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.9);
xline(2*N_fixed, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.9);
hold off;
xlabel('Frequency $k$', 'Interpreter', 'latex');
ylabel('$|I_N(\cos k\theta) - I(\cos k\theta)|$', 'Interpreter', 'latex');
title('(b) single modes, fixed $N = 16$', 'Interpreter', 'latex');
set(gca, 'YScale', 'log');
ylim([1e-17, 1e1]); xlim([-0.5, 32.5]);
grid on; set(gca, 'GridAlpha', 0.3);

exportgraphics(gcf, fullfile(output_dir, 'band_limited.pdf'), 'ContentType', 'vector');
exportgraphics(gcf, fullfile(output_dir, 'band_limited.png'), 'Resolution', 300);
fprintf('Figure saved to %s\n', fullfile(output_dir, 'band_limited.pdf'));
