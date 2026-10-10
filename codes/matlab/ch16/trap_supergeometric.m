% trap_supergeometric.m
% Chapter 16: Integration of Periodic Functions
%
% Computational Etude 16.5: Supergeometric Decay on e^cos(x).
%
% f_5(x) = exp(cos x) is entire (analytic in the whole complex plane),
% so the strip-analyticity bound applies for any width a > 0.
% Optimising a as a function of N gives the supergeometric rate
% (e/(2N))^N.  Each new sample point gives more than one extra digit.
%
% Generates Figure 16.5: supergeometric.pdf
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
% Setup
% -------------------------------------------------------------------------
f5 = @(x) exp(cos(x));
I_exact = 2 * pi * besseli(0, 1);
fprintf('Exact integral: 2*pi*I_0(1) = %.16f\n\n', I_exact);

N_values = 1:16;
errors = zeros(size(N_values));
for i = 1:length(N_values)
    N = N_values(i);
    theta = 2*pi*(0:N-1)/N;
    errors(i) = abs((2*pi/N) * sum(f5(theta)) - I_exact);
end
errors = max(errors, 1e-17);

% Asymptotic from Weideman: 2*sqrt(2*pi/N) * (e/(2N))^N
ee = exp(1);
theory = 2 * sqrt(2 * pi ./ N_values) .* (ee ./ (2 * N_values)).^N_values;
theory = max(theory, 1e-17);

% Print table
fprintf('%4s  %22s  %14s  %14s\n', 'N', 'I_N', '|I_N - I|', 'theory');
for i = 1:length(N_values)
    N = N_values(i);
    theta = 2*pi*(0:N-1)/N;
    I_N = (2*pi/N) * sum(f5(theta));
    fprintf('%4d  %22.16f  %14.4e  %14.4e\n', N, I_N, errors(i), theory(i));
end

% -------------------------------------------------------------------------
% Plot
% -------------------------------------------------------------------------
fig = figure('Position', [100, 100, 700, 500]);
semilogy(N_values, errors, 'o-', 'Color', CORAL, 'MarkerSize', 6, ...
    'LineWidth', 1.2, 'MarkerFaceColor', CORAL, ...
    'DisplayName', 'Trapezoidal error');
hold on;
semilogy(N_values, theory, '--', 'Color', NAVY, 'LineWidth', 1.0, ...
    'DisplayName', '$2\sqrt{2\pi/N}\,(e/2N)^N$');
hold off;
xlabel('Number of nodes $N$', 'Interpreter', 'latex');
ylabel('Absolute error $|I_N - I|$', 'Interpreter', 'latex');
title('Supergeometric decay for $f(x) = e^{\cos x}$', 'Interpreter', 'latex');
legend('Location', 'northeast', 'Interpreter', 'latex', ...
    'EdgeColor', [0.5 0.5 0.5]);
grid on; set(gca, 'GridAlpha', 0.3);
xlim([0, 17]); ylim([1e-18, 1e2]);

exportgraphics(gcf, fullfile(output_dir, 'supergeometric.pdf'), 'ContentType', 'vector');
exportgraphics(gcf, fullfile(output_dir, 'supergeometric.png'), 'Resolution', 300);
fprintf('\nFigure saved to %s\n', fullfile(output_dir, 'supergeometric.pdf'));
