% trap_real_line.m
% Chapter 16: Integration of Periodic Functions
%
% Computational Etude 16.7: Trapezoidal Rule on the Real Line.
%
% Compute (1/sqrt(pi)) * int_{-infty}^{infty} exp(-x^2) dx = 1 by the
% truncated trapezoidal rule with step h = 2*pi/N.  Theorem 5.1 of
% Trefethen-Weideman (2014) predicts O(exp(-pi^2/h)) decay; the actual
% error reaches machine precision around N = 12.
%
% Generates Figure 16.7: real_line_gaussian.pdf
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
% The Goodwin example: w(x) = exp(-x^2)/sqrt(pi), I = 1
% -------------------------------------------------------------------------
w = @(x) exp(-x.^2) / sqrt(pi);
I_exact = 1;

N_values = 1:12;
h_values = 2 * pi ./ N_values;
errors = zeros(size(N_values));
for i = 1:length(N_values)
    h = h_values(i);
    n_max = max(ceil(28 / h), 30);
    k = -n_max:n_max;
    I_h = h * sum(w(k * h));
    errors(i) = abs(I_h - I_exact);
end
errors = max(errors, 1e-17);

theory = exp(-pi^2 ./ h_values);

% Print table
fprintf('%4s  %10s  %22s  %14s\n', 'N', 'h', 'I_h', '|I_h - I|');
for i = 1:length(N_values)
    h = h_values(i);
    n_max = max(ceil(28 / h), 30);
    k = -n_max:n_max;
    I_h = h * sum(w(k * h));
    fprintf('%4d  %10.6f  %22.16f  %14.4e\n', N_values(i), h, I_h, errors(i));
end

% -------------------------------------------------------------------------
% Plot
% -------------------------------------------------------------------------
fig = figure('Position', [100, 100, 700, 500]);
semilogy(N_values, errors, 'o-', 'Color', CORAL, 'MarkerSize', 6, ...
    'LineWidth', 1.2, 'MarkerFaceColor', CORAL, ...
    'DisplayName', '$|I_h - I|$, $h = 2\pi/N$');
hold on;
semilogy(N_values, theory, '--', 'Color', NAVY, 'LineWidth', 1.0, ...
    'DisplayName', '$e^{-\pi^2/h}$');
hold off;
xlabel('$N$ (with step $h = 2\pi/N$)', 'Interpreter', 'latex');
ylabel('Absolute error', 'Interpreter', 'latex');
title('Real-line trapezoidal rule on $e^{-x^2}/\sqrt{\pi}$', ...
    'Interpreter', 'latex');
legend('Location', 'northeast', 'Interpreter', 'latex', ...
    'EdgeColor', [0.5 0.5 0.5]);
grid on; set(gca, 'GridAlpha', 0.3);
xlim([0, 13]); ylim([1e-18, 1e1]);

exportgraphics(gcf, fullfile(output_dir, 'real_line_gaussian.pdf'), 'ContentType', 'vector');
exportgraphics(gcf, fullfile(output_dir, 'real_line_gaussian.png'), 'Resolution', 300);
fprintf('\nFigure saved to %s\n', fullfile(output_dir, 'real_line_gaussian.pdf'));
