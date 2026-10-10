% trap_subgeometric.m
% Chapter 16: Integration of Periodic Functions
%
% Computational Etude 16.6: Subgeometric Decay on Weideman's f_6.
%
% f_6(x) = exp((cos x - 1)/(cos x + 1)) is C^infty-periodic but has an
% essential singularity at x = pi.  It therefore lies in NO strip of
% analyticity, and the trapezoidal-rule error decays at the
% subgeometric rate exp(-(3/2) N^{2/3}).
%
% The diagnostic visual: replot the error against N^{2/3}; the curve
% should become a straight line.
%
% Generates Figure 16.6: subgeometric.pdf  (two-panel figure)
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
% f_6 with safe handling at x = pi
% -------------------------------------------------------------------------
f6 = @(x) (cos(x) + 1 > 1e-15) .* exp((cos(x) - 1) ./ max(cos(x) + 1, 1e-15));

% Reference value: 2 * e * pi * (1 - erf(1))
I_exact = 2 * exp(1) * pi * (1 - erf(1));
fprintf('Exact integral: 2*e*pi*(1 - erf(1)) = %.16f\n\n', I_exact);

N_values = [4 6 8 10 12 16 20 24 30 40 50 60 80 100 120 160 200];
errors = zeros(size(N_values));
for i = 1:length(N_values)
    N = N_values(i);
    theta = 2*pi*(0:N-1)/N;
    errors(i) = abs((2*pi/N) * sum(f6(theta)) - I_exact);
end
errors = max(errors, 1e-17);

% Theoretical envelope: 8*sqrt(pi/3) * exp(-(3/2)*N^{2/3})
theory = 8 * sqrt(pi / 3) * exp(-1.5 * double(N_values).^(2/3));

% -------------------------------------------------------------------------
% Two-panel plot
% -------------------------------------------------------------------------
fig = figure('Position', [100, 100, 1100, 450]);

subplot(1, 2, 1);
semilogy(N_values, errors, 'o-', 'Color', CORAL, 'MarkerSize', 5, ...
    'LineWidth', 1.2, 'MarkerFaceColor', CORAL, ...
    'DisplayName', 'Trapezoidal error');
hold on;
semilogy(N_values, theory, '--', 'Color', NAVY, 'LineWidth', 1.0, ...
    'DisplayName', '$8\sqrt{\pi/3}\,e^{-3 N^{2/3}/2}$');
hold off;
xlabel('Number of nodes $N$', 'Interpreter', 'latex');
ylabel('Absolute error $|I_N - I|$', 'Interpreter', 'latex');
title('(a) error vs $N$: looks slower than geometric', 'Interpreter', 'latex');
legend('Location', 'northeast', 'Interpreter', 'latex', ...
    'EdgeColor', [0.5 0.5 0.5]);
grid on; set(gca, 'GridAlpha', 0.3);
ylim([1e-17, 1e0]);

subplot(1, 2, 2);
N23 = double(N_values).^(2/3);
semilogy(N23, errors, 'o-', 'Color', CORAL, 'MarkerSize', 5, ...
    'LineWidth', 1.2, 'MarkerFaceColor', CORAL, ...
    'DisplayName', 'Trapezoidal error');
hold on;
semilogy(N23, theory, '--', 'Color', NAVY, 'LineWidth', 1.0, ...
    'DisplayName', '$8\sqrt{\pi/3}\,e^{-3 N^{2/3}/2}$');
hold off;
xlabel('$N^{2/3}$', 'Interpreter', 'latex');
ylabel('Absolute error $|I_N - I|$', 'Interpreter', 'latex');
title('(b) error vs $N^{2/3}$: linear on semilog', 'Interpreter', 'latex');
legend('Location', 'northeast', 'Interpreter', 'latex', ...
    'EdgeColor', [0.5 0.5 0.5]);
grid on; set(gca, 'GridAlpha', 0.3);
ylim([1e-17, 1e0]);

exportgraphics(gcf, fullfile(output_dir, 'subgeometric.pdf'), 'ContentType', 'vector');
exportgraphics(gcf, fullfile(output_dir, 'subgeometric.png'), 'Resolution', 300);
fprintf('Figure saved to %s\n', fullfile(output_dir, 'subgeometric.pdf'));
