% trap_poisson_kernel.m
% Chapter 16: Integration of Periodic Functions
%
% Computational Etude 16.4: Geometric Convergence on the Poisson Kernel.
%
% f_4(x) = 1 / (a - cos x), a > 1.  Predicted geometric rate r^N with
% r = a - sqrt(a^2 - 1), and an explicit closed-form trapezoidal error
%
%   I(f_4) - T_N(f_4) = -8*pi * r/(1 - r^2) * r^N / (1 - r^N)
%
% derived in Weideman (2002), eq. (21).
%
% Generates Figure 16.4: poisson_kernel.pdf  (two-panel figure)
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
a = 2;
f4 = @(x) 1 ./ (a - cos(x));
I_exact = 2 * pi / sqrt(a^2 - 1);
r = a - sqrt(a^2 - 1);
alpha = log(a + sqrt(a^2 - 1));

fprintf('Poisson kernel f_4 = 1 / (a - cos x), a = %g\n', a);
fprintf('Exact integral: 2*pi/sqrt(a^2 - 1) = %.10f\n', I_exact);
fprintf('Geometric rate r = %.10f\n\n', r);

N_values = 2:2:50;
errors = zeros(size(N_values));
for i = 1:length(N_values)
    N = N_values(i);
    theta = 2*pi*(0:N-1)/N;
    errors(i) = abs((2*pi/N) * sum(f4(theta)) - I_exact);
end
errors = max(errors, 1e-17);

theory = 8 * pi * r / (1 - r^2) * r.^N_values ./ (1 - r.^N_values);
theory = max(theory, 1e-17);

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
    'DisplayName', 'closed-form (Weideman 2002)');
hold off;
xlabel('Number of nodes $N$', 'Interpreter', 'latex');
ylabel('Absolute error $|I_N - I|$', 'Interpreter', 'latex');
title('Geometric decay at rate $r^N$ for $a = 2$', 'Interpreter', 'latex');
legend('Location', 'northeast', 'Interpreter', 'latex', ...
    'EdgeColor', [0.5 0.5 0.5]);
grid on; set(gca, 'GridAlpha', 0.3);
ylim([1e-18, 1e1]);

subplot(1, 2, 2);
% Strip of analyticity and pole picture
yline(0, 'k', 'LineWidth', 1.0);
hold on;
% Strip
fill([-0.5, 2*pi+0.5, 2*pi+0.5, -0.5], [-alpha, -alpha, alpha, alpha], ...
    NAVY, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
yline(alpha, '--', 'Color', NAVY, 'LineWidth', 0.8);
yline(-alpha, '--', 'Color', NAVY, 'LineWidth', 0.8);
% Trapezoidal nodes for N = 12
N_demo = 12;
theta_nodes = 2*pi*(0:N_demo-1)/N_demo;
plot(theta_nodes, zeros(size(theta_nodes)), 'o', 'Color', CORAL, ...
    'MarkerSize', 6, 'MarkerFaceColor', CORAL);
% Poles
plot(0, alpha, 'x', 'Color', NAVY, 'MarkerSize', 10, 'LineWidth', 2);
plot(0, -alpha, 'x', 'Color', NAVY, 'MarkerSize', 10, 'LineWidth', 2);
plot(2*pi, alpha, 'x', 'Color', NAVY, 'MarkerSize', 10, 'LineWidth', 2);
plot(2*pi, -alpha, 'x', 'Color', NAVY, 'MarkerSize', 10, 'LineWidth', 2);
hold off;
xlim([-0.5, 2*pi + 0.5]); ylim([-2, 2]);
xticks([0, pi/2, pi, 3*pi/2, 2*pi]);
xticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'});
xlabel('$\mathrm{Re}\,\theta$', 'Interpreter', 'latex');
ylabel('$\mathrm{Im}\,\theta$', 'Interpreter', 'latex');
title('Poles in the complex $\theta$-plane', 'Interpreter', 'latex');
grid on; set(gca, 'GridAlpha', 0.3);

exportgraphics(gcf, fullfile(output_dir, 'poisson_kernel.pdf'), 'ContentType', 'vector');
exportgraphics(gcf, fullfile(output_dir, 'poisson_kernel.png'), 'Resolution', 300);
fprintf('Figure saved to %s\n', fullfile(output_dir, 'poisson_kernel.pdf'));
