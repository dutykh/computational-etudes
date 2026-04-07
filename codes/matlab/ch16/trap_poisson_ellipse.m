% trap_poisson_ellipse.m
% Chapter 16: Integration of Periodic Functions
%
% Computational Etude 16.1: Poisson's Ellipse, the Original Paradox.
%
% Reproduces Figure 1.1 of Trefethen and Weideman, SIAM Rev. 2014:
% the trapezoidal-rule perimeter of an ellipse with semi-axes 1/(2*pi)
% and 0.6/(2*pi).  The integral is (2/pi) * E(0.36), and the
% trapezoidal error decays as 3^(-N) because the integrand has branch
% points at theta = +/- i*log(3) in the complex plane.
%
% Reference: L. N. Trefethen and J. A. C. Weideman, "The exponentially
%            convergent trapezoidal rule", SIAM Rev. 56(3):385-458 (2014).
%
% Generates Figure 16.1: poisson_ellipse.pdf
%
% Author: Dr. Denys Dutykh
%         Mathematics Department
%         Khalifa University of Science and Technology
%         Abu Dhabi, UAE
%
% Part of the book "Computational Etudes: A Spectral Approach"
% https://github.com/dutykh/computational-etudes

clear; close all; clc;

% Publication-quality figure settings
set(0, 'DefaultAxesFontSize', 11);
set(0, 'DefaultTextFontSize', 11);
set(0, 'DefaultLineLineWidth', 1.5);
set(0, 'DefaultAxesLineWidth', 0.8);

% Book color scheme
NAVY  = [20, 45, 110] / 255;
CORAL = [231, 76, 60] / 255;

% Output path
output_dir = fullfile(fileparts(mfilename('fullpath')), ...
    '..', '..', '..', 'textbook', 'figures', 'ch16', 'matlab');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

% -------------------------------------------------------------------------
% Poisson's integrand and the exact value
% -------------------------------------------------------------------------
integrand = @(t) sqrt(1 - 0.36 * sin(t).^2);

% MATLAB ellipke uses the parameter m = k^2; with m = 0.36 we get
% [K(0.36), E(0.36)].
[~, E] = ellipke(0.36);
I_exact = (2 / pi) * E;
fprintf('Exact value: I = %.16f\n\n', I_exact);

% Trapezoidal rule for N = 4, 8, ..., 200
N_values = 4:4:200;
errors = zeros(size(N_values));
for i = 1:length(N_values)
    N = N_values(i);
    theta = 2*pi*(0:N-1)/N;
    I_N = (2*pi/N) * sum(integrand(theta)) / (2*pi);
    errors(i) = abs(I_N - I_exact);
end
errors = max(errors, 1e-17);

% Theoretical envelope: 3^{-N}
theory = max(3.^(-double(N_values)), 1e-17);

% Print Poisson's table
fprintf('%4s  %22s  %14s\n', 'N/4', 'I_N', '|I_N - I|');
for i = 1:10
    N = N_values(i);
    theta = 2*pi*(0:N-1)/N;
    I_N = (2*pi/N) * sum(integrand(theta)) / (2*pi);
    fprintf('%4d  %22.18f  %14.4e\n', N/4, I_N, errors(i));
end

% -------------------------------------------------------------------------
% Plot
% -------------------------------------------------------------------------
fig = figure('Position', [100, 100, 700, 480]);
semilogy(N_values/4, errors, 'o-', 'Color', CORAL, 'MarkerSize', 4, ...
    'LineWidth', 1.1, 'MarkerFaceColor', CORAL, 'DisplayName', 'Trapezoidal rule');
hold on;
semilogy(N_values/4, theory, '--', 'Color', NAVY, 'LineWidth', 1.0, ...
    'DisplayName', '3^{-N}');
hold off;
xlabel('$N/4$ (number of independent samples)', 'Interpreter', 'latex');
ylabel('Absolute error $|I_N - I|$', 'Interpreter', 'latex');
title("Poisson's ellipse: trapezoidal convergence is geometric");
legend('Location', 'northeast', 'EdgeColor', [0.5 0.5 0.5]);
grid on;
set(gca, 'GridAlpha', 0.3);
xlim([0, 50]);
ylim([1e-18, 1e0]);

exportgraphics(gcf, fullfile(output_dir, 'poisson_ellipse.pdf'), 'ContentType', 'vector');
exportgraphics(gcf, fullfile(output_dir, 'poisson_ellipse.png'), 'Resolution', 300);
fprintf('\nFigure saved to %s\n', fullfile(output_dir, 'poisson_ellipse.pdf'));
