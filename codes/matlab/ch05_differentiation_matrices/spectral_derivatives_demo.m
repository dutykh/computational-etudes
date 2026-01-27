%% spectral_derivatives_demo.m
%
% Demonstrates the periodic spectral differentiation matrix by computing
% first and second derivatives of a smooth periodic function.
%
% Test function: u(x) = exp(sin^2(x))
% - First derivative:  u'(x)  = sin(2x) * exp(sin^2(x))
% - Second derivative: u''(x) = [sin^2(2x) + 2*cos(2x)] * exp(sin^2(x))
%
% With N = 64 grid points, the spectral method achieves near machine precision!
%
% Author: Dr. Denys Dutykh
%         Mathematics Department
%         Khalifa University of Science and Technology
%         Abu Dhabi, UAE
%
% Part of the book "Computational Études: A Spectral Approach"
% https://github.com/dutykh/computational-etudes

clear; close all; clc;

%% Configuration
% Book color scheme
NAVY = [20, 45, 110] / 255;
SKY = [120, 150, 210] / 255;
CORAL = [231, 76, 60] / 255;
TEAL = [22, 160, 133] / 255;
PURPLE = [142, 68, 173] / 255;

% Output path
script_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch05', 'matlab');
output_file = fullfile(output_dir, 'spectral_derivatives_demo.pdf');

%% Main computation
% Number of grid points
N = 64;

% Construct the spectral differentiation matrix
[D, x] = spectral_matrix_periodic(N);

% Fine grid for plotting exact solutions
x_fine = linspace(0, 2*pi, 200);

% Test function: u(x) = exp(sin^2(x))
u = exp(sin(x).^2);
u_fine = exp(sin(x_fine).^2);

% Exact first derivative: u'(x) = sin(2x) * exp(sin^2(x))
u1_exact = sin(2*x) .* exp(sin(x).^2);
u1_fine = sin(2*x_fine) .* exp(sin(x_fine).^2);

% Exact second derivative: u''(x) = [sin^2(2x) + 2*cos(2x)] * exp(sin^2(x))
u2_exact = (sin(2*x).^2 + 2*cos(2*x)) .* exp(sin(x).^2);
u2_fine = (sin(2*x_fine).^2 + 2*cos(2*x_fine)) .* exp(sin(x_fine).^2);

% Numerical derivatives using spectral differentiation matrix
u1_num = D * u(:);       % First derivative: D * u
u2_num = D * (D * u(:)); % Second derivative: D^2 * u

% Compute errors
error_u1 = max(abs(u1_num - u1_exact(:)));
error_u2 = max(abs(u2_num - u2_exact(:)));

%% Create figure
fig = figure('Units', 'inches', 'Position', [1, 1, 11, 4.5]);

% Left panel: u(x) and u'(x)
subplot(1, 2, 1);
hold on;

plot(x_fine, u_fine, '-', 'Color', NAVY, 'LineWidth', 2, ...
     'DisplayName', '$u(x) = e^{\sin^2 x}$');
plot(x_fine, u1_fine, '-', 'Color', CORAL, 'LineWidth', 2, ...
     'DisplayName', '$u''(x)$ exact');
plot(x, u1_num, 'o', 'Color', TEAL, 'MarkerSize', 7, ...
     'MarkerFaceColor', 'white', 'LineWidth', 2, ...
     'DisplayName', '$u''(x)$ spectral');

xlabel('$x$', 'Interpreter', 'latex');
ylabel('$u(x)$, $u''(x)$', 'Interpreter', 'latex');
title(sprintf('First Derivative ($N = %d$)', N), 'Interpreter', 'latex');
xlim([0, 2*pi]);
xticks([0, pi/2, pi, 3*pi/2, 2*pi]);
xticklabels({'$0$', '$\pi/2$', '$\pi$', '$3\pi/2$', '$2\pi$'});
set(gca, 'TickLabelInterpreter', 'latex');
legend('Location', 'northeast', 'Interpreter', 'latex');
grid on;
set(gca, 'GridAlpha', 0.3);

% Add error annotation
text(0.05, 0.05, sprintf('Max error: %.2e', error_u1), ...
     'Units', 'normalized', 'FontSize', 10, 'Color', TEAL, ...
     'BackgroundColor', 'white', 'EdgeColor', TEAL, ...
     'Margin', 3);

hold off;

% Right panel: u''(x)
subplot(1, 2, 2);
hold on;

plot(x_fine, u2_fine, '-', 'Color', PURPLE, 'LineWidth', 2, ...
     'DisplayName', '$u''''(x)$ exact');
plot(x, u2_num, 's', 'Color', TEAL, 'MarkerSize', 7, ...
     'MarkerFaceColor', 'white', 'LineWidth', 2, ...
     'DisplayName', '$u''''(x)$ spectral');

xlabel('$x$', 'Interpreter', 'latex');
ylabel('$u''''(x)$', 'Interpreter', 'latex');
title(sprintf('Second Derivative ($N = %d$)', N), 'Interpreter', 'latex');
xlim([0, 2*pi]);
xticks([0, pi/2, pi, 3*pi/2, 2*pi]);
xticklabels({'$0$', '$\pi/2$', '$\pi$', '$3\pi/2$', '$2\pi$'});
set(gca, 'TickLabelInterpreter', 'latex');
legend('Location', 'northeast', 'Interpreter', 'latex');
grid on;
set(gca, 'GridAlpha', 0.3);

% Add error annotation
text(0.05, 0.05, sprintf('Max error: %.2e', error_u2), ...
     'Units', 'normalized', 'FontSize', 10, 'Color', TEAL, ...
     'BackgroundColor', 'white', 'EdgeColor', TEAL, ...
     'Margin', 3);

hold off;

%% Save figure
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

exportgraphics(fig, output_file, 'ContentType', 'vector');
fprintf('Figure saved to: %s\n', output_file);

png_file = strrep(output_file, '.pdf', '.png');
exportgraphics(fig, png_file, 'Resolution', 300);
fprintf('PNG saved to: %s\n', png_file);

%% Print results
fprintf('\n');
fprintf('============================================================\n');
fprintf('Spectral Differentiation Demo\n');
fprintf('============================================================\n');
fprintf('Test function: u(x) = exp(sin²(x))\n');
fprintf('Number of grid points: N = %d\n', N);
fprintf('------------------------------------------------------------\n');
fprintf('First derivative  max error: %.4e\n', error_u1);
fprintf('Second derivative max error: %.4e\n', error_u2);
fprintf('------------------------------------------------------------\n');
fprintf('With only 64 points, we achieve near machine precision!\n');
fprintf('============================================================\n');

%% Convergence study for table
fprintf('\nConvergence Study (for table):\n');
fprintf('--------------------------------------------------\n');
fprintf('%6s  %14s  %14s\n', 'N', 'Error u''', 'Error u''''');
fprintf('--------------------------------------------------\n');

for N_test = [8, 16, 32, 64]
    [D_test, x_test] = spectral_matrix_periodic(N_test);
    u_test = exp(sin(x_test).^2);
    u1_exact_test = sin(2*x_test) .* exp(sin(x_test).^2);
    u2_exact_test = (sin(2*x_test).^2 + 2*cos(2*x_test)) .* exp(sin(x_test).^2);
    u1_num_test = D_test * u_test(:);
    u2_num_test = D_test * (D_test * u_test(:));
    err1 = max(abs(u1_num_test - u1_exact_test(:)));
    err2 = max(abs(u2_num_test - u2_exact_test(:)));
    fprintf('%6d  %14.4e  %14.4e\n', N_test, err1, err2);
end

fprintf('--------------------------------------------------\n');
