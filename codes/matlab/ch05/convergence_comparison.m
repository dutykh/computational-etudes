%% convergence_comparison.m
%
% The "Computational Étude" of Chapter 5: Comparing the accuracy of finite
% difference methods (orders 2, 4, 6) against the spectral method.
%
% Test function: u(x) = 1 / (2 + sin(x)) on [0, 2π)
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
% Part of "Computational Études: A Spectral Approach"
% https://github.com/dutykh/computational-etudes

clear; close all; clc;

%% Configuration
N_values = [4, 6, 8, 10, 12, 14, 16, 20, 24, 28, 32, 40, 48, 56, 64];
orders = [2, 4, 6];

% Colors
NAVY = [0.078 0.176 0.431];
SKY = [0.471 0.588 0.824];
CORAL = [0.906 0.298 0.235];
TEAL = [0.086 0.627 0.522];

% Output path
script_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch05', 'matlab');
output_file = fullfile(output_dir, 'convergence_comparison.pdf');

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% Test function and derivative
u_func = @(x) 1 ./ (2 + sin(x));
u_deriv = @(x) -cos(x) ./ (2 + sin(x)).^2;

%% Compute errors
errors_fd = zeros(length(orders), length(N_values));
errors_spectral = zeros(1, length(N_values));

for i = 1:length(N_values)
    N = N_values(i);
    h = 2 * pi / N;
    x = h * (0:N-1)';

    u = u_func(x);
    du_exact = u_deriv(x);

    % Finite difference methods
    for j = 1:length(orders)
        order = orders(j);
        if order < N
            [D_fd, ~] = fd_matrix_periodic(N, order);
            du_fd = D_fd * u;
            errors_fd(j, i) = max(abs(du_fd - du_exact));
        else
            errors_fd(j, i) = NaN;
        end
    end

    % Spectral method
    [D_spec, ~] = spectral_matrix_periodic(N);
    du_spec = D_spec * u;
    errors_spectral(i) = max(abs(du_spec - du_exact));
end

%% Create figure
fig = figure('Units', 'inches', 'Position', [1, 1, 8, 5.5]);

% Plot finite difference errors
colors = {CORAL, TEAL, SKY};
markers = {'o', 's', '^'};
for j = 1:length(orders)
    semilogy(N_values, errors_fd(j, :), ['-' markers{j}], ...
             'Color', colors{j}, 'LineWidth', 1.5, 'MarkerSize', 5, ...
             'DisplayName', sprintf('FD%d (order %d)', orders(j), orders(j)));
    hold on;
end

% Plot spectral errors
semilogy(N_values, errors_spectral, '-d', 'Color', NAVY, ...
         'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'Spectral');

% Reference lines
N_ref = linspace(8, 64, 100);
C2 = errors_fd(1, 5) * N_values(5)^2;
C4 = errors_fd(2, 5) * N_values(5)^4;
C6 = errors_fd(3, 5) * N_values(5)^6;

% Use lighter colors to simulate transparency (blend with white at 40% opacity)
CORAL_LIGHT = CORAL * 0.4 + [1 1 1] * 0.6;
TEAL_LIGHT = TEAL * 0.4 + [1 1 1] * 0.6;
SKY_LIGHT = SKY * 0.4 + [1 1 1] * 0.6;

semilogy(N_ref, C2 ./ N_ref.^2, '--', 'Color', CORAL_LIGHT, 'LineWidth', 1, 'HandleVisibility', 'off');
semilogy(N_ref, C4 ./ N_ref.^4, '--', 'Color', TEAL_LIGHT, 'LineWidth', 1, 'HandleVisibility', 'off');
semilogy(N_ref, C6 ./ N_ref.^6, '--', 'Color', SKY_LIGHT, 'LineWidth', 1, 'HandleVisibility', 'off');

% Machine precision line
yline(2.2e-16, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'HandleVisibility', 'off');
text(66, 4e-16, 'Machine precision', 'FontSize', 9, 'Color', [0.5 0.5 0.5], ...
     'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

% Slope annotations - positioned at end of reference lines
text(68, C2 / 68^2 * 1.8, 'O(N^{-2})', 'FontSize', 10, 'Color', CORAL, ...
     'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
text(68, C4 / 68^4 * 2.5, 'O(N^{-4})', 'FontSize', 10, 'Color', TEAL, ...
     'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
text(68, C6 / 68^6 * 3.0, 'O(N^{-6})', 'FontSize', 10, 'Color', SKY, ...
     'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

% Spectral annotation
text(12, 5e-10, {'Spectral', 'accuracy!'}, 'FontSize', 10, 'Color', NAVY, ...
     'FontWeight', 'bold', 'HorizontalAlignment', 'center');
annotation('arrow', [0.28, 0.38], [0.38, 0.22], 'Color', NAVY);

% Labels
xlabel('Number of grid points N', 'FontSize', 11);
ylabel('Maximum error ||u'' - Du||_\infty', 'FontSize', 11);
title('Differentiation Error: Finite Differences vs Spectral Method', 'FontSize', 12);

% Test function info
annotation('textbox', [0.55, 0.75, 0.25, 0.15], ...
           'String', {'Test function:', 'u(x) = 1/(2 + sin(x))', 'Domain: [0, 2\pi) (periodic)'}, ...
           'FontSize', 9, 'EdgeColor', [0.5 0.5 0.5], 'BackgroundColor', 'white');

xlim([0, 70]);
ylim([1e-16, 1e1]);
legend('Location', 'northeast', 'FontSize', 9);
grid on;
box on;

%% Save figure
exportgraphics(fig, output_file, 'ContentType', 'vector');
exportgraphics(fig, strrep(output_file, '.pdf', '.png'), 'Resolution', 300);

fprintf('Figure saved to: %s\n', output_file);

%% Print results table
fprintf('\nConvergence Study Results:\n');
fprintf('======================================================================\n');
fprintf('%4s  %12s  %12s  %12s  %12s\n', 'N', 'FD2', 'FD4', 'FD6', 'Spectral');
fprintf('----------------------------------------------------------------------\n');
for i = 1:length(N_values)
    fprintf('%4d  %12.4e  %12.4e  %12.4e  %12.4e\n', ...
            N_values(i), errors_fd(1, i), errors_fd(2, i), errors_fd(3, i), errors_spectral(i));
end
fprintf('======================================================================\n');

close(fig);
