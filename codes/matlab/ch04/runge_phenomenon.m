%% runge_phenomenon.m
%
% Demonstrates the Runge phenomenon: the failure of polynomial interpolation
% on equispaced grids for certain smooth functions. Named after Carl Runge
% who discovered this counterintuitive behavior in 1901.
%
% The Runge function f(x) = 1 / (1 + 25x²) is smooth and infinitely
% differentiable, yet polynomial interpolation on equispaced nodes diverges
% as the degree increases.
%
% Author: Dr. Denys Dutykh
%         Mathematics Department
%         Khalifa University of Science and Technology
%         Abu Dhabi, UAE
%
% Part of the book "Computational Études: A Spectral Approach"

clear; close all; clc;

%% Configuration
N_FINE = 1000;

% Book color scheme
NAVY = [20, 45, 110] / 255;
SKY = [120, 150, 210] / 255;
CORAL = [231, 76, 60] / 255;
TEAL = [22, 160, 133] / 255;

% Output path
script_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch04', 'matlab');
output_file = fullfile(output_dir, 'runge_phenomenon.pdf');

%% The Runge function
runge = @(x) 1 ./ (1 + 25*x.^2);

%% Fine grid for plotting
x_fine = linspace(-1, 1, N_FINE)';
f_exact = runge(x_fine);

%% Interpolation for different N values
N_values = [6, 10, 14];
colors = {SKY, TEAL, CORAL};

%% Create figure
fig = figure('Units', 'inches', 'Position', [1, 1, 7, 5]);
hold on; box on;

% Plot exact function
plot(x_fine, f_exact, 'Color', NAVY, 'LineWidth', 2, 'DisplayName', 'Runge function');

% Interpolate for each N
for i = 1:length(N_values)
    N = N_values(i);
    color = colors{i};

    % Equispaced nodes
    x_nodes = linspace(-1, 1, N+1)';
    f_nodes = runge(x_nodes);

    % Lagrange interpolation
    p_interp = lagrange_interp(x_nodes, f_nodes, x_fine);

    % Plot interpolant
    plot(x_fine, p_interp, '--', 'Color', color, 'LineWidth', 1.2, ...
         'DisplayName', sprintf('$p_{%d}(x)$ (equispaced)', N));

    % Plot nodes
    plot(x_nodes, f_nodes, 'o', 'Color', color, 'MarkerSize', 4, ...
         'HandleVisibility', 'off');
end

xlabel('$x$', 'Interpreter', 'latex');
ylabel('$f(x)$', 'Interpreter', 'latex');
xlim([-1, 1]);
ylim([-0.5, 1.5]);
title('Runge Phenomenon: Failure of Equispaced Interpolation', 'FontSize', 11);
legend('Location', 'northeast', 'Interpreter', 'latex');

% Add annotation
text(-0.6, -0.35, 'Oscillations grow near boundaries', 'FontSize', 9, 'Color', [0.5, 0.5, 0.5]);

% Add formula
text(-0.9, 0.9, '$f(x) = \frac{1}{1 + 25x^2}$', 'Interpreter', 'latex', 'FontSize', 11);

set(gca, 'FontSize', 10, 'TickLabelInterpreter', 'latex');
grid on;
set(gca, 'GridAlpha', 0.2);

%% Save figure
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

exportgraphics(fig, output_file, 'ContentType', 'vector');
fprintf('Figure saved to: %s\n', output_file);

png_file = strrep(output_file, '.pdf', '.png');
exportgraphics(fig, png_file, 'Resolution', 300);

%% Print error statistics
fprintf('\nMaximum interpolation error at boundaries:\n');
fprintf('%s\n', repmat('-', 1, 50));
for N = N_values
    x_nodes = linspace(-1, 1, N+1)';
    f_nodes = runge(x_nodes);
    p_interp = lagrange_interp(x_nodes, f_nodes, x_fine);
    max_error = max(abs(runge(x_fine) - p_interp));
    fprintf('N = %2d: max |f - p_N| = %.4f\n', N, max_error);
end
fprintf('%s\n', repmat('-', 1, 50));
fprintf('Note: Error INCREASES with N for equispaced nodes!\n');

%% Lagrange interpolation function
function p = lagrange_interp(x_nodes, f_nodes, x_eval)
    n = length(x_nodes);
    p = zeros(size(x_eval));

    for k = 1:n
        L_k = ones(size(x_eval));
        for j = 1:n
            if j ~= k
                L_k = L_k .* (x_eval - x_nodes(j)) / (x_nodes(k) - x_nodes(j));
            end
        end
        p = p + f_nodes(k) * L_k;
    end
end
