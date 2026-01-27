%% chebyshev_success.m
%
% Demonstrates successful interpolation of the Runge function using Chebyshev
% nodes. In contrast to the disastrous equispaced interpolation, Chebyshev
% points provide rapid and uniform convergence.
%
% Chebyshev-Gauss-Lobatto points:
%     x_j = cos(jπ/N),  j = 0, 1, ..., N
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
output_file = fullfile(output_dir, 'chebyshev_success.pdf');

%% The Runge function
runge = @(x) 1 ./ (1 + 25*x.^2);

%% Fine grid for plotting
x_fine = linspace(-1, 1, N_FINE)';
f_exact = runge(x_fine);

%% Different polynomial degrees
N_values = [6, 10, 14];
colors = {SKY, TEAL, CORAL};

%% Create figure with two panels
fig = figure('Units', 'inches', 'Position', [1, 1, 10, 4]);

% Left panel: Function and interpolants
subplot(1, 2, 1);
hold on; box on;

plot(x_fine, f_exact, 'Color', NAVY, 'LineWidth', 2, 'DisplayName', 'Runge function');

for i = 1:length(N_values)
    N = N_values(i);
    color = colors{i};

    % Chebyshev nodes
    j = 0:N;
    x_nodes = cos(j * pi / N)';
    f_nodes = runge(x_nodes);

    % Interpolate
    p_interp = lagrange_interp(x_nodes, f_nodes, x_fine);

    plot(x_fine, p_interp, '--', 'Color', color, 'LineWidth', 1.2, ...
         'DisplayName', sprintf('$p_{%d}(x)$ (Chebyshev)', N));
    plot(x_nodes, f_nodes, 'o', 'Color', color, 'MarkerSize', 5, ...
         'HandleVisibility', 'off');
end

xlabel('$x$', 'Interpreter', 'latex');
ylabel('$f(x)$', 'Interpreter', 'latex');
xlim([-1, 1]);
ylim([-0.1, 1.1]);
title('Chebyshev Interpolation: Success', 'FontSize', 11);
legend('Location', 'northeast', 'Interpreter', 'latex');

text(-0.9, 0.95, '$f(x) = \frac{1}{1 + 25x^2}$', 'Interpreter', 'latex', 'FontSize', 11);

set(gca, 'FontSize', 10, 'TickLabelInterpreter', 'latex');
grid on;
set(gca, 'GridAlpha', 0.2);

% Right panel: Error
subplot(1, 2, 2);
hold on; box on;

for i = 1:length(N_values)
    N = N_values(i);
    color = colors{i};

    j = 0:N;
    x_nodes = cos(j * pi / N)';
    f_nodes = runge(x_nodes);
    p_interp = lagrange_interp(x_nodes, f_nodes, x_fine);
    error_val = abs(runge(x_fine) - p_interp);

    semilogy(x_fine, error_val + 1e-16, 'Color', color, 'LineWidth', 1.5, ...
             'DisplayName', sprintf('$N = %d$', N));
end

xlabel('$x$', 'Interpreter', 'latex');
ylabel('$|f(x) - p_N(x)|$', 'Interpreter', 'latex');
xlim([-1, 1]);
ylim([1e-10, 1]);
title('Interpolation Error (log scale)', 'FontSize', 11);
legend('Location', 'northeast', 'Interpreter', 'latex');

text(0.4, 1e-3, 'Error decreases with $N$', 'Interpreter', 'latex', ...
     'FontSize', 9, 'Color', [0.5, 0.5, 0.5]);

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

%% Print comparison
fprintf('\nMaximum interpolation error (Chebyshev vs Equispaced):\n');
fprintf('%s\n', repmat('-', 1, 60));
fprintf('%4s %15s %15s %12s\n', 'N', 'Chebyshev', 'Equispaced', 'Ratio');
fprintf('%s\n', repmat('-', 1, 60));
for N = N_values
    % Chebyshev
    j = 0:N;
    x_cheb = cos(j * pi / N)';
    f_cheb = runge(x_cheb);
    p_cheb = lagrange_interp(x_cheb, f_cheb, x_fine);
    err_cheb = max(abs(runge(x_fine) - p_cheb));

    % Equispaced
    x_equi = linspace(-1, 1, N+1)';
    f_equi = runge(x_equi);
    p_equi = lagrange_interp(x_equi, f_equi, x_fine);
    err_equi = max(abs(runge(x_fine) - p_equi));

    fprintf('%4d %15.2e %15.2e %12.1fx\n', N, err_cheb, err_equi, err_equi/err_cheb);
end
fprintf('%s\n', repmat('-', 1, 60));

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
