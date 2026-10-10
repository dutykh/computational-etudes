%% lagrange_basis.m
%
% Visualizes Lagrange basis polynomials for equispaced vs Chebyshev nodes.
%
% The Lagrange basis polynomials L_k(x) satisfy L_k(x_j) = delta_{kj}.
% For equispaced nodes, they develop large oscillations near boundaries.
% For Chebyshev nodes, they remain well-behaved.
%
% Author: Dr. Denys Dutykh
%         Mathematics Department
%         Khalifa University of Science and Technology
%         Abu Dhabi, UAE
%
% Part of the book "Computational Études: A Spectral Approach"

clear; close all; clc;

%% Configuration
N = 10;
N_FINE = 1000;

% Book color scheme
NAVY = [20, 45, 110] / 255;
SKY = [120, 150, 210] / 255;
CORAL = [231, 76, 60] / 255;
TEAL = [22, 160, 133] / 255;
PURPLE = [142, 68, 173] / 255;
ORANGE = [230, 126, 34] / 255;

% Output path
script_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch04', 'matlab');
output_file = fullfile(output_dir, 'lagrange_basis.pdf');

%% Generate nodes
x_equi = linspace(-1, 1, N+1)';
j = 0:N;
x_cheb = cos(j * pi / N)';

% Fine grid
x_fine = linspace(-1, 1, N_FINE)';

% Indices to plot
indices = [1, 2, N/2+1, N, N+1];  % 1-indexed
colors = {CORAL, ORANGE, TEAL, SKY, PURPLE};

%% Create figure
fig = figure('Units', 'inches', 'Position', [1, 1, 11, 4.5]);

% Left panel: Equispaced basis functions
subplot(1, 2, 1);
hold on; box on;

for i = 1:length(indices)
    k = indices(i);
    color = colors{i};
    L_k = lagrange_basis_k(x_equi, k, x_fine);
    plot(x_fine, L_k, 'Color', color, 'LineWidth', 1.2, ...
         'DisplayName', sprintf('$L_{%d}(x)$', k-1));
end

% Mark nodes
plot(x_equi, zeros(size(x_equi)), 'o', 'Color', NAVY, 'MarkerSize', 5, ...
     'MarkerFaceColor', NAVY, 'HandleVisibility', 'off');

yline(0, '-', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5, 'HandleVisibility', 'off');
yline(1, '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5, 'Alpha', 0.5, 'HandleVisibility', 'off');

xlabel('$x$', 'Interpreter', 'latex');
ylabel('$L_k(x)$', 'Interpreter', 'latex');
title(sprintf('Equispaced Nodes ($N = %d$)', N), 'FontSize', 11);
xlim([-1, 1]);
ylim([-2.5, 2.5]);
legend('Location', 'northwest', 'Interpreter', 'latex', 'NumColumns', 2);

text(-0.5, -2.2, 'Large oscillations near boundaries', 'FontSize', 9, 'Color', [0.5, 0.5, 0.5]);

set(gca, 'FontSize', 10, 'TickLabelInterpreter', 'latex');
grid on;
set(gca, 'GridAlpha', 0.2);

% Right panel: Chebyshev basis functions
subplot(1, 2, 2);
hold on; box on;

for i = 1:length(indices)
    k = indices(i);
    color = colors{i};
    L_k = lagrange_basis_k(x_cheb, k, x_fine);
    plot(x_fine, L_k, 'Color', color, 'LineWidth', 1.2, ...
         'DisplayName', sprintf('$L_{%d}(x)$', k-1));
end

% Mark nodes
plot(x_cheb, zeros(size(x_cheb)), 'o', 'Color', NAVY, 'MarkerSize', 5, ...
     'MarkerFaceColor', NAVY, 'HandleVisibility', 'off');

yline(0, '-', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5, 'HandleVisibility', 'off');
yline(1, '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5, 'Alpha', 0.5, 'HandleVisibility', 'off');

xlabel('$x$', 'Interpreter', 'latex');
ylabel('$L_k(x)$', 'Interpreter', 'latex');
title(sprintf('Chebyshev Nodes ($N = %d$)', N), 'FontSize', 11);
xlim([-1, 1]);
ylim([-2.5, 2.5]);
legend('Location', 'northwest', 'Interpreter', 'latex', 'NumColumns', 2);

text(0.3, 1.5, 'Bounded oscillations', 'FontSize', 9, 'Color', [0.5, 0.5, 0.5]);

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

%% Print maximum values
fprintf('\nMaximum |L_k(x)| for N = %d:\n', N);
fprintf('%s\n', repmat('-', 1, 50));
fprintf('%4s %15s %15s\n', 'k', 'Equispaced', 'Chebyshev');
fprintf('%s\n', repmat('-', 1, 50));
for k = 1:N+1
    L_equi = lagrange_basis_k(x_equi, k, x_fine);
    L_cheb = lagrange_basis_k(x_cheb, k, x_fine);
    fprintf('%4d %15.4f %15.4f\n', k-1, max(abs(L_equi)), max(abs(L_cheb)));
end
fprintf('%s\n', repmat('-', 1, 50));

% Compute Lebesgue constants
lebesgue_equi = zeros(size(x_fine));
lebesgue_cheb = zeros(size(x_fine));
for k = 1:N+1
    lebesgue_equi = lebesgue_equi + abs(lagrange_basis_k(x_equi, k, x_fine));
    lebesgue_cheb = lebesgue_cheb + abs(lagrange_basis_k(x_cheb, k, x_fine));
end
fprintf('\nLebesgue constants:\n');
fprintf('  Equispaced: Lambda_%d = %.2f\n', N, max(lebesgue_equi));
fprintf('  Chebyshev:  Lambda_%d = %.2f\n', N, max(lebesgue_cheb));
fprintf('  Ratio: %.1fx\n', max(lebesgue_equi)/max(lebesgue_cheb));

%% Lagrange basis function
function L_k = lagrange_basis_k(x_nodes, k, x_eval)
    n = length(x_nodes);
    L_k = ones(size(x_eval));

    for j = 1:n
        if j ~= k
            L_k = L_k .* (x_eval - x_nodes(j)) / (x_nodes(k) - x_nodes(j));
        end
    end
end
