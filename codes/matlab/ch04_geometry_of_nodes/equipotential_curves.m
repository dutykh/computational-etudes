%% equipotential_curves.m
%
% Visualizes the potential theory explanation for polynomial interpolation
% convergence. Based on Trefethen's Spectral Methods in MATLAB, Program 10.
%
% The potential function and its equipotential curves determine the
% region of convergence for polynomial interpolation.
%
% For Chebyshev density, equipotentials are Bernstein ellipses.
%
% Author: Dr. Denys Dutykh
%         Mathematics Department
%         Khalifa University of Science and Technology
%         Abu Dhabi, UAE
%
% Part of the book "Computational Études: A Spectral Approach"

clear; close all; clc;

%% Configuration
N_GRID = 500;

% Book color scheme
NAVY = [20, 45, 110] / 255;
SKY = [120, 150, 210] / 255;
CORAL = [231, 76, 60] / 255;
TEAL = [22, 160, 133] / 255;

% Output path
script_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch04', 'matlab');
output_file = fullfile(output_dir, 'equipotential_curves.pdf');

%% Create complex grid
x = linspace(-1.5, 1.5, N_GRID);
y = linspace(-1.0, 1.0, N_GRID);
[X, Y] = meshgrid(x, y);
Z = X + 1i * Y;

%% Compute potentials
% Uniform (equispaced) potential
N_nodes = 64;
x_nodes = linspace(-1, 1, N_nodes);
phi_uniform = zeros(size(Z));
for k = 1:N_nodes
    phi_uniform = phi_uniform + log(abs(Z - x_nodes(k)) + 1e-15);
end
phi_uniform = phi_uniform / N_nodes;

% Chebyshev potential: phi = log|z + sqrt(z^2 - 1)| - log(2)
sqrt_term = sqrt(Z.^2 - 1);
result = Z + sqrt_term;
mask = abs(result) < 1;
result(mask) = Z(mask) - sqrt_term(mask);
phi_chebyshev = log(abs(result)) - log(2);

%% Contour levels
levels_uniform = linspace(-0.6, 0.6, 13);
levels_chebyshev = log([1.1, 1.2, 1.4, 1.6, 2.0, 2.5, 3.0]) - log(2);

% Runge function singularity location
pole_y = 0.2;

%% Create figure
fig = figure('Units', 'inches', 'Position', [1, 1, 11, 4.5]);

% Left panel: Uniform potential
subplot(1, 2, 1);
hold on; box on;

contour(X, Y, phi_uniform, levels_uniform, 'Color', SKY, 'LineWidth', 0.8);

% Draw interval [-1,1]
plot([-1, 1], [0, 0], 'Color', NAVY, 'LineWidth', 2);
plot([-1, 1], [0, 0], 'o', 'Color', NAVY, 'MarkerSize', 6, 'MarkerFaceColor', NAVY);

% Mark singularities
plot(0, pole_y, 'x', 'Color', CORAL, 'MarkerSize', 10, 'LineWidth', 2);
plot(0, -pole_y, 'x', 'Color', CORAL, 'MarkerSize', 10, 'LineWidth', 2);

xlabel('$\mathrm{Re}(z)$', 'Interpreter', 'latex');
ylabel('$\mathrm{Im}(z)$', 'Interpreter', 'latex');
title('Equispaced Nodes (Uniform Density)', 'FontSize', 11);
xlim([-1.5, 1.5]);
ylim([-1, 1]);
axis equal;

h = plot(nan, nan, 'x', 'Color', CORAL, 'MarkerSize', 10, 'LineWidth', 2);
legend(h, sprintf('Poles at $\\pm%.1fi$', pole_y), 'Interpreter', 'latex', 'Location', 'northeast');

set(gca, 'FontSize', 10, 'TickLabelInterpreter', 'latex');
grid on;
set(gca, 'GridAlpha', 0.2);

% Right panel: Chebyshev potential
subplot(1, 2, 2);
hold on; box on;

contour(X, Y, real(phi_chebyshev), levels_chebyshev, 'Color', TEAL, 'LineWidth', 0.8);

% Draw interval [-1,1]
plot([-1, 1], [0, 0], 'Color', NAVY, 'LineWidth', 2);
plot([-1, 1], [0, 0], 'o', 'Color', NAVY, 'MarkerSize', 6, 'MarkerFaceColor', NAVY);

% Mark singularities
plot(0, pole_y, 'x', 'Color', CORAL, 'MarkerSize', 10, 'LineWidth', 2);
plot(0, -pole_y, 'x', 'Color', CORAL, 'MarkerSize', 10, 'LineWidth', 2);

% Critical Bernstein ellipse
rho_critical = abs(0.2i + sqrt(-1.04));
theta = linspace(0, 2*pi, 200);
a_crit = (rho_critical + 1/rho_critical) / 2;
b_crit = (rho_critical - 1/rho_critical) / 2;
plot(a_crit * cos(theta), b_crit * sin(theta), '--', 'Color', CORAL, 'LineWidth', 1.5);

xlabel('$\mathrm{Re}(z)$', 'Interpreter', 'latex');
ylabel('$\mathrm{Im}(z)$', 'Interpreter', 'latex');
title('Chebyshev Nodes (Bernstein Ellipses)', 'FontSize', 11);
xlim([-1.5, 1.5]);
ylim([-1, 1]);
axis equal;

h1 = plot(nan, nan, 'x', 'Color', CORAL, 'MarkerSize', 10, 'LineWidth', 2);
h2 = plot(nan, nan, '--', 'Color', CORAL, 'LineWidth', 1.5);
legend([h1, h2], {sprintf('Poles at $\\pm%.1fi$', pole_y), ...
        sprintf('Critical ellipse ($\\rho\\approx%.2f$)', rho_critical)}, ...
       'Interpreter', 'latex', 'Location', 'northeast');

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

%% Print analysis
fprintf('\nPotential Theory Analysis:\n');
fprintf('%s\n', repmat('-', 1, 60));
fprintf('Runge function poles at z = +/- %.1fi\n', pole_y);
fprintf('\nChebyshev (Bernstein ellipse) analysis:\n');
fprintf('  Critical rho = %.4f\n', rho_critical);
fprintf('  Convergence rate: O(rho^{-N}) = O(%.4f^N)\n', 1/rho_critical);
fprintf('%s\n', repmat('-', 1, 60));
