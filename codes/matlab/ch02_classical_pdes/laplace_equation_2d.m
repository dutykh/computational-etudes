%% laplace_equation_2d.m
%
% Visualizes the solution of the Laplace equation in a periodic strip using the
% truncated Fourier series derived in Chapter 2. The boundary condition is a
% smooth function demonstrating how different Fourier modes decay toward the
% interior.
%
% Laplace equation in strip D = [0, 2*pi] x [0, 1]:
%     u_xx + u_yy = 0
%
% Boundary conditions:
%     u(x, 0) = f(x) = sin(x) + 0.5*sin(3x)   (prescribed at bottom)
%     u(x, 1) = 0                              (zero at top)
%     u periodic in x with period 2*pi
%
% Analytical solution:
%     u(x,y) = a_0*(1-y) + sum [a_n*cos(nx) + b_n*sin(nx)] * sinh(n*(1-y))/sinh(n)
%
% where a_n, b_n are Fourier coefficients of f(x).
%
% Author: Dr. Denys Dutykh
%         Mathematics Department
%         Khalifa University of Science and Technology
%         Abu Dhabi, UAE
%
% Part of the book "Computational Études: A Spectral Approach"
% https://github.com/dutykh/computational-etudes

clear; clc; close all;

%% Configuration
N_MODES = 50;           % Number of Fourier modes in truncation
NX = 300;               % Grid points in x direction
NY = 150;               % Grid points in y direction

% Output path (relative to this script's location)
script_dir = fileparts(mfilename('fullpath'));
output_file = fullfile(script_dir, '..', '..', '..', 'book', 'figures', 'ch02', 'matlab', 'laplace_solution.pdf');

%% Compute Fourier coefficients for f(x) = sin(x) + 0.5*sin(3x)
% a_0 = 0 (integral of sin over [0, 2*pi] is zero)
% a_n = 0 for all n (no cosine terms)
% b_1 = 1, b_3 = 0.5, b_n = 0 otherwise

a0 = 0.0;
a_n = zeros(1, N_MODES + 1);  % 1-indexed: a_n(n+1) = coefficient of cos(nx)
b_n = zeros(1, N_MODES + 1);  % 1-indexed: b_n(n+1) = coefficient of sin(nx)

b_n(2) = 1.0;   % b_1 = 1
b_n(4) = 0.5;   % b_3 = 0.5

%% Create 2D grid
x = linspace(0, 2*pi, NX);
y = linspace(0, 1, NY);
[X, Y] = meshgrid(x, y);

%% Evaluate solution on grid
% u(x,y) = a_0*(1-y) + sum [a_n*cos(nx) + b_n*sin(nx)] * sinh(n*(1-y))/sinh(n)

U = a0 * (1 - Y);

for n = 1:N_MODES
    % Skip if both coefficients are zero
    if abs(a_n(n + 1)) < 1e-15 && abs(b_n(n + 1)) < 1e-15
        continue;
    end

    % y-dependent factor: sinh(n*(1-y)) / sinh(n)
    sinh_n = sinh(n);
    if sinh_n > 1e-10
        y_factor = sinh(n * (1 - Y)) / sinh_n;
    else
        y_factor = zeros(size(Y));
    end

    % Add contribution from this mode
    U = U + (a_n(n + 1) * cos(n * X) + b_n(n + 1) * sin(n * X)) .* y_factor;
end

%% Create figure
fig = figure('Units', 'inches', 'Position', [1, 1, 6, 3.5]);

% Determine symmetric color limits
vmax = max(abs(U(:)));
vmin = -vmax;

% Create filled contour plot
levels = linspace(vmin, vmax, 25);
contourf(X, Y, U, levels, 'LineStyle', 'none');
hold on;

% Add contour lines for emphasis
contour(X, Y, U, linspace(vmin, vmax, 13), 'Color', [0.2, 0.2, 0.2], 'LineWidth', 0.3);

% Custom colormap: navy -> white -> sky blue
n_colors = 256;
navy_rgb = [20, 45, 110] / 255;
white_rgb = [1, 1, 1];
sky_rgb = [120, 150, 210] / 255;

cmap = zeros(n_colors, 3);
for i = 1:n_colors/2
    t = (i - 1) / (n_colors/2 - 1);
    cmap(i, :) = navy_rgb * (1 - t) + white_rgb * t;
end
for i = 1:n_colors/2
    t = (i - 1) / (n_colors/2 - 1);
    cmap(n_colors/2 + i, :) = white_rgb * (1 - t) + sky_rgb * t;
end
colormap(cmap);

% Colorbar
cbar = colorbar;
ylabel(cbar, '$u(x,y)$', 'Interpreter', 'latex', 'FontSize', 11);

% Axis labels
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 11);
ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 11);

% Custom x-ticks with pi notation
set(gca, 'XTick', [0, pi/2, pi, 3*pi/2, 2*pi]);
set(gca, 'XTickLabel', {'$0$', '$\pi/2$', '$\pi$', '$3\pi/2$', '$2\pi$'}, ...
    'TickLabelInterpreter', 'latex');

% y-ticks
set(gca, 'YTick', [0, 0.25, 0.5, 0.75, 1.0]);

% Color limits
caxis([vmin, vmax]);

box on;

%% Save figure
% Ensure output directory exists
[output_dir, ~, ~] = fileparts(output_file);
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Export to PDF
exportgraphics(fig, output_file, 'ContentType', 'vector', 'Resolution', 300);
fprintf('Figure saved to: %s\n', output_file);

close(fig);
