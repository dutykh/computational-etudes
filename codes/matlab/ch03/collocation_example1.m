%% collocation_example1.m
%
% Demonstrates the collocation (pseudospectral) method for solving a boundary
% value problem with a three-coefficient polynomial approximation.
%
% Problem:
%     u''(x) - (4x² + 2)u(x) = 0,   -1 ≤ x ≤ 1
%     u(-1) = 1,  u(1) = 1
%
% Exact solution:
%     u_exact(x) = exp(x² - 1)
%
% Trial function (automatically satisfies BCs):
%     u₂(x) = 1 + (1 - x²)(a₀ + a₁x + a₂x²)
%
% Collocation points:
%     x = -1/2, 0, 1/2
%
% Result:
%     a₀ = -73/118,  a₁ = 0,  a₂ = -14/59
%
% Author: Dr. Denys Dutykh
%         Mathematics Department
%         Khalifa University of Science and Technology
%         Abu Dhabi, UAE
%
% Part of the book "Computational Études: A Spectral Approach"

clear; close all; clc;

%% Configuration
N_POINTS = 500;

% Book color scheme
NAVY = [20, 45, 110] / 255;
SKY = [120, 150, 210] / 255;
CORAL = [231, 76, 60] / 255;

% Output path
script_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch03', 'matlab');
output_file = fullfile(output_dir, 'collocation_example1.pdf');

%% Exact solution
u_exact = @(x) exp(x.^2 - 1);

%% Approximate solution coefficients
% These are found by solving the collocation system (see derivation in chapter)
a0 = -73/118;
a1 = 0;
a2 = -14/59;

% Approximate solution using trial function form
u_approx = @(x) 1 + (1 - x.^2) .* (a0 + a1*x + a2*x.^2);

%% Compute solutions
x = linspace(-1, 1, N_POINTS)';
u_ex = u_exact(x);
u_ap = u_approx(x);
error_val = u_ex - u_ap;

% Collocation points
x_coll = [-0.5; 0; 0.5];
u_coll = u_approx(x_coll);

%% Create two-panel figure
fig = figure('Units', 'inches', 'Position', [1, 1, 9, 3.5]);

% Left panel: Exact vs Approximate
subplot(1, 2, 1);
hold on; box on;

plot(x, u_ex, 'Color', NAVY, 'LineWidth', 1.5, 'DisplayName', 'Exact');

% Plot approximate with markers at intervals
idx = 1:25:N_POINTS;
plot(x(idx), u_ap(idx), 'o', 'Color', SKY, 'MarkerSize', 4, ...
     'DisplayName', 'Collocation ($N=3$)');

% Mark collocation points
plot(x_coll, u_coll, 's', 'Color', CORAL, 'MarkerSize', 8, ...
     'MarkerFaceColor', 'none', 'LineWidth', 1.5, ...
     'DisplayName', 'Collocation points');

xlabel('$x$', 'Interpreter', 'latex');
ylabel('$u(x)$', 'Interpreter', 'latex');
xlim([-1, 1]);
title('Exact vs Approximate', 'FontSize', 11);
legend('Location', 'north', 'Interpreter', 'latex');

set(gca, 'FontSize', 10, 'TickLabelInterpreter', 'latex');
grid on;
set(gca, 'GridAlpha', 0.2);

% Right panel: Error
subplot(1, 2, 2);
hold on; box on;

plot(x, error_val, 'Color', NAVY, 'LineWidth', 1.5);
yline(0, '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5);

xlabel('$x$', 'Interpreter', 'latex');
ylabel('$u_{\mathrm{exact}} - u_{\mathrm{approx}}$', 'Interpreter', 'latex');
xlim([-1, 1]);
title('Error', 'FontSize', 11);

% Annotate max error
[max_err, max_idx] = max(abs(error_val));
text(0.3, error_val(max_idx) - 0.005, sprintf('Max error: %.4f', max_err), ...
     'FontSize', 9, 'Color', NAVY);

set(gca, 'FontSize', 10, 'TickLabelInterpreter', 'latex');
grid on;
set(gca, 'GridAlpha', 0.2);

%% Save figure
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Save as PDF
exportgraphics(fig, output_file, 'ContentType', 'vector');
fprintf('Figure saved to: %s\n', output_file);

% Save as PNG
png_file = strrep(output_file, '.pdf', '.png');
exportgraphics(fig, png_file, 'Resolution', 300);

%% Print error table
fprintf('\nError comparison at key points:\n');
fprintf('%s\n', repmat('-', 1, 60));
fprintf('%8s %12s %12s %12s\n', 'x', 'u_exact', 'u_approx', 'error');
fprintf('%s\n', repmat('-', 1, 60));
test_pts = [-1.0, -0.5, 0.0, 0.5, 1.0];
for xi = test_pts
    ue = u_exact(xi);
    ua = u_approx(xi);
    err = ue - ua;
    fprintf('%8.2f %12.5f %12.5f %12.5f\n', xi, ue, ua, err);
end
fprintf('%s\n', repmat('-', 1, 60));
fprintf('Maximum absolute error: %.5f\n', max(abs(error_val)));
