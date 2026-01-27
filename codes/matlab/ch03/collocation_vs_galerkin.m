%% collocation_vs_galerkin.m
%
% Compares the Collocation (pseudospectral) and Galerkin methods for solving
% a reaction-diffusion boundary value problem.
%
% Problem:
%     u''(x) - 4u(x) = -1,   -1 ≤ x ≤ 1
%     u(-1) = 0,  u(1) = 0
%
% Exact solution:
%     u_exact(x) = (1/4)(1 - cosh(2x)/cosh(2))
%
% Basis functions (satisfy homogeneous BCs):
%     φ₀(x) = (1 - x²)
%     φ₁(x) = (1 - x²)x² = x² - x⁴
%
% Trial function:
%     u₁(x) = a₀φ₀(x) + a₁φ₁(x)
%
% Methods compared:
%     1. Collocation: Force residual to zero at x = 0 and x = 0.5
%     2. Galerkin: Force residual orthogonal to basis functions
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
GREEN = [39, 174, 96] / 255;

% Output path
script_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch03', 'matlab');
output_file = fullfile(output_dir, 'collocation_vs_galerkin.pdf');

%% Define functions
% Exact solution
u_exact = @(x) 0.25 * (1 - cosh(2*x) / cosh(2));

% Basis functions
phi0 = @(x) 1 - x.^2;
phi1 = @(x) x.^2 - x.^4;

% Second derivatives
phi0_xx = @(x) -2 * ones(size(x));
phi1_xx = @(x) 2 - 12*x.^2;

% Operator L = d²/dx² - 4
L_phi0 = @(x) phi0_xx(x) - 4*phi0(x);
L_phi1 = @(x) phi1_xx(x) - 4*phi1(x);

%% Solve Collocation system
% Collocation points: x = 0 and x = 0.5
x1 = 0.0;
x2 = 0.5;

% Build matrix A where A*[a0; a1] = b
A_coll = [L_phi0(x1), L_phi1(x1);
          L_phi0(x2), L_phi1(x2)];

% RHS: f = -1, so R = Lu - f = Lu + 1 = 0 means Lu = -1
b_coll = [-1; -1];

% Solve
coeffs_coll = A_coll \ b_coll;
a0_coll = coeffs_coll(1);
a1_coll = coeffs_coll(2);

fprintf('Collocation coefficients:\n');
fprintf('  a₀ = %.6f\n', a0_coll);
fprintf('  a₁ = %.6f\n', a1_coll);

%% Solve Galerkin system
% Matrix entries: A_{ij} = ∫_{-1}^{1} L[φⱼ](x) φᵢ(x) dx

A00 = integral(@(x) L_phi0(x) .* phi0(x), -1, 1);
A01 = integral(@(x) L_phi1(x) .* phi0(x), -1, 1);
A10 = integral(@(x) L_phi0(x) .* phi1(x), -1, 1);
A11 = integral(@(x) L_phi1(x) .* phi1(x), -1, 1);

A_gal = [A00, A01; A10, A11];

% RHS: b_i = ∫_{-1}^{1} f(x) φᵢ(x) dx = ∫_{-1}^{1} (-1) φᵢ(x) dx
b0 = integral(@(x) -1 * phi0(x), -1, 1);
b1 = integral(@(x) -1 * phi1(x), -1, 1);

b_gal = [b0; b1];

% Solve
coeffs_gal = A_gal \ b_gal;
a0_gal = coeffs_gal(1);
a1_gal = coeffs_gal(2);

fprintf('\nGalerkin coefficients:\n');
fprintf('  a₀ = %.6f\n', a0_gal);
fprintf('  a₁ = %.6f\n', a1_gal);

%% Define approximate solutions
u_coll_fn = @(x) a0_coll * phi0(x) + a1_coll * phi1(x);
u_gal_fn = @(x) a0_gal * phi0(x) + a1_gal * phi1(x);

%% Compute solutions on grid
x = linspace(-1, 1, N_POINTS)';
u_ex = u_exact(x);
u_coll_val = u_coll_fn(x);
u_gal_val = u_gal_fn(x);

error_coll = u_ex - u_coll_val;
error_gal = u_ex - u_gal_val;

% Collocation points for visualization
x_coll_pts = [0.0; 0.5];
u_coll_pts = u_coll_fn(x_coll_pts);

%% Create two-panel figure
fig = figure('Units', 'inches', 'Position', [1, 1, 9, 3.5]);

% Left panel: All solutions
subplot(1, 2, 1);
hold on; box on;

plot(x, u_ex, 'Color', NAVY, 'LineWidth', 1.8, 'DisplayName', 'Exact');
plot(x, u_coll_val, '--', 'Color', CORAL, 'LineWidth', 1.5, ...
     'DisplayName', 'Collocation');
plot(x, u_gal_val, '-.', 'Color', GREEN, 'LineWidth', 1.5, ...
     'DisplayName', 'Galerkin');
plot(x_coll_pts, u_coll_pts, 's', 'Color', CORAL, 'MarkerSize', 8, ...
     'MarkerFaceColor', 'none', 'LineWidth', 1.5, ...
     'DisplayName', 'Collocation pts');

xlabel('$x$', 'Interpreter', 'latex');
ylabel('$u(x)$', 'Interpreter', 'latex');
xlim([-1, 1]);
title('Solutions Comparison', 'FontSize', 11);
legend('Location', 'south', 'Interpreter', 'latex', 'NumColumns', 2);

set(gca, 'FontSize', 10, 'TickLabelInterpreter', 'latex');
grid on;
set(gca, 'GridAlpha', 0.2);

% Right panel: Errors
subplot(1, 2, 2);
hold on; box on;

plot(x, error_coll, 'Color', CORAL, 'LineWidth', 1.5, 'DisplayName', 'Collocation');
plot(x, error_gal, 'Color', GREEN, 'LineWidth', 1.5, 'DisplayName', 'Galerkin');
yline(0, '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5);

xlabel('$x$', 'Interpreter', 'latex');
ylabel('$u_{\mathrm{exact}} - u_{\mathrm{approx}}$', 'Interpreter', 'latex');
xlim([-1, 1]);
title('Error Comparison', 'FontSize', 11);
legend('Location', 'northeast', 'Interpreter', 'latex');

set(gca, 'FontSize', 10, 'TickLabelInterpreter', 'latex');
grid on;
set(gca, 'GridAlpha', 0.2);

%% Save figure
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

exportgraphics(fig, output_file, 'ContentType', 'vector');
fprintf('\nFigure saved to: %s\n', output_file);

png_file = strrep(output_file, '.pdf', '.png');
exportgraphics(fig, png_file, 'Resolution', 300);

%% Print comparison table
fprintf('\n%s\n', repmat('=', 1, 65));
fprintf('Comparison at u(0) - the central maximum\n');
fprintf('%s\n', repmat('=', 1, 65));
fprintf('%-20s %-15s %-15s\n', 'Method', 'u(0)', 'Abs. Error');
fprintf('%s\n', repmat('-', 1, 65));
fprintf('%-20s %-15.6f %-15.6f\n', 'Exact', u_exact(0), 0);
fprintf('%-20s %-15.6f %-15.6f\n', 'Collocation (N=2)', ...
        u_coll_fn(0), abs(u_exact(0) - u_coll_fn(0)));
fprintf('%-20s %-15.6f %-15.6f\n', 'Galerkin (N=2)', ...
        u_gal_fn(0), abs(u_exact(0) - u_gal_fn(0)));
fprintf('%s\n', repmat('=', 1, 65));

% RMS errors
rms_coll = sqrt(mean(error_coll.^2));
rms_gal = sqrt(mean(error_gal.^2));
fprintf('\nRMS Error (Collocation): %.6f\n', rms_coll);
fprintf('RMS Error (Galerkin):    %.6f\n', rms_gal);

% Max errors
max_coll = max(abs(error_coll));
max_gal = max(abs(error_gal));
fprintf('\nMax Error (Collocation): %.6f\n', max_coll);
fprintf('Max Error (Galerkin):    %.6f\n', max_gal);
