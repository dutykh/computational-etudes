%% ho_clamped_beam.m - Clamped beam under exponential load
%
% Chapter 14: Higher-Order Boundary Value Problems
%
% Computational Etude 14.1: Clamped beam under exponential load.
%
% Solves the fourth-order BVP (biharmonic-type equation):
%
%   u^(4)(x) = e^x,   -1 < x < 1,   u(+-1) = u'(+-1) = 0
%
% using a Chebyshev spectral method with the polynomial trick from
% Trefethen's "Spectral Methods in MATLAB" (Program 38).
%
% The key idea is to write u(x) = (1 - x^2) * q(x), which automatically
% satisfies both Dirichlet conditions u(+-1) = 0.  The factor (1 - x^2)
% also enforces u'(+-1) = 0 because u'(x) = -2x*q(x) + (1 - x^2)*q'(x),
% so u'(+-1) = -2*(+-1)*q(+-1) + 0 = 0 only if q(+-1) = 0.  Instead of
% enforcing this directly, the fourth derivative operator is reformulated
% via the product rule as:
%
%   D4full = [diag(1 - x^2) * D^4 - 8 * diag(x) * D^3 - 12 * D^2] * S
%
% where S = diag(0, 1/(1 - x_1^2), ..., 1/(1 - x_{N-1}^2), 0) maps v
% back to q on interior points and annihilates boundary values.
%
% The exact solution is u(x) = e^x + c_3*x^3 + c_2*x^2 + c_1*x + c_0,
% where the constants are determined by the four clamped boundary conditions.
%
% Generates Figure: clamped_beam.pdf
%   Left panel:  solution (circles) vs exact (line) for N = 15
%   Right panel: convergence (max error) and condition number vs N
%
% Author: Dr. Denys Dutykh
%         Mathematics Department
%         Khalifa University of Science and Technology
%         Abu Dhabi, UAE
%
% Part of the book "Computational Etudes: A Spectral Approach"
% https://github.com/dutykh/computational-etudes
%
% Last modified: February 2026

clear;
close all;
clc;

% Add path to Chapter 7 functions
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'ch07'));

%% Publication-quality figure settings
set(groot, 'DefaultAxesFontSize', 10);
set(groot, 'DefaultAxesFontName', 'CMU Serif');
set(groot, 'DefaultTextFontSize', 10);
set(groot, 'DefaultTextFontName', 'CMU Serif');
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');
set(groot, 'DefaultAxesLineWidth', 0.8);
set(groot, 'DefaultLineLineWidth', 1.5);

%% Color scheme
NAVY   = [20,  45, 110] / 255;
SKY    = [120, 150, 210] / 255;
CORAL  = [231,  76,  60] / 255;
TEAL   = [22,  160, 133] / 255;
PURPLE = [142,  68, 173] / 255;
ORANGE = [230, 126,  34] / 255;

%% Output directory
script_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch14', 'matlab');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% Solve for N = 15 (demonstration)
N_demo = 15;
[D4int, x] = biharmonic_clamped_matrix(N_demo);

% Right-hand side at interior points
f = exp(x(2:N_demo));

% Solve the linear system
v = D4int \ f;

% Reconstruct full solution (boundary values are zero)
u_num = zeros(N_demo + 1, 1);
u_num(2:N_demo) = v;

% Exact solution on the full grid
u_ex = exact_solution(x);

demo_error = max(abs(u_num - u_ex));
fprintf('N = %d: max |u_num - u_exact| = %.4e\n', N_demo, demo_error);

% Fine grid for smooth exact solution curve
x_fine = linspace(-1, 1, 500)';
u_fine = exact_solution(x_fine);

%% Convergence study
N_values = 4:2:32;
errors = zeros(size(N_values));
cond_numbers = zeros(size(N_values));

for k = 1:length(N_values)
    N_val = N_values(k);
    [D4int_n, x_n] = biharmonic_clamped_matrix(N_val);
    f_n = exp(x_n(2:N_val));

    v_n = D4int_n \ f_n;

    u_n = zeros(N_val + 1, 1);
    u_n(2:N_val) = v_n;

    u_ex_n = exact_solution(x_n);
    err = max(abs(u_n - u_ex_n));
    errors(k) = max(err, 1e-16);
    cond_numbers(k) = cond(D4int_n);
end

%% Print convergence table
fprintf('\n%4s  %12s  %14s\n', 'N', 'Max Error', 'cond(D4int)');
fprintf('%s\n', repmat('-', 1, 36));
for k = 1:length(N_values)
    fprintf('%4d  %12.4e  %14.4e\n', N_values(k), errors(k), cond_numbers(k));
end
fprintf('%s\n', repmat('-', 1, 36));

%% Create figure
fig = figure('Units', 'inches', 'Position', [1, 1, 11, 4.5]);

% ---- Left panel: solution for N = 15 ----
subplot(1, 2, 1);
plot(x_fine, u_fine, '-', 'Color', NAVY, 'LineWidth', 1.5, ...
     'DisplayName', 'Exact solution');
hold on;
plot(x, u_num, 'o', 'Color', CORAL, 'MarkerSize', 6, ...
     'MarkerFaceColor', CORAL, 'MarkerEdgeColor', 'w', ...
     'LineWidth', 0.5, 'DisplayName', sprintf('Spectral ($N = %d$)', N_demo));
hold off;
xlabel('$x$');
ylabel('$u(x)$');
title('Clamped beam: $u^{(4)} = e^x$', 'FontSize', 11);
legend('Location', 'northwest', 'FontSize', 9);
xlim([-1.05, 1.05]);
grid on;
set(gca, 'GridAlpha', 0.3);
box on;

% ---- Right panel: convergence + condition number ----
% Primary y-axis: max error (semilog)
ax2 = subplot(1, 2, 2);
semilogy(N_values, errors, 'o-', 'Color', CORAL, 'LineWidth', 1.5, ...
         'MarkerSize', 5, 'MarkerFaceColor', CORAL, 'MarkerEdgeColor', 'w');
hold on;
yline(2.2e-16, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
text(33, 4e-16, '$\epsilon_{\mathrm{mach}}$', 'FontSize', 9, ...
     'Color', [0.5 0.5 0.5], 'HorizontalAlignment', 'right', ...
     'VerticalAlignment', 'bottom');
hold off;
xlabel('$N$');
ylabel('Max error', 'Color', CORAL);
ax2.YColor = CORAL;
xlim([2, 34]);
ylim([1e-16, 1e0]);
title('Convergence and conditioning', 'FontSize', 11);
grid on;
set(gca, 'GridAlpha', 0.3);
box on;

% Secondary y-axis: condition number
ax2b = axes('Position', ax2.Position, 'YAxisLocation', 'right', ...
            'Color', 'none', 'XTick', [], 'Box', 'off');
hold on;
semilogy(N_values, cond_numbers, 's--', 'Color', TEAL, 'LineWidth', 1.2, ...
         'MarkerSize', 4, 'MarkerFaceColor', TEAL, 'MarkerEdgeColor', 'w', ...
         'Parent', ax2b);
hold off;
ylabel('$\mathrm{cond}(D_4)$', 'Color', TEAL);
ax2b.YColor = TEAL;
xlim([2, 34]);

% Add combined legend
legend(ax2, {'Max error', '$\kappa(D_4)$'}, 'Location', 'center left', ...
       'FontSize', 9);

sgtitle('Computational Etude 14.1: Clamped Beam Under Exponential Load', ...
        'FontSize', 12);

%% Save figure
exportgraphics(fig, fullfile(output_dir, 'clamped_beam.pdf'), ...
               'ContentType', 'vector');
fprintf('\nFigure saved to: %s\n', fullfile(output_dir, 'clamped_beam.pdf'));

close(fig);


%% =====================================================================
%  Local functions
%  =====================================================================

function [D4int, x] = biharmonic_clamped_matrix(N)
% BIHARMONIC_CLAMPED_MATRIX  Build fourth-order operator for clamped BCs.
%
%   The polynomial trick: write u(x) = (1 - x^2) * q(x), which
%   automatically satisfies u(+-1) = 0 and u'(+-1) = 0.  The fourth
%   derivative operator is reformulated using the product rule.
%
%   Inputs:
%     N - Number of Chebyshev intervals (grid has N+1 points)
%
%   Outputs:
%     D4int - Interior block of the fourth-order operator, size (N-1) x (N-1)
%     x     - Full Chebyshev grid (N+1 points)

    [D, x] = cheb_matrix(N);

    % Powers of the differentiation matrix
    D2 = D * D;
    D3 = D2 * D;
    D4 = D2 * D2;

    % S maps v -> q on interior points, zeros at boundaries
    s = zeros(N + 1, 1);
    s(2:N) = 1.0 ./ (1.0 - x(2:N).^2);
    S = diag(s);

    % Full fourth-order operator with the polynomial trick
    % Derived from the Leibniz rule for d^4/dx^4 [(1 - x^2) * q(x)]
    D4full = (diag(1.0 - x.^2) * D4 - 8.0 * diag(x) * D3 - 12.0 * D2) * S;

    % Extract interior block (strip first and last rows/columns)
    D4int = D4full(2:N, 2:N);
end

function u = exact_solution(x)
% EXACT_SOLUTION  Exact solution of u^(4) = e^x with clamped BCs.
%
%   The general solution is:
%     u(x) = e^x + c3*x^3 + c2*x^2 + c1*x + c0
%
%   since e^x is a particular solution and {x^3, x^2, x, 1} span the
%   nullspace of d^4/dx^4.  The four constants are determined by the
%   boundary conditions u(+-1) = u'(+-1) = 0.

    e_plus  = exp(1.0);
    e_minus = exp(-1.0);

    % 4x4 linear system: A * [c3; c2; c1; c0] = b
    A = [ 1.0,  1.0,  1.0, 1.0;    % u(1) = 0
         -1.0,  1.0, -1.0, 1.0;    % u(-1) = 0
          3.0,  2.0,  1.0, 0.0;    % u'(1) = 0
          3.0, -2.0,  1.0, 0.0];   % u'(-1) = 0
    b = [-e_plus; -e_minus; -e_plus; -e_minus];

    c = A \ b;
    c3 = c(1); c2 = c(2); c1 = c(3); c0 = c(4);

    u = exp(x) + c3 * x.^3 + c2 * x.^2 + c1 * x + c0;
end
