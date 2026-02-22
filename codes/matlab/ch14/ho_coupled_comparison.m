%% ho_coupled_comparison.m - Direct D^4 vs coupled D^2 system comparison
%
% Chapter 14: Higher-Order Boundary Value Problems
%
% Computational Etude 14.3: Direct vs. coupled system comparison.
%
% Solves the clamped beam problem:
%
%   u^(4)(x) = e^x,   -1 < x < 1,   u(+-1) = u'(+-1) = 0
%
% using TWO approaches and compares their accuracy and conditioning.
%
% Approach 1 (Direct D^4):
%   The polynomial trick from Trefethen's "Spectral Methods in MATLAB"
%   (Program 38). Write u(x) = (1 - x^2) * q(x), which automatically
%   satisfies u(+-1) = 0 and (via the product rule) u'(+-1) = 0.
%   The fourth-order operator is reformulated as:
%     D4full = [diag(1 - x^2) * D^4 - 8*diag(x) * D^3 - 12*D^2] * S
%   and the interior system D4int * v = f is solved directly.
%
% Approach 2 (Coupled second-order system):
%   Introduce w = u'' to split the fourth-order ODE into two coupled
%   second-order equations:
%     u'' = w,      with u(+-1) = 0
%     w'' = e^x,    with u'(+-1) = 0
%   Assembled as a 2(N+1) x 2(N+1) block system.
%
% The direct approach yields a condition number scaling as O(N^8),
% while the coupled system, involving only D^2, scales as O(N^4).
%
% Generates Figure: coupled_comparison.pdf
%   Left panel:  convergence comparison (max error vs N)
%   Right panel: condition number comparison (loglog)
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

%% Convergence and conditioning study
N_values = 6:2:40;
errors_direct  = zeros(size(N_values));
errors_coupled = zeros(size(N_values));
cond_direct    = zeros(size(N_values));
cond_coupled   = zeros(size(N_values));

for k = 1:length(N_values)
    N_val = N_values(k);

    % --- Direct D^4 approach ---
    [u_d, x_d, D4int] = solve_direct_d4(N_val);
    u_exact_d = exact_solution(x_d);
    err_d = max(abs(u_d - u_exact_d));
    errors_direct(k) = max(err_d, 1e-16);
    cond_direct(k) = cond(D4int);

    % --- Coupled D^2 approach ---
    [u_c, x_c, A_block] = solve_coupled_d2(N_val);
    u_exact_c = exact_solution(x_c);
    err_c = max(abs(u_c - u_exact_c));
    errors_coupled(k) = max(err_c, 1e-16);
    cond_coupled(k) = cond(A_block);
end

%% Print convergence table
fprintf('%4s  %13s  %13s  %14s  %14s\n', 'N', 'Err (Direct)', ...
        'Err (Coupled)', 'cond(D4int)', 'cond(Block)');
fprintf('%s\n', repmat('-', 1, 66));
for k = 1:length(N_values)
    fprintf('%4d  %13.4e  %13.4e  %14.4e  %14.4e\n', N_values(k), ...
            errors_direct(k), errors_coupled(k), ...
            cond_direct(k), cond_coupled(k));
end
fprintf('%s\n', repmat('-', 1, 66));

%% Compute fitted slopes
N_arr = double(N_values);
n_fit = floor(length(N_values) / 2);
log_N = log(N_arr(n_fit:end));

log_cond_d = log(cond_direct(n_fit:end));
p_d = polyfit(log_N, log_cond_d, 1);
slope_d = p_d(1);

log_cond_c = log(cond_coupled(n_fit:end));
p_c = polyfit(log_N, log_cond_c, 1);
slope_c = p_c(1);

fprintf('\nFitted condition number slopes (log-log):\n');
fprintf('  Direct D^4:  %.2f  (expected ~8)\n', slope_d);
fprintf('  Coupled D^2: %.2f  (expected ~4)\n', slope_c);

%% Create figure
fig = figure('Units', 'inches', 'Position', [1, 1, 11, 4.5]);

% ---- Left panel: max error vs N (semilog) ----
subplot(1, 2, 1);
semilogy(N_values, errors_direct, 'o-', 'Color', CORAL, 'LineWidth', 1.5, ...
         'MarkerSize', 5, 'MarkerFaceColor', CORAL, 'MarkerEdgeColor', 'w', ...
         'DisplayName', 'Direct $D^4$ (polynomial trick)');
hold on;
semilogy(N_values, errors_coupled, 's-', 'Color', TEAL, 'LineWidth', 1.5, ...
         'MarkerSize', 5, 'MarkerFaceColor', TEAL, 'MarkerEdgeColor', 'w', ...
         'DisplayName', 'Coupled $D^2$ system');
yline(2.2e-16, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
text(43, 4e-16, '$\epsilon_{\mathrm{mach}}$', 'FontSize', 9, ...
     'Color', [0.5 0.5 0.5], 'HorizontalAlignment', 'right', ...
     'VerticalAlignment', 'bottom');
hold off;

xlabel('$N$');
ylabel('Max error');
title('Convergence comparison', 'FontSize', 11);
legend('Location', 'northeast', 'FontSize', 9);
xlim([4, 44]);
ylim([1e-16, 1e2]);
grid on;
set(gca, 'GridAlpha', 0.3);
box on;

% ---- Right panel: condition number vs N (loglog) ----
subplot(1, 2, 2);
loglog(N_values, cond_direct, 'o-', 'Color', CORAL, 'LineWidth', 1.5, ...
       'MarkerSize', 5, 'MarkerFaceColor', CORAL, 'MarkerEdgeColor', 'w', ...
       'DisplayName', 'Direct $D^4$: $\mathcal{O}(N^8)$');
hold on;
loglog(N_values, cond_coupled, 's-', 'Color', TEAL, 'LineWidth', 1.5, ...
       'MarkerSize', 5, 'MarkerFaceColor', TEAL, 'MarkerEdgeColor', 'w', ...
       'DisplayName', 'Coupled $D^2$: $\mathcal{O}(N^4)$');

% Reference lines (visual guides)
N_ref = [10, 40];
ref_N_mid = N_arr(floor(length(N_arr) / 2));

% Direct reference: N^8
cond_d_mid = cond_direct(floor(length(N_values) / 2));
loglog(N_ref, cond_d_mid * (N_ref / ref_N_mid).^8, '--', ...
       'Color', CORAL, 'LineWidth', 0.8, 'HandleVisibility', 'off');

% Coupled reference: N^4
cond_c_mid = cond_coupled(floor(length(N_values) / 2));
loglog(N_ref, cond_c_mid * (N_ref / ref_N_mid).^4, '--', ...
       'Color', TEAL, 'LineWidth', 0.8, 'HandleVisibility', 'off');
hold off;

xlabel('$N$');
ylabel('Condition number');
title('Conditioning comparison', 'FontSize', 11);
legend('Location', 'northwest', 'FontSize', 9);
grid on;
set(gca, 'GridAlpha', 0.3);
box on;

sgtitle('Computational Etude 14.3: Direct $D^4$ vs.\ Coupled $D^2$ System', ...
        'FontSize', 12);

%% Save figure
exportgraphics(fig, fullfile(output_dir, 'coupled_comparison.pdf'), ...
               'ContentType', 'vector');
fprintf('\nFigure saved to: %s\n', fullfile(output_dir, 'coupled_comparison.pdf'));

close(fig);


%% =====================================================================
%  Local functions
%  =====================================================================

function u = exact_solution(x)
% EXACT_SOLUTION  Exact solution of u^(4) = e^x with clamped BCs.
%
%   u(x) = e^x + c3*x^3 + c2*x^2 + c1*x + c0
%   with constants determined by u(+-1) = u'(+-1) = 0.

    e_plus  = exp(1.0);
    e_minus = exp(-1.0);

    A = [ 1.0,  1.0,  1.0, 1.0;
         -1.0,  1.0, -1.0, 1.0;
          3.0,  2.0,  1.0, 0.0;
          3.0, -2.0,  1.0, 0.0];
    b = [-e_plus; -e_minus; -e_plus; -e_minus];

    c = A \ b;
    c3 = c(1); c2 = c(2); c1 = c(3); c0 = c(4);

    u = exp(x) + c3 * x.^3 + c2 * x.^2 + c1 * x + c0;
end

function [u, x, D4int] = solve_direct_d4(N)
% SOLVE_DIRECT_D4  Solve u^(4) = e^x with clamped BCs via polynomial trick.
%
%   Writes u(x) = (1 - x^2) * q(x) so that all four BCs are built in.
%   The operator is reformulated using the Leibniz rule and restricted
%   to interior points.

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
    D4full = (diag(1.0 - x.^2) * D4 - 8.0 * diag(x) * D3 - 12.0 * D2) * S;

    % Extract interior block
    D4int = D4full(2:N, 2:N);

    % Right-hand side at interior points
    f = exp(x(2:N));

    % Solve
    v = D4int \ f;

    % Reconstruct full solution (boundary values are zero)
    u = zeros(N + 1, 1);
    u(2:N) = v;
end

function [u, x, A_full] = solve_coupled_d2(N)
% SOLVE_COUPLED_D2  Solve u^(4) = e^x via coupled second-order system.
%
%   Introduce w = u'' to obtain:
%     u'' = w        (with u(+-1) = 0)
%     w'' = e^x      (with u'(+-1) = 0 encoded via D)
%
%   The full 2(N+1) x 2(N+1) block system is:
%     [ A11  A12 ] [ u ]   [ b1 ]
%     [ A21  A22 ] [ w ] = [ b2 ]

    [D, x] = cheb_matrix(N);
    D2 = D * D;
    n = N + 1;
    I_mat = eye(n);
    Z = zeros(n, n);

    % Block matrices
    A11 = D2;           % u'' part
    A12 = -I_mat;       % -w part  (u'' - w = 0)
    A21 = Z;            % coupling for w equation
    A22 = D2;           % w'' part

    % Right-hand sides
    b1 = zeros(n, 1);
    b2 = exp(x);

    % ----- BCs for the u-block (Dirichlet: u(+-1) = 0) -----
    % Row 1 (x(1) = 1): u(1) = 0
    A11(1, :) = 0.0;
    A11(1, 1) = 1.0;
    A12(1, :) = 0.0;
    b1(1) = 0.0;

    % Row N+1 (x(N+1) = -1): u(-1) = 0
    A11(N + 1, :) = 0.0;
    A11(N + 1, N + 1) = 1.0;
    A12(N + 1, :) = 0.0;
    b1(N + 1) = 0.0;

    % ----- BCs for the w-block (Neumann: u'(+-1) = 0) -----
    % Row 1 of w-block: u'(x(1)) = u'(1) = D(1,:) * u = 0
    A21(1, :) = D(1, :);
    A22(1, :) = 0.0;
    b2(1) = 0.0;

    % Row N+1 of w-block: u'(x(N+1)) = u'(-1) = D(N+1,:) * u = 0
    A21(N + 1, :) = D(N + 1, :);
    A22(N + 1, :) = 0.0;
    b2(N + 1) = 0.0;

    % Assemble the full block system
    A_full = [A11, A12; A21, A22];
    rhs = [b1; b2];

    % Solve
    sol = A_full \ rhs;
    u = sol(1:n);
end
