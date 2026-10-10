% quad_aliasing_chebyshev.m
% Chapter 15: Quadrature in Spectral Methods
%
% Computational Etude 15.6: Aliasing and Chebyshev Quadrature.
%
% Illustrates the aliasing phenomenon in Clenshaw--Curtis quadrature by
% examining the quadrature error E_n(T_k) = |I(T_k) - I_n(T_k)| for
% individual Chebyshev polynomials T_k, and its impact on the integration
% of |x|^3.
%
% Left panel: For n = 30, the CC quadrature error |E_n(T_k)| is plotted
% against even k from 0 to 200. Three regimes are visible:
%     (1) k <= n: errors are at machine zero (exact integration),
%     (2) n < k <= 2n: errors are small, O(n^{-2}),
%     (3) k > 2n: aliasing errors of O(1) at multiples of 2n.
%
% Right panel: For f(x) = |x|^3, the Chebyshev expansion coefficients
% |a_k| are plotted alongside the aliased error contributions |a_k * E_n(T_k)|.
%
% Reference: L. N. Trefethen, "Is Gauss Quadrature Better than
%            Clenshaw--Curtis?", SIAM Review 50(1):67--87, 2008, Figure 3.
%
% Generates Figure 15.6: aliasing_chebyshev.pdf
%
% Author: Dr. Denys Dutykh
%         Mathematics Department
%         Khalifa University of Science and Technology
%         Abu Dhabi, UAE
%
% Part of the book "Computational Etudes: A Spectral Approach"
% https://github.com/dutykh/computational-etudes

clear; close all; clc;

% Publication-quality figure settings
set(0, 'DefaultAxesFontSize', 11);
set(0, 'DefaultTextFontSize', 11);
set(0, 'DefaultLineLineWidth', 1.5);
set(0, 'DefaultAxesLineWidth', 0.8);

% Book color scheme
NAVY   = [20, 45, 110] / 255;
SKY    = [120, 150, 210] / 255;
CORAL  = [231, 76, 60] / 255;
TEAL   = [22, 160, 133] / 255;
PURPLE = [142, 68, 173] / 255;
ORANGE = [230, 126, 34] / 255;

% Output path
output_dir = fullfile(fileparts(mfilename('fullpath')), ...
    '..', '..', '..', 'textbook', 'figures', 'ch15', 'matlab');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

% -------------------------------------------------------------------------
% Parameters
% -------------------------------------------------------------------------
n = 30;       % CC quadrature with n+1 = 31 points
k_max = 200;  % Maximum Chebyshev degree to examine

% Compute CC nodes and weights
[x_cc, w_cc] = clenshaw_curtis_weights(n);

% -------------------------------------------------------------------------
% Left panel: |E_n(T_k)| for even k = 0, 2, 4, ..., k_max
% -------------------------------------------------------------------------
k_even = 0:2:k_max;
errors = zeros(size(k_even));

for i = 1:length(k_even)
    k = k_even(i);

    % Evaluate T_k at CC nodes: T_k(x) = cos(k * arccos(x))
    Tk_at_nodes = cos(k * acos(x_cc));

    % Numerical integral via CC quadrature
    I_n = w_cc' * Tk_at_nodes;

    % Exact integral
    I_exact = exact_integral_Tk(k);

    errors(i) = abs(I_n - I_exact);
end

% -------------------------------------------------------------------------
% Right panel: Chebyshev coefficients of |x|^3 and aliased errors
% -------------------------------------------------------------------------
f_abs3 = @(x) abs(x).^3;

% Compute many Chebyshev coefficients using a large number of points
N_coeffs = 512;
a_k = chebyshev_coefficients(f_abs3, N_coeffs);

% Only even coefficients are nonzero for this even function
k_coeffs = 0:2:min(k_max, N_coeffs - 1);
a_even = abs(a_k(k_coeffs + 1));  % +1 for 1-based indexing

% Compute aliased error contributions: |a_k * E_n(T_k)|
n_common = min(length(k_coeffs), length(k_even));
k_common = k_even(1:n_common);
a_common = zeros(1, n_common);
e_common = zeros(1, n_common);

for i = 1:n_common
    k = k_common(i);
    if k < N_coeffs
        a_common(i) = abs(a_k(k + 1));  % +1 for 1-based indexing
    end
    % Compute E_n(T_k)
    Tk_at_nodes = cos(k * acos(x_cc));
    I_n = w_cc' * Tk_at_nodes;
    I_exact = exact_integral_Tk(k);
    e_common(i) = abs(I_n - I_exact);
end

aliased_product = a_common .* e_common;

% -------------------------------------------------------------------------
% Create figure
% -------------------------------------------------------------------------
fl = 1e-17;

fig = figure('Position', [100, 100, 1100, 450]);

% --- Left panel: CC quadrature errors for T_k ---
subplot(1, 2, 1);

errors_plot = max(errors, fl);
semilogy(k_even, errors_plot, 'o', 'Color', NAVY, 'MarkerSize', 3, ...
    'MarkerEdgeColor', NAVY, 'MarkerFaceColor', NAVY);
hold on;

% Shade the region n <= k <= 2n
patch([n, 2*n, 2*n, n], [1e-17, 1e-17, 1e1, 1e1], SKY, ...
    'FaceAlpha', 0.12, 'EdgeColor', 'none', ...
    'DisplayName', '$n \leq k \leq 2n$');

% Vertical dashed lines at k = n and k = 2n
xline(n, '--', 'Color', CORAL, 'LineWidth', 0.8, 'Alpha', 0.7);
xline(2*n, '--', 'Color', CORAL, 'LineWidth', 0.8, 'Alpha', 0.7);

% Annotate the regions
text(n/2, 2e-16, 'exact', 'HorizontalAlignment', 'center', ...
    'FontSize', 9, 'Color', TEAL, 'FontAngle', 'italic');
text(1.5*n, 1e-3, {'small', 'errors'}, 'HorizontalAlignment', 'center', ...
    'FontSize', 9, 'Color', NAVY, 'FontAngle', 'italic');
text(2.5*n + 20, 0.5, 'aliasing', 'HorizontalAlignment', 'center', ...
    'FontSize', 9, 'Color', CORAL, 'FontAngle', 'italic');

hold off;

xlabel('Even degree $k$', 'Interpreter', 'latex');
ylabel('$|E_n(T_k)|$', 'Interpreter', 'latex');
title('(a) CC error for Chebyshev polynomials ($n = 30$)', ...
    'Interpreter', 'latex');
xlim([-2, k_max + 5]);
ylim([1e-17, 1e1]);
legend('', '$n \leq k \leq 2n$', 'Location', 'southeast', ...
    'Interpreter', 'latex', 'EdgeColor', [0.5 0.5 0.5], 'FontSize', 9);
grid on; set(gca, 'GridAlpha', 0.3);

% --- Right panel: |a_k| and |a_k * E_n(T_k)| for |x|^3 ---
subplot(1, 2, 2);

a_plot = max(a_even, fl);
alias_plot = max(aliased_product(1:length(k_coeffs)), fl);

semilogy(k_coeffs, a_plot, 'o', 'Color', NAVY, 'MarkerSize', 3, ...
    'MarkerEdgeColor', NAVY, 'MarkerFaceColor', NAVY, ...
    'DisplayName', '$|a_k|$ (Chebyshev coefficients)');
hold on;
semilogy(k_common(1:length(k_coeffs)), alias_plot, 's', 'Color', CORAL, ...
    'MarkerSize', 3, 'MarkerFaceColor', 'w', ...
    'MarkerEdgeColor', CORAL, 'LineWidth', 0.7, ...
    'DisplayName', '$|a_k \cdot E_n(T_k)|$ (aliased error)');

% Shade the region n <= k <= 2n
patch([n, 2*n, 2*n, n], [1e-17, 1e-17, 1e1, 1e1], SKY, ...
    'FaceAlpha', 0.12, 'EdgeColor', 'none', 'HandleVisibility', 'off');

hold off;

xlabel('Even degree $k$', 'Interpreter', 'latex');
ylabel('Magnitude', 'Interpreter', 'latex');
title('(b) Aliasing for $f(x) = |x|^3$', 'Interpreter', 'latex');
xlim([-2, k_max + 5]);
ylim([1e-17, 1e1]);
legend('Location', 'northeast', 'Interpreter', 'latex', ...
    'EdgeColor', [0.5 0.5 0.5], 'FontSize', 9);
grid on; set(gca, 'GridAlpha', 0.3);

% Save output
exportgraphics(gcf, fullfile(output_dir, 'aliasing_chebyshev.pdf'), 'ContentType', 'vector');
exportgraphics(gcf, fullfile(output_dir, 'aliasing_chebyshev.png'), 'Resolution', 300);
fprintf('Figure saved to %s\n', fullfile(output_dir, 'aliasing_chebyshev.pdf'));

% =========================================================================
% Local functions
% =========================================================================

function I = exact_integral_Tk(k)
    % Exact integral of T_k(x) over [-1, 1].
    % I(T_k) = 2 / (1 - k^2) for even k, 0 for odd k.
    if mod(k, 2) == 1
        I = 0.0;
    else
        I = 2.0 / (1.0 - k^2);
    end
end

function a = chebyshev_coefficients(f, N)
    % Compute Chebyshev expansion coefficients a_k of f(x) using
    % the DCT via evaluation at Chebyshev points.
    j = (0:N-1)';
    x = cos(pi * j / (N - 1));  % Chebyshev extremal points
    fx = f(x);

    % DCT-I via FFT
    v = [fx; fx(end-1:-1:2)];
    coeffs = real(fft(v));
    a = coeffs(1:N) / (N - 1);
    a(1) = a(1) / 2;
    a(N) = a(N) / 2;
end

function [x, w] = clenshaw_curtis_weights(n)
    % Clenshaw-Curtis nodes and weights via DCT-I / FFT.
    if n == 0
        x = 0; w = 2; return;
    end
    if n == 1
        x = [1; -1]; w = [1; 1]; return;
    end
    x = cos(pi * (0:n)' / n);
    c = zeros(n+1, 1);
    c(1:2:end) = 2 ./ (1 - (0:2:n)'.^2);
    v = [c; c(n:-1:2)];  % length 2n, mirror
    f = real(fft(v));
    w = f(1:n+1) / n;
    w(1) = w(1) / 2;
    w(end) = w(end) / 2;
end
