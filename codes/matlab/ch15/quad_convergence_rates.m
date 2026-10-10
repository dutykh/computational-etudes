% quad_convergence_rates.m
% Chapter 15: Quadrature in Spectral Methods
%
% Computational Etude 15.10: Asymptotic Convergence Rate Verification.
%
% Verifies the theoretical asymptotic convergence rates for two quadrature
% strategies applied to the integral
%
%     I = int_{-inf}^{inf} e^{-x^2} cos(x^3) dx
%
% using the rescaling technique of Trefethen and Weideman (2014).
%
% Left panel: Gauss--Hermite errors plotted against n^{1/2}. Theory
% predicts O(exp(-C * n^{1/2})) convergence. On a semilog plot with
% horizontal axis sqrt(n), the data should fall on a straight line.
%
% Right panel: Truncated Gauss--Legendre errors plotted against n^{2/3}.
% For a Gauss--Legendre rule on [-L * n^{1/3}, L * n^{1/3}], theory
% predicts O(exp(-C * n^{2/3})) convergence.
%
% Reference: L. N. Trefethen and J. A. C. Weideman, "The exponentially
%            convergent trapezoidal rule", SIAM Review 56(3):385--458, 2014.
%
% Generates Figure 15.10: convergence_rates.pdf
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
% Reference value via high-accuracy adaptive quadrature
% -------------------------------------------------------------------------
integrand = @(x) exp(-x.^2) .* cos(x.^3);
I_exact = integral(integrand, -Inf, Inf, 'AbsTol', 1e-15, 'RelTol', 1e-15);
fprintf('Reference integral value: %.16e\n', I_exact);

% -------------------------------------------------------------------------
% Gauss--Hermite errors for n = 10, 15, ..., 500
% -------------------------------------------------------------------------
n_gh = 10:5:500;
err_gh = zeros(size(n_gh));

for i = 1:length(n_gh)
    n = n_gh(i);
    [x_h, w_h] = gauss_hermite(n);
    I_gh = w_h' * cos(x_h.^3);
    err_gh(i) = abs(I_gh - I_exact);
end

% -------------------------------------------------------------------------
% Truncated Gauss--Legendre on [-L * n^{1/3}, L * n^{1/3}]
% -------------------------------------------------------------------------
L = 3.0;  % Scaling parameter
n_gl = 10:5:500;
err_gl = zeros(size(n_gl));

for i = 1:length(n_gl)
    n = n_gl(i);
    % Truncation half-width grows as n^{1/3}
    half_width = L * n^(1.0 / 3.0);

    % Gauss--Legendre nodes on [-1, 1], mapped to [-half_width, half_width]
    [s, ws] = gauss_legendre(n);
    x_mapped = half_width * s;
    w_mapped = half_width * ws;

    % Evaluate the full integrand
    f_vals = exp(-x_mapped.^2) .* cos(x_mapped.^3);
    I_gl_val = w_mapped' * f_vals;
    err_gl(i) = abs(I_gl_val - I_exact);
end

% -------------------------------------------------------------------------
% Prepare data for plotting
% -------------------------------------------------------------------------
fl = 1e-16;

% Gauss--Hermite: plot vs sqrt(n)
sqrt_n_gh = sqrt(n_gh);
log_err_gh = log10(max(err_gh, fl));

% Truncated GL: plot vs n^{2/3}
n23_gl = n_gl.^(2.0 / 3.0);
log_err_gl = log10(max(err_gl, fl));

% Linear fits on the valid (above-floor) data
mask_gh = err_gh > fl;
if sum(mask_gh) > 2
    p_gh = polyfit(sqrt_n_gh(mask_gh), log_err_gh(mask_gh), 1);
    fit_gh = polyval(p_gh, sqrt_n_gh);
    slope_gh = p_gh(1);
else
    fit_gh = [];
    slope_gh = NaN;
end

mask_gl = err_gl > fl;
if sum(mask_gl) > 2
    p_gl = polyfit(n23_gl(mask_gl), log_err_gl(mask_gl), 1);
    fit_gl = polyval(p_gl, n23_gl);
    slope_gl = p_gl(1);
else
    fit_gl = [];
    slope_gl = NaN;
end

% -------------------------------------------------------------------------
% Create figure
% -------------------------------------------------------------------------
fig = figure('Position', [100, 100, 1100, 450]);

% --- Left panel: Gauss--Hermite vs sqrt(n) ---
subplot(1, 2, 1);

semilogy(sqrt_n_gh, max(err_gh, fl), 'o', 'Color', NAVY, ...
    'MarkerSize', 3, 'MarkerEdgeColor', NAVY, 'MarkerFaceColor', NAVY, ...
    'DisplayName', 'Gauss-Hermite error');
hold on;

if ~isempty(fit_gh)
    semilogy(sqrt_n_gh, 10.^fit_gh, '--', 'Color', CORAL, 'LineWidth', 1.2, ...
        'DisplayName', sprintf('Linear fit (slope = %.2f)', slope_gh));
end

yline(eps, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8, ...
    'HandleVisibility', 'off');

hold off;

xlabel('$n^{1/2}$', 'Interpreter', 'latex');
ylabel('Absolute error $|I - I_n|$', 'Interpreter', 'latex');
title('(a) Gauss-Hermite: $O(e^{-C\sqrt{n}})$', 'Interpreter', 'latex');
legend('Location', 'northeast', 'EdgeColor', [0.5 0.5 0.5], 'FontSize', 9);
grid on; set(gca, 'GridAlpha', 0.3);
ylim([1e-16, 1e1]);

% --- Right panel: Truncated GL vs n^{2/3} ---
subplot(1, 2, 2);

semilogy(n23_gl, max(err_gl, fl), 'o', 'Color', TEAL, ...
    'MarkerSize', 3, 'MarkerEdgeColor', TEAL, 'MarkerFaceColor', TEAL, ...
    'DisplayName', 'Truncated GL error');
hold on;

if ~isempty(fit_gl)
    semilogy(n23_gl, 10.^fit_gl, '--', 'Color', CORAL, 'LineWidth', 1.2, ...
        'DisplayName', sprintf('Linear fit (slope = %.2f)', slope_gl));
end

yline(eps, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8, ...
    'HandleVisibility', 'off');

hold off;

xlabel('$n^{2/3}$', 'Interpreter', 'latex');
ylabel('Absolute error $|I - I_n|$', 'Interpreter', 'latex');
title('(b) Truncated GL: $O(e^{-C\, n^{2/3}})$', 'Interpreter', 'latex');
legend('Location', 'northeast', 'EdgeColor', [0.5 0.5 0.5], 'FontSize', 9);
grid on; set(gca, 'GridAlpha', 0.3);
ylim([1e-16, 1e1]);

% Save output
exportgraphics(gcf, fullfile(output_dir, 'convergence_rates.pdf'), 'ContentType', 'vector');
exportgraphics(gcf, fullfile(output_dir, 'convergence_rates.png'), 'Resolution', 300);
fprintf('Figure saved to %s\n', fullfile(output_dir, 'convergence_rates.pdf'));

% =========================================================================
% Local functions
% =========================================================================

function [x, w] = gauss_hermite(n)
    % Gauss-Hermite nodes and weights via Golub-Welsch.
    j = (1:n-1)';
    beta = sqrt(j / 2);
    J = diag(beta, 1) + diag(beta, -1);
    [V, D] = eig(J);
    [x, idx] = sort(diag(D));
    w = sqrt(pi) * V(1, idx)'.^2;
end

function [x, w] = gauss_legendre(n)
    % n-point Gauss-Legendre nodes and weights via Golub-Welsch.
    if n == 1
        x = 0; w = 2; return;
    end
    j = (1:n-1)';
    beta = j ./ sqrt(4*j.^2 - 1);
    J = diag(beta, 1) + diag(beta, -1);
    [V, D] = eig(J);
    [x, idx] = sort(diag(D));
    w = 2 * V(1, idx)'.^2;
end
