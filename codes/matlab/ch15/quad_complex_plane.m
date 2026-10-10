% quad_complex_plane.m
% Chapter 15: Quadrature in Spectral Methods
%
% Computational Etude 15.7: Quadrature in the Complex Plane.
%
% Reproduces the essential content of Figure 4 from Trefethen (2008),
% comparing Newton--Cotes, Gauss--Legendre, and Clenshaw--Curtis
% quadrature rules through rational approximation in the complex plane.
%
% A quadrature rule with nodes x_k and weights w_k defines a rational
% approximant r_n(z) = sum_{k} w_k / (z - x_k) to the function
% phi(z) = log((z + 1) / (z - 1)).
%
% The contour plot of |phi(z) - r_n(z)| reveals the accuracy of each
% quadrature rule throughout the complex plane.
%
% The figure shows a 1 x 3 panel layout for n = 16 (17 points):
%     (a) Newton--Cotes
%     (b) Gauss--Legendre
%     (c) Clenshaw--Curtis
%
% Reference: L. N. Trefethen, "Is Gauss Quadrature Better than
%            Clenshaw--Curtis?", SIAM Review 50(1):67--87, 2008, Figure 4.
%
% Generates Figure 15.7: complex_plane.pdf
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
n = 16;  % 17-point rules

% Compute nodes and weights for all three rules
[x_nc, w_nc] = newton_cotes_weights(n);
[x_gl, w_gl] = gauss_legendre(n + 1);
[x_cc, w_cc] = clenshaw_curtis_weights(n);

rules_x = {x_nc, x_gl, x_cc};
rules_w = {w_nc, w_gl, w_cc};
rule_titles = {'(a) Newton-Cotes', '(b) Gauss-Legendre', '(c) Clenshaw-Curtis'};

% Complex grid
nx = 400; ny = 300;
x_re = linspace(-2.5, 2.5, nx);
y_im = linspace(-2.0, 2.0, ny);
[X, Y] = meshgrid(x_re, y_im);
Z = X + 1i * Y;

% Mask points too close to the branch cut [-1, 1] on the real axis
mask_distance = 0.03;
on_cut = (abs(Y) < mask_distance) & (X >= -1.0) & (X <= 1.0);

% Contour levels: 10^{-13}, 10^{-12}, ..., 10^0
levels = logspace(-13, 0, 14);

% Compute phi(z) = log((z+1)/(z-1)) once
phi_vals = log((Z + 1.0) ./ (Z - 1.0));

% -------------------------------------------------------------------------
% NEW (panel d): convergence on a benchmark integral
%   f(x) = 1/(1 + 16 x^2),  exact integral = (1/2) atan(4) ≈ 0.66291...
% -------------------------------------------------------------------------
f_bench  = @(x) 1.0 ./ (1.0 + 16.0 .* x.^2);
I_exact  = 0.5 * atan(4.0);
ns_conv  = [4 6 8 12 16 24 32 48 64 96 128];
err_nc   = zeros(size(ns_conv));
err_gl   = zeros(size(ns_conv));
err_cc   = zeros(size(ns_conv));
for k = 1:numel(ns_conv)
    nv = ns_conv(k);
    [xn, wn] = newton_cotes_weights(nv);
    [xg, wg] = gauss_legendre(nv + 1);
    [xc, wc] = clenshaw_curtis_weights(nv);
    err_nc(k) = abs(sum(wn .* f_bench(xn)) - I_exact);
    err_gl(k) = abs(sum(wg .* f_bench(xg)) - I_exact);
    err_cc(k) = abs(sum(wc .* f_bench(xc)) - I_exact);
end

% -------------------------------------------------------------------------
% Create 2x2 figure
% -------------------------------------------------------------------------
fig = figure('Position', [50, 100, 1100, 900]);

for p = 1:3
    if p == 1
        subplot(2, 2, 1);
    elseif p == 2
        subplot(2, 2, 2);
    else
        subplot(2, 2, 3);
    end

    x_q = rules_x{p};
    w_q = rules_w{p};

    % Compute rational approximation r_n(z) = sum_k w_k / (z - x_k)
    r_vals = zeros(size(Z));
    for k = 1:length(x_q)
        r_vals = r_vals + w_q(k) ./ (Z - x_q(k));
    end

    % Approximation error
    err = abs(phi_vals - r_vals);

    % Mask branch cut region
    err(on_cut) = NaN;

    % Contour plot
    contour(X, Y, err, levels, 'LineColor', 'k', 'LineWidth', 0.5);
    hold on;

    % Mark the interval [-1, 1]
    plot([-1, 1], [0, 0], '-', 'Color', NAVY, 'LineWidth', 2.0);

    % Mark quadrature nodes
    plot(x_q, zeros(size(x_q)), 'o', 'Color', CORAL, 'MarkerSize', 3, ...
        'MarkerFaceColor', CORAL, 'MarkerEdgeColor', CORAL);

    hold off;

    title(rule_titles{p}, 'FontSize', 11);
    xlabel('Re$(z)$', 'Interpreter', 'latex');
    ylabel('Im$(z)$', 'Interpreter', 'latex');
    xlim([-2.5, 2.5]);
    ylim([-2.0, 2.0]);
    axis equal;
    xlim([-2.5, 2.5]);
    ylim([-2.0, 2.0]);
    grid on; set(gca, 'GridAlpha', 0.15);
end

% Panel (d): NEW convergence on benchmark integral
subplot(2, 2, 4);
semilogy(ns_conv, err_nc + 1e-18, '-o', 'Color', CORAL, 'LineWidth', 1.0, ...
         'MarkerFaceColor', 'w', 'MarkerSize', 5, ...
         'DisplayName', '(a) Newton-Cotes (Runge: diverges)'); hold on;
semilogy(ns_conv, err_gl + 1e-18, '-s', 'Color', NAVY, 'LineWidth', 1.0, ...
         'MarkerFaceColor', 'w', 'MarkerSize', 5, ...
         'DisplayName', '(b) Gauss-Legendre (geometric)');
semilogy(ns_conv, err_cc + 1e-18, '-^', 'Color', TEAL, 'LineWidth', 1.0, ...
         'MarkerFaceColor', 'w', 'MarkerSize', 5, ...
         'DisplayName', '(c) Clenshaw-Curtis (geometric)');
hold off;
xlabel('n (rule order, n+1 points)');
ylabel('|I_n - I|');
title('(d) convergence on int_{-1}^1 dx / (1 + 16 x^2)', 'FontSize', 11);
grid on; box on; set(gca, 'GridAlpha', 0.3);
legend('Location', 'northeast', 'FontSize', 8);

% Save output
exportgraphics(gcf, fullfile(output_dir, 'complex_plane.pdf'), 'ContentType', 'vector');
exportgraphics(gcf, fullfile(output_dir, 'complex_plane.png'), 'Resolution', 300);
fprintf('Figure saved to %s\n', fullfile(output_dir, 'complex_plane.pdf'));
fprintf('  exact integral = %.10f\n', I_exact);
fprintf('  NC at n=64: error = %.3e\n', err_nc(9));
fprintf('  GL at n=64: error = %.3e\n', err_gl(9));
fprintf('  CC at n=64: error = %.3e\n', err_cc(9));

% =========================================================================
% Local functions
% =========================================================================

function [x, w] = newton_cotes_weights(n)
    x = linspace(-1, 1, n+1)';
    V = zeros(n+1);
    for k = 0:n
        V(k+1,:) = x.^k;
    end
    rhs = arrayfun(@(k) (1 - (-1)^(k+1))/(k+1), (0:n)');
    w = V \ rhs;
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
