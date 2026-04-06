% quad_gauss_cc_construction.m
% Chapter 15: Quadrature in Spectral Methods
%
% Computational Etude 15.4: Constructing Gauss and Clenshaw--Curtis Rules.
%
% Implements two fundamental quadrature construction algorithms from scratch
% and validates them against known results.
%
% Gauss--Legendre via the Golub--Welsch Algorithm
% ------------------------------------------------
% The n-point Gauss--Legendre nodes are eigenvalues of the n x n symmetric
% tridiagonal Jacobi matrix J with:
%
%     J(j, j)   = 0                          (diagonal)
%     J(j, j+1) = beta_j = j / sqrt(4j^2 - 1)  (off-diagonal, j = 1, ..., n-1)
%
% The weights are w_j = 2 * (v_j(1))^2, where v_j is the eigenvector
% corresponding to the j-th eigenvalue.
%
% Clenshaw--Curtis via FFT
% -------------------------
% The n+1 Clenshaw--Curtis nodes are x_j = cos(j pi / n), j = 0, ..., n.
% The weights are computed from the Chebyshev moments using the DCT-I.
%
% Both methods are validated by integrating e^x (analytic: e - 1/e) and
% 1/(1+x^2) (analytic: pi/2).
%
% Generates Figure 15.4: gauss_cc_construction.pdf (two-panel figure)
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
% Test functions and their exact integrals over [-1, 1]
% -------------------------------------------------------------------------
f1 = @(x) exp(x);
I1_exact = exp(1) - exp(-1);  % e - e^{-1}

f2 = @(x) 1.0 ./ (1.0 + x.^2);
I2_exact = pi / 2.0;  % 2 * arctan(1) = pi/2

n_values = 2:2:50;
nn = length(n_values);

err_gw_f1  = zeros(nn, 1);
err_cc_f1  = zeros(nn, 1);

err_gw_f2  = zeros(nn, 1);
err_cc_f2  = zeros(nn, 1);

for i = 1:nn
    n = n_values(i);

    % Golub--Welsch Gauss--Legendre (n+1 points)
    [x_gw, w_gw] = gauss_legendre(n + 1);
    err_gw_f1(i) = abs(w_gw' * f1(x_gw) - I1_exact);
    err_gw_f2(i) = abs(w_gw' * f2(x_gw) - I2_exact);

    % Clenshaw--Curtis (n+1 points)
    [x_cc, w_cc] = clenshaw_curtis_weights(n);
    err_cc_f1(i) = abs(w_cc' * f1(x_cc) - I1_exact);
    err_cc_f2(i) = abs(w_cc' * f2(x_cc) - I2_exact);
end

% Clamp zeros
fl = 1e-16;
err_gw_f1  = max(err_gw_f1, fl);
err_cc_f1  = max(err_cc_f1, fl);
err_gw_f2  = max(err_gw_f2, fl);
err_cc_f2  = max(err_cc_f2, fl);

% -------------------------------------------------------------------------
% Create two-panel figure
% -------------------------------------------------------------------------
fig = figure('Position', [100, 100, 1100, 450]);

% --- Panel (a): e^x ---
subplot(1, 2, 1);
semilogy(n_values, err_gw_f1, 'o-', 'Color', NAVY, 'MarkerSize', 4, ...
    'LineWidth', 1.2, 'MarkerFaceColor', NAVY, ...
    'DisplayName', 'Golub--Welsch (Gauss)');
hold on;
semilogy(n_values, err_cc_f1, 'd-', 'Color', TEAL, 'MarkerSize', 4, ...
    'LineWidth', 1.2, 'MarkerFaceColor', TEAL, ...
    'DisplayName', 'CC (FFT)');
hold off;

xlabel('Number of points $n$', 'Interpreter', 'latex');
ylabel('Absolute error $|I - I_n|$', 'Interpreter', 'latex');
title('(a) $f(x) = e^x$', 'Interpreter', 'latex');
legend('Location', 'best', 'EdgeColor', [0.5 0.5 0.5], 'FontSize', 9);
grid on; set(gca, 'GridAlpha', 0.3);
xlim([n_values(1) - 1, n_values(end) + 1]);

% --- Panel (b): 1/(1+x^2) ---
subplot(1, 2, 2);
semilogy(n_values, err_gw_f2, 'o-', 'Color', NAVY, 'MarkerSize', 4, ...
    'LineWidth', 1.2, 'MarkerFaceColor', NAVY, ...
    'DisplayName', 'Golub--Welsch (Gauss)');
hold on;
semilogy(n_values, err_cc_f2, 'd-', 'Color', TEAL, 'MarkerSize', 4, ...
    'LineWidth', 1.2, 'MarkerFaceColor', TEAL, ...
    'DisplayName', 'CC (FFT)');
hold off;

xlabel('Number of points $n$', 'Interpreter', 'latex');
ylabel('Absolute error $|I - I_n|$', 'Interpreter', 'latex');
title('(b) $f(x) = 1/(1 + x^2)$', 'Interpreter', 'latex');
legend('Location', 'best', 'EdgeColor', [0.5 0.5 0.5], 'FontSize', 9);
grid on; set(gca, 'GridAlpha', 0.3);
xlim([n_values(1) - 1, n_values(end) + 1]);

% Save output
exportgraphics(gcf, fullfile(output_dir, 'gauss_cc_construction.pdf'), 'ContentType', 'vector');
exportgraphics(gcf, fullfile(output_dir, 'gauss_cc_construction.png'), 'Resolution', 300);
fprintf('Figure saved to %s\n', fullfile(output_dir, 'gauss_cc_construction.pdf'));

% =========================================================================
% Local functions
% =========================================================================

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
