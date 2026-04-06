% quad_convergence_race.m
% Chapter 15: Quadrature in Spectral Methods
%
% Computational Etude 15.5: The Convergence Race.
%
% Reproduces the essential content of Figure 2 from Trefethen (2008),
% comparing the convergence of Gauss--Legendre and Clenshaw--Curtis
% quadrature for six functions of varying regularity on [-1, 1]:
%
%     (a) x^{20}        -- polynomial, finite exactness threshold
%     (b) e^x            -- entire function, supergeometric convergence
%     (c) e^{-x^2}       -- entire, supergeometric convergence
%     (d) 1/(1+16x^2)    -- analytic with poles at x = +- i/4, geometric rate
%     (e) e^{-1/x^2}     -- C^infinity but not analytic at x = 0
%     (f) |x|^3          -- C^2 but not C^3, algebraic convergence
%
% The six-panel (3 x 2) figure shows the absolute quadrature error
% |I - I_n| vs the number of nodes n on a semilogarithmic scale for
% n = 1, 2, ..., 40.
%
% Reference: L. N. Trefethen, "Is Gauss Quadrature Better than
%            Clenshaw--Curtis?", SIAM Review 50(1):67--87, 2008, Figure 2.
%
% Generates Figure 15.5: convergence_race.pdf
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
% Test functions
% -------------------------------------------------------------------------
f_poly  = @(x) x.^20;
f_exp   = @(x) exp(x);
f_gauss = @(x) exp(-x.^2);
f_runge = @(x) 1.0 ./ (1.0 + 16.0 * x.^2);
f_flat  = @(x) exp(-1.0 ./ x.^2) .* (x ~= 0);
f_abs3  = @(x) abs(x).^3;

functions = {f_poly, f_exp, f_gauss, f_runge, f_flat, f_abs3};

% -------------------------------------------------------------------------
% Exact integrals over [-1, 1]
% -------------------------------------------------------------------------
I_poly  = 2.0 / 21.0;                      % integral x^20 dx
I_exp   = exp(1) - exp(-1);                % e - e^{-1}
I_gauss = sqrt(pi) * erf(1.0);            % sqrt(pi) * erf(1)
I_runge = atan(4.0) / 2.0;                % arctan(4)/2
I_flat  = integral(@(x) exp(-1./x.^2).*(x~=0), -1, 1, ...
    'AbsTol', 1e-15, 'RelTol', 1e-15);    % numerical high-accuracy reference
I_abs3  = 0.5;                             % 2 * integral_0^1 x^3 dx = 1/2

integrals = [I_poly, I_exp, I_gauss, I_runge, I_flat, I_abs3];

titles = {'$x^{20}$', '$e^x$', '$e^{-x^2}$', ...
          '$1/(1 + 16x^2)$', '$e^{-1/x^2}$', '$|x|^3$'};

% -------------------------------------------------------------------------
% Compute errors
% -------------------------------------------------------------------------
n_max = 40;
n_values = 1:n_max;

err_gauss = zeros(6, n_max);
err_cc    = zeros(6, n_max);

for idx = 1:6
    f = functions{idx};
    I_exact = integrals(idx);

    for i = 1:n_max
        n = n_values(i);

        % Gauss--Legendre with n+1 points
        [x_gl, w_gl] = gauss_legendre(n + 1);
        I_gl = w_gl' * f(x_gl);
        err_gauss(idx, i) = abs(I_gl - I_exact);

        % Clenshaw--Curtis with n+1 points
        [x_cc, w_cc] = clenshaw_curtis_weights(n);
        I_cc = w_cc' * f(x_cc);
        err_cc(idx, i) = abs(I_cc - I_exact);
    end
end

% Clamp zeros for log plot
fl = 1e-16;
err_gauss = max(err_gauss, fl);
err_cc    = max(err_cc, fl);

% -------------------------------------------------------------------------
% Create 3x2 figure
% -------------------------------------------------------------------------
fig = figure('Position', [100, 100, 1000, 1000]);

for idx = 1:6
    subplot(3, 2, idx);

    % Gauss--Legendre: filled circles
    semilogy(n_values, err_gauss(idx, :), 'o', 'Color', NAVY, ...
        'MarkerSize', 3.5, 'MarkerFaceColor', NAVY, ...
        'MarkerEdgeColor', NAVY, 'DisplayName', 'Gauss-Legendre');
    hold on;
    % Clenshaw--Curtis: open circles
    semilogy(n_values, err_cc(idx, :), 'o', 'Color', TEAL, ...
        'MarkerSize', 3.5, 'MarkerFaceColor', 'w', ...
        'MarkerEdgeColor', TEAL, 'LineWidth', 0.8, ...
        'DisplayName', 'Clenshaw-Curtis');
    hold off;

    title(titles{idx}, 'Interpreter', 'latex', 'FontSize', 11);
    xlim([0, n_max + 1]);
    ylim([fl * 0.5, inf]);
    grid on; set(gca, 'GridAlpha', 0.3);

    % x-label only on bottom row
    if idx >= 5
        xlabel('$n$', 'Interpreter', 'latex');
    end

    % y-label only on left column
    if mod(idx, 2) == 1
        ylabel('$|I - I_n|$', 'Interpreter', 'latex');
    end

    % Legend only on the first panel
    if idx == 1
        legend('Location', 'northeast', 'EdgeColor', [0.5 0.5 0.5], ...
            'FontSize', 8);
    end
end

% Save output
exportgraphics(gcf, fullfile(output_dir, 'convergence_race.pdf'), 'ContentType', 'vector');
exportgraphics(gcf, fullfile(output_dir, 'convergence_race.png'), 'Resolution', 300);
fprintf('Figure saved to %s\n', fullfile(output_dir, 'convergence_race.pdf'));

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
