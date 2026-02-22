%% ho_orr_sommerfeld.m - Orr-Sommerfeld eigenvalue problem
%
% Chapter 14: Higher-Order Boundary Value Problems
%
% Computational Etude 14.6: Orr-Sommerfeld spectrum at R = 5772, alpha = 1.
%
% The Orr-Sommerfeld equation governs linear stability of parallel shear
% flows. For plane Poiseuille flow U(x) = 1 - x^2 (parabolic profile
% between walls at x = +-1), it reads:
%
%   (1/(i*alpha*R))*(D4 - 2*alpha^2*D2 + alpha^4*I)*v
%     - U*(D2 - alpha^2*I)*v + U''*v = lambda*(-(D2 - alpha^2*I))*v
%
% This is a generalized eigenvalue problem A v = lambda B v where:
%   A = (D4 - 2*alpha^2*D2 + alpha^4*I)/(i*alpha*R)
%       - diag(1 - x^2)*(D2 - alpha^2*I) + U''*I
%   B = -(D2 - alpha^2*I)
%
% The clamped D4 is built via the polynomial trick (Trefethen, Program 38):
% writing u(x) = (1 - x^2)*q(x) enforces u(+-1) = u'(+-1) = 0.
%
% At R = 5772.22, the rightmost eigenvalue lambda_max ~ -0.00007819...,
% placing the flow just barely stable (critical Reynolds number).
%
% Based on Program 40 from Trefethen's "Spectral Methods in MATLAB".
%
% Generates Figure: orr_sommerfeld.pdf - 2x2 grid of spectra
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

%% Parameters
R = 5772.0;          % Critical Reynolds number for plane Poiseuille flow
alpha = 1.0;         % Streamwise wavenumber
N_values = [40, 60, 80, 100];

fprintf('Orr-Sommerfeld Eigenvalue Problem\n');
fprintf('%s\n', repmat('=', 1, 60));
fprintf('  Flow:             Plane Poiseuille  U(x) = 1 - x^2\n');
fprintf('  Reynolds number:  R = %.0f\n', R);
fprintf('  Wavenumber:       alpha = %g\n', alpha);
fprintf('\n');

%% Compute eigenvalues for each N
results = cell(length(N_values), 1);
lam_maxes = zeros(length(N_values), 1);

for idx = 1:length(N_values)
    N = N_values(idx);
    [eigenvalues, lam_max] = orr_sommerfeld_eigenvalues(N, R, alpha);
    results{idx} = eigenvalues;
    lam_maxes(idx) = lam_max;
    fprintf('  N = %3d:  lambda_max = %+.10f %+.10fi\n', ...
            N, real(lam_max), imag(lam_max));
end

fprintf('\n  Expected:  lambda_max ~ -0.00007819...\n');
fprintf('%s\n', repmat('=', 1, 60));

%% Create figure: 2x2 grid of spectra
fig = figure('Units', 'inches', 'Position', [1, 1, 10, 8]);

for idx = 1:length(N_values)
    N = N_values(idx);
    eigenvalues = results{idx};
    lam_max = lam_maxes(idx);

    subplot(2, 2, idx);

    % Plot eigenvalues as dots in the complex plane
    plot(real(eigenvalues), imag(eigenvalues), '.', ...
         'Color', NAVY, 'MarkerSize', 5);
    hold on;

    % Mark the imaginary axis (stability boundary)
    xline(0, ':', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.6);

    hold off;

    % Title with N and lambda_max
    title(sprintf('$N = %d$ \\quad $\\lambda_{\\max} = %.8f$', ...
          N, real(lam_max)), 'FontSize', 11);

    % Axis labels
    xlabel('$\mathrm{Re}\,\lambda$');
    ylabel('$\mathrm{Im}\,\lambda$');

    % Axis limits
    xlim([-0.8, 0.2]);
    ylim([-1.0, 0.0]);

    % Grid
    grid on;
    set(gca, 'GridAlpha', 0.2);

    % Spine styling
    box on;
end

sgtitle(sprintf(['Computational \\''{E}tude 14.6: Orr--Sommerfeld Spectrum ', ...
        '($R = %.0f$, $\\alpha = %g$)'], R, alpha), ...
        'FontSize', 13, 'Interpreter', 'latex');

%% Save figure
exportgraphics(fig, fullfile(output_dir, 'orr_sommerfeld.pdf'), ...
               'ContentType', 'vector');
fprintf('\nFigure saved to: %s\n', fullfile(output_dir, 'orr_sommerfeld.pdf'));

close(fig);


%% =====================================================================
%  Local functions
%  =====================================================================

function [eigenvalues, lam_max] = orr_sommerfeld_eigenvalues(N, R, alpha)
% ORR_SOMMERFELD_EIGENVALUES  Compute Orr-Sommerfeld eigenvalues.
%
%   The clamped D4 is built via the polynomial trick (Trefethen, Program 38):
%   writing v(x) = (1-x^2)*q(x) enforces v(+-1) = v'(+-1) = 0 automatically.
%
%   The Dirichlet D2 is simply the interior block of D^2.
%
%   Boundary conditions v'(+-1) = 0 are imposed via boundary bordering:
%   replacing the first and last rows of A with D(1, 2:N), D(N+1, 2:N).

    [D_cheb, x] = cheb_matrix(N);

    % Powers of the differentiation matrix
    D2_full = D_cheb * D_cheb;
    D3_full = D2_full * D_cheb;
    D4_full = D3_full * D_cheb;

    % --- Clamped D4 via polynomial trick ---
    s = zeros(N + 1, 1);
    s(2:N) = 1.0 ./ (1.0 - x(2:N).^2);
    S = diag(s);

    D4mat = (diag(1.0 - x.^2) * D4_full ...
             - 8.0 * diag(x) * D3_full ...
             - 12.0 * D2_full) * S;

    % Restrict to interior points
    D4_int = D4mat(2:N, 2:N);

    % --- Dirichlet D2: strip boundary rows/columns ---
    D2_int = D2_full(2:N, 2:N);

    % Interior grid points
    x_int = x(2:N);

    M = N - 1;  % number of interior points
    I_int = eye(M);

    % Base flow: U(x) = 1 - x^2, U''(x) = -2
    U_diag = diag(1.0 - x_int.^2);

    % Shift operator: S_op = D2 - alpha^2 I
    S_op = D2_int - alpha^2 * I_int;

    % S^2 = D4 - 2*alpha^2*D2 + alpha^4*I
    S2 = D4_int - 2.0 * alpha^2 * D2_int + alpha^4 * I_int;

    % Assemble A and B matrices
    % A = S^2 / (i*alpha*R) - U*S_op + U'' (U'' = -2)
    A = S2 / (1i * alpha * R) - U_diag * S_op + (-2.0) * I_int;

    % B = -S_op = -(D2 - alpha^2*I)
    B = -S_op;

    % Impose Neumann boundary conditions via row replacement:
    % v'(+1) = D(1, 2:N) . v_int = 0  => row 1
    % v'(-1) = D(N+1, 2:N) . v_int = 0  => row M
    A(1, :)   = D_cheb(1, 2:N);
    B(1, :)   = 0.0;

    A(M, :)   = D_cheb(N + 1, 2:N);
    B(M, :)   = 0.0;

    % Solve generalized eigenvalue problem
    eigenvalues = eig(A, B);

    % Filter out infinite and NaN eigenvalues
    finite_mask = isfinite(eigenvalues) & (abs(eigenvalues) < 1e8);
    eigenvalues = eigenvalues(finite_mask);

    % Find the rightmost eigenvalue (largest real part)
    [~, idx_max] = max(real(eigenvalues));
    lam_max = eigenvalues(idx_max);
end
