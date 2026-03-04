%% bc_qnm_poschl_teller.m - Quasinormal modes of the Poschl-Teller potential
%
% Chapter 13: Boundary Conditions in Spectral Methods
%
% Crown Jewel: Quasinormal mode (QNM) computation for the wave equation
% with the Poschl-Teller potential barrier.
%
% Eigenvalue problem:
%   psi'' + (omega^2 - V(x)) psi = 0,   V(x) = V0 * sech^2(x)
%
% Outgoing radiation BCs:
%   psi'(-L) + i*omega*psi(-L) = 0
%   psi'(L)  - i*omega*psi(L)  = 0
%
% Substitution lambda = -i*omega converts to a quadratic eigenvalue problem:
%   psi'' - V*psi - lambda^2 * psi = 0
%
% QEP: (M0 + lambda*M1 + lambda^2*M2) * psi = 0
%   Interior:   M0 = D^2 - diag(V),  M1 = 0,   M2 = -I
%   Row 1 (x=+L):    M0(1,:) = D(1,:),   M1(1,1) = 1,    M2(1,:) = 0
%   Row N+1 (x=-L):  M0(N+1,:) = D(N+1,:), M1(N+1,N+1) = -1, M2(N+1,:) = 0
%
% Solved via polyeig (equivalent to companion linearisation A*v=lambda*B*v,
% but numerically robust when M2 has zero boundary rows making B singular).
%
% Parameters: V0 = 2, L = 8, N = 80
% Exact fundamental QNM: omega_0 = sqrt(7)/2 - i/2
%
% Generates:
%   Figure 1: qnm_spectrum.pdf       - complex omega-plane scatter plot
%   Figure 2: qnm_eigenfunctions.pdf - Re/Im of first two QNMs (2 panels)
%   Figure 3: qnm_convergence.pdf    - error vs N and error vs L (2 panels)
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
output_dir = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch13', 'matlab');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% Parameters
V0    = 2.0;       % Potential depth
L     = 8.0;       % Domain half-width
N     = 80;        % Chebyshev intervals

% Exact fundamental QNM
omega_exact = sqrt(7) / 2 - 1i / 2;

%% Solve the QEP
[omega_all, eigvecs, x_phys] = solve_qnm(N, L, V0);

% Filter physical modes: Im(omega) < 0 (damped), finite, not too large
phys_mask = isfinite(omega_all) & ...
            (imag(omega_all) < -0.01) & ...
            (abs(imag(omega_all)) < 5.0) & ...
            (abs(real(omega_all)) < 10.0);
omega_phys = omega_all(phys_mask);
eigvecs_phys = eigvecs(:, phys_mask);

% Sort by |Im(omega)| (least damped first)
[~, sort_idx] = sort(abs(imag(omega_phys)));
omega_phys = omega_phys(sort_idx);
eigvecs_phys = eigvecs_phys(:, sort_idx);

% Find fundamental: closest to known exact value (robust against spurious
% modes that may have smaller |Im(omega)| than the true fundamental)
[~, fund_idx] = min(abs(omega_phys - omega_exact));

% Reorder so that the fundamental is first, rest sorted by |Im|
rest_idx = [1:fund_idx-1, fund_idx+1:length(omega_phys)];
omega_phys  = omega_phys([fund_idx, rest_idx]);
eigvecs_phys = eigvecs_phys(:, [fund_idx, rest_idx]);

fprintf('Fundamental QNM:\n');
fprintf('  Computed:  omega = %.10f %+.10fi\n', real(omega_phys(1)), imag(omega_phys(1)));
fprintf('  Exact:     omega = %.10f %+.10fi\n', real(omega_exact), imag(omega_exact));
fprintf('  Error:     %.2e\n', abs(omega_phys(1) - omega_exact));

%% =====================================================================
%  Figure 1: QNM spectrum in the complex omega-plane
%  =====================================================================
fig1 = figure('Units', 'inches', 'Position', [1, 1, 7, 5]);

% Plot all eigenvalues (light gray)
plot(real(omega_all), imag(omega_all), '.', 'Color', [0.8 0.8 0.8], ...
     'MarkerSize', 4, 'DisplayName', 'All eigenvalues');
hold on;

% Plot physical modes (Im(omega) < 0)
plot(real(omega_phys), imag(omega_phys), 'o', 'Color', NAVY, ...
     'MarkerSize', 6, 'MarkerFaceColor', NAVY, ...
     'DisplayName', 'Damped modes ($\mathrm{Im}\,\omega < 0$)');

% Highlight exact fundamental mode
plot(real(omega_exact), imag(omega_exact), 'p', 'Color', CORAL, ...
     'MarkerSize', 14, 'MarkerFaceColor', CORAL, ...
     'DisplayName', 'Exact $\omega_0$');

hold off;

xlabel('$\mathrm{Re}(\omega)$');
ylabel('$\mathrm{Im}(\omega)$');
title('Quasinormal Mode Spectrum (P\"oschl--Teller)', 'FontSize', 12);
legend('Location', 'best', 'FontSize', 9);
xlim([-4, 4]);
ylim([-6, 1]);
yline(0, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
grid on;
set(gca, 'GridAlpha', 0.3);
box on;

exportgraphics(fig1, fullfile(output_dir, 'qnm_spectrum.pdf'), ...
               'ContentType', 'vector');
fprintf('Figure 1 saved to: %s\n', fullfile(output_dir, 'qnm_spectrum.pdf'));
close(fig1);

%% =====================================================================
%  Figure 2: Eigenfunctions of first two QNMs
%  =====================================================================
fig2 = figure('Units', 'inches', 'Position', [1, 1, 11, 4.5]);

for mode = 1:2
    subplot(1, 2, mode);

    psi = eigvecs_phys(:, mode);
    % Normalize by max absolute value
    psi = psi / max(abs(psi));

    plot(x_phys, real(psi), '-', 'Color', NAVY, 'LineWidth', 1.5, ...
         'DisplayName', '$\mathrm{Re}(\psi)$');
    hold on;
    plot(x_phys, imag(psi), '--', 'Color', CORAL, 'LineWidth', 1.5, ...
         'DisplayName', '$\mathrm{Im}(\psi)$');
    plot(x_phys, abs(psi), ':', 'Color', TEAL, 'LineWidth', 1.2, ...
         'DisplayName', '$|\psi|$');
    hold off;

    xlabel('$x$');
    ylabel('$\psi(x)$');
    omega_k = omega_phys(mode);
    title(sprintf('Mode %d: $\\omega = %.4f %+.4fi$', mode, real(omega_k), imag(omega_k)), ...
          'FontSize', 11);
    legend('Location', 'best', 'FontSize', 9);
    xlim([-L, L]);
    grid on;
    set(gca, 'GridAlpha', 0.3);
    box on;
end

sgtitle('QNM Eigenfunctions', 'FontSize', 12);

exportgraphics(fig2, fullfile(output_dir, 'qnm_eigenfunctions.pdf'), ...
               'ContentType', 'vector');
fprintf('Figure 2 saved to: %s\n', fullfile(output_dir, 'qnm_eigenfunctions.pdf'));
close(fig2);

%% =====================================================================
%  Figure 3: Convergence studies
%  =====================================================================
fig3 = figure('Units', 'inches', 'Position', [1, 1, 11, 4.5]);

% ---- Left panel: error vs N (fixed L = 8) ----
subplot(1, 2, 1);

N_values = 20:10:120;
errors_N = zeros(size(N_values));

for j = 1:length(N_values)
    N_val = N_values(j);
    [omega_test, ~, ~] = solve_qnm(N_val, L, V0);
    phys = omega_test(imag(omega_test) < 0);
    if ~isempty(phys)
        % Find closest to exact
        [~, idx] = min(abs(phys - omega_exact));
        errors_N(j) = abs(phys(idx) - omega_exact);
    else
        errors_N(j) = NaN;
    end
end

semilogy(N_values, errors_N, 'o-', 'Color', NAVY, 'LineWidth', 1.5, ...
         'MarkerSize', 6, 'MarkerFaceColor', NAVY, 'MarkerEdgeColor', 'w');
xlabel('$N$');
ylabel('$|\omega_0 - \omega_{\mathrm{exact}}|$');
title(sprintf('Convergence vs $N$ ($L = %g$)', L), 'FontSize', 11);
grid on;
set(gca, 'GridAlpha', 0.3);
box on;

% ---- Right panel: error vs L (fixed N = 80) ----
subplot(1, 2, 2);

L_values = 2:1:14;
errors_L = zeros(size(L_values));

for j = 1:length(L_values)
    L_val = L_values(j);
    [omega_test, ~, ~] = solve_qnm(N, L_val, V0);
    phys = omega_test(imag(omega_test) < 0);
    if ~isempty(phys)
        [~, idx] = min(abs(phys - omega_exact));
        errors_L(j) = abs(phys(idx) - omega_exact);
    else
        errors_L(j) = NaN;
    end
end

semilogy(L_values, errors_L, 's-', 'Color', CORAL, 'LineWidth', 1.5, ...
         'MarkerSize', 6, 'MarkerFaceColor', CORAL, 'MarkerEdgeColor', 'w');
xlabel('$L$');
ylabel('$|\omega_0 - \omega_{\mathrm{exact}}|$');
title(sprintf('Convergence vs $L$ ($N = %d$)', N), 'FontSize', 11);
grid on;
set(gca, 'GridAlpha', 0.3);
box on;

sgtitle('QNM Convergence (P\"oschl--Teller, $V_0 = 2$)', 'FontSize', 12);

exportgraphics(fig3, fullfile(output_dir, 'qnm_convergence.pdf'), ...
               'ContentType', 'vector');
fprintf('Figure 3 saved to: %s\n', fullfile(output_dir, 'qnm_convergence.pdf'));
close(fig3);


%% =====================================================================
%  Local functions
%  =====================================================================

function [omega_all, eigvecs_orig, x_phys] = solve_qnm(N, L, V0)
% SOLVE_QNM  Solve the quasinormal mode problem via quadratic eigenvalue problem.
%
%   QEP: (M0 + lambda*M1 + lambda^2*M2) * psi = 0
%   where lambda = -i*omega, so omega = i*lambda.
%
%   Uses polyeig, which performs companion linearisation internally but is
%   numerically robust when M2 has zero boundary rows (making a naive
%   companion B matrix singular and causing eig(A,B) to fail).

    [D_cheb, xi] = cheb_matrix(N);

    % Map from [-1, 1] to [-L, L]: D_phys = D_cheb / L
    D = D_cheb / L;
    D2 = D * D;
    x_phys = L * xi;

    % Potential on the grid
    V = V0 * sech(x_phys).^2;

    I_mat = eye(N + 1);

    % Interior QEP matrices
    M0 = D2 - diag(V);
    M1 = zeros(N + 1);
    M2 = -I_mat;

    % Boundary condition at x = +L (row 1, xi = 1):
    %   Outgoing wave: psi'(+L) - i*omega*psi(+L) = 0
    %   With lambda = -i*omega (so i*omega = -lambda):
    %     psi'(+L) + lambda*psi(+L) = 0
    %   => M0(1,:) = D(1,:), M1(1,1) = 1, M2(1,:) = 0
    M0(1, :) = D(1, :);
    M1(1, :) = 0;
    M1(1, 1) = 1;
    M2(1, :) = 0;

    % Boundary condition at x = -L (row N+1, xi = -1):
    %   Outgoing wave: psi'(-L) + i*omega*psi(-L) = 0
    %   With lambda = -i*omega:
    %     psi'(-L) - lambda*psi(-L) = 0
    %   => M0(N+1,:) = D(N+1,:), M1(N+1,N+1) = -1, M2(N+1,:) = 0
    M0(N+1, :)     = D(N+1, :);
    M1(N+1, :)     = 0;
    M1(N+1, N+1)   = -1;
    M2(N+1, :)     = 0;

    % Solve QEP via polyeig (companion linearisation done internally)
    [X, lambda_all] = polyeig(M0, M1, M2);

    % Convert: omega = i * lambda
    omega_all = 1i * lambda_all;

    % Eigenvectors are already in the original (N+1)-dimensional space
    eigvecs_orig = X;
end
