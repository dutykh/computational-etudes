%% ho_beam_eigenmodes.m - Vibration modes of a clamped beam
%
% Chapter 14: Higher-Order Boundary Value Problems
%
% Computational Etude 14.2: Vibration Modes of a Clamped Beam.
%
% Solves the eigenvalue problem for a clamped-clamped Euler-Bernoulli beam:
%
%   u''''(x) = lambda * u(x),   -1 < x < 1,
%   u(+-1) = u'(+-1) = 0.
%
% The four boundary conditions (two Dirichlet and two Neumann) are imposed
% using the boundary-bordering technique:
%
%   1. Build the Chebyshev differentiation matrix D on N+1 points.
%   2. Compute D4 = D^4 (fourth-derivative matrix).
%   3. Extract the interior block D4_int = D4(2:N, 2:N) (Dirichlet u=0).
%   4. Replace the first and last rows of D4_int with the derivative
%      conditions D(1, 2:N) and D(N+1, 2:N) to enforce u'(+-1) = 0.
%   5. Solve the resulting generalised eigenvalue problem A v = lambda B v
%      where B has zero rows at the replaced positions.
%
% Reference eigenvalues are computed from the transcendental equations:
%   Symmetric modes:       tan(beta) + tanh(beta) = 0
%   Antisymmetric modes:   tan(beta) - tanh(beta) = 0
% In both cases lambda = beta^4.
%
% Generates Figure: beam_eigenmodes.pdf (2 panels)
%   Left:  first 6 eigenmodes of the clamped beam
%   Right: relative eigenvalue error vs mode number
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
N = 30;
n_display = 6;

%% Compute reference eigenvalues
n_ref_modes = 20;
lam_exact = compute_reference_eigenvalues(n_ref_modes);

fprintf('Vibration Modes of a Clamped Beam\n');
fprintf('%s\n', repmat('=', 1, 60));
fprintf('  Equation: u''''(x) = lambda u(x),  x in (-1, 1)\n');
fprintf('  BCs: u(+-1) = u''(+-1) = 0\n');
fprintf('  N = %d Chebyshev points\n\n', N);

fprintf('  Reference eigenvalues (from transcendental equation):\n');
fprintf('  %s\n', repmat('-', 1, 50));
fprintf('  %6s %20s %14s\n', 'Mode', 'lambda (exact)', 'beta');
fprintf('  %s\n', repmat('-', 1, 50));
for k = 1:min(10, length(lam_exact))
    beta_k = lam_exact(k)^0.25;
    fprintf('  %6d %20.6f %14.6f\n', k, lam_exact(k), beta_k);
end
fprintf('  %s\n\n', repmat('-', 1, 50));

%% Solve eigenvalue problem
[eigenvalues, eigenvectors, x_int] = solve_beam_eigenmodes(N);

% Full Chebyshev grid (for interpolation)
[~, x_full] = cheb_matrix(N);

n_compare = min(length(eigenvalues), length(lam_exact));
fprintf('  Number of computed positive eigenvalues: %d\n\n', length(eigenvalues));

fprintf('  Comparison of numerical vs exact eigenvalues:\n');
fprintf('  %s\n', repmat('-', 1, 70));
fprintf('  %6s %18s %18s %14s\n', 'Mode', 'lambda (num)', 'lambda (exact)', 'Rel. Error');
fprintf('  %s\n', repmat('-', 1, 70));
for k = 1:n_compare
    rel_err = abs(eigenvalues(k) - lam_exact(k)) / lam_exact(k);
    fprintf('  %6d %18.6f %18.6f %14.4e\n', k, eigenvalues(k), lam_exact(k), rel_err);
end
fprintf('  %s\n\n', repmat('-', 1, 70));

%% Fine grid for smooth plotting
x_fine = linspace(-1, 1, 500)';

%% Create figure
fig = figure('Units', 'inches', 'Position', [1, 1, 12, 5]);

% ---- Left panel: first 6 eigenmodes ----
subplot(1, 2, 1);
colors_mode = {NAVY, SKY, CORAL, TEAL, PURPLE, ORANGE};

hold on;
for k = 1:n_display
    if k > length(eigenvalues)
        break;
    end

    % Extract eigenmode at interior points
    u_int = real(eigenvectors(:, k));

    % Build full solution including boundary zeros
    u_full = zeros(N + 1, 1);
    u_full(2:N) = u_int;

    % Normalise so max|u| = 1
    u_full = u_full / max(abs(u_full));

    % Fix sign convention: make the first peak (from the left) positive
    % x_full is descending (1 to -1), so reverse for left-to-right
    u_asc = flipud(u_full);
    [~, idx_max] = max(abs(u_asc));
    first_extremum = u_asc(idx_max);
    if first_extremum < 0
        u_full = -u_full;
    end

    % Interpolate onto fine grid using barycentric formula
    u_fine_k = bary_interp(x_full, u_full, x_fine);

    lam_k = eigenvalues(k);
    plot(x_fine, u_fine_k, '-', 'Color', colors_mode{k}, 'LineWidth', 1.5, ...
         'DisplayName', sprintf('Mode %d ($\\lambda_{%d} \\approx %.1f$)', k, k, lam_k));
end
hold off;

xlabel('$x$');
ylabel('$u(x)$');
title(sprintf('First %d eigenmodes ($N = %d$)', n_display, N));
legend('Location', 'south', 'FontSize', 7.5, 'NumColumns', 1);
xlim([-1.05, 1.05]);
grid on;
set(gca, 'GridAlpha', 0.3);
box on;

% ---- Right panel: eigenvalue convergence ----
subplot(1, 2, 2);

modes = 1:n_compare;
rel_errors = abs(eigenvalues(1:n_compare) - lam_exact(1:n_compare)) ...
             ./ lam_exact(1:n_compare);
rel_errors = max(rel_errors, 1e-16);

semilogy(modes, rel_errors, 'o-', 'Color', NAVY, 'LineWidth', 1.2, ...
         'MarkerSize', 6, 'MarkerFaceColor', NAVY, 'MarkerEdgeColor', 'w', ...
         'LineWidth', 1.2);

xlabel('Mode number $n$');
ylabel('$|\lambda_n^{\mathrm{num}} - \lambda_n^{\mathrm{exact}}| / \lambda_n^{\mathrm{exact}}$');
title(sprintf('Relative eigenvalue error ($N = %d$)', N));
grid on;
set(gca, 'GridAlpha', 0.3);
box on;

sgtitle('Clamped Beam Eigenmodes: $u'''''''' = \lambda\, u$, $u(\pm 1) = u''(\pm 1) = 0$', ...
        'FontSize', 12);

%% Save figure
exportgraphics(fig, fullfile(output_dir, 'beam_eigenmodes.pdf'), ...
               'ContentType', 'vector');
fprintf('Figure saved to: %s\n', fullfile(output_dir, 'beam_eigenmodes.pdf'));

close(fig);


%% =====================================================================
%  Local functions
%  =====================================================================

function [A, B, x_int] = biharmonic_clamped_eig_matrix(N)
% BIHARMONIC_CLAMPED_EIG_MATRIX  Build matrices for clamped beam eigenvalues.
%
%   Uses the boundary-bordering technique:
%     1. Build D^4 on the (N+1)-point Chebyshev grid.
%     2. Extract interior rows/columns [2:N, 2:N] for Dirichlet BCs.
%     3. Replace the first row with D(1, 2:N)   (enforces u'(1) = 0).
%     4. Replace the last row with D(N+1, 2:N)  (enforces u'(-1) = 0).
%
%   Outputs:
%     A     - Modified fourth-derivative matrix, size (N-1) x (N-1)
%     B     - Mass matrix (identity with zeroed BC rows), size (N-1) x (N-1)
%     x_int - Interior Chebyshev grid points

    [D, x] = cheb_matrix(N);

    % Fourth-derivative matrix
    D2 = D * D;
    D4 = D2 * D2;

    % Extract interior block (Dirichlet u(+-1) = 0)
    A = D4(2:N, 2:N);

    % Replace first row with derivative condition at x = +1 (index 1)
    % u'(+1) = D(1,:) * u = D(1, 2:N) * u_int = 0
    A(1, :) = D(1, 2:N);

    % Replace last row with derivative condition at x = -1 (index N+1)
    % u'(-1) = D(N+1,:) * u = D(N+1, 2:N) * u_int = 0
    A(end, :) = D(N + 1, 2:N);

    % Mass matrix: identity except zero rows where BCs were imposed
    B = eye(N - 1);
    B(1, :) = 0.0;
    B(end, :) = 0.0;

    x_int = x(2:N);
end

function [eigenvalues, eigenvectors, x_int] = solve_beam_eigenmodes(N)
% SOLVE_BEAM_EIGENMODES  Solve the clamped beam eigenvalue problem.
%
%   Solves u'''' = lambda u with u(+-1) = u'(+-1) = 0.
%
%   Outputs:
%     eigenvalues  - Computed eigenvalues sorted ascending
%     eigenvectors - Corresponding eigenvectors (columns)
%     x_int        - Interior grid points

    [A, B, x_int] = biharmonic_clamped_eig_matrix(N);

    % Solve generalised eigenvalue problem A v = lambda B v
    [V, D_eig] = eig(A, B);
    eigenvalues = diag(D_eig);

    % Filter out infinite/NaN eigenvalues (from the singular rows of B)
    valid = isfinite(eigenvalues) & (abs(eigenvalues) < 1e15);
    eigenvalues = eigenvalues(valid);
    V = V(:, valid);

    % Take real parts (eigenvalues should be real and positive)
    eigenvalues = real(eigenvalues);

    % Sort by ascending real part
    [eigenvalues, idx] = sort(eigenvalues);
    V = real(V(:, idx));

    % Keep only positive eigenvalues (physical modes)
    pos_mask = eigenvalues > 0;
    eigenvalues = eigenvalues(pos_mask);
    eigenvectors = V(:, pos_mask);
end

function lambdas = compute_reference_eigenvalues(n_modes)
% COMPUTE_REFERENCE_EIGENVALUES  Reference eigenvalues for clamped beam.
%
%   Computes eigenvalues from the transcendental equations:
%     Symmetric modes:       tan(beta) + tanh(beta) = 0
%     Antisymmetric modes:   tan(beta) - tanh(beta) = 0
%   with lambda = beta^4.

    betas = [];
    eps_search = 0.01;
    n_search = n_modes + 10;

    for k = 0:n_search - 1
        % Continuity interval for tan: (k*pi - pi/2, k*pi + pi/2)
        lo = k * pi - pi / 2 + eps_search;
        hi = k * pi + pi / 2 - eps_search;
        if lo < eps_search
            lo = eps_search;
        end

        % Check for symmetric mode root: tan(beta) + tanh(beta) = 0
        sym_lo = tan(lo) + tanh(lo);
        sym_hi = tan(hi) + tanh(hi);
        if sym_lo * sym_hi < 0
            root = fzero(@(b) tan(b) + tanh(b), [lo, hi]);
            if root > 0.5 && all(abs(root - betas) > 0.01)
                betas = [betas; root]; %#ok<AGROW>
            end
        end

        % Check for antisymmetric mode root: tan(beta) - tanh(beta) = 0
        antisym_lo = tan(lo) - tanh(lo);
        antisym_hi = tan(hi) - tanh(hi);
        if antisym_lo * antisym_hi < 0
            root = fzero(@(b) tan(b) - tanh(b), [lo, hi]);
            if root > 0.5 && all(abs(root - betas) > 0.01)
                betas = [betas; root]; %#ok<AGROW>
            end
        end
    end

    % Sort betas and compute eigenvalues lambda = beta^4
    betas = sort(betas);
    lambdas = betas(1:min(n_modes, length(betas))).^4;
end

function p = bary_interp(x_nodes, f_nodes, x_eval)
% BARY_INTERP  Barycentric interpolation at Chebyshev nodes.

    N = length(x_nodes) - 1;

    % Barycentric weights for Chebyshev-Gauss-Lobatto points
    w = ones(N + 1, 1);
    w(1)   = 0.5;
    w(N+1) = 0.5;
    w = w .* (-1).^(0:N)';

    p = zeros(size(x_eval));
    for j = 1:length(x_eval)
        xj = x_eval(j);
        d = xj - x_nodes;

        [min_d, idx] = min(abs(d));
        if min_d < 1e-15
            p(j) = f_nodes(idx);
        else
            terms = w ./ d;
            p(j) = sum(terms .* f_nodes) / sum(terms);
        end
    end
end
