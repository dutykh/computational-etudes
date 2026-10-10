%% harmonic_oscillator.m
%
% Solving the quantum harmonic oscillator eigenvalue problem:
%   -u'' + x^2 u = lambda u,  x in R
%
% using Chebyshev spectral methods on a truncated domain [-L, L],
% with a comparison to the periodic (Fourier) spectral method.
%
% Exact eigenvalues: lambda_n = 2n + 1 for n = 0, 1, 2, ...
% Exact eigenfunctions: Hermite functions u_n(x) = H_n(x) * exp(-x^2/2)
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
% Part of "Computational Etudes: A Spectral Approach"
% https://github.com/dutykh/computational-etudes

clear; close all; clc;

%% Import Chebyshev matrix from Chapter 7
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'ch07'));

%% Configuration
L = 8;  % Domain half-width

% Colors
NAVY = [0.078 0.176 0.431];
CORAL = [0.906 0.298 0.235];
TEAL = [0.086 0.627 0.522];
PURPLE = [0.608 0.349 0.714];
ORANGE = [0.902 0.494 0.133];

colors = {NAVY, CORAL, TEAL, PURPLE, ORANGE};

% Output path
script_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch08', 'matlab');
output_file = fullfile(output_dir, 'harmonic_oscillator.pdf');

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% Create figure with 2 panels
fig = figure('Units', 'inches', 'Position', [1, 1, 12, 5]);

%% Panel 1: Chebyshev eigenfunctions
subplot(1, 2, 1);

N = 64;
[eigenvalues, eigenvectors, x_full] = solve_chebyshev(N, L);
x_int = x_full(2:N);

% Fine uniform grid for smooth exact curves
x_fine = linspace(-6, 6, 500)';

for n = 0:4
    % Numerical eigenfunction (interior only)
    u_int = eigenvectors(:, n+1);

    % Embed in full grid with boundary zeros
    u_num = zeros(N+1, 1);
    u_num(2:N) = u_int;

    % Exact eigenfunction on Chebyshev grid (for sign/scale matching)
    u_exact = hermite_function(n, x_full);

    % Exact eigenfunction on fine grid (for plotting)
    u_exact_fine = hermite_function(n, x_fine);

    % Fix sign to match exact
    if u_num' * u_exact < 0
        u_num = -u_num;
    end

    % Scale to match
    scale = max(abs(u_exact)) / max(abs(u_num));
    u_num = u_num * scale;

    % Plot exact on fine grid and numerical on Chebyshev grid
    plot(x_fine, u_exact_fine + 2*n, '-', 'Color', colors{n+1}, 'LineWidth', 1.5, ...
         'HandleVisibility', 'off');
    hold on;
    plot(x_full(1:4:end), u_num(1:4:end) + 2*n, 'o', 'Color', colors{n+1}, ...
         'MarkerSize', 4, 'HandleVisibility', 'off');
end

xlabel('x');
ylabel('u_n(x) (shifted)');
title(sprintf('Chebyshev Eigenfunctions (N=%d, L=%d)', N, L));
xlim([-6, 6]);

% Generic legend entries (compact, 2 items only)
plot(nan, nan, '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.5, 'DisplayName', 'Exact');
plot(nan, nan, 'o', 'Color', [0.4 0.4 0.4], 'MarkerSize', 4, 'DisplayName', 'Chebyshev');
legend('Location', 'northeast', 'FontSize', 9);
box on;

%% Panel 2: Convergence comparison (Chebyshev vs Fourier)
subplot(1, 2, 2);

N_values = [12, 18, 24, 30, 36, 42, 48];
exact_eigenvalues = [1, 3, 5, 7];  % First four: 2n+1

errors_cheb = zeros(length(N_values), 4);
errors_four = zeros(length(N_values), 4);

for i = 1:length(N_values)
    Nval = N_values(i);

    % Chebyshev
    [eigs_c, ~, ~] = solve_chebyshev(Nval, L);
    for j = 1:4
        errors_cheb(i, j) = max(abs(eigs_c(j) - exact_eigenvalues(j)), 1e-16);
    end

    % Fourier
    [eigs_f, ~, ~] = solve_fourier(Nval, L);
    for j = 1:4
        errors_four(i, j) = max(abs(eigs_f(j) - exact_eigenvalues(j)), 1e-16);
    end
end

% Plot lambda_0 and lambda_1 convergence for both methods
semilogy(N_values, errors_four(:, 1), 'o-', 'Color', NAVY, ...
         'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'Fourier \lambda_0');
hold on;
semilogy(N_values, errors_cheb(:, 1), 's--', 'Color', CORAL, ...
         'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'Chebyshev \lambda_0');
semilogy(N_values, errors_four(:, 2), '^-', 'Color', TEAL, ...
         'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'Fourier \lambda_1');
semilogy(N_values, errors_cheb(:, 2), 'v--', 'Color', PURPLE, ...
         'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'Chebyshev \lambda_1');

% Machine precision line
yline(2.2e-16, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'HandleVisibility', 'off');
text(50, 4e-16, 'Machine precision', 'FontSize', 9, 'Color', [0.5 0.5 0.5], ...
     'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

xlabel('Number of grid points N');
ylabel('Eigenvalue error');
title(sprintf('Fourier vs Chebyshev (L=%d)', L));
legend('Location', 'northeast', 'FontSize', 9);
xlim([8, 52]);
ylim([1e-16, 1e2]);
grid on;
box on;

%% Save figure
exportgraphics(fig, output_file, 'ContentType', 'vector');
exportgraphics(fig, strrep(output_file, '.pdf', '.png'), 'Resolution', 300);

fprintf('Figure saved to: %s\n', output_file);

%% Print comparison table
fprintf('\n==============================================================================\n');
fprintf('Fourier vs Chebyshev Eigenvalue Comparison (L = 8)\n');
fprintf('Exact eigenvalues: lambda_n = 2n + 1\n');
fprintf('==============================================================================\n');
fprintf('%4s  %16s  %16s  %16s  %16s\n', ...
        'N', 'Fourier err(l0)', 'Cheb err(l0)', 'Fourier err(l1)', 'Cheb err(l1)');
fprintf('------------------------------------------------------------------------------\n');

for i = 1:length(N_values)
    fprintf('%4d  %16.2e  %16.2e  %16.2e  %16.2e\n', ...
            N_values(i), errors_four(i, 1), errors_cheb(i, 1), ...
            errors_four(i, 2), errors_cheb(i, 2));
end

fprintf('==============================================================================\n');

% Chebyshev eigenvalue table (N=36)
fprintf('\nChebyshev eigenvalues at N = 36, L = 8:\n');
[eigs_c, ~, ~] = solve_chebyshev(36, L);
for j = 0:3
    fprintf('  lambda_%d = %.14f  (error: %.2e)\n', j, eigs_c(j+1), abs(eigs_c(j+1) - (2*j+1)));
end

fprintf('\nFourier eigenvalues at N = 36, L = 8:\n');
[eigs_f, ~, ~] = solve_fourier(36, L);
for j = 0:3
    fprintf('  lambda_%d = %.14f  (error: %.2e)\n', j, eigs_f(j+1), abs(eigs_f(j+1) - (2*j+1)));
end

close(fig);


%% Helper function: solve using Chebyshev spectral method
function [eigenvalues, eigenvectors, x_full] = solve_chebyshev(N, L)
    % Chebyshev points on [-1, 1]
    [D, t] = cheb_matrix(N);

    % Map to physical domain [-L, L]
    x_full = L * t;

    % Second derivative matrix on physical domain
    D2 = (D * D) / L^2;

    % Extract interior system (strip boundary rows and columns)
    D2_int = D2(2:N, 2:N);
    x_int = x_full(2:N);

    % Build operator: A = -D2 + diag(x^2)
    A = -D2_int + diag(x_int.^2);

    % Solve eigenvalue problem
    [eigenvectors, D_eig] = eig(A);
    eigenvalues = diag(D_eig);

    % Sort by eigenvalue
    [eigenvalues, idx] = sort(real(eigenvalues));
    eigenvectors = eigenvectors(:, idx);
end


%% Helper function: solve using periodic (Fourier) spectral method
function [eigenvalues, eigenvectors, x] = solve_fourier(N, L)
    % Set up grid
    h = 2 * pi / N;
    x_periodic = h * (0:N-1)';
    x = L * (x_periodic - pi) / pi;  % Map to [-L, L]

    % Second derivative matrix
    column = zeros(N, 1);
    column(1) = -pi^2 / (3*h^2) - 1/6;
    for k = 1:N-1
        column(k+1) = -0.5 * ((-1)^k) / sin(h*k/2)^2;
    end
    D2 = toeplitz(column) * (pi/L)^2;

    % Potential matrix
    V = diag(x.^2);

    % Full operator: -D2 + V
    A = -D2 + V;

    % Solve eigenvalue problem
    [eigenvectors, D_eig] = eig(A);
    eigenvalues = diag(D_eig);

    % Sort by eigenvalue
    [eigenvalues, idx] = sort(real(eigenvalues));
    eigenvectors = eigenvectors(:, idx);
end


%% Helper function: Hermite function (normalized eigenfunction)
function psi = hermite_function(n, x)
    % psi_n(x) = (1/sqrt(2^n n! sqrt(pi))) * H_n(x) * exp(-x^2/2)

    % Normalization constant
    norm_const = 1 / sqrt(2^n * factorial(n) * sqrt(pi));

    % Physicist's Hermite polynomial via recurrence
    H_n = hermite_poly(n, x);

    psi = norm_const * H_n .* exp(-x.^2 / 2);
end


%% Helper function: Hermite polynomial
function H = hermite_poly(n, x)
    % Physicist's Hermite polynomial via recurrence relation
    % H_0 = 1, H_1 = 2x, H_{n+1} = 2x*H_n - 2n*H_{n-1}

    if n == 0
        H = ones(size(x));
    elseif n == 1
        H = 2 * x;
    else
        H_prev2 = ones(size(x));  % H_0
        H_prev1 = 2 * x;          % H_1
        for k = 2:n
            H = 2 * x .* H_prev1 - 2 * (k-1) * H_prev2;
            H_prev2 = H_prev1;
            H_prev1 = H;
        end
    end
end
