%% harmonic_oscillator.m
%
% Solving the quantum harmonic oscillator eigenvalue problem:
%   -u'' + x^2 u = lambda u,  x in R
%
% using the periodic spectral method on a truncated domain [-L, L].
%
% Exact eigenvalues: lambda_n = 2n + 1 for n = 0, 1, 2, ...
% Exact eigenfunctions: Hermite functions u_n(x) = H_n(x) * exp(-x^2/2)
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
% Part of "Computational Etudes: A Spectral Approach"
% https://github.com/dutykh/computational-etudes

clear; close all; clc;

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
output_dir = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch06_smoothness', 'matlab');
output_file = fullfile(output_dir, 'harmonic_oscillator.pdf');

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% Create figure with 2 panels
fig = figure('Units', 'inches', 'Position', [1, 1, 12, 5]);

%% Panel 1: Eigenfunctions
subplot(1, 2, 1);

N = 64;
[eigenvalues, eigenvectors, x] = solve_harmonic_oscillator(N, L);

for n = 0:4
    % Numerical eigenfunction
    u_num = eigenvectors(:, n+1);

    % Exact eigenfunction (Hermite function)
    u_exact = hermite_function(n, x);

    % Fix sign to match exact
    if u_num' * u_exact < 0
        u_num = -u_num;
    end

    % Scale to match
    scale = max(abs(u_exact)) / max(abs(u_num));
    u_num = u_num * scale;

    % Plot exact and numerical
    plot(x, u_exact + 2*n, '-', 'Color', colors{n+1}, 'LineWidth', 1.5);
    hold on;
    plot(x(1:4:end), u_num(1:4:end) + 2*n, 'o', 'Color', colors{n+1}, ...
         'MarkerSize', 4);
end

xlabel('x');
ylabel('u_n(x) (shifted)');
title(sprintf('Harmonic Oscillator Eigenfunctions (N=%d, L=%d)', N, L));
xlim([-6, 6]);

% Add legend entries manually
plot(nan, nan, '-', 'Color', NAVY, 'LineWidth', 1.5, 'DisplayName', 'Exact');
plot(nan, nan, 'o', 'Color', NAVY, 'MarkerSize', 4, 'DisplayName', 'Spectral');
legend('Location', 'northeast', 'FontSize', 9);

% Add annotation
text(4.5, 7.5, {'Lines: exact', 'Dots: spectral'}, 'FontSize', 9, ...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
     'BackgroundColor', 'white', 'EdgeColor', [0.5 0.5 0.5]);
box on;

%% Panel 2: Eigenvalue convergence
subplot(1, 2, 2);

N_values = [6, 12, 18, 24, 30, 36, 42, 48];
exact_eigenvalues = [1, 3, 5, 7];  % First four: 2n+1
markers = {'o', 's', '^', 'd'};

errors = zeros(length(N_values), 4);

for i = 1:length(N_values)
    [eigs, ~, ~] = solve_harmonic_oscillator(N_values(i), L);
    for j = 1:4
        errors(i, j) = abs(eigs(j) - exact_eigenvalues(j));
    end
end

for j = 1:4
    semilogy(N_values, errors(:, j), ['-' markers{j}], 'Color', colors{j}, ...
             'LineWidth', 1.5, 'MarkerSize', 6, ...
             'DisplayName', sprintf('\\lambda_%d = %d', j-1, exact_eigenvalues(j)));
    hold on;
end

% Machine precision line
yline(2.2e-16, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'HandleVisibility', 'off');
text(50, 4e-16, 'Machine precision', 'FontSize', 9, 'Color', [0.5 0.5 0.5], ...
     'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

xlabel('Number of grid points N');
ylabel('Eigenvalue error |\lambda_{computed} - \lambda_{exact}|');
title(sprintf('Eigenvalue Convergence (L=%d)', L));
legend('Location', 'northeast', 'FontSize', 9);
xlim([0, 52]);
ylim([1e-16, 1e2]);
grid on;
box on;

%% Save figure
exportgraphics(fig, output_file, 'ContentType', 'vector');
exportgraphics(fig, strrep(output_file, '.pdf', '.png'), 'Resolution', 300);

fprintf('Figure saved to: %s\n', output_file);

%% Print eigenvalue table (replicating Trefethen's Output 8)
fprintf('\n======================================================================\n');
fprintf('Harmonic Oscillator Eigenvalue Convergence\n');
fprintf('Exact eigenvalues: lambda_n = 2n + 1\n');
fprintf('======================================================================\n');

for N = [6, 12, 18, 24, 30, 36]
    [eigs, ~, ~] = solve_harmonic_oscillator(N, L);
    fprintf('\nN = %d\n', N);
    for j = 0:3
        fprintf('  lambda_%d = %.14f  (error: %.2e)\n', j, eigs(j+1), abs(eigs(j+1) - (2*j+1)));
    end
end

fprintf('\n======================================================================\n');

close(fig);


%% Helper function: solve harmonic oscillator
function [eigenvalues, eigenvectors, x] = solve_harmonic_oscillator(N, L)
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
