%% fourier_diff_comparison.m
%
% Chapter 11: Fourier Pseudospectral Methods for Periodic PDEs
%
% Compares Fourier differentiation via dense matrix multiplication (O(N^2))
% with FFT-based differentiation (O(N log N)).
%
% Test function: u(x) = exp(sin(x)) on [0, 2*pi)
% Exact derivative: u'(x) = cos(x) * exp(sin(x))
%
% Generates two figures:
%   - fourier_diff_accuracy.pdf: spectral convergence of both methods
%   - fourier_diff_timing.pdf: timing comparison (matrix vs FFT)
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
% Part of "Computational Etudes: A Spectral Approach"
% https://github.com/dutykh/computational-etudes
% Date: February 2026

clear; close all; clc;

%% Configuration
NAVY  = [0.078 0.176 0.431];
SKY   = [0.471 0.588 0.824];
CORAL = [0.906 0.298 0.235];
TEAL  = [0.086 0.627 0.522];

% Output path
script_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch11', 'matlab');

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% Figure 1: Accuracy comparison
N_values = 4:2:50;
err_matrix = zeros(size(N_values));
err_fft    = zeros(size(N_values));

for idx = 1:length(N_values)
    N = N_values(idx);
    x = 2 * pi * (0:N-1)' / N;
    u = exp(sin(x));
    du_exact = cos(x) .* exp(sin(x));

    % Matrix method
    D = fourier_diff_matrix(N);
    du_mat = D * u;
    err_matrix(idx) = max(abs(du_mat - du_exact));

    % FFT method
    du_fft = fourier_diff_fft(u);
    err_fft(idx) = max(abs(du_fft - du_exact));
end

fig1 = figure('Units', 'inches', 'Position', [1, 1, 6, 4.5]);

semilogy(N_values, err_matrix, 'o-', 'Color', NAVY, 'MarkerSize', 5, ...
         'LineWidth', 1.5, 'DisplayName', 'Matrix $Du$');
hold on;
semilogy(N_values, err_fft, 's--', 'Color', CORAL, 'MarkerSize', 5, ...
         'LineWidth', 1.5, 'DisplayName', '$\mathcal{F}^{-1}(ik\,\hat{u})$');

xlabel('$N$ (number of grid points)', 'Interpreter', 'latex', 'FontSize', 11);
ylabel('$\|u_x^{\mathrm{num}} - u_x^{\mathrm{exact}}\|_\infty$', ...
       'Interpreter', 'latex', 'FontSize', 11);
title('Spectral Convergence: $u(x) = e^{\sin x}$', ...
      'Interpreter', 'latex', 'FontSize', 12);
legend('Interpreter', 'latex', 'FontSize', 10, 'Location', 'northeast');
grid on;
set(gca, 'GridAlpha', 0.3);
ylim([1e-16, inf]);
box on;

exportgraphics(fig1, fullfile(output_dir, 'fourier_diff_accuracy.pdf'), ...
               'ContentType', 'vector');
exportgraphics(fig1, fullfile(output_dir, 'fourier_diff_accuracy.png'), ...
               'Resolution', 300);
fprintf('Figure saved: %s\n', fullfile(output_dir, 'fourier_diff_accuracy.pdf'));
close(fig1);

%% Figure 2: Timing comparison
N_timing = 2.^(3:13);  % 8 to 8192
time_matrix = NaN(size(N_timing));
time_fft    = zeros(size(N_timing));
n_repeats   = 5;

for idx = 1:length(N_timing)
    N = N_timing(idx);
    x = 2 * pi * (0:N-1)' / N;
    u = exp(sin(x));

    % Time matrix method (only up to N = 2048)
    if N <= 2048
        D = fourier_diff_matrix(N);
        tic;
        for r = 1:n_repeats
            du = D * u; %#ok<NASGU>
        end
        time_matrix(idx) = toc / n_repeats;
    end

    % Time FFT method
    tic;
    for r = 1:n_repeats
        du = fourier_diff_fft(u); %#ok<NASGU>
    end
    time_fft(idx) = toc / n_repeats;
end

fig2 = figure('Units', 'inches', 'Position', [1, 1, 6, 4.5]);

% Matrix timing (only where valid)
valid = ~isnan(time_matrix);
loglog(N_timing(valid), time_matrix(valid), 'o-', 'Color', NAVY, ...
       'MarkerSize', 5, 'LineWidth', 1.5, 'DisplayName', 'Matrix $O(N^2)$');
hold on;
loglog(N_timing, time_fft, 's-', 'Color', CORAL, ...
       'MarkerSize', 5, 'LineWidth', 1.5, 'DisplayName', 'FFT $O(N \log N)$');

% Reference slopes
N_ref = [64, 2048];
idx_last = find(valid, 1, 'last');
scale_n2 = time_matrix(idx_last) * (N_ref / N_timing(idx_last)).^2;
loglog(N_ref, scale_n2, ':', 'Color', SKY, 'LineWidth', 1.0, ...
       'DisplayName', '$\sim N^2$');

N_ref2 = [64, 8192];
scale_nlogn = time_fft(end) * (N_ref2 .* log(N_ref2)) / ...
              (N_timing(end) * log(N_timing(end)));
loglog(N_ref2, scale_nlogn, ':', 'Color', [0.902 0.494 0.133], ...
       'LineWidth', 1.0, 'DisplayName', '$\sim N \log N$');

xlabel('$N$', 'Interpreter', 'latex', 'FontSize', 11);
ylabel('Time (seconds)', 'FontSize', 11);
title('Computational Cost: Matrix vs FFT Differentiation', 'FontSize', 12);
legend('Interpreter', 'latex', 'FontSize', 10, 'Location', 'northwest');
grid on;
set(gca, 'GridAlpha', 0.3);
box on;

exportgraphics(fig2, fullfile(output_dir, 'fourier_diff_timing.pdf'), ...
               'ContentType', 'vector');
exportgraphics(fig2, fullfile(output_dir, 'fourier_diff_timing.png'), ...
               'Resolution', 300);
fprintf('Figure saved: %s\n', fullfile(output_dir, 'fourier_diff_timing.pdf'));
close(fig2);


%% -----------------------------------------------------------------------
%  Nested functions
%  -----------------------------------------------------------------------

function D = fourier_diff_matrix(N)
    % Build the N x N Fourier differentiation matrix for the periodic grid
    % x_j = 2*pi*j/N, j = 0, ..., N-1.
    D = zeros(N);
    h = 2 * pi / N;
    for i = 1:N
        for j = 1:N
            if i ~= j
                D(i, j) = 0.5 * (-1)^(i + j) / tan((i - j) * h / 2);
            end
        end
    end
end

function du = fourier_diff_fft(u)
    % Compute the derivative of u on a periodic grid using FFT.
    % u : column vector of length N (values at x_j = 2*pi*j/N)
    N = length(u);
    u_hat = fft(u);

    % Wavenumbers in FFT order
    if mod(N, 2) == 0
        k = [0:N/2-1, 0, -N/2+1:-1]';
    else
        k = [0:(N-1)/2, -(N-1)/2:-1]';
    end

    du_hat = 1i * k .* u_hat;
    du = real(ifft(du_hat));
end
