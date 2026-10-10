% trap_fft_coefficients.m
% Chapter 16: Integration of Periodic Functions
%
% Computational Etude 16.8: FFT Computation of Fourier Coefficients.
%
% The FFT of a sample vector, divided by N, is exactly the periodic
% trapezoidal rule applied to the Fourier coefficient integrals.
% Test on f(x) = 1/(2 - cos x), whose Fourier coefficients are known
% in closed form: c_k = (1/sqrt(3)) * (2 - sqrt(3))^|k|.
%
% Generates Figure 16.8: fft_coefficients.pdf  (two-panel figure)
%
% Author: Dr. Denys Dutykh
%         Mathematics Department
%         Khalifa University of Science and Technology
%         Abu Dhabi, UAE
%
% Part of the book "Computational Etudes: A Spectral Approach"
% https://github.com/dutykh/computational-etudes

clear; close all; clc;

set(0, 'DefaultAxesFontSize', 11);
set(0, 'DefaultTextFontSize', 11);
set(0, 'DefaultLineLineWidth', 1.5);
set(0, 'DefaultAxesLineWidth', 0.8);

NAVY  = [20, 45, 110] / 255;
CORAL = [231, 76, 60] / 255;

output_dir = fullfile(fileparts(mfilename('fullpath')), ...
    '..', '..', '..', 'textbook', 'figures', 'ch16', 'matlab');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

% -------------------------------------------------------------------------
% f(x) = 1 / (2 - cos x), Fourier coefficients c_k = r^|k| / sqrt(3)
% -------------------------------------------------------------------------
a = 2;
f = @(x) 1 ./ (a - cos(x));
r = a - sqrt(a^2 - 1);
N = 32;
theta = 2 * pi * (0:N-1) / N;
c_fft = fft(f(theta)) / N;
c_fft_sym = fftshift(c_fft);
k_sym = -N/2 : N/2 - 1;
c_exact = (1 / sqrt(a^2 - 1)) * r .^ abs(k_sym);
err = abs(c_fft_sym - c_exact);
err = max(err, 1e-17);

fprintf('FFT vs exact Fourier coefficients (N = %d)\n', N);
fprintf('%5s  %14s  %14s  %14s\n', 'k', '|c_k| FFT', '|c_k| exact', 'error');
for kk = 0:N/2-1
    idx = N/2 + kk + 1;
    fprintf('%5d  %14.6e  %14.6e  %14.4e\n', ...
        kk, abs(c_fft_sym(idx)), c_exact(idx), err(idx));
end

% -------------------------------------------------------------------------
% Two-panel plot
% -------------------------------------------------------------------------
fig = figure('Position', [100, 100, 1100, 450]);

subplot(1, 2, 1);
semilogy(k_sym, abs(c_exact), '-', 'Color', NAVY, 'LineWidth', 1.2, ...
    'DisplayName', 'exact');
hold on;
semilogy(k_sym, abs(c_fft_sym), 'o', 'Color', CORAL, 'MarkerSize', 5, ...
    'MarkerFaceColor', CORAL, 'DisplayName', sprintf('FFT, N = %d', N));
hold off;
xlabel('Mode index $k$', 'Interpreter', 'latex');
ylabel('$|\hat f_k|$', 'Interpreter', 'latex');
title('(a) Fourier-coefficient magnitudes for $1/(2 - \cos x)$', ...
    'Interpreter', 'latex');
legend('Location', 'northeast', 'Interpreter', 'latex', ...
    'EdgeColor', [0.5 0.5 0.5]);
grid on; set(gca, 'GridAlpha', 0.3);
ylim([1e-12, 1e1]);

subplot(1, 2, 2);
semilogy(k_sym, err, 'o-', 'Color', CORAL, 'MarkerSize', 5, ...
    'LineWidth', 1.0, 'MarkerFaceColor', CORAL, ...
    'DisplayName', '$|\hat f_k^{\mathrm{FFT}} - \hat f_k|$');
xlabel('Mode index $k$', 'Interpreter', 'latex');
ylabel('Absolute error in $\hat f_k$', 'Interpreter', 'latex');
title('(b) FFT error: machine precision in the resolved band', ...
    'Interpreter', 'latex');
legend('Location', 'northeast', 'Interpreter', 'latex', ...
    'EdgeColor', [0.5 0.5 0.5]);
grid on; set(gca, 'GridAlpha', 0.3);
ylim([1e-18, 1e0]);

exportgraphics(gcf, fullfile(output_dir, 'fft_coefficients.pdf'), 'ContentType', 'vector');
exportgraphics(gcf, fullfile(output_dir, 'fft_coefficients.png'), 'Resolution', 300);
fprintf('\nFigure saved to %s\n', fullfile(output_dir, 'fft_coefficients.pdf'));
