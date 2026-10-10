% smoothness_spectra.m
% Chapter 9: Physical and Fourier Space on Grids
%
% Demonstrates the relationship between function smoothness and spectral decay.
% Shows spectra of: exp(sin(x)), periodic hat, and square wave.
%
% Author: Dr. Denys Dutykh
% Date: February 2026

function smoothness_spectra()
    N = 256;
    x = 2*pi*(0:N-1)/N;

    % Create figure with 3 panels
    figure('Position', [100, 100, 1000, 300]);

    % 1. Analytic: exp(sin(x))
    subplot(1, 3, 1);
    f1 = exp(sin(x));
    [k1, spec1] = compute_spectrum(f1, N);
    semilogy(k1, spec1, 'b.-', 'MarkerSize', 4, 'LineWidth', 0.8);
    xlabel('k', 'FontSize', 11);
    ylabel('|\^f_k|', 'FontSize', 11);
    title('exp(sin x) (analytic)', 'FontSize', 11);
    xlim([0, N/2]);
    ylim([1e-16, 10]);
    grid on;
    set(gca, 'GridAlpha', 0.3);

    % 2. C^0 but not C^1: periodic hat function
    subplot(1, 3, 2);
    f2 = max(0, 1 - abs(x - pi) / (pi/2));
    [k2, spec2] = compute_spectrum(f2, N);
    loglog(k2, spec2, 'b.-', 'MarkerSize', 4, 'LineWidth', 0.8); hold on;
    % Reference line for k^{-2}
    k_ref = [5, 100];
    loglog(k_ref, 0.5 * k_ref.^(-2), 'r--', 'LineWidth', 1.5);
    xlabel('k', 'FontSize', 11);
    ylabel('|\^f_k|', 'FontSize', 11);
    title('Hat function (C^0)', 'FontSize', 11);
    xlim([1, N/2]);
    legend('Spectrum', 'k^{-2}', 'Location', 'northeast');
    grid on;
    set(gca, 'GridAlpha', 0.3);
    hold off;

    % 3. Discontinuous: square wave
    subplot(1, 3, 3);
    f3 = sign(sin(x));
    [k3, spec3] = compute_spectrum(f3, N);
    loglog(k3, spec3, 'b.-', 'MarkerSize', 4, 'LineWidth', 0.8); hold on;
    % Reference line for k^{-1}
    k_ref = [1, 100];
    loglog(k_ref, 1.5 * k_ref.^(-1), 'r--', 'LineWidth', 1.5);
    xlabel('k', 'FontSize', 11);
    ylabel('|\^f_k|', 'FontSize', 11);
    title('Square wave (discontinuous)', 'FontSize', 11);
    xlim([1, N/2]);
    legend('Spectrum', 'k^{-1}', 'Location', 'northeast');
    grid on;
    set(gca, 'GridAlpha', 0.3);
    hold off;

    % Save figure
    set(gcf, 'PaperPositionMode', 'auto');
    print('-dpdf', '../../../textbook/figures/ch09/matlab/smoothness_spectra.pdf');
    fprintf('Figure saved to textbook/figures/ch09/matlab/smoothness_spectra.pdf\n');
end

function [k, f_hat_mag] = compute_spectrum(f, N)
    % Compute Fourier spectrum of a periodic function
    f_hat = fft(f) / N;
    % Take positive wavenumbers only (excluding k=0)
    k = 1:N/2-1;
    f_hat_mag = abs(f_hat(2:N/2));
end
