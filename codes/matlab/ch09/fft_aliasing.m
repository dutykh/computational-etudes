% fft_aliasing.m
% Chapter 9: Physical and Fourier Space on Grids
%
% Demonstrates aliasing in FFT of high-frequency data.
% sin(17x) sampled at N=32 points shows spike at k=-15 due to aliasing.
% This is Etude 4 from the chapter.
%
% Author: Dr. Denys Dutykh
% Date: February 2026

function fft_aliasing()
    % Grid parameters
    N = 32;
    j = 0:N-1;
    x = 2*pi*j/N;

    % High-frequency function: k=17 > Nyquist limit N/2=16
    u = sin(17*x);

    % Compute FFT
    u_hat = fft(u) / N;

    % Wavenumber vector
    k = [0:N/2-1, -N/2:-1];

    % Shift for plotting
    k_shifted = fftshift(k);
    u_hat_shifted = fftshift(u_hat);

    % Find the aliased location
    fprintf('N = %d, Nyquist limit = %d\n', N, N/2);
    fprintf('True wavenumber: 17\n');
    fprintf('Alias: 17 = -15 + 32, so appears at k = -15\n');
    [~, max_idx] = max(abs(u_hat_shifted));
    fprintf('Maximum |u_hat_k| at k = %d\n', k_shifted(max_idx));

    % Create figure
    figure('Position', [100, 100, 700, 350]);

    % Stem plot of spectrum
    stem(k_shifted, abs(u_hat_shifted), 'b-', 'LineWidth', 1.5, 'MarkerSize', 6);
    hold on;

    % Highlight the aliased peak
    aliased_k = -15;
    aliased_idx = find(k_shifted == aliased_k);
    plot(aliased_k, abs(u_hat_shifted(aliased_idx)), 'ro', 'MarkerSize', 10, ...
         'MarkerFaceColor', 'r');

    % Add annotation
    text(aliased_k + 3, 0.4, sprintf('k = 17 aliases to\nk = %d', aliased_k), ...
         'FontSize', 10);
    annotation('arrow', [0.35, 0.28], [0.68, 0.82]);

    xlabel('k (wavenumber)', 'FontSize', 12);
    ylabel('|\^u_k|', 'FontSize', 12);
    title('FFT of sin(17x) with N = 32 points: aliasing to k = -15', 'FontSize', 12);
    xlim([-N/2 - 1, N/2 + 1]);
    ylim([0, 0.6]);
    legend('Spectrum', 'Aliased peak at k = -15', 'Location', 'northeast');
    grid on;
    set(gca, 'GridAlpha', 0.3);
    hold off;

    % Save figure
    set(gcf, 'PaperPositionMode', 'auto');
    print('-dpdf', '../../../textbook/figures/ch09/matlab/fft_aliasing.pdf');
    fprintf('Figure saved to textbook/figures/ch09/matlab/fft_aliasing.pdf\n');
end
