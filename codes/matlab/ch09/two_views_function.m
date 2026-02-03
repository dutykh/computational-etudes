% two_views_function.m
% Chapter 9: Physical and Fourier Space on Grids
%
% Demonstrates the two views of a function: physical space and Fourier space.
% Shows a smooth function and its rapidly decaying Fourier transform.
%
% Author: Dr. Denys Dutykh
% Date: February 2026

function two_views_function()
    % Parameters
    L = 20.0;
    N = 2^10;

    % Physical grid
    x = linspace(-L/2, L/2, N+1);
    x(end) = [];  % Remove endpoint for periodicity
    dx = x(2) - x(1);

    % Test function: Gaussian-modulated cosine
    u = exp(-x.^2) .* cos(3*x);

    % Wavenumber vector
    k = 2*pi * [0:N/2-1, -N/2:-1] / (N*dx);

    % Approximate Fourier transform via FFT
    u_hat = dx * fft(u);

    % Shift for plotting (center k=0)
    k_shifted = fftshift(k);
    u_hat_shifted = fftshift(u_hat);

    % Create figure
    figure('Position', [100, 100, 900, 320]);

    % Left: Physical space
    subplot(1, 2, 1);
    plot(x, u, 'b-', 'LineWidth', 1.2);
    xlabel('x', 'FontSize', 12);
    ylabel('u(x)', 'FontSize', 12);
    title('Physical space', 'FontSize', 12);
    xlim([-8, 8]);
    grid on;
    set(gca, 'GridAlpha', 0.3);

    % Right: Fourier space (log scale)
    subplot(1, 2, 2);
    semilogy(k_shifted, abs(u_hat_shifted), 'b-', 'LineWidth', 1.2);
    xlabel('k', 'FontSize', 12);
    ylabel('|\^u(k)|', 'FontSize', 12);
    title('Fourier space (approximate spectrum)', 'FontSize', 12);
    xlim([-15, 15]);
    ylim([1e-16, 10]);
    grid on;
    set(gca, 'GridAlpha', 0.3);

    % Save figure
    set(gcf, 'PaperPositionMode', 'auto');
    print('-dpdf', '../../../textbook/figures/ch09/matlab/two_views_function.pdf');
    fprintf('Figure saved to textbook/figures/ch09/matlab/two_views_function.pdf\n');
end
