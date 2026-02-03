% zero_padding_interpolation.m
% Chapter 9: Physical and Fourier Space on Grids
%
% Demonstrates band-limited interpolation via zero-padding in Fourier space.
% This is Etude 4 from the chapter.
%
% Author: Dr. Denys Dutykh
% Date: February 2026

function zero_padding_interpolation()
    % Coarse grid parameters
    N = 32;
    j = 0:N-1;
    x_coarse = 2*pi*j/N;

    % Test function: exp(sin(x))
    v = exp(sin(x_coarse));

    % Interpolate via zero-padding
    q = 4;
    M = q * N;
    x_fine = 2*pi*(0:M-1)/M;
    v_fine = zero_pad_interpolate(v, q);

    % Exact function for comparison
    v_exact = exp(sin(x_fine));

    % Compute interpolation error
    error = max(abs(v_fine - v_exact));
    fprintf('N = %d coarse points, M = %d fine points\n', N, M);
    fprintf('Maximum interpolation error: %.2e\n', error);

    % Create figure
    figure('Position', [100, 100, 700, 350]);

    % Plot interpolant
    plot(x_fine, v_fine, 'b-', 'LineWidth', 1.2); hold on;

    % Plot coarse samples
    plot(x_coarse, v, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'w', ...
         'LineWidth', 2);

    xlabel('x', 'FontSize', 12);
    ylabel('exp(sin x)', 'FontSize', 12);
    title('Periodic band-limited interpolation via FFT zero-padding', 'FontSize', 12);
    legend('Band-limited interpolant', sprintf('Coarse samples (N = %d)', N), ...
           'Location', 'northeast');
    xlim([0, 2*pi]);
    grid on;
    set(gca, 'GridAlpha', 0.3);
    hold off;

    % Save figure
    set(gcf, 'PaperPositionMode', 'auto');
    print('-dpdf', '../../../textbook/figures/ch09/matlab/zero_padding_interpolation.pdf');
    fprintf('Figure saved to textbook/figures/ch09/matlab/zero_padding_interpolation.pdf\n');
end

function v_fine = zero_pad_interpolate(v, q)
    % Interpolate a periodic function by zero-padding in Fourier space.
    N = length(v);
    M = q * N;

    % Forward FFT
    v_hat = fft(v);

    % Zero-pad in Fourier space
    v_hat_padded = zeros(1, M);

    % Copy low frequencies
    v_hat_padded(1:N/2) = v_hat(1:N/2);
    v_hat_padded(end-N/2+2:end) = v_hat(N/2+2:end);

    % Handle Nyquist mode: split between +N/2 and -N/2
    v_hat_padded(N/2+1) = v_hat(N/2+1) / 2;
    v_hat_padded(M-N/2+1) = v_hat(N/2+1) / 2;

    % Inverse FFT and scale
    v_fine = real(ifft(v_hat_padded)) * q;
end
