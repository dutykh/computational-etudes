% eye_aliasing.m
% Chapter 9: Physical and Fourier Space on Grids
%
% Demonstrates visual aliasing: sin(n) for integer n appears to oscillate
% slowly when plotted, because the eye aliases the high frequency to a
% low-frequency pattern.
%
% Author: Dr. Denys Dutykh
% Date: February 2026

function eye_aliasing()
    % Integer arguments: sin(n) for n = 1, 2, ..., 3000
    n = 1:3000;
    y = sin(n);

    % Create figure
    figure('Position', [100, 100, 700, 300]);

    plot(n, y, '.', 'Color', [0.08, 0.18, 0.43], 'MarkerSize', 1);

    xlabel('n', 'FontSize', 12);
    ylabel('sin(n)', 'FontSize', 12);
    title('Visual aliasing: sin(n) for integer n appears to oscillate slowly', ...
          'FontSize', 12);
    xlim([0, 3000]);
    ylim([-1.3, 1.3]);
    grid on;
    set(gca, 'GridAlpha', 0.3);

    % Save figure
    set(gcf, 'PaperPositionMode', 'auto');
    print('-dpdf', '../../../textbook/figures/ch09/matlab/eye_aliasing.pdf');
    fprintf('Figure saved to textbook/figures/ch09/matlab/eye_aliasing.pdf\n');
end
