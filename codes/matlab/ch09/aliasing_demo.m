% aliasing_demo.m
% Chapter 9: Physical and Fourier Space on Grids
%
% Demonstrates aliasing: sin(pi*x) and sin(9*pi*x) coincide on grid h=1/4.
% This is Etude 1 from the chapter.
%
% Author: Dr. Denys Dutykh
% Date: February 2026

function aliasing_demo()
    % Grid spacing
    h = 0.25;

    % Continuous plotting grid
    x_cont = linspace(0, 2, 500);

    % Discrete grid points
    x_grid = 0:h:2;

    % Two functions with different frequencies
    f_cont = sin(pi * x_cont);      % Wavenumber k = pi
    g_cont = sin(9 * pi * x_cont);  % Wavenumber k = 9*pi

    % Samples at grid points
    f_samples = sin(pi * x_grid);
    g_samples = sin(9 * pi * x_grid);

    % Verify aliasing
    max_diff = max(abs(f_samples - g_samples));
    fprintf('Maximum difference at grid points: %.2e\n', max_diff);
    fprintf('Aliasing relation: 9*pi = pi + 2*pi*4 = pi + 2*pi/h\n');

    % Create figure
    figure('Position', [100, 100, 700, 350]);

    % Plot continuous functions
    plot(x_cont, f_cont, 'b-', 'LineWidth', 1.5); hold on;
    plot(x_cont, g_cont, 'r--', 'LineWidth', 1.5);

    % Plot grid samples
    plot(x_grid, f_samples, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'w', ...
         'LineWidth', 2);

    % Add vertical lines at grid points
    for i = 1:length(x_grid)
        xline(x_grid(i), 'Color', [0.5, 0.5, 0.5], 'Alpha', 0.3);
    end

    xlabel('x', 'FontSize', 12);
    ylabel('Value', 'FontSize', 12);
    title('Aliasing on grid h = 1/4: sin(\pi x) and sin(9\pi x) coincide', 'FontSize', 12);
    legend('sin(\pi x)', 'sin(9\pi x)', 'Grid samples', 'Location', 'northeast');
    xlim([-0.05, 2.05]);
    ylim([-1.3, 1.3]);
    grid on;
    set(gca, 'GridAlpha', 0.3);
    hold off;

    % Save figure
    set(gcf, 'PaperPositionMode', 'auto');
    print('-dpdf', '../../../textbook/figures/ch09/matlab/aliasing_demo.pdf');
    fprintf('Figure saved to textbook/figures/ch09/matlab/aliasing_demo.pdf\n');
end
