% sinc_interpolation.m
% Chapter 9: Physical and Fourier Space on Grids
%
% Demonstrates band-limited interpolation using the sinc function.
% Shows interpolation of delta, square wave, and hat functions.
% This is Etude 2 from the chapter.
%
% Author: Dr. Denys Dutykh
% Date: February 2026

function sinc_interpolation()
    % Grid parameters
    h = 1;
    xmax = 10;
    x_grid = -xmax:h:xmax;  % Grid points
    x_fine = (-xmax - h/20):h/10:(xmax + h/20);  % Fine plotting grid

    % Create figure with 3 panels
    figure('Position', [100, 100, 700, 700]);

    % Three grid functions
    titles = {'Discrete delta function', 'Discrete square wave', 'Discrete hat function'};

    for plt = 1:3
        subplot(3, 1, plt);

        switch plt
            case 1
                v = (x_grid == 0);              % Delta function
            case 2
                v = (abs(x_grid) <= 3);         % Square wave
            case 3
                v = max(0, 1 - abs(x_grid)/3);  % Hat function
        end

        % Plot discrete samples
        plot(x_grid, v, 'ko', 'MarkerSize', 8); hold on;

        % Compute and plot interpolant
        p = zeros(size(x_fine));
        for i = 1:length(x_grid)
            p = p + v(i) * sinc_kernel((x_fine - x_grid(i)) / h);
        end
        plot(x_fine, p, 'b-', 'LineWidth', 1.0);

        xlim([-xmax, xmax]);
        ylim([-0.5, 1.5]);
        yticks([0, 0.5, 1]);
        title(titles{plt}, 'FontSize', 11);
        grid on;
        set(gca, 'GridAlpha', 0.3);
        legend('Grid values', 'Band-limited interpolant', 'Location', 'northeast');
        hold off;
    end

    xlabel('x', 'FontSize', 12);

    % Save figure
    set(gcf, 'PaperPositionMode', 'auto');
    print('-dpdf', '../../../textbook/figures/ch09/matlab/sinc_interpolation.pdf');
    fprintf('Figure saved to textbook/figures/ch09/matlab/sinc_interpolation.pdf\n');
end

function s = sinc_kernel(z)
    % Compute sinc(z) = sin(pi*z)/(pi*z), with sinc(0) = 1
    s = ones(size(z));
    mask = z ~= 0;
    s(mask) = sin(pi*z(mask)) ./ (pi*z(mask));
end
