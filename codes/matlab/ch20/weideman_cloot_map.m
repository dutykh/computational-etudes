function weideman_cloot_map()
%% weideman_cloot_map - Chapter 20: illustrative figure for the
%% Weideman--Cloot (1990) sinh map.
%
% Three panels:
%   (a) the map y = sinh(ell t) for ell in {0.5, 1, 1.5, 2}
%   (b) sinh-mapped grids at fixed N = 64
%   (c) sech(y) shown in y-space and t-space (rescaled) at ell = 1
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
% Part of "Computational Etudes: A Spectral Approach"

    script_dir = fileparts(mfilename('fullpath'));
    out_dir = fullfile(script_dir, '..', '..', '..', ...
                       'textbook', 'figures', 'ch20', 'matlab');
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end
    configure_style();
    [NAVY, CORAL, TEAL, ORANGE] = colours();

    ell_values = [0.5, 1.0, 1.5, 2.0];
    palette_a  = {NAVY, CORAL, TEAL, ORANGE};
    N_grid     = 64;
    ell_demo   = 1.0;

    fig = figure('Units', 'inches', 'Position', [1, 1, 13.0, 4.0], 'Color', 'w');
    tl = tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

    % Panel (a)
    nexttile(tl); hold on;
    t_dense = linspace(-pi, pi, 401);
    for k = 1:numel(ell_values)
        ell = ell_values(k);
        y_dense = sinh(ell .* t_dense);
        y_max = sinh(ell * pi);
        plot(t_dense, y_dense, 'Color', palette_a{k}, 'LineWidth', 1.4, ...
             'DisplayName', sprintf('$\\ell = %g$, $y_{\\max} \\approx %.1f$', ell, y_max));
    end
    yline(0, '-k', 'LineWidth', 0.4, 'Alpha', 0.5, 'HandleVisibility', 'off');
    xline(-pi, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.6, 'HandleVisibility', 'off');
    xline(+pi, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.6, 'HandleVisibility', 'off');
    xlim([-pi - 0.1, pi + 0.1]); ylim([-30, 30]);
    xlabel('computational coordinate $t \in [-\pi, \pi]$', 'Interpreter', 'latex');
    ylabel('physical coordinate $y$', 'Interpreter', 'latex');
    title('(a) the map $y = \sinh(\ell\, t)$', 'Interpreter', 'latex');
    grid on; box on;
    legend('Location', 'northwest', 'FontSize', 8, 'Interpreter', 'latex');

    % Panel (b)
    nexttile(tl); hold on;
    y_axis = linspace(-30, 30, 601);
    profile = 0.4 ./ cosh(y_axis);
    plot(y_axis, profile + 0.45, 'Color', [0.5 0.5 0.5 0.4], ...
         'LineWidth', 0.8, 'HandleVisibility', 'off');
    text(15, 0.88, 'target sech(y)', 'FontSize', 8, ...
         'Color', [0.5 0.5 0.5]);
    row_y = linspace(0.32, 0.04, numel(ell_values));
    j = 0:N_grid-1;
    t_nodes = -pi + 2 * pi * j / N_grid;
    for k = 1:numel(ell_values)
        ell = ell_values(k);
        y_nodes = sinh(ell .* t_nodes);
        y_show = y_nodes(abs(y_nodes) <= 30);
        plot(y_show, row_y(k) * ones(size(y_show)), '|', ...
             'Color', palette_a{k}, 'MarkerSize', 12, 'LineWidth', 1.1, ...
             'DisplayName', sprintf('$\\ell = %g$', ell));
    end
    xlim([-30, 30]); ylim([0, 1]);
    xlabel('physical coordinate $y$', 'Interpreter', 'latex');
    set(gca, 'YTick', []);
    title(sprintf('(b) sinh-mapped grids at $N = %d$ for several $\\ell$', N_grid), ...
          'Interpreter', 'latex');
    grid on; box on;
    legend('Location', 'southeast', 'FontSize', 9, 'NumColumns', 2, ...
           'Interpreter', 'latex');

    % Panel (c)
    nexttile(tl); hold on;
    y_max = sinh(ell_demo * pi);
    y_phys = linspace(-y_max, y_max, 801);
    sech_in_y = 1.0 ./ cosh(y_phys);
    plot(y_phys, sech_in_y, 'Color', NAVY, 'LineWidth', 1.3, ...
         'DisplayName', '$\mathrm{sech}(y)$ in $y$-space');
    t_axis = linspace(-pi, pi, 801);
    sech_in_t = 1.0 ./ cosh(sinh(ell_demo .* t_axis));
    t_scaled = t_axis * (y_max / pi);
    plot(t_scaled, sech_in_t, '--', 'Color', CORAL, 'LineWidth', 1.3, ...
         'DisplayName', 'same function in $t$-space (scaled)');
    yline(0, '-k', 'LineWidth', 0.4, 'Alpha', 0.5, 'HandleVisibility', 'off');
    xlim([-y_max, y_max]); ylim([-0.1, 1.1]);
    xlabel('physical coordinate $y$ (or scaled $t$)', 'Interpreter', 'latex');
    ylabel('function value', 'Interpreter', 'latex');
    title(sprintf('(c) $\\mathrm{sech}(y)$ in two coordinate frames at $\\ell = %g$', ell_demo), ...
          'Interpreter', 'latex');
    grid on; box on;
    legend('Location', 'northeast', 'FontSize', 9, 'Interpreter', 'latex');

    exportgraphics(fig, fullfile(out_dir, 'weideman_cloot_map.pdf'), ...
                   'ContentType', 'vector');
    exportgraphics(fig, fullfile(out_dir, 'weideman_cloot_map.png'), ...
                   'Resolution', 300);
    close(fig);

    fprintf('[Section 20.x / illustrative figure]  Weideman--Cloot sinh map\n');
    fprintf('  ell values shown: %s\n', mat2str(ell_values));
    fprintf('  y_max = sinh(ell pi) per ell:');
    for k = 1:numel(ell_values)
        fprintf(' %.2f', sinh(ell_values(k) * pi));
    end
    fprintf('\n');
    fprintf('  figure: %s\n', fullfile(out_dir, 'weideman_cloot_map.pdf'));
end

function configure_style()
    set(groot, 'defaultAxesFontName', 'CMU Serif');
    set(groot, 'defaultTextFontName', 'CMU Serif');
    set(groot, 'defaultAxesFontSize', 11);
    set(groot, 'defaultAxesLineWidth', 0.8);
    set(groot, 'defaultLineLineWidth', 1.0);
    set(groot, 'defaultAxesBox', 'on');
end

function [NAVY, CORAL, TEAL, ORANGE] = colours()
    NAVY   = [0x14, 0x2D, 0x6E] / 255;
    CORAL  = [0xE7, 0x4C, 0x3C] / 255;
    TEAL   = [0x16, 0xA0, 0x85] / 255;
    ORANGE = [0xE6, 0x7E, 0x22] / 255;
end
