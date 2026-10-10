function rational_chebyshev_map()
%% rational_chebyshev_map - Illustrative figure for the rational Chebyshev
%% TB_n basis (chapter 18, etude 18.4 introductory section).
%
% Three panels:
%   (a) the map x = ell * t / sqrt(1 - t^2) for ell in {1, 2, 4, 8}
%   (b) collocation grids on the real line at fixed N = 24 for the same ell
%   (c) first five basis functions TB_n(x) at ell = 2 over x in [-10, 10]
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
% Part of "Computational Etudes: A Spectral Approach"

    script_dir = fileparts(mfilename('fullpath'));
    addpath(fullfile(script_dir, '..', 'ch07'));
    addpath(script_dir);                     % rational_chebyshev
    out_dir = fullfile(script_dir, '..', '..', '..', ...
                       'textbook', 'figures', 'ch18', 'matlab');
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end
    configure_style();
    [NAVY, SKY, CORAL, TEAL, ORANGE, PURPLE, GOLD] = colours();  %#ok<ASGLU>

    ell_values = [1.0, 2.0, 4.0, 8.0];
    palette_a  = {NAVY, CORAL, TEAL, ORANGE};
    N_grid     = 24;
    ell_basis  = 2.0;

    fig = figure('Units', 'inches', 'Position', [1, 1, 13.0, 4.0], 'Color', 'w');
    tl = tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

    % ----- Panel (a): the map x(t) -----
    nexttile(tl); hold on;
    t_dense = linspace(-0.99, 0.99, 401);
    for k = 1:numel(ell_values)
        ell = ell_values(k);
        x_dense = ell .* t_dense ./ sqrt(1 - t_dense.^2);
        plot(t_dense, x_dense, 'Color', palette_a{k}, 'LineWidth', 1.4, ...
             'DisplayName', sprintf('$\\ell = %g$', ell));
    end
    yline(0, '-k', 'LineWidth', 0.4, 'Alpha', 0.5, 'HandleVisibility', 'off');
    xline(-1, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.6, 'HandleVisibility', 'off');
    xline( 1, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.6, 'HandleVisibility', 'off');
    xlim([-1.05, 1.05]); ylim([-15, 15]);
    xlabel('computational coordinate $t$', 'Interpreter', 'latex');
    ylabel('physical coordinate $x$', 'Interpreter', 'latex');
    title('(a) the map $x = \ell\, t / \sqrt{1 - t^2}$', 'Interpreter', 'latex');
    grid on; box on;
    legend('Location', 'northwest', 'FontSize', 9, 'Interpreter', 'latex');

    % ----- Panel (b): grid clustering -----
    nexttile(tl); hold on;
    x_ref = linspace(-15, 15, 601);
    profile = 0.4 * exp(-(x_ref / 4).^2);
    plot(x_ref, profile + 0.45, 'Color', [0.5 0.5 0.5 0.4], 'LineWidth', 0.8, ...
         'HandleVisibility', 'off');
    text(8.5, 0.88, 'target $e^{-(x/4)^2}$', 'FontSize', 8, ...
         'Color', [0.5 0.5 0.5], 'Interpreter', 'latex');
    row_y = linspace(0.32, 0.04, numel(ell_values));
    for k = 1:numel(ell_values)
        ell = ell_values(k);
        x_int = rational_chebyshev(N_grid, ell);   % returns x_int on positional output 1...
        % rational_chebyshev returns [D1_x, D2_x, x_int]; we need the third
        [~, ~, x_int] = rational_chebyshev(N_grid, ell);
        x_show = x_int(abs(x_int) <= 15);
        plot(x_show, row_y(k) * ones(size(x_show)), '|', ...
             'Color', palette_a{k}, 'MarkerSize', 12, 'LineWidth', 1.1, ...
             'DisplayName', sprintf('$\\ell = %g$', ell));
    end
    xlim([-15, 15]); ylim([0, 1]);
    xlabel('physical coordinate $x$', 'Interpreter', 'latex');
    set(gca, 'YTick', []);
    title(sprintf('(b) collocation grids at $N = %d$ for several $\\ell$', ...
                  N_grid), 'Interpreter', 'latex');
    grid on; box on;
    legend('Location', 'southeast', 'FontSize', 9, 'NumColumns', 2, ...
           'Interpreter', 'latex');

    % ----- Panel (c): first five basis functions at ell = ell_basis -----
    nexttile(tl); hold on;
    x_basis = linspace(-10, 10, 801);
    palette_c = {NAVY, CORAL, TEAL, PURPLE, GOLD};
    for n = 0:4
        plot(x_basis, TB_n(n, x_basis, ell_basis), ...
             'Color', palette_c{n + 1}, 'LineWidth', 1.2, ...
             'DisplayName', sprintf('$\\mathrm{TB}_{%d}$', n));
    end
    yline(0, '-k', 'LineWidth', 0.4, 'Alpha', 0.5, 'HandleVisibility', 'off');
    yline(+1, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5, 'HandleVisibility', 'off');
    yline(-1, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5, 'HandleVisibility', 'off');
    xlim([-10, 10]); ylim([-1.2, 1.25]);
    xlabel('physical coordinate $x$', 'Interpreter', 'latex');
    ylabel('$\mathrm{TB}_n(x)$', 'Interpreter', 'latex');
    title(sprintf('(c) basis functions $\\mathrm{TB}_0,\\dots,\\mathrm{TB}_4$ at $\\ell = %g$', ...
                  ell_basis), 'Interpreter', 'latex');
    grid on; box on;
    legend('Location', 'southeast', 'FontSize', 9, 'NumColumns', 5, ...
           'Interpreter', 'latex');

    exportgraphics(fig, fullfile(out_dir, 'rational_chebyshev_map.pdf'), ...
                   'ContentType', 'vector');
    exportgraphics(fig, fullfile(out_dir, 'rational_chebyshev_map.png'), ...
                   'Resolution', 300);
    close(fig);

    fprintf('[Etude 18.4 / illustrative figure]  rational Chebyshev map\n');
    fprintf('  ell values shown: %s\n', mat2str(ell_values));
    fprintf('  N = %d for grid panel; ell = %g for basis panel\n', N_grid, ell_basis);
    fprintf('  figure: %s\n', fullfile(out_dir, 'rational_chebyshev_map.pdf'));
end

function v = TB_n(n, x, ell)
    % Rational Chebyshev basis: TB_n(x) = T_n( x / sqrt(ell^2 + x^2) ).
    t = x ./ sqrt(ell^2 + x.^2);
    if n == 0
        v = ones(size(t)); return;
    elseif n == 1
        v = t; return;
    end
    Tk_minus2 = ones(size(t));
    Tk_minus1 = t;
    Tk = Tk_minus1;
    for kk = 2:n
        Tk = 2 .* t .* Tk_minus1 - Tk_minus2;
        Tk_minus2 = Tk_minus1;
        Tk_minus1 = Tk;
    end
    v = Tk;
end

function configure_style()
    set(groot, 'defaultAxesFontName', 'CMU Serif');
    set(groot, 'defaultTextFontName', 'CMU Serif');
    set(groot, 'defaultAxesFontSize', 11);
    set(groot, 'defaultAxesLineWidth', 0.8);
    set(groot, 'defaultLineLineWidth', 1.0);
    set(groot, 'defaultAxesBox', 'on');
end

function [NAVY, SKY, CORAL, TEAL, ORANGE, PURPLE, GOLD] = colours()
    NAVY   = [0x14, 0x2D, 0x6E] / 255;
    SKY    = [0x78, 0x96, 0xD2] / 255;
    CORAL  = [0xE7, 0x4C, 0x3C] / 255;
    TEAL   = [0x16, 0xA0, 0x85] / 255;
    ORANGE = [0xE6, 0x7E, 0x22] / 255;
    PURPLE = [0x8E, 0x44, 0xAD] / 255;
    GOLD   = [0xD4, 0xA0, 0x17] / 255;
end
