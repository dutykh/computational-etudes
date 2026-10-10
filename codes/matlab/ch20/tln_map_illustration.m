function tln_map_illustration()
%% tln_map_illustration - Chapter 20: illustrative figure for the TL_n map.
%
% Three panels:
%   (a) the map y = ell (1+x)/(1-x) for ell in {1, 2, 4, 8}
%   (b) collocation grids on [0, infty) at fixed N = 24 for the same ell
%   (c) first five basis functions TL_n(y) at ell = 2 over y in [0, 20]
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
% Part of "Computational Etudes: A Spectral Approach"

    script_dir = fileparts(mfilename('fullpath'));
    addpath(fullfile(script_dir, '..', 'ch07'));
    out_dir = fullfile(script_dir, '..', '..', '..', ...
                       'textbook', 'figures', 'ch20', 'matlab');
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end
    configure_style();
    [NAVY, CORAL, TEAL, ORANGE, PURPLE, GOLD] = colours();

    ell_values = [1.0, 2.0, 4.0, 8.0];
    palette_a  = {NAVY, CORAL, TEAL, ORANGE};
    N_grid     = 24;
    ell_basis  = 2.0;

    fig = figure('Units', 'inches', 'Position', [1, 1, 11.0, 8.0], 'Color', 'w');
    tl = tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    % ----- Panel (a) -----
    nexttile(tl); hold on;
    x_dense = linspace(-1.0, 0.985, 401);
    for k = 1:numel(ell_values)
        ell = ell_values(k);
        y_dense = ell .* (1 + x_dense) ./ (1 - x_dense);
        plot(x_dense, y_dense, 'Color', palette_a{k}, 'LineWidth', 1.4, ...
             'DisplayName', sprintf('$\\ell = %g$', ell));
    end
    yline(0, '-k', 'LineWidth', 0.4, 'Alpha', 0.5, 'HandleVisibility', 'off');
    xline(+1, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.6, 'HandleVisibility', 'off');
    xlim([-1.05, 1.05]); ylim([0, 60]);
    xlabel('computational coordinate $x$', 'Interpreter', 'latex');
    ylabel('physical coordinate $y$', 'Interpreter', 'latex');
    title('(a) the map $y = \ell\, (1 + x) / (1 - x)$', 'Interpreter', 'latex');
    grid on; box on;
    legend('Location', 'northwest', 'FontSize', 9, 'Interpreter', 'latex');

    % ----- Panel (b) -----
    nexttile(tl); hold on;
    y_ref = linspace(0, 30, 601);
    profile = 0.4 ./ cosh(y_ref / 4);
    plot(y_ref, profile + 0.45, 'Color', [0.5 0.5 0.5 0.4], 'LineWidth', 0.8, ...
         'HandleVisibility', 'off');
    text(18, 0.88, 'target sech(y/4)', 'FontSize', 8, ...
         'Color', [0.5 0.5 0.5]);
    row_y = linspace(0.32, 0.04, numel(ell_values));
    [~, x_int_full] = cheb_matrix(N_grid);
    x_int = x_int_full(2:end-1);
    for k = 1:numel(ell_values)
        ell = ell_values(k);
        y_int = ell .* (1 + x_int) ./ (1 - x_int);
        y_show = y_int(y_int <= 30);
        plot(y_show, row_y(k) * ones(size(y_show)), '|', ...
             'Color', palette_a{k}, 'MarkerSize', 12, 'LineWidth', 1.1, ...
             'DisplayName', sprintf('$\\ell = %g$', ell));
    end
    xlim([0, 30]); ylim([0, 1]);
    xlabel('physical coordinate $y$', 'Interpreter', 'latex');
    set(gca, 'YTick', []);
    title(sprintf('(b) collocation grids at $N = %d$ for several $\\ell$', N_grid), ...
          'Interpreter', 'latex');
    grid on; box on;
    legend('Location', 'northeast', 'FontSize', 9, 'NumColumns', 2, ...
           'Interpreter', 'latex');

    % ----- Panel (c) -----
    nexttile(tl); hold on;
    y_basis = linspace(0.001, 20, 801);
    palette_c = {NAVY, CORAL, TEAL, PURPLE, GOLD};
    for n = 0:4
        plot(y_basis, TL_n(n, y_basis, ell_basis), ...
             'Color', palette_c{n + 1}, 'LineWidth', 1.2, ...
             'DisplayName', sprintf('$\\mathrm{TL}_{%d}$', n));
    end
    yline(0, '-k', 'LineWidth', 0.4, 'Alpha', 0.5, 'HandleVisibility', 'off');
    yline(+1, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5, 'HandleVisibility', 'off');
    yline(-1, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5, 'HandleVisibility', 'off');
    xlim([0, 20]); ylim([-1.2, 1.25]);
    xlabel('physical coordinate $y$', 'Interpreter', 'latex');
    ylabel('$\mathrm{TL}_n(y)$', 'Interpreter', 'latex');
    title(sprintf('(c) basis functions $\\mathrm{TL}_0,\\dots,\\mathrm{TL}_4$ at $\\ell = %g$', ell_basis), ...
          'Interpreter', 'latex');
    grid on; box on;
    legend('Location', 'southeast', 'FontSize', 9, 'NumColumns', 5, ...
           'Interpreter', 'latex');

    % ----- Panel (d): coefficient diagnostic for sech(y) -----
    nexttile(tl); hold on; set(gca, 'YScale', 'log');
    N_diag = 48;
    diag_pal = {CORAL, NAVY, TEAL};
    ell_diag = [0.5, 2.0, 8.0];
    for k = 1:numel(ell_diag)
        a = abs(tln_coeffs_sech(N_diag, ell_diag(k)));
        plot(0:N_diag, a + 1e-18, '-o', 'Color', diag_pal{k}, ...
             'MarkerSize', 3, 'LineWidth', 1.0, ...
             'DisplayName', sprintf('$\\ell = %g$', ell_diag(k)));
    end
    grid on; box on;
    xlabel('degree $n$', 'Interpreter', 'latex');
    ylabel('$|a_n|$ of $\mathrm{TL}_n$ expansion', 'Interpreter', 'latex');
    title('(d) coefficient diagnostic: $\mathrm{TL}_n$ expansion of $\mathrm{sech}(y)$', ...
          'Interpreter', 'latex');
    ylim([1e-17, 10]);
    legend('Location', 'northeast', 'FontSize', 9, 'Interpreter', 'latex');

    exportgraphics(fig, fullfile(out_dir, 'tln_map_illustration.pdf'), ...
                   'ContentType', 'vector');
    exportgraphics(fig, fullfile(out_dir, 'tln_map_illustration.png'), ...
                   'Resolution', 300);
    close(fig);

    fprintf('[Etude 20.x / illustrative figure]  rational Chebyshev TL_n map\n');
    fprintf('  ell values shown: %s\n', mat2str(ell_values));
    fprintf('  N = %d for grid panel; ell = %g for basis panel\n', N_grid, ell_basis);
    fprintf('  figure: %s\n', fullfile(out_dir, 'tln_map_illustration.pdf'));
end

function v = TL_n(n, y, ell)
    t = (y - ell) ./ (y + ell);
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

function a = tln_coeffs_sech(N, ell)
    j = (0:N).';
    x = cos(pi * j / N);
    fv = zeros(size(x));
    interior = abs(x) < 1 - 1e-12;
    y = ell * (1 + x(interior)) ./ (1 - x(interior));
    fv(interior) = 1 ./ cosh(y);
    fv(x < -1 + 1e-12) = 1.0;
    ext = [fv(:); fv(N:-1:2)];
    A = real(fft(ext)) / N;
    A(1) = A(1) * 0.5;
    A(N+1) = A(N+1) * 0.5;
    a = A(1:N+1);
end

function [NAVY, CORAL, TEAL, ORANGE, PURPLE, GOLD] = colours()
    NAVY   = [0x14, 0x2D, 0x6E] / 255;
    CORAL  = [0xE7, 0x4C, 0x3C] / 255;
    TEAL   = [0x16, 0xA0, 0x85] / 255;
    ORANGE = [0xE6, 0x7E, 0x22] / 255;
    PURPLE = [0x8E, 0x44, 0xAD] / 255;
    GOLD   = [0xD4, 0xA0, 0x17] / 255;
end
