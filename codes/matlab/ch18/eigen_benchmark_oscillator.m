function eigen_benchmark_oscillator(varargin)
%% eigen_benchmark_oscillator - Etude 18.4: the infinite-interval tax.
%
% Quantum harmonic oscillator u_xx + (lambda - x^2) u = 0 on the real line,
% via the rational_chebyshev algebraic map x = ell t / sqrt(1 - t^2). Exact
% spectrum lambda_j = 2j+1. Reproduces Boyd (2000) Figs 7.4 and 7.6.
%
% The map parameter `ell` (Greek script ell) replaces what earlier drafts
% of this textbook called `L`, so it does not collide with the truncation
% half-width `L` of chapter 20.
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
% Part of "Computational Etudes: A Spectral Approach"

    dump_path = parse_args(varargin{:});
    script_dir = fileparts(mfilename('fullpath'));
    addpath(fullfile(script_dir, '..', 'ch07'));
    addpath(script_dir);                         % for rational_chebyshev
    out_dir = fullfile(script_dir, '..', '..', '..', ...
                       'textbook', 'figures', 'ch18', 'matlab');
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end
    configure_style();
    [NAVY, SKY, CORAL, PURPLE] = colours();

    Ns       = [16, 32];
    ell_scan = [2.0, 4.0, 8.0];
    ell_best = 4.0;
    N_scan   = 32;

    fig = figure('Units', 'inches', 'Position', [1, 1, 11.0, 4.5], 'Color', 'w');
    tl = tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    % --- Panel A: N-scan at ell = ell_best ---
    nexttile(tl); hold on;
    cols_A = {NAVY, CORAL}; marks_A = {'o', 's'};
    counts_N = zeros(1, numel(Ns));
    for k = 1:numel(Ns)
        lam = solve_oscillator(Ns(k), ell_best);
        j = (0:numel(lam)-1)';
        lam_exact = 2 * j + 1;
        err = abs(lam - lam_exact);
        counts_N(k) = count_good(err, 1e-2);
        plot(j + 1, max(err, 1e-17), [marks_A{k} '-'], ...
            'Color', cols_A{k}, 'MarkerFaceColor', 'w', ...
            'MarkerSize', 5, 'LineWidth', 0.6, ...
            'DisplayName', sprintf('$N = %d$, $\\ell = %.1f$', Ns(k), ell_best));
    end
    yline(1e-2, '--k', 'LineWidth', 0.6, 'Alpha', 0.5, 'HandleVisibility', 'off');
    set(gca, 'YScale', 'log');
    hold off;
    xlabel('mode number $j$', 'Interpreter', 'latex');
    ylabel('$|\lambda^{\mathrm{num}}_j - \lambda^{\mathrm{exact}}_j|$', 'Interpreter', 'latex');
    title('infinite-interval tax: oscillator spectrum', 'Interpreter', 'latex');
    ylim([1e-17, 1e4]);
    grid on; box on;
    legend('Location', 'southeast', 'FontSize', 9, 'Interpreter', 'latex');

    % --- Panel B: ell-scan at N = N_scan ---
    nexttile(tl); hold on;
    cols_B = {SKY, NAVY, PURPLE}; marks_B = {'^', 'o', 'v'};
    counts_ell = zeros(1, numel(ell_scan));
    for k = 1:numel(ell_scan)
        lam = solve_oscillator(N_scan, ell_scan(k));
        j = (0:numel(lam)-1)';
        lam_exact = 2 * j + 1;
        err = abs(lam - lam_exact);
        counts_ell(k) = count_good(err, 1e-2);
        plot(j + 1, max(err, 1e-17), [marks_B{k} '-'], ...
            'Color', cols_B{k}, 'MarkerFaceColor', 'w', ...
            'MarkerSize', 5, 'LineWidth', 0.6, ...
            'DisplayName', sprintf('$\\ell = %.1f$', ell_scan(k)));
    end
    yline(1e-2, '--k', 'LineWidth', 0.6, 'Alpha', 0.5, 'HandleVisibility', 'off');
    set(gca, 'YScale', 'log');
    hold off;
    xlabel('mode number $j$', 'Interpreter', 'latex');
    ylabel('$|\lambda^{\mathrm{num}}_j - \lambda^{\mathrm{exact}}_j|$', 'Interpreter', 'latex');
    title(sprintf('$\\ell$-scan at $N = %d$: map-parameter sensitivity', N_scan), ...
        'Interpreter', 'latex');
    ylim([1e-17, 1e4]);
    grid on; box on;
    legend('Location', 'southeast', 'FontSize', 9, 'Interpreter', 'latex');

    exportgraphics(fig, fullfile(out_dir, 'eigen_benchmark_oscillator.pdf'), ...
                   'ContentType', 'vector');
    exportgraphics(fig, fullfile(out_dir, 'eigen_benchmark_oscillator.png'), ...
                   'Resolution', 300);
    close(fig);

    fprintf('[Etude 18.4]  harmonic oscillator on (-inf, +inf)\n');
    fprintf('  ell = %.1f scan over N:\n', ell_best);
    for k = 1:numel(Ns)
        fprintf('    N = %2d:  %d good eigenvalues (|err| < 1e-2)\n', Ns(k), counts_N(k));
    end
    fprintf('  N = %d scan over ell:\n', N_scan);
    for k = 1:numel(ell_scan)
        fprintf('    ell = %.1f:  %d good eigenvalues\n', ell_scan(k), counts_ell(k));
    end
    fprintf('  figure: %s\n', fullfile(out_dir, 'eigen_benchmark_oscillator.pdf'));

    if ~isempty(dump_path)
        r = struct();
        r.Ns = Ns; r.ell_best = ell_best; r.ell_scan = ell_scan;
        for k = 1:numel(Ns)
            r.counts_by_N.(sprintf('N%d', Ns(k))) = counts_N(k);
            r.spectra.(sprintf('N%d_ell%.0f', Ns(k), ell_best)) = solve_oscillator(Ns(k), ell_best)';
        end
        for k = 1:numel(ell_scan)
            r.counts_by_ell.(sprintf('ell%.0f', ell_scan(k))) = counts_ell(k);
            r.spectra.(sprintf('N%d_ell%.0f', N_scan, ell_scan(k))) = solve_oscillator(N_scan, ell_scan(k))';
        end
        fid = fopen(dump_path, 'w');
        fprintf(fid, '%s', jsonencode(r));
        fclose(fid);
        fprintf('  dumped: %s\n', dump_path);
    end
end

function dump_path = parse_args(varargin)
    dump_path = '';
    k = 1;
    while k <= numel(varargin)
        if strcmp(varargin{k}, '--dump')
            dump_path = char(varargin{k+1}); k = k + 2;
        else
            k = k + 1;
        end
    end
end

function lam = solve_oscillator(N, ell)
    [~, D2_x, x] = rational_chebyshev(N, ell);
    H = -D2_x + diag(x.^2);
    lam_all = eig(H);
    keep = abs(imag(lam_all)) < 1e-6;
    lam = real(lam_all(keep));
    lam = sort(lam(lam > 0));
end

function k = count_good(err, tol)
    idx = find(err >= tol, 1, 'first');
    if isempty(idx); k = numel(err); else; k = idx - 1; end
end

function configure_style()
    set(groot, 'defaultAxesFontName', 'CMU Serif');
    set(groot, 'defaultTextFontName', 'CMU Serif');
    set(groot, 'defaultAxesFontSize', 11);
    set(groot, 'defaultAxesLineWidth', 0.8);
    set(groot, 'defaultLineLineWidth', 1.0);
    set(groot, 'defaultAxesBox', 'on');
end

function [NAVY, SKY, CORAL, PURPLE] = colours()
    NAVY   = [0x14, 0x2D, 0x6E] / 255;
    SKY    = [0x78, 0x96, 0xD2] / 255;
    CORAL  = [0xE7, 0x4C, 0x3C] / 255;
    PURPLE = [0x8E, 0x44, 0xAD] / 255;
end
