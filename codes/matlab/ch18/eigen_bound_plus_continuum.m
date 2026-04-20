function eigen_bound_plus_continuum(varargin)
%% eigen_bound_plus_continuum - Etude 18.6: bound states + continuum.
%
% Pöschl-Teller potential V(x) = -nu(nu+1) sech^2(x) on the real line.
% For integer nu, exactly nu bound states at E = -(nu-j)^2 plus continuum.
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
% Part of "Computational Etudes: A Spectral Approach"

    [N1, N2, nu, tol, dump_path] = parse_args(varargin{:});
    script_dir = fileparts(mfilename('fullpath'));
    addpath(fullfile(script_dir, '..', 'ch07'));
    addpath(script_dir);
    out_dir = fullfile(script_dir, '..', '..', '..', ...
                       'textbook', 'figures', 'ch18', 'matlab');
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end
    configure_style();
    [NAVY, CORAL, TEAL] = colours();

    lam1 = solve_poschl_teller(N1, nu);
    lam2 = solve_poschl_teller(N2, nu);
    expected = -((nu:-1:1).^2)';
    rep = spectrum_verify(lam1, lam2, tol);

    fig = figure('Units', 'inches', 'Position', [1, 1, 11.0, 4.5], 'Color', 'w');
    tl = tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    % Panel A: spectra overlay
    nexttile(tl); hold on;
    j = (1:25)';
    plot(j, lam1(1:25), 'o', 'Color', NAVY, 'MarkerFaceColor', 'w', ...
        'MarkerSize', 6, 'LineWidth', 1.1, ...
        'DisplayName', sprintf('$N_1 = %d$', N1));
    plot(j, lam2(1:25), 'x', 'Color', CORAL, 'MarkerSize', 7, ...
        'LineWidth', 1.1, 'DisplayName', sprintf('$N_2 = %d$', N2));
    for k = 1:numel(expected)
        yline(expected(k), '--', 'Color', TEAL, 'LineWidth', 0.6, 'Alpha', 0.6, ...
              'HandleVisibility', 'off');
    end
    yline(0.0, '-k', 'LineWidth', 0.6, 'Alpha', 0.4, 'HandleVisibility', 'off');
    hold off;
    xlabel('mode number $j$', 'Interpreter', 'latex');
    ylabel('eigenvalue $E$', 'Interpreter', 'latex');
    title(sprintf('Pöschl-Teller $V = -%d \\cdot %d\\,\\mathrm{sech}^2 x$ (exact: %d bound states)', ...
          nu, nu+1, numel(expected)), 'Interpreter', 'latex');
    ylim([min(-20, min(lam1(1:25))-1), max(30, max(lam1(1:25))+2)]);
    grid on; box on;
    legend('Location', 'northwest', 'FontSize', 9, 'Interpreter', 'latex');

    % Panel B: drift
    nexttile(tl); hold on;
    j = (1:numel(rep.lam1))';
    inv_ord = 1 ./ max(rep.delta_ordinal, 1e-18);
    inv_nst = 1 ./ max(rep.delta_nearest, 1e-18);
    plot(j, inv_ord, 'o', 'Color', NAVY, 'MarkerFaceColor', 'w', ...
        'MarkerSize', 5, 'LineWidth', 1.0, ...
        'DisplayName', '$1/\delta_\mathrm{ordinal}$');
    plot(j, inv_nst, 'x', 'Color', CORAL, 'MarkerSize', 7, ...
        'LineWidth', 1.0, 'DisplayName', '$1/\delta_\mathrm{nearest}$');
    yline(1/tol, '--k', 'LineWidth', 0.6, 'Alpha', 0.5, 'HandleVisibility', 'off');
    set(gca, 'YScale', 'log');
    hold off;
    xlabel('mode number $j$', 'Interpreter', 'latex');
    ylabel('$1/\delta_j$  (trusted = high)', 'Interpreter', 'latex');
    title(sprintf('drift diagnostic: trusted = %d', rep.n_trusted), 'Interpreter', 'latex');
    xlim([0, 30]);
    grid on; box on;
    legend('Location', 'northeast', 'FontSize', 9, 'Interpreter', 'latex');

    exportgraphics(fig, fullfile(out_dir, 'eigen_bound_plus_continuum.pdf'), 'ContentType', 'vector');
    exportgraphics(fig, fullfile(out_dir, 'eigen_bound_plus_continuum.png'), 'Resolution', 300);
    close(fig);

    fprintf('[Etude 18.6]  Pöschl-Teller with nu = %g\n', nu);
    fprintf('  expected bound states: %s\n', mat2str(expected'));
    fprintf('  trusted (drift < %g): %d\n', tol, rep.n_trusted);
    fprintf('  figure: %s\n', fullfile(out_dir, 'eigen_bound_plus_continuum.pdf'));

    if ~isempty(dump_path)
        r = struct('N1', N1, 'N2', N2, 'nu', nu, 'tol', tol, ...
            'expected_bound', expected', ...
            'n_trusted', rep.n_trusted, ...
            'lam1_low', lam1(1:25)', 'lam2_low', lam2(1:25)', ...
            'trusted_first_k', double(rep.trusted(1:10))');
        fid = fopen(dump_path, 'w');
        fprintf(fid, '%s', jsonencode(r));
        fclose(fid);
        fprintf('  dumped: %s\n', dump_path);
    end
end

function [N1, N2, nu, tol, dump_path] = parse_args(varargin)
    N1 = 60; N2 = 96; nu = 4.0; tol = 1e-3; dump_path = '';
    k = 1;
    while k <= numel(varargin)
        s = varargin{k};
        switch s
            case '--N1';   N1 = as_num(varargin{k+1}); k = k + 2;
            case '--N2';   N2 = as_num(varargin{k+1}); k = k + 2;
            case '--nu';   nu = as_num(varargin{k+1}); k = k + 2;
            case '--tol';  tol = as_num(varargin{k+1}); k = k + 2;
            case '--dump'; dump_path = char(varargin{k+1}); k = k + 2;
            otherwise; k = k + 1;
        end
    end
end

function v = as_num(x)
    if ischar(x) || isstring(x); v = str2double(x); else; v = x; end
end

function lam = solve_poschl_teller(N, nu)
    [~, D2, x] = rational_chebyshev(N, 6.0);
    V = -nu * (nu + 1) ./ cosh(x).^2;
    H = -D2 + diag(V);
    lam_all = eig(H);
    keep = abs(imag(lam_all)) < 1e-6;
    lam = sort(real(lam_all(keep)));
end

function configure_style()
    set(groot, 'defaultAxesFontName', 'CMU Serif');
    set(groot, 'defaultTextFontName', 'CMU Serif');
    set(groot, 'defaultAxesFontSize', 11);
    set(groot, 'defaultAxesLineWidth', 0.8);
    set(groot, 'defaultLineLineWidth', 1.0);
    set(groot, 'defaultAxesBox', 'on');
end

function [NAVY, CORAL, TEAL] = colours()
    NAVY  = [0x14, 0x2D, 0x6E] / 255;
    CORAL = [0xE7, 0x4C, 0x3C] / 255;
    TEAL  = [0x16, 0xA0, 0x85] / 255;
end
