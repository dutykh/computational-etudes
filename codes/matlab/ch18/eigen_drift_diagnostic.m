function eigen_drift_diagnostic(varargin)
%% eigen_drift_diagnostic - Etude 18.5: build a spectrum lie detector.
%
% Applies the spectrum_verify utility to the Dirichlet Laplacian and the
% harmonic oscillator (on the real line via rational_chebyshev), plotting
% reciprocal drift vs. mode number.
%
% Reproduces Boyd (2000) Fig 7.7.
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
% Part of "Computational Etudes: A Spectral Approach"

    [N1, N2, tol, dump_path] = parse_args(varargin{:});
    script_dir = fileparts(mfilename('fullpath'));
    addpath(fullfile(script_dir, '..', 'ch07'));
    addpath(script_dir);   % for rational_chebyshev and spectrum_verify
    out_dir = fullfile(script_dir, '..', '..', '..', ...
                       'textbook', 'figures', 'ch18', 'matlab');
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end
    configure_style();
    [NAVY, CORAL] = colours();

    % --- Problem A: Dirichlet Laplacian ---
    lam1_A = solve_laplacian(N1);
    lam2_A = solve_laplacian(N2);
    rep_A  = spectrum_verify(lam1_A, lam2_A, tol);

    % --- Problem B: Harmonic oscillator ---
    lam1_B = solve_oscillator(N1, 4.0);
    lam2_B = solve_oscillator(N2, 4.0);
    rep_B  = spectrum_verify(lam1_B, lam2_B, tol);

    fig = figure('Units', 'inches', 'Position', [1, 1, 11.0, 4.5], 'Color', 'w');
    tl = tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile(tl);
    plot_panel(rep_A, sprintf('Dirichlet Laplacian, $N_1 = %d$, $N_2 = %d$', N1, N2), tol, NAVY, CORAL);

    nexttile(tl);
    plot_panel(rep_B, sprintf('Harmonic oscillator (real line), $N_1 = %d$, $N_2 = %d$', N1, N2), tol, NAVY, CORAL);

    exportgraphics(fig, fullfile(out_dir, 'eigen_drift_diagnostic.pdf'), ...
                   'ContentType', 'vector');
    exportgraphics(fig, fullfile(out_dir, 'eigen_drift_diagnostic.png'), ...
                   'Resolution', 300);
    close(fig);

    fprintf('[Etude 18.5]  drift diagnostic at N1 = %d, N2 = %d, tol = %g\n', N1, N2, tol);
    fprintf('  Dirichlet Laplacian : trusted = %d of %d\n', rep_A.n_trusted, numel(rep_A.lam1));
    fprintf('  Harmonic oscillator  : trusted = %d of %d\n', rep_B.n_trusted, numel(rep_B.lam1));
    fprintf('  figure: %s\n', fullfile(out_dir, 'eigen_drift_diagnostic.pdf'));

    if ~isempty(dump_path)
        r = struct('N1', N1, 'N2', N2, 'tol', tol);
        r.laplacian = struct( ...
            'n_trusted',    rep_A.n_trusted, ...
            'lam1',         rep_A.lam1', ...
            'lam2',         rep_A.lam2', ...
            'delta_ordinal', rep_A.delta_ordinal', ...
            'delta_nearest', rep_A.delta_nearest', ...
            'trusted',      double(rep_A.trusted'));
        r.oscillator = struct( ...
            'n_trusted',    rep_B.n_trusted, ...
            'lam1',         rep_B.lam1', ...
            'lam2',         rep_B.lam2', ...
            'delta_ordinal', rep_B.delta_ordinal', ...
            'delta_nearest', rep_B.delta_nearest', ...
            'trusted',      double(rep_B.trusted'));
        fid = fopen(dump_path, 'w');
        fprintf(fid, '%s', jsonencode(r));
        fclose(fid);
        fprintf('  dumped: %s\n', dump_path);
    end
end


function plot_panel(rep, title_str, tol, NAVY, CORAL)
    j = (1:numel(rep.lam1))';
    inv_ord = 1.0 ./ max(rep.delta_ordinal, 1e-18);
    inv_nst = 1.0 ./ max(rep.delta_nearest, 1e-18);
    hold on;
    plot(j, inv_ord, 'o', 'Color', NAVY, 'MarkerFaceColor', 'w', ...
        'MarkerSize', 5, 'LineWidth', 1.0, ...
        'DisplayName', '$1/\delta_{j,\mathrm{ordinal}}$');
    plot(j, inv_nst, 'x', 'Color', CORAL, 'MarkerSize', 7, ...
        'LineWidth', 1.0, ...
        'DisplayName', '$1/\delta_{j,\mathrm{nearest}}$');
    yline(1 / tol, '--k', 'LineWidth', 0.6, 'Alpha', 0.5, 'HandleVisibility', 'off');
    hold off;
    set(gca, 'YScale', 'log');
    xlabel('mode number $j$', 'Interpreter', 'latex');
    ylabel('$1 / \delta_j$  (trusted = high)', 'Interpreter', 'latex');
    title(title_str, 'Interpreter', 'latex');
    grid on; box on;
    legend('Location', 'northeast', 'FontSize', 9, 'Interpreter', 'latex');
end


function [N1, N2, tol, dump_path] = parse_args(varargin)
    N1 = 32; N2 = 48; tol = 1e-3; dump_path = '';
    k = 1;
    while k <= numel(varargin)
        switch varargin{k}
            case '--N1'
                N1 = varargin{k+1};
                if ischar(N1) || isstring(N1); N1 = str2double(N1); end
                k = k + 2;
            case '--N2'
                N2 = varargin{k+1};
                if ischar(N2) || isstring(N2); N2 = str2double(N2); end
                k = k + 2;
            case '--tol'
                tol = varargin{k+1};
                if ischar(tol) || isstring(tol); tol = str2double(tol); end
                k = k + 2;
            case '--dump'
                dump_path = char(varargin{k+1}); k = k + 2;
            otherwise
                k = k + 1;
        end
    end
end


function lam = solve_laplacian(N)
    [D, ~] = cheb_matrix(N);
    D2 = D * D;
    A = -D2(2:N, 2:N);
    lam = sort(real(eig(A)));
end


function lam = solve_oscillator(N, ell)
    [~, D2, x] = rational_chebyshev(N, ell);
    H = -D2 + diag(x.^2);
    lam_all = eig(H);
    keep = abs(imag(lam_all)) < 1e-6;
    lam = real(lam_all(keep));
    lam = sort(lam(lam > 0));
end


function configure_style()
    set(groot, 'defaultAxesFontName', 'CMU Serif');
    set(groot, 'defaultTextFontName', 'CMU Serif');
    set(groot, 'defaultAxesFontSize', 11);
    set(groot, 'defaultAxesLineWidth', 0.8);
    set(groot, 'defaultLineLineWidth', 1.0);
    set(groot, 'defaultAxesBox', 'on');
end


function [NAVY, CORAL] = colours()
    NAVY  = [0x14, 0x2D, 0x6E] / 255;
    CORAL = [0xE7, 0x4C, 0x3C] / 255;
end
