function eigen_benchmark_finite(varargin)
%% eigen_benchmark_finite - Etude 18.3: finite-interval benchmark + N/2 rule.
%
% Dirichlet Laplacian at N = 16, 32, 64; absolute eigenvalue error vs.
% mode number, semilog. Reproduces Boyd (2000) Figs 7.1 and 7.3.
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
% Part of "Computational Etudes: A Spectral Approach"

    [dump_path] = parse_args(varargin{:});
    script_dir = fileparts(mfilename('fullpath'));
    addpath(fullfile(script_dir, '..', 'ch07'));
    out_dir = fullfile(script_dir, '..', '..', '..', ...
                       'textbook', 'figures', 'ch18', 'matlab');
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end
    configure_style();
    [NAVY, CORAL, TEAL] = colours();

    Ns = [16, 32, 64];
    cols = {NAVY, CORAL, TEAL};
    marks = {'o', 's', '^'};
    counts = zeros(1, numel(Ns));
    spectra = cell(1, numel(Ns));

    fig = figure('Units', 'inches', 'Position', [1, 1, 7.6, 4.6], 'Color', 'w');
    hold on;
    for k = 1:numel(Ns)
        N = Ns(k);
        lam = solve_dirichlet(N);
        spectra{k} = lam;
        j = (1:numel(lam))';
        lam_exact = (j * pi / 2).^2;
        err = abs(lam - lam_exact);
        counts(k) = count_good(err, 1e-2);
        err_plot = max(err, 1e-17);
        plot(j, err_plot, [marks{k} '-'], 'Color', cols{k}, ...
            'MarkerFaceColor', 'w', 'MarkerSize', 5, ...
            'LineWidth', 0.6, 'DisplayName', sprintf('$N = %d$', N));
        xline(N/2, ':', 'Color', cols{k}, 'LineWidth', 0.7, 'Alpha', 0.5, ...
              'HandleVisibility', 'off');
    end
    yline(1e-2, '--k', 'LineWidth', 0.6, 'Alpha', 0.5, ...
          'HandleVisibility', 'off');
    text(1.5, 2e-2, '0.01 tolerance', 'FontSize', 8, 'Color', [0 0 0 0.7]);
    set(gca, 'YScale', 'log');
    hold off;
    xlabel('mode number $j$', 'Interpreter', 'latex');
    ylabel('$|\lambda^{\mathrm{num}}_j - \lambda^{\mathrm{exact}}_j|$', ...
        'Interpreter', 'latex');
    title('Dirichlet Laplacian: errors at three resolutions', ...
        'Interpreter', 'latex');
    xlim([0, max(Ns) + 2]); ylim([1e-17, 1e6]);
    grid on; box on;
    legend('Location', 'southeast', 'FontSize', 9, 'Interpreter', 'latex');

    exportgraphics(fig, fullfile(out_dir, 'eigen_benchmark_finite.pdf'), ...
                   'ContentType', 'vector');
    exportgraphics(fig, fullfile(out_dir, 'eigen_benchmark_finite.png'), ...
                   'Resolution', 300);
    close(fig);

    fprintf('[Etude 18.3]\n');
    for k = 1:numel(Ns)
        fprintf('  N = %2d:  %d eigenvalues with |err| < 1e-2   (heuristic N/2 = %d)\n', ...
                Ns(k), counts(k), floor(Ns(k)/2));
    end
    fprintf('  figure: %s\n', fullfile(out_dir, 'eigen_benchmark_finite.pdf'));

    if ~isempty(dump_path)
        r = struct();
        r.Ns = Ns;
        for k = 1:numel(Ns)
            r.spectra.(sprintf('N%d', Ns(k))) = spectra{k}(:)';
            r.good_counts.(sprintf('N%d', Ns(k))) = counts(k);
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

function lam = solve_dirichlet(N)
    [D, ~] = cheb_matrix(N);
    D2 = D * D;
    A = -D2(2:N, 2:N);
    lam = sort(real(eig(A)));
end

function k = count_good(err, tol)
    idx = find(err >= tol, 1, 'first');
    if isempty(idx)
        k = numel(err);
    else
        k = idx - 1;
    end
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
