function eigen_physically_spurious(varargin)
%% eigen_physically_spurious - Etude 18.7: manufacturing fake instability.
%
% Gottlieb-Orszag streamfunction: nu u_xxxx = lambda u_xx, clamped BCs.
% Naive formulation gives spurious large POSITIVE eigenvalues; cured
% formulation does not.
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
% Part of "Computational Etudes: A Spectral Approach"

    dump_path = parse_args(varargin{:});
    script_dir = fileparts(mfilename('fullpath'));
    addpath(fullfile(script_dir, '..', 'ch07'));
    out_dir = fullfile(script_dir, '..', '..', '..', ...
                       'textbook', 'figures', 'ch18', 'matlab');
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end
    configure_style();
    [NAVY, CORAL, TEAL] = colours();

    Ns_scan = [16, 24, 32, 48, 64, 96];
    nu = 1.0;

    [A1, B1] = naive_go(32, nu); lam1 = sfilt(eig(A1, B1));
    [A2, B2] = naive_go(48, nu); lam2 = sfilt(eig(A2, B2));
    pos1 = lam1(lam1 > 1);  pos2 = lam2(lam2 > 1);

    max_positive = nan(1, numel(Ns_scan));
    for k = 1:numel(Ns_scan)
        [A, B] = naive_go(Ns_scan(k), nu);
        lam = sfilt(eig(A, B));
        pos = lam(lam > 1);
        if ~isempty(pos); max_positive(k) = max(pos); end
    end

    [Ac1, Bc1] = cured_go(32, nu); lam_c1 = sfilt(eig(Ac1, Bc1));
    [Ac2, Bc2] = cured_go(48, nu); lam_c2 = sfilt(eig(Ac2, Bc2));
    pc1 = lam_c1(lam_c1 > 1); pc2 = lam_c2(lam_c2 > 1);

    fig = figure('Units', 'inches', 'Position', [1, 1, 14.0, 4.2], 'Color', 'w');
    tl = tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile(tl); hold on;
    scatter(1:numel(lam1), lam1, 24, 'MarkerEdgeColor', NAVY, ...
        'MarkerFaceColor', 'w', 'LineWidth', 1.0, ...
        'DisplayName', sprintf('$N = 32$'));
    scatter(1:numel(lam2), lam2, 18, 'x', 'MarkerEdgeColor', CORAL, ...
        'LineWidth', 1.0, 'DisplayName', sprintf('$N = 48$'));
    yline(0, '-k', 'LineWidth', 0.6, 'Alpha', 0.5, 'HandleVisibility', 'off');
    hold off;
    xlabel('mode number $j$', 'Interpreter', 'latex');
    ylabel('$\lambda$ (naive)', 'Interpreter', 'latex');
    title('naive: spurious $\lambda > 0$ visible', 'Interpreter', 'latex');
    grid on; box on;
    legend('Location', 'northeast', 'FontSize', 9, 'Interpreter', 'latex');

    nexttile(tl);
    mask = isfinite(max_positive);
    loglog(Ns_scan(mask), max_positive(mask), 'o-', 'Color', NAVY, ...
        'MarkerFaceColor', 'w', 'MarkerSize', 6, 'LineWidth', 0.8, ...
        'DisplayName', 'largest positive (naive)');
    hold on;
    if any(mask)
        Nm = Ns_scan(mask); pm = max_positive(mask);
        loglog(Nm, pm(1) * (Nm / Nm(1)).^4, '--', 'Color', TEAL, ...
            'LineWidth', 0.8, 'DisplayName', '$N^4$ reference');
    end
    hold off;
    xlabel('$N$', 'Interpreter', 'latex');
    ylabel('largest positive $\lambda$', 'Interpreter', 'latex');
    title('scaling of spurious mode magnitude', 'Interpreter', 'latex');
    grid on; box on;
    legend('Location', 'northwest', 'FontSize', 9, 'Interpreter', 'latex');

    nexttile(tl); hold on;
    scatter(1:numel(lam_c1), lam_c1, 24, 'MarkerEdgeColor', NAVY, ...
        'MarkerFaceColor', 'w', 'LineWidth', 1.0, ...
        'DisplayName', '$N = 32$ (cured)');
    scatter(1:numel(lam_c2), lam_c2, 18, 'x', 'MarkerEdgeColor', CORAL, ...
        'LineWidth', 1.0, 'DisplayName', '$N = 48$ (cured)');
    yline(0, '-k', 'LineWidth', 0.6, 'Alpha', 0.5, 'HandleVisibility', 'off');
    hold off;
    xlabel('mode number $j$', 'Interpreter', 'latex');
    ylabel('$\lambda$ (cured)', 'Interpreter', 'latex');
    title(sprintf('cured: %d spurious at $N=32$, %d at $N=48$', ...
        numel(pc1), numel(pc2)), 'Interpreter', 'latex');
    grid on; box on;
    legend('Location', 'northeast', 'FontSize', 9, 'Interpreter', 'latex');

    exportgraphics(fig, fullfile(out_dir, 'eigen_physically_spurious.pdf'), 'ContentType', 'vector');
    exportgraphics(fig, fullfile(out_dir, 'eigen_physically_spurious.png'), 'Resolution', 300);
    close(fig);

    fprintf('[Etude 18.7]  Gottlieb-Orszag streamfunction (nu = %g)\n', nu);
    fprintf('  Naive positives at N=32: %d  at N=48: %d\n', numel(pos1), numel(pos2));
    fprintf('  Cured positives at N=32: %d  at N=48: %d\n', numel(pc1), numel(pc2));
    fprintf('  figure: %s\n', fullfile(out_dir, 'eigen_physically_spurious.pdf'));

    if ~isempty(dump_path)
        r = struct('nu', nu, 'Ns_scan', Ns_scan, 'max_positive', max_positive);
        r.naive_n_positive = struct('N32', numel(pos1), 'N48', numel(pos2));
        r.cured_n_positive = struct('N32', numel(pc1), 'N48', numel(pc2));
        r.naive_positives_32 = pos1';
        r.naive_positives_48 = pos2';
        r.cured_spectrum_32_first5 = lam_c1(1:min(5, numel(lam_c1)))';
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
        else; k = k + 1; end
    end
end

function [A, B] = naive_go(N, nu)
    [D, ~] = cheb_matrix(N);
    D2 = D * D; D4 = D2 * D2;
    A = nu * D4; B = D2;
    ID = eye(N + 1);
    A(1,   :) = ID(1,   :); B(1,   :) = 0;
    A(2,   :) = D(1,    :); B(2,   :) = 0;
    A(N,   :) = D(N+1,  :); B(N,   :) = 0;
    A(N+1, :) = ID(N+1, :); B(N+1, :) = 0;
end

function [A, B] = cured_go(N, nu)
    [D, ~] = cheb_matrix(N);
    D2 = D * D;
    D2i = D2(2:N, 2:N);
    A = nu * (D2i * D2i);
    B = D2i;
    A(1, :)   = D(1,   2:N); B(1, :)   = 0;
    A(end, :) = D(N+1, 2:N); B(end, :) = 0;
end

function s = sfilt(lam)
    lam = lam(isfinite(lam));
    lam = real(lam(abs(imag(lam)) < 1e-6));
    s = sort(lam);
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
