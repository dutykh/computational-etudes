function eigen_heinrichs_condition(varargin)
%% eigen_heinrichs_condition - Etude 18.8: condition-number surgery.
% Fourth-order clamped eigenproblem: naive D^4 bordering vs. Heinrichs.
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
% Part of "Computational Etudes: A Spectral Approach"

    dump_path = parse_args(varargin{:});
    script_dir = fileparts(mfilename('fullpath'));
    addpath(fullfile(script_dir, '..', 'ch07'));
    addpath(script_dir);
    out_dir = fullfile(script_dir, '..', '..', '..', ...
                       'textbook', 'figures', 'ch18', 'matlab');
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end
    configure_style();
    [NAVY, CORAL] = colours();

    Ns = [12, 16, 24, 32, 48, 64, 96];
    betas = exact_clamped_beta(6);
    lam_exact = betas.^4;

    kappa_naive = zeros(size(Ns));
    kappa_hein  = zeros(size(Ns));
    drift_naive = zeros(size(Ns));
    drift_hein  = zeros(size(Ns));

    for k = 1:numel(Ns)
        [kn, lam_n] = cond_naive(Ns(k));
        [kh, lam_h] = cond_heinrichs(Ns(k));
        kappa_naive(k) = kn; kappa_hein(k) = kh;
        drift_naive(k) = abs(lam_n(1) - lam_exact(1));
        drift_hein(k)  = abs(lam_h(1) - lam_exact(1));
    end

    fig = figure('Units', 'inches', 'Position', [1, 1, 11.0, 4.5], 'Color', 'w');
    tl = tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile(tl);
    loglog(Ns, kappa_naive, 'o-', 'Color', NAVY, 'MarkerFaceColor', 'w', ...
        'MarkerSize', 6, 'LineWidth', 0.8, 'DisplayName', 'naive ($D^4$ bordered)');
    hold on;
    loglog(Ns, kappa_hein, 's-', 'Color', CORAL, 'MarkerFaceColor', 'w', ...
        'MarkerSize', 5, 'LineWidth', 0.8, 'DisplayName', 'Heinrichs $(1{-}x^2)^2 T_j$');
    loglog(Ns, kappa_naive(1) * (Ns / Ns(1)).^8, '--', 'Color', NAVY, ...
        'LineWidth', 0.6, 'DisplayName', '$N^8$');
    loglog(Ns, kappa_hein(1) * (Ns / Ns(1)).^4, '--', 'Color', CORAL, ...
        'LineWidth', 0.6, 'DisplayName', '$N^4$');
    hold off;
    xlabel('$N$', 'Interpreter', 'latex');
    ylabel('$\kappa$ (condition number)', 'Interpreter', 'latex');
    title('conditioning of the fourth-derivative matrix', 'Interpreter', 'latex');
    grid on; box on;
    legend('Location', 'northwest', 'FontSize', 9, 'Interpreter', 'latex');

    nexttile(tl);
    loglog(Ns, max(drift_naive, 1e-17), 'o-', 'Color', NAVY, ...
        'MarkerFaceColor', 'w', 'MarkerSize', 6, 'LineWidth', 0.8, ...
        'DisplayName', 'naive: $|\lambda_1 -$ exact$|$');
    hold on;
    loglog(Ns, max(drift_hein, 1e-17), 's-', 'Color', CORAL, ...
        'MarkerFaceColor', 'w', 'MarkerSize', 5, 'LineWidth', 0.8, ...
        'DisplayName', 'Heinrichs: $|\lambda_1 -$ exact$|$');
    hold off;
    xlabel('$N$', 'Interpreter', 'latex');
    ylabel('error in first eigenvalue', 'Interpreter', 'latex');
    title('first-eigenvalue error vs. $N$', 'Interpreter', 'latex');
    grid on; box on;
    legend('Location', 'northeast', 'FontSize', 9, 'Interpreter', 'latex');

    exportgraphics(fig, fullfile(out_dir, 'eigen_heinrichs_condition.pdf'), 'ContentType', 'vector');
    exportgraphics(fig, fullfile(out_dir, 'eigen_heinrichs_condition.png'), 'Resolution', 300);
    close(fig);

    fprintf('[Etude 18.8]  Heinrichs condition-number surgery\n');
    fprintf('  exact lambda_1..4 = %s\n', mat2str(lam_exact(1:4)'));
    fprintf('  N, kappa(naive), kappa(Heinrichs):\n');
    for k = 1:numel(Ns)
        fprintf('    N=%3d   kappa_naive=%.2e   kappa_Heinrichs=%.2e\n', ...
                Ns(k), kappa_naive(k), kappa_hein(k));
    end
    fprintf('  figure: %s\n', fullfile(out_dir, 'eigen_heinrichs_condition.pdf'));

    if ~isempty(dump_path)
        r = struct('Ns', Ns, 'kappa_naive', kappa_naive, ...
            'kappa_heinrichs', kappa_hein, ...
            'drift_naive', drift_naive, 'drift_heinrichs', drift_hein, ...
            'lam_exact_first4', lam_exact(1:4)');
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

function betas = exact_clamped_beta(n)
    f = @(b) cos(2*b) .* cosh(2*b) - 1;
    betas = [];
    for j = 1:n
        b_asym = (2*j + 1) * pi / 4;
        lo = b_asym - 0.5; hi = b_asym + 0.5;
        flo = f(lo); fhi = f(hi);
        if flo * fhi < 0
            for iter = 1:200
                mid = 0.5*(lo+hi); fmid = f(mid);
                if flo * fmid < 0; hi = mid; fhi = fmid;
                else; lo = mid; flo = fmid; end
                if hi - lo < 1e-14; break; end
            end
            betas(end+1) = 0.5*(lo+hi);
        end
    end
end

function [k, lam] = cond_naive(N)
    [A, B] = heinrichs_basis('clamped_naive', N);
    lam_all = eig(A, B);
    lam_all = lam_all(isfinite(lam_all));
    lam = sort(real(lam_all(abs(imag(lam_all)) < 1e-6)));
    lam = lam(lam > 1);
    if numel(lam) < 1; lam = nan(4,1); else; lam = lam(1:min(4,numel(lam))); end
    k = cond(A);
end

function [k, lam] = cond_heinrichs(N)
    [A, M, ~] = heinrichs_basis('clamped_heinrichs', N);
    lam_all = eig(A, M);
    lam_all = lam_all(isfinite(lam_all));
    lam = sort(real(lam_all(abs(imag(lam_all)) < 1e-6)));
    lam = lam(lam > 1);
    if numel(lam) < 1; lam = nan(4,1); else; lam = lam(1:min(4,numel(lam))); end
    A_std = M \ A;
    k = cond(A_std);
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
