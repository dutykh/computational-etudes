function lanczos_economization(varargin)
%% lanczos_economization - Etude 21.6: Chebyshev surrogate for an
% expensive determinant.  Compare naive bracketing+bisection to a 17-pt
% Chebyshev surrogate root-found via the colleague matrix.
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
% Part of "Computational Etudes: A Spectral Approach"

    dump_path = parse_args(varargin{:});
    script_dir = fileparts(mfilename('fullpath'));
    addpath(script_dir);
    cm = tricks_common(); cm.configure_style();
    out_dir = cm.output_dir(script_dir);

    M_size = 20;
    SHIFTS = (1:M_size).' .^ 2;
    rng(0);
    Q = qr(randn(M_size));
    diag_v = 0.05 * randn(M_size, 1);
    T2 = Q * diag(diag_v) * Q.';

    expensive_D = @(lam) det(T2 + diag(lam - SHIFTS));

    a = 0.5; b = 30.0;

    K_dense = 60;
    [roots_A, ~, ~] = roots_naive(K_dense, a, b, expensive_D);
    n_LUs_A = K_dense + 30 * numel(roots_A);

    K = 17;
    [roots_B, lam_B, D_B] = roots_chebyshev(K, a, b, expensive_D);
    n_LUs_B = K;

    seeds = SHIFTS(SHIFTS >= a & SHIFTS <= b);
    exact_roots = zeros(size(seeds));
    for s = 1:numel(seeds)
        x = double(seeds(s));
        for it = 1:80
            d  = expensive_D(x);
            dp = (expensive_D(x + 1e-7) - expensive_D(x - 1e-7)) / 2e-7;
            if abs(dp) < 1e-30; break; end
            x = x - d / dp;
        end
        exact_roots(s) = x;
    end
    exact_roots = sort(exact_roots);

    err_A = arrayfun(@(r) min(abs(r - exact_roots)), sort(roots_A));
    err_B = arrayfun(@(r) min(abs(r - exact_roots)), sort(roots_B));

    fig = figure('Units', 'inches', 'Position', [1, 1, 11.5, 4.4], 'Color', 'w');
    tl = tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile(tl); hold on;
    lam_dense = linspace(a, b, 600);
    D_dense = arrayfun(expensive_D, lam_dense);
    plot(lam_dense, D_dense, 'Color', cm.NAVY, 'LineWidth', 1.0, ...
         'DisplayName', '$D(\lambda)=\det A(\lambda)$');
    yline(0, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.4, 'Alpha', 0.5);
    scatter(roots_A, zeros(size(roots_A)), 80, cm.CORAL, 'x', 'LineWidth', 1.4, ...
            'DisplayName', sprintf('naive: %d roots (%d LUs)', numel(roots_A), n_LUs_A));
    scatter(roots_B, zeros(size(roots_B)) + 0.05 * max(abs(D_dense)), 80, ...
            cm.TEAL, 'o', 'LineWidth', 1.4, ...
            'DisplayName', sprintf('Lanczos: %d roots (%d LUs)', numel(roots_B), n_LUs_B));
    scatter(lam_B, D_B, 20, cm.TEAL, 'filled', ...
            'DisplayName', sprintf('Lanczos sample points ($K=%d$)', K));
    hold off;
    xlabel('$\lambda$', 'Interpreter', 'latex');
    ylabel('$D(\lambda)$', 'Interpreter', 'latex');
    title('the determinant and the two strategies'' detected roots', ...
          'Interpreter', 'latex');
    legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 8);
    grid on; box on;
    xlim([a - 1, b + 1]);

    nexttile(tl);
    counts = [n_LUs_A, n_LUs_B];
    accs   = [max(err_A), max(err_B)];
    bar([1 2], counts, 0.5, 'FaceColor', 'flat', 'EdgeColor', 'k', 'LineWidth', 0.6, ...
        'CData', [cm.CORAL; cm.TEAL]);
    set(gca, 'XTick', [1 2], 'XTickLabel', {sprintf('naive\n+ 30 bisects'), ...
                                             sprintf('Lanczos\nK=%d', K)});
    ylabel('# expensive LU calls');
    title(sprintf('cost: %d vs %d LU factorisations', n_LUs_B, n_LUs_A));
    text(1, counts(1) + 5, sprintf('%d LUs', counts(1)), 'HorizontalAlignment', 'center', 'FontSize', 9);
    text(2, counts(2) + 5, sprintf('%d LUs', counts(2)), 'HorizontalAlignment', 'center', 'FontSize', 9);
    text(1, counts(1) * 0.5, sprintf('max root err\n%.1e', accs(1)), ...
         'HorizontalAlignment', 'center', 'FontSize', 8.5, 'Color', 'w', 'FontWeight', 'bold');
    text(2, counts(2) * 0.5, sprintf('max root err\n%.1e', accs(2)), ...
         'HorizontalAlignment', 'center', 'FontSize', 8.5, 'Color', 'w', 'FontWeight', 'bold');
    ylim([0, max(counts) * 1.15]);

    exportgraphics(fig, fullfile(out_dir, 'lanczos_economization.pdf'), 'ContentType', 'vector');
    exportgraphics(fig, fullfile(out_dir, 'lanczos_economization.png'), 'Resolution', 300);
    close(fig);

    fprintf('[Etude 21.6]  Lanczos economization for an expensive determinant\n');
    fprintf('  exact roots in [%.1f, %.1f]: ', a, b);
    fprintf('%.6f ', exact_roots); fprintf('\n');
    fprintf('  naive  : %d roots, %d LUs, max err = %.3e\n', numel(roots_A), n_LUs_A, max(err_A));
    fprintf('  Lanczos: %d roots, %d LUs, max err = %.3e\n', numel(roots_B), n_LUs_B, max(err_B));
    fprintf('  figure: %s\n', fullfile(out_dir, 'lanczos_economization.pdf'));

    if ~isempty(dump_path)
        r = struct('exact_roots', exact_roots.', 'roots_naive', sort(roots_A).', ...
                   'roots_lanczos', sort(roots_B).', ...
                   'n_LUs_naive', n_LUs_A, 'n_LUs_lanczos', n_LUs_B, ...
                   'max_err_naive', max(err_A), 'max_err_lanczos', max(err_B));
        fid = fopen(dump_path, 'w'); fprintf(fid, '%s', jsonencode(r)); fclose(fid);
    end
end

function dump_path = parse_args(varargin)
    dump_path = ''; k = 1;
    while k <= numel(varargin)
        if strcmp(varargin{k}, '--dump'); dump_path = char(varargin{k+1}); k = k + 2;
        else; k = k + 1; end
    end
end

function [roots_out, lam, D_samples] = roots_chebyshev(K, a, b, fn)
    j = 0:K-1;
    t = cos(pi * j / (K - 1));
    lam = 0.5 * (a + b) + 0.5 * (b - a) * t;
    D_samples = arrayfun(fn, lam);
    ext = [D_samples, fliplr(D_samples(2:K-1))];
    A = real(fft(ext)) / (K - 1);
    A(1) = 0.5 * A(1);
    A(K) = 0.5 * A(K);
    coeffs = A(1:K);
    n = K - 1;
    C = zeros(n, n);
    for i = 1:n-1
        C(i, i + 1) = 0.5;
        C(i + 1, i) = 0.5;
    end
    C(1, 2) = 1.0;
    last = -coeffs(1:n) / (2.0 * coeffs(n + 1));
    last(n - 1) = last(n - 1) + 0.5;
    C(n, :) = last;
    eigs_C = eig(C);
    real_in = sort(real(eigs_C(abs(imag(eigs_C)) < 1e-6 & ...
                              -1.0 <= real(eigs_C) & real(eigs_C) <= 1.0)));
    roots_out = 0.5 * (a + b) + 0.5 * (b - a) * real_in;
end

function [roots_out, lam, D] = roots_naive(K_dense, a, b, fn)
    lam = linspace(a, b, K_dense);
    D = arrayfun(fn, lam);
    roots_out = [];
    for k = 1:K_dense - 1
        if D(k) * D(k + 1) < 0
            l = lam(k); r = lam(k + 1); fl = D(k);
            for it = 1:30
                mid = 0.5 * (l + r);
                fm = fn(mid);
                if fl * fm <= 0
                    r = mid;
                else
                    l = mid; fl = fm;
                end
            end
            roots_out = [roots_out; 0.5 * (l + r)];
        end
    end
end
