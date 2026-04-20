function eigen_power_inverse(varargin)
%% eigen_power_inverse - Etude 18.9: power and inverse iteration.
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
% Part of "Computational Etudes: A Spectral Approach"

    dump_path = parse_args(varargin{:});
    script_dir = fileparts(mfilename('fullpath'));
    addpath(fullfile(script_dir, '..', 'ch07'));
    out_dir = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch18', 'matlab');
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end
    configure_style(); [NAVY, CORAL, TEAL] = colours();

    N = 32;
    A = build_matrix(N);
    lam_ref = sort(real(eig(A)));
    [hist_pow, ~] = power_method(A, 80, 0);
    lam_max = lam_ref(end);

    shifts = [5, 90, 250];
    hists = cell(1, numel(shifts)); targets = zeros(1, numel(shifts));
    for k = 1:numel(shifts)
        [hk, ~] = inverse_iteration(A, shifts(k), 30, 1);
        hists{k} = hk;
        [~, idx] = min(abs(lam_ref - hk(end)));
        targets(k) = lam_ref(idx);
    end

    lam_a = lam_ref(4); lam_b = lam_ref(5);
    mu_bad = 0.5 * (lam_a + lam_b);
    [hist_bad, ~] = inverse_iteration(A, mu_bad, 30, 1);

    fig = figure('Units', 'inches', 'Position', [1, 1, 13.5, 4.2], 'Color', 'w');
    tl = tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile(tl);
    semilogy(1:numel(hist_pow), max(abs(hist_pow - lam_max), 1e-18), 'o-', ...
        'Color', NAVY, 'MarkerFaceColor', 'w', 'MarkerSize', 4, 'LineWidth', 0.8);
    xlabel('iteration $k$', 'Interpreter', 'latex');
    ylabel('$|\mu^{(k)} - \lambda_{\max}|$', 'Interpreter', 'latex');
    title(sprintf('power method $\\to$ $\\lambda_{\\max} = %.3f$', lam_max), 'Interpreter', 'latex');
    ylim([1e-18, 1e4]); grid on; box on;

    nexttile(tl); hold on;
    cols = {NAVY, CORAL, TEAL};
    for k = 1:numel(shifts)
        err = max(abs(hists{k} - targets(k)), 1e-18);
        semilogy(1:numel(err), err, 'o-', 'Color', cols{k}, ...
            'MarkerFaceColor', 'w', 'MarkerSize', 4, 'LineWidth', 0.8, ...
            'DisplayName', sprintf('shift $\\mu=%g$ $\\to$ $\\lambda=%.3f$', shifts(k), targets(k)));
    end
    set(gca, 'YScale', 'log'); hold off;
    xlabel('iteration $k$', 'Interpreter', 'latex');
    ylabel('$|\mu^{(k)} - \lambda_{\rm target}|$', 'Interpreter', 'latex');
    title('inverse iteration: three shifts, three modes', 'Interpreter', 'latex');
    ylim([1e-18, 1e4]); grid on; box on;
    legend('Location', 'northeast', 'FontSize', 8, 'Interpreter', 'latex');

    nexttile(tl); hold on;
    err_a = max(abs(hist_bad - lam_a), 1e-18);
    err_b = max(abs(hist_bad - lam_b), 1e-18);
    semilogy(1:numel(hist_bad), err_a, 'o-', 'Color', NAVY, ...
        'MarkerFaceColor', 'w', 'MarkerSize', 4, 'LineWidth', 0.8, ...
        'DisplayName', sprintf('dist to $\\lambda=%.3f$', lam_a));
    semilogy(1:numel(hist_bad), err_b, 's-', 'Color', CORAL, ...
        'MarkerFaceColor', 'w', 'MarkerSize', 4, 'LineWidth', 0.8, ...
        'DisplayName', sprintf('dist to $\\lambda=%.3f$', lam_b));
    set(gca, 'YScale', 'log'); hold off;
    xlabel('iteration $k$', 'Interpreter', 'latex');
    ylabel('$|\mu^{(k)} - \lambda|$', 'Interpreter', 'latex');
    title(sprintf('cautionary: $\\mu=%.3f$ between two modes', mu_bad), 'Interpreter', 'latex');
    grid on; box on;
    legend('Location', 'northeast', 'FontSize', 8, 'Interpreter', 'latex');

    exportgraphics(fig, fullfile(out_dir, 'eigen_power_inverse.pdf'), 'ContentType', 'vector');
    exportgraphics(fig, fullfile(out_dir, 'eigen_power_inverse.png'), 'Resolution', 300);
    close(fig);

    fprintf('[Etude 18.9]  power and inverse iteration, N = %d\n', N);
    fprintf('  power method final -> %.6f (exact lam_max = %.6f)\n', hist_pow(end), lam_max);
    for k = 1:numel(shifts)
        fprintf('  inverse iter, shift %g  ->  lambda = %.6f\n', shifts(k), targets(k));
    end
    fprintf('  figure: %s\n', fullfile(out_dir, 'eigen_power_inverse.pdf'));

    if ~isempty(dump_path)
        r = struct('N', N, 'lam_max', lam_max, 'final_power_method', hist_pow(end), ...
            'shifts', shifts, 'targets', targets, 'bad_shift', mu_bad, ...
            'bad_pair', [lam_a, lam_b], 'final_bad', hist_bad(end));
        fid = fopen(dump_path, 'w'); fprintf(fid, '%s', jsonencode(r)); fclose(fid);
        fprintf('  dumped: %s\n', dump_path);
    end
end

function dump_path = parse_args(varargin)
    dump_path = ''; k = 1;
    while k <= numel(varargin)
        if strcmp(varargin{k}, '--dump'); dump_path = char(varargin{k+1}); k = k+2;
        else; k = k + 1; end
    end
end

function A = build_matrix(N)
    [D, ~] = cheb_matrix(N);
    D2 = D * D;
    A = -D2(2:N, 2:N);
end

function [hist, v] = power_method(A, max_iter, seed)
    rng(seed);
    n = size(A, 1);
    v = randn(n, 1); v = v / norm(v);
    hist = zeros(max_iter, 1);
    for k = 1:max_iter
        w = A * v;
        hist(k) = v' * w;
        v = w / norm(w);
    end
end

function [hist, v] = inverse_iteration(A, shift, max_iter, seed)
    rng(seed);
    n = size(A, 1);
    v = randn(n, 1); v = v / norm(v);
    hist = zeros(max_iter + 1, 1);
    [L, U, P] = lu(A - shift * eye(n));
    for k = 1:max_iter
        w = U \ (L \ (P * v));
        hist(k) = v' * (A * v);
        v = w / norm(w);
    end
    hist(end) = v' * (A * v);
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
