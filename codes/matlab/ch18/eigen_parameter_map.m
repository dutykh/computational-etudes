function eigen_parameter_map(varargin)
%% eigen_parameter_map - Etude 18.10: map the neutral curve safely.
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
% Part of "Computational Etudes: A Spectral Approach"

    dump_path = parse_args(varargin{:});
    script_dir = fileparts(mfilename('fullpath'));
    addpath(fullfile(script_dir, '..', 'ch07'));
    out_dir = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch18', 'matlab');
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end
    configure_style(); [NAVY, CORAL] = colours();

    N = 40;
    alpha_coarse = linspace(-2.0, 2.0, 9);
    beta_coarse  = linspace(0.0, 8.0, 9);
    lambda1 = zeros(numel(alpha_coarse), numel(beta_coarse));
    for i = 1:numel(alpha_coarse)
        for j = 1:numel(beta_coarse)
            [lambda1(i,j), ~] = qr_lowest_two(N, alpha_coarse(i), beta_coarse(j));
        end
    end

    beta_fine = linspace(0.0, 8.0, 40);
    lam1_safe = zeros(size(beta_fine));
    for k = 1:numel(beta_fine)
        [lam1_safe(k), ~] = qr_lowest_two(N, 0.0, beta_fine(k));
    end

    [~, l2_0] = qr_lowest_two(N, 0.0, beta_fine(1));
    shift = l2_0;
    lam1_naive = zeros(size(beta_fine));
    lam1_naive(1) = l2_0;
    for k = 2:numel(beta_fine)
        lam = inv_iter_shift(N, 0.0, beta_fine(k), shift, 20);
        lam1_naive(k) = lam; shift = lam;
    end

    fig = figure('Units', 'inches', 'Position', [1, 1, 11.5, 4.5], 'Color', 'w');
    tl = tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile(tl);
    [AA, BB] = meshgrid(alpha_coarse, beta_coarse);
    % lambda1 is indexed (alpha, beta), meshgrid returns (beta, alpha), so transpose
    contourf(AA, BB, lambda1', 18, 'LineColor', 'none');
    hold on;
    contour(AA, BB, lambda1', 8, 'LineColor', 'w', 'LineWidth', 0.5);
    scatter(AA(:), BB(:), 14, 'w', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.6);
    hold off;
    colormap(gca, viridis_fallback()); colorbar;
    xlabel('$\alpha$', 'Interpreter', 'latex');
    ylabel('$\beta$', 'Interpreter', 'latex');
    title('coarse QR grid (9$\times$9): lowest eigenvalue', 'Interpreter', 'latex');

    nexttile(tl); hold on;
    plot(beta_fine, lam1_safe, 'o-', 'Color', NAVY, 'MarkerFaceColor', 'w', ...
        'MarkerSize', 4, 'LineWidth', 0.8, ...
        'DisplayName', 'SAFE: QR re-anchor each step');
    plot(beta_fine, lam1_naive, 'x-', 'Color', CORAL, 'MarkerSize', 6, ...
        'LineWidth', 0.8, ...
        'DisplayName', 'NAIVE: inverse iter seeded from wrong branch');
    hold off;
    xlabel('$\beta$ (with $\alpha = 0$)', 'Interpreter', 'latex');
    ylabel('eigenvalue', 'Interpreter', 'latex');
    title('eigenvalue tracking: safe vs. naive local iteration', 'Interpreter', 'latex');
    grid on; box on;
    legend('Location', 'southeast', 'FontSize', 9, 'Interpreter', 'latex');

    exportgraphics(fig, fullfile(out_dir, 'eigen_parameter_map.pdf'), 'ContentType', 'vector');
    exportgraphics(fig, fullfile(out_dir, 'eigen_parameter_map.png'), 'Resolution', 300);
    close(fig);

    fprintf('[Etude 18.10]  parameter map\n');
    fprintf('  coarse grid: %d x %d = %d QR solves\n', ...
            numel(alpha_coarse), numel(beta_coarse), numel(alpha_coarse)*numel(beta_coarse));
    fprintf('  safe tracker: starts at %.4f, ends at %.4f\n', lam1_safe(1), lam1_safe(end));
    fprintf('  naive tracker: starts at %.4f, ends at %.4f\n', lam1_naive(1), lam1_naive(end));
    fprintf('  figure: %s\n', fullfile(out_dir, 'eigen_parameter_map.pdf'));

    if ~isempty(dump_path)
        r = struct('N', N, 'alpha_grid', alpha_coarse, 'beta_grid', beta_coarse, ...
            'lambda1_coarse', lambda1, 'beta_fine', beta_fine, ...
            'lam1_fine_safe', lam1_safe, 'lam1_fine_naive', lam1_naive);
        fid = fopen(dump_path, 'w'); fprintf(fid, '%s', jsonencode(r)); fclose(fid);
        fprintf('  dumped: %s\n', dump_path);
    end
end

function dump_path = parse_args(varargin)
    dump_path = ''; k = 1;
    while k <= numel(varargin)
        if strcmp(varargin{k}, '--dump'); dump_path = char(varargin{k+1}); k = k + 2;
        else; k = k + 1; end
    end
end

function [l1, l2] = qr_lowest_two(N, alpha, beta)
    [D, x] = cheb_matrix(N);
    D2 = D * D;
    A = -D2(2:N, 2:N) + diag(alpha + beta * x(2:N).^2);
    lam = sort(real(eig(A)));
    l1 = lam(1); l2 = lam(2);
end

function lam = inv_iter_shift(N, alpha, beta, shift, max_iter)
    [D, x] = cheb_matrix(N);
    D2 = D * D;
    A = -D2(2:N, 2:N) + diag(alpha + beta * x(2:N).^2);
    n = size(A, 1);
    rng(0);
    v = randn(n, 1); v = v / norm(v);
    [L, U, P] = lu(A - shift * eye(n));
    lam = 0.0;
    for k = 1:max_iter
        w = U \ (L \ (P * v));
        v = w / norm(w);
        lam = v' * (A * v);
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

function [NAVY, CORAL] = colours()
    NAVY  = [0x14, 0x2D, 0x6E] / 255;
    CORAL = [0xE7, 0x4C, 0x3C] / 255;
end

function cmap = viridis_fallback()
%% Simple approximation of viridis for older MATLAB releases without it.
    try
        cmap = viridis(256);
    catch
        cmap = parula(256);
    end
end
