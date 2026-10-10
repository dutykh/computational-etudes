function eigen_smooth_lie(varargin)
%% eigen_smooth_lie - Computational Etude 18.1: The smooth lie.
%
% Chapter 18: Linear Spectral Eigenproblems.
%
% Solves the benchmark eigenvalue problem
%
%     u_xx + lambda u = 0,   u(+-1) = 0
%
% with Chebyshev pseudospectral collocation at N = 16. Exact eigenvalues
% are lambda_j = (j pi / 2)^2, eigenfunctions cos(j pi x/2) (j odd) and
% sin(j pi x / 2) (j even).
%
% Reproduces Boyd (2000) Fig 7.2: low modes recovered to near machine
% precision, high modes SMOOTH but quantitatively wrong.
%
% Usage:
%   eigen_smooth_lie                       % default N=16, no JSON dump
%   eigen_smooth_lie('--N', 32)            % override N
%   eigen_smooth_lie('--dump', '/tmp/x.json')  % also write JSON
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
% Part of "Computational Etudes: A Spectral Approach"
% https://github.com/dutykh/computational-etudes

    [N, dump_path] = parse_args(varargin{:});

    script_dir = fileparts(mfilename('fullpath'));
    addpath(fullfile(script_dir, '..', 'ch07'));
    out_dir = fullfile(script_dir, '..', '..', '..', ...
                       'textbook', 'figures', 'ch18', 'matlab');
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end

    configure_style();
    [NAVY, ~, CORAL] = colours();  %#ok<ASGLU>

    low_mode  = 3;
    high_mode = 15;
    [lam, V, x] = solve_problem(N);

    idx_low  = low_mode;
    idx_high = high_mode;

    u_low_num  = normalise_max([0; V(:, idx_low); 0]);
    u_high_num = normalise_max([0; V(:, idx_high); 0]);

    xfine = linspace(-1, 1, 801)';
    u_low_ex  = normalise_max(exact_eigenfunction(xfine, low_mode));
    u_high_ex = normalise_max(exact_eigenfunction(xfine, high_mode));

    lam_num_low   = lam(idx_low);
    lam_num_high  = lam(idx_high);
    lam_ex_low    = exact_eigenvalue(low_mode);
    lam_ex_high   = exact_eigenvalue(high_mode);

    fig = figure('Units', 'inches', 'Position', [1, 1, 10.5, 4.2], 'Color', 'w');
    tl = tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile(tl);
    hold on;
    plot(xfine, u_low_ex, 'Color', NAVY, 'LineWidth', 1.3, ...
         'DisplayName', sprintf('exact mode %d', low_mode));
    plot(x, u_low_num, 'o', 'Color', CORAL, 'MarkerFaceColor', 'w', ...
         'MarkerSize', 6, 'LineWidth', 1.1, ...
         'DisplayName', sprintf('numerical (N=%d)', N));
    hold off;
    title(sprintf('mode $j=%d$: $\\lambda_{\\mathrm{num}}=%.4f$, $\\lambda_{\\mathrm{exact}}=%.4f$', ...
                  low_mode, lam_num_low, lam_ex_low), ...
          'Interpreter', 'latex', 'FontSize', 10);
    xlabel('$x$', 'Interpreter', 'latex');
    ylabel('$u(x)$', 'Interpreter', 'latex');
    xlim([-1.05, 1.05]); ylim([-1.15, 1.15]);
    legend('Location', 'south', 'FontSize', 9, 'Interpreter', 'latex');
    grid on; box on;

    nexttile(tl);
    hold on;
    plot(xfine, u_high_ex, 'Color', NAVY, 'LineWidth', 1.3, ...
         'DisplayName', sprintf('exact mode %d', high_mode));
    plot(x, u_high_num, 'o-', 'Color', CORAL, 'MarkerFaceColor', 'w', ...
         'MarkerSize', 6, 'LineWidth', 0.8, ...
         'DisplayName', sprintf('numerical (N=%d)', N));
    hold off;
    title(sprintf('mode $j=%d$: $\\lambda_{\\mathrm{num}}=%.4f$, $\\lambda_{\\mathrm{exact}}=%.4f$', ...
                  high_mode, lam_num_high, lam_ex_high), ...
          'Interpreter', 'latex', 'FontSize', 10);
    xlabel('$x$', 'Interpreter', 'latex');
    ylabel('$u(x)$', 'Interpreter', 'latex');
    xlim([-1.05, 1.05]); ylim([-1.15, 1.15]);
    legend('Location', 'south', 'FontSize', 9, 'Interpreter', 'latex');
    grid on; box on;

    exportgraphics(fig, fullfile(out_dir, 'eigen_smooth_lie.pdf'), ...
                   'ContentType', 'vector');
    exportgraphics(fig, fullfile(out_dir, 'eigen_smooth_lie.png'), ...
                   'Resolution', 300);
    close(fig);

    fprintf('[Etude 18.1]  N = %d\n', N);
    fprintf('  mode %2d:  lambda_num = %.6e   exact = %.6e   |err| = %.2e\n', ...
            low_mode,  lam_num_low,  lam_ex_low,  abs(lam_num_low  - lam_ex_low));
    fprintf('  mode %2d:  lambda_num = %.6e   exact = %.6e   |err| = %.2e\n', ...
            high_mode, lam_num_high, lam_ex_high, abs(lam_num_high - lam_ex_high));
    fprintf('  figure: %s\n', fullfile(out_dir, 'eigen_smooth_lie.pdf'));

    if ~isempty(dump_path)
        results = struct( ...
            'N',               N, ...
            'low_mode',        low_mode, ...
            'high_mode',       high_mode, ...
            'eigvals_sorted',  lam(:)', ...
            'lam_num_low',     lam_num_low, ...
            'lam_num_high',    lam_num_high, ...
            'lam_exact_low',   lam_ex_low, ...
            'lam_exact_high',  lam_ex_high, ...
            'abs_err_low',     abs(lam_num_low  - lam_ex_low), ...
            'abs_err_high',    abs(lam_num_high - lam_ex_high));
        fid = fopen(dump_path, 'w');
        fprintf(fid, '%s', jsonencode(results));
        fclose(fid);
        fprintf('  dumped: %s\n', dump_path);
    end
end


function [N, dump_path] = parse_args(varargin)
    N = 16;
    dump_path = '';
    k = 1;
    while k <= numel(varargin)
        switch varargin{k}
            case '--N'
                N = varargin{k+1};
                if ischar(N) || isstring(N); N = str2double(N); end
                k = k + 2;
            case '--dump'
                dump_path = char(varargin{k+1});
                k = k + 2;
            otherwise
                k = k + 1;
        end
    end
end


function [lam, V, x] = solve_problem(N)
    [D, x] = cheb_matrix(N);
    D2 = D * D;
    A = -D2(2:N, 2:N);
    [V, Lam] = eig(A);
    lam_all = diag(Lam);
    [~, order] = sort(real(lam_all));
    lam = real(lam_all(order));
    V = real(V(:, order));
    % sign-fix: first interior component positive
    for k = 1:size(V, 2)
        if V(1, k) < 0
            V(:, k) = -V(:, k);
        end
    end
end


function y = exact_eigenvalue(j)
    y = (j * pi / 2)^2;
end


function y = exact_eigenfunction(x, j)
    if mod(j, 2) == 1
        y = cos(j * pi * x / 2);
    else
        y = sin(j * pi * x / 2);
    end
end


function u = normalise_max(u)
    m = max(abs(u));
    if m == 0; return; end
    u = u / m;
    [~, idx] = max(abs(u));
    if u(idx) < 0
        u = -u;
    end
end


function configure_style()
    set(groot, 'defaultAxesFontName', 'CMU Serif');
    set(groot, 'defaultTextFontName', 'CMU Serif');
    set(groot, 'defaultAxesFontSize', 11);
    set(groot, 'defaultAxesLineWidth', 0.8);
    set(groot, 'defaultLineLineWidth', 1.0);
    set(groot, 'defaultAxesBox', 'on');
    set(groot, 'defaultAxesTickDir', 'in');
end


function [NAVY, SKY, CORAL] = colours()
    NAVY  = [0x14, 0x2D, 0x6E] / 255;
    SKY   = [0x78, 0x96, 0xD2] / 255;
    CORAL = [0xE7, 0x4C, 0x3C] / 255;
end
