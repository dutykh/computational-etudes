function eigen_two_formulations(varargin)
%% eigen_two_formulations - Computational Etude 18.2: Assemble the pencil.
%
% Same problem as Etude 18.1 (u_xx + lambda u = 0, u(+-1) = 0) solved two
% ways:
%   (A) Boundary bordering: (N+1) x (N+1) pencil A v = lambda B v;
%   (B) Basis recombination: interior -D^2_int u = lambda u.
%
% Both recover the same finite spectrum.
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
    [NAVY, CORAL, TEAL] = colours();

    [A_pen, B_pen, ~] = boundary_bordering(N);
    A_reg = basis_recombination(N);

    lam_pen = solve_pencil(A_pen, B_pen);
    lam_reg = solve_regular(A_reg);
    m = min(numel(lam_pen), numel(lam_reg));
    lam_pen = lam_pen(1:m); lam_reg = lam_reg(1:m);

    j = (1:m)';
    lam_exact = (j * pi / 2).^2;
    err_pen = abs(lam_pen - lam_exact);
    err_reg = abs(lam_reg - lam_exact);
    diff_pr = abs(lam_pen - lam_reg);

    fig = figure('Units', 'inches', 'Position', [1, 1, 10.5, 4.2], 'Color', 'w');
    tl = tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile(tl);
    hold on;
    semilogy(j, max(err_pen, 1e-17), 'o', 'Color', NAVY, ...
        'MarkerFaceColor', 'w', 'MarkerSize', 6, 'LineWidth', 1.1, ...
        'DisplayName', 'boundary bordering (pencil)');
    semilogy(j, max(err_reg, 1e-17), 's', 'Color', CORAL, ...
        'MarkerFaceColor', 'w', 'MarkerSize', 5, 'LineWidth', 1.1, ...
        'DisplayName', 'basis recombination (regular)');
    semilogy(j, max(diff_pr, 1e-17), 'x', 'Color', TEAL, ...
        'MarkerSize', 7, 'LineWidth', 1.2, ...
        'DisplayName', '|pencil $-$ regular|');
    yline(1e-2, '--k', 'LineWidth', 0.6);
    set(gca, 'YScale', 'log');
    hold off;
    xlabel('mode number $j$', 'Interpreter', 'latex');
    ylabel('$|\lambda^{\mathrm{num}}_j - \lambda^{\mathrm{exact}}_j|$', ...
        'Interpreter', 'latex');
    title(sprintf('both formulations, $N = %d$', N), 'Interpreter', 'latex');
    ylim([1e-17, 1e5]);
    grid on; box on;
    legend('Location', 'southeast', 'FontSize', 8, 'Interpreter', 'latex');

    nexttile(tl);
    logA = log10(max(abs(A_pen), 1e-20));
    imagesc(logA); colormap(gca, flipud(bone));
    set(gca, 'YDir', 'reverse');
    axis equal tight;
    title('$\log_{10}|A|$ for boundary-bordered pencil', 'Interpreter', 'latex');
    xlabel('column index'); ylabel('row index');
    colorbar;

    exportgraphics(fig, fullfile(out_dir, 'eigen_two_formulations.pdf'), ...
                   'ContentType', 'vector');
    exportgraphics(fig, fullfile(out_dir, 'eigen_two_formulations.png'), ...
                   'Resolution', 300);
    close(fig);

    fprintf('[Etude 18.2]  N = %d\n', N);
    fprintf('  max |lam_pencil - lam_regular| = %.3e\n', max(diff_pr));
    fprintf('  max |lam_pencil  - exact|      = %.3e\n', max(err_pen));
    fprintf('  max |lam_regular - exact|      = %.3e\n', max(err_reg));
    fprintf('  figure: %s\n', fullfile(out_dir, 'eigen_two_formulations.pdf'));

    if ~isempty(dump_path)
        r = struct('N', N, ...
            'lam_pencil',        lam_pen(:)', ...
            'lam_regular',       lam_reg(:)', ...
            'lam_exact',         lam_exact(:)', ...
            'max_abs_agreement', max(diff_pr), ...
            'max_err_pencil',    max(err_pen), ...
            'max_err_regular',   max(err_reg));
        fid = fopen(dump_path, 'w');
        fprintf(fid, '%s', jsonencode(r));
        fclose(fid);
        fprintf('  dumped: %s\n', dump_path);
    end
end

function [N, dump_path] = parse_args(varargin)
    N = 16; dump_path = '';
    k = 1;
    while k <= numel(varargin)
        switch varargin{k}
            case '--N'
                N = varargin{k+1};
                if ischar(N) || isstring(N); N = str2double(N); end
                k = k + 2;
            case '--dump'
                dump_path = char(varargin{k+1}); k = k + 2;
            otherwise
                k = k + 1;
        end
    end
end

function [A, B, x] = boundary_bordering(N)
    [D, x] = cheb_matrix(N);
    D2 = D * D;
    A = -D2;
    B = eye(N + 1);
    A(1,   :) = 0; A(1,   1)   = 1; B(1,   :) = 0;
    A(N+1, :) = 0; A(N+1, N+1) = 1; B(N+1, :) = 0;
end

function A_int = basis_recombination(N)
    [D, ~] = cheb_matrix(N);
    D2 = D * D;
    A_int = -D2(2:N, 2:N);
end

function lam = solve_pencil(A, B)
    lam_all = eig(A, B);
    keep = isfinite(lam_all) & abs(imag(lam_all)) < 1e-6;
    lam = sort(real(lam_all(keep)));
end

function lam = solve_regular(A)
    lam_all = eig(A);
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
