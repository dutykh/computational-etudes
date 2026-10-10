function rescue_naive_vs_tailored(varargin)
%% rescue_naive_vs_tailored - Etude 21.1: a rescue story in miniature.
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
% Part of "Computational Etudes: A Spectral Approach"
%
% f(y) = sech(y) is approximated on [-L, L] by Chebyshev collocation.
% Same N = 48, but two truncation lengths:
%   naive    L = 30 -- 'wide enough to be safe', wastes resolution
%   tailored L =  8 -- just enough room for the support of f
% The naive choice yields ~10% error; the tailored choice ~1e-4.

    dump_path = parse_args(varargin{:});
    script_dir = fileparts(mfilename('fullpath'));
    addpath(script_dir);
    cm = tricks_common(); cm.configure_style();
    out_dir = cm.output_dir(script_dir);

    N = 48;
    L_naive = 30.0;
    L_tail  = 8.0;

    y_ref = linspace(-35.0, 35.0, 6001);
    [a_naive, err_naive] = trunc_cheb_err(N, L_naive, y_ref);
    [a_tail,  err_tail]  = trunc_cheb_err(N, L_tail,  y_ref);

    inside_naive = abs(y_ref) <= L_naive;
    inside_tail  = abs(y_ref) <= L_tail;
    e_naive = max(err_naive(inside_naive), [], 'omitnan');
    e_tail  = max(err_tail(inside_tail),   [], 'omitnan');

    fig = figure('Units', 'inches', 'Position', [1, 1, 11.5, 4.4], 'Color', 'w');
    tl = tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    % --- Left panel ----------------------------------------------------
    nexttile(tl); hold on;
    plot(y_ref, sech(y_ref), 'Color', cm.NAVY, 'LineWidth', 1.4, ...
         'DisplayName', '$f(y) = \mathrm{sech}\,y$');
    y_naive_grid = L_naive * cgl(N);
    scatter(y_naive_grid, sech(y_naive_grid), 22, 'o', ...
            'MarkerEdgeColor', cm.CORAL, 'MarkerFaceColor', 'w', ...
            'LineWidth', 1.0, ...
            'DisplayName', sprintf('naive ($L = %d$)', round(L_naive)));
    y_tail_grid = L_tail * cgl(N);
    scatter(y_tail_grid, sech(y_tail_grid), 28, '^', ...
            'MarkerEdgeColor', cm.TEAL, 'MarkerFaceColor', 'w', ...
            'LineWidth', 1.0, ...
            'DisplayName', sprintf('tailored ($L = %d$)', round(L_tail)));
    xline(0, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.4, 'Alpha', 0.4);
    xlabel('$y$', 'Interpreter', 'latex');
    ylabel('$f(y)$', 'Interpreter', 'latex');
    title(sprintf('same $f$, same $N=%d$, two truncation lengths $L$', N), ...
          'Interpreter', 'latex');
    legend('Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 9);
    grid on; box on;
    xlim([-22 22]); ylim([-0.06 1.12]); hold off;

    % --- Right panel ---------------------------------------------------
    nexttile(tl); hold on;
    semilogy(y_ref(inside_naive), err_naive(inside_naive) + 1e-18, ...
             'Color', cm.CORAL, 'LineWidth', 1.0, ...
             'DisplayName', sprintf('naive ($L=%d$): max err $\\approx$%.1e', round(L_naive), e_naive));
    semilogy(y_ref(inside_tail),  err_tail(inside_tail)  + 1e-18, ...
             'Color', cm.TEAL, 'LineWidth', 1.0, ...
             'DisplayName', sprintf('tailored ($L=%d$): max err $\\approx$%.1e', round(L_tail), e_tail));
    set(gca, 'YScale', 'log');
    yline(1e-15, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.4, 'Alpha', 0.5, ...
          'HandleVisibility', 'off');
    xlabel('$y$', 'Interpreter', 'latex');
    ylabel('$|f(y) - f_N(y)|$', 'Interpreter', 'latex');
    title('pointwise error: a single tuned $L$ buys decimals', 'Interpreter', 'latex');
    legend('Location', 'south', 'Interpreter', 'latex', 'FontSize', 9);
    grid on; box on;
    xlim([-22 22]); ylim([1e-16 1e0]); hold off;

    exportgraphics(fig, fullfile(out_dir, 'rescue_naive_vs_tailored.pdf'), ...
                   'ContentType', 'vector');
    exportgraphics(fig, fullfile(out_dir, 'rescue_naive_vs_tailored.png'), ...
                   'Resolution', 300);
    close(fig);

    fprintf('[Etude 21.1]  rescue story\n');
    fprintf('  N = %d, f(y) = sech(y), |f(L)| = sech(L) ~ truncation noise floor\n', N);
    fprintf('  naive    L = %5.2f: max err = %.3e\n', L_naive, e_naive);
    fprintf('  tailored L = %5.2f: max err = %.3e\n', L_tail,  e_tail);
    fprintf('  ratio (tailored / naive) = %.2e\n', e_tail / e_naive);
    fprintf('  figure: %s\n', fullfile(out_dir, 'rescue_naive_vs_tailored.pdf'));

    if ~isempty(dump_path)
        r = struct('N', N, 'L_naive', L_naive, 'L_tailored', L_tail, ...
                   'err_naive', e_naive, 'err_tailored', e_tail);
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

function x = cgl(N)
    x = cos(pi * (0:N) / N);
end

function a = cheb_coeffs(v)
    N = numel(v) - 1;
    ext = [v(:); v(N:-1:2)']';
    A = real(fft(ext)) / N;
    A(1) = 0.5 * A(1);
    A(N+1) = 0.5 * A(N+1);
    a = A(1:N+1);
end

function val = clenshaw(a, t)
    N = numel(a) - 1;
    T0 = ones(size(t));
    if N == 0; val = a(1) * T0; return; end
    T1 = t;
    val = a(1) * T0 + a(2) * T1;
    for n = 2:N
        Tk = 2.0 * t .* T1 - T0;
        val = val + a(n+1) * Tk;
        T0 = T1; T1 = Tk;
    end
end

function [a, err] = trunc_cheb_err(N, L, y_ref)
    t_grid = cgl(N);
    y_grid = L * t_grid;
    a = cheb_coeffs(sech(y_grid));
    t_ref = y_ref / L;
    inside = abs(t_ref) <= 1.0;
    approx = nan(size(y_ref));
    approx(inside) = clenshaw(a, t_ref(inside));
    err = abs(approx - sech(y_ref));
    err(~inside) = NaN;
end
