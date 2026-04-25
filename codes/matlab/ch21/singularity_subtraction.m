function singularity_subtraction(varargin)
%% singularity_subtraction - Etude 21.3: subtract a known singular part.
% f(x) = e^x + 0.1 (1+x)^(2/3) on [-1, 1].  Direct Chebyshev expansion of
% f decays algebraically; subtracting the singular part leaves e^x, which
% has geometric decay.
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
% Part of "Computational Etudes: A Spectral Approach"

    dump_path = parse_args(varargin{:});
    script_dir = fileparts(mfilename('fullpath'));
    addpath(script_dir);
    cm = tricks_common(); cm.configure_style();
    out_dir = cm.output_dir(script_dir);

    csing = 0.1;
    f_full     = @(x) exp(x) + csing * (1.0 + x).^(2/3);
    f_singular = @(x) csing * (1.0 + x).^(2/3);
    f_smooth   = @(x) exp(x);

    N_show = 80;
    x_show = cgl(N_show);
    a_naive = cheb_coeffs(f_full(x_show));
    a_trick = cheb_coeffs(f_smooth(x_show));

    Ns = [8 12 16 24 32 48 64 96 128 192 256];
    err_naive = zeros(size(Ns));
    err_trick = zeros(size(Ns));
    x_eval = linspace(-1.0, 1.0, 5001);
    for k = 1:numel(Ns)
        N = Ns(k);
        xg = cgl(N);
        a_n = cheb_coeffs(f_full(xg));
        a_t = cheb_coeffs(f_smooth(xg));
        approx_n = clenshaw(a_n, x_eval);
        approx_t = clenshaw(a_t, x_eval) + f_singular(x_eval);
        err_naive(k) = max(abs(approx_n - f_full(x_eval)));
        err_trick(k) = max(abs(approx_t - f_full(x_eval)));
    end

    fig = figure('Units', 'inches', 'Position', [1, 1, 11.5, 4.4], 'Color', 'w');
    tl = tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile(tl); hold on;
    n_axis = 0:N_show;
    semilogy(n_axis, abs(a_naive) + 1e-20, 'o-', 'Color', cm.CORAL, ...
             'MarkerFaceColor', 'w', 'MarkerSize', 4, 'LineWidth', 0.8, ...
             'DisplayName', 'naive: $|a_n|$ for $f(x)$');
    semilogy(n_axis, abs(a_trick) + 1e-20, '^-', 'Color', cm.TEAL, ...
             'MarkerFaceColor', 'w', 'MarkerSize', 4, 'LineWidth', 0.8, ...
             'DisplayName', 'trick: $|a_n|$ for $f - c(1+x)^{2/3}$');
    set(gca, 'YScale', 'log');
    hold off;
    xlabel('Chebyshev degree $n$', 'Interpreter', 'latex');
    ylabel('$|a_n|$', 'Interpreter', 'latex');
    title(sprintf('coefficient decay at $N = %d$', N_show), 'Interpreter', 'latex');
    ylim([1e-18, 5]);
    legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 9);
    grid on; box on;

    nexttile(tl); hold on;
    loglog(Ns, err_naive, 'o-', 'Color', cm.CORAL, 'MarkerFaceColor', 'w', ...
           'MarkerSize', 6, 'LineWidth', 1.0, 'DisplayName', 'naive: max $|f - f_N|$');
    loglog(Ns, err_trick, '^-', 'Color', cm.TEAL, 'MarkerFaceColor', 'w', ...
           'MarkerSize', 6, 'LineWidth', 1.0, ...
           'DisplayName', 'trick: max $|f - (\mathrm{smooth}_N + \mathrm{singular})|$');
    yline(1e-15, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.4, 'Alpha', 0.5);
    set(gca, 'XScale', 'log', 'YScale', 'log');
    hold off;
    xlabel('$N$', 'Interpreter', 'latex');
    ylabel('max pointwise error', 'Interpreter', 'latex');
    title('subtraction reaches machine $\epsilon$ at modest $N$', 'Interpreter', 'latex');
    ylim([1e-17, 1]);
    legend('Location', 'southwest', 'Interpreter', 'latex', 'FontSize', 9);
    grid on; box on;

    exportgraphics(fig, fullfile(out_dir, 'singularity_subtraction.pdf'), 'ContentType', 'vector');
    exportgraphics(fig, fullfile(out_dir, 'singularity_subtraction.png'), 'Resolution', 300);
    close(fig);

    fprintf('[Etude 21.3]  singularity subtraction (1-D surrogate)\n');
    for k = 1:numel(Ns)
        fprintf('  N = %4d: naive err = %.3e, trick err = %.3e\n', ...
                Ns(k), err_naive(k), err_trick(k));
    end
    fprintf('  figure: %s\n', fullfile(out_dir, 'singularity_subtraction.pdf'));

    if ~isempty(dump_path)
        r = struct('Ns', Ns, 'err_naive', err_naive, 'err_trick', err_trick);
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
    ext = [v(:); flipud(v(2:N).')]';
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
