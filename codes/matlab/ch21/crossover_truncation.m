function crossover_truncation(varargin)
%% crossover_truncation - Etude 21.4: cross-over truncation.
%
% Boyd's cartoon: a_n = 10 exp(-n/3) + 1e-6 / n^5.
% Real example: f(x) = exp(-((x-0.3)/0.1)^2) + 1e-7 (x+1)^(1/3).
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
% Part of "Computational Etudes: A Spectral Approach"

    dump_path = parse_args(varargin{:});
    script_dir = fileparts(mfilename('fullpath'));
    addpath(script_dir);
    cm = tricks_common(); cm.configure_style();
    out_dir = cm.output_dir(script_dir);

    n_axis = 1:400;
    head_geom = 10.0 * exp(-n_axis / 3.0);
    tail_alg  = 1.0e-6 ./ n_axis.^5;
    a_total = head_geom + tail_alg;
    log_diff = log10(tail_alg + 1e-300) - log10(head_geom + 1e-300);
    [~, idx] = min(abs(log_diff));
    n_cross = double(n_axis(idx));

    f_real = @(x) exp(-((x - 0.3) / 0.1).^2) + 1e-7 * (x + 1.0).^(1/3);
    N = 400;
    xg = cgl(N);
    a = cheb_coeffs(f_real(xg));

    fig = figure('Units', 'inches', 'Position', [1, 1, 11.5, 4.4], 'Color', 'w');
    tl = tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile(tl); hold on;
    semilogy(n_axis, head_geom, 'Color', cm.NAVY, 'LineWidth', 1.2, ...
             'DisplayName', '$10\,e^{-n/3}$ (geometric head)');
    semilogy(n_axis, tail_alg, 'Color', cm.CORAL, 'LineWidth', 1.2, ...
             'DisplayName', '$10^{-6}/n^5$ (algebraic tail)');
    semilogy(n_axis, a_total, 'Color', cm.PURPLE, 'LineWidth', 1.6, ...
             'DisplayName', '$a_n = $ head $+$ tail');
    xline(n_cross, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8, 'Alpha', 0.6, ...
          'DisplayName', sprintf('$n_\\times \\approx %d$', n_cross));
    set(gca, 'YScale', 'log');
    hold off;
    xlabel('index $n$', 'Interpreter', 'latex');
    ylabel('$|a_n|$', 'Interpreter', 'latex');
    title('Panel A. Boyd''s cartoon: a slope change at $n_\times$', ...
          'Interpreter', 'latex');
    ylim([1e-22, 50]);
    legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 9);
    grid on; box on;

    nexttile(tl); hold on;
    n_axis2 = 0:N;
    semilogy(n_axis2, abs(a) + 1e-300, 'Color', cm.NAVY, 'LineWidth', 1.0, ...
             'Marker', '.', 'MarkerSize', 4, 'MarkerFaceColor', cm.NAVY, ...
             'DisplayName', 'computed $|a_n|$ for $f(x)$');
    set(gca, 'YScale', 'log');
    hold off;
    xlabel('Chebyshev degree $n$', 'Interpreter', 'latex');
    ylabel('$|a_n|$', 'Interpreter', 'latex');
    title('Panel B. $f(x) = e^{-((x-0.3)/0.1)^2} + 10^{-7}(x+1)^{1/3}$', ...
          'Interpreter', 'latex');
    ylim([1e-18, 5]);
    legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 9);
    grid on; box on;

    exportgraphics(fig, fullfile(out_dir, 'crossover_truncation.pdf'), 'ContentType', 'vector');
    exportgraphics(fig, fullfile(out_dir, 'crossover_truncation.png'), 'Resolution', 300);
    close(fig);

    fprintf('[Etude 21.4]  cross-over truncation\n');
    fprintf('  cartoon n_cross approx %d\n', n_cross);
    fprintf('  real example: |a_50|  = %.3e\n', abs(a(51)));
    fprintf('  real example: |a_200| = %.3e\n', abs(a(201)));
    fprintf('  real example: |a_400| = %.3e\n', abs(a(401)));
    fprintf('  figure: %s\n', fullfile(out_dir, 'crossover_truncation.pdf'));

    if ~isempty(dump_path)
        r = struct('n_cross_cartoon', n_cross, 'a_real', abs(a));
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
