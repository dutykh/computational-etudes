function ioakimidis_root(varargin)
%% ioakimidis_root - Etude 21.7: non-iterative root via Chebyshev quadrature.
% f(x) = sin(x - pi/4) / sqrt(1 + 10 x^2) on [-1, 1] has the root x = pi/4.
% Ioakimidis recovers it geometrically; bisection on the same number of
% evaluations is much slower.
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
% Part of "Computational Etudes: A Spectral Approach"

    dump_path = parse_args(varargin{:});
    script_dir = fileparts(mfilename('fullpath'));
    addpath(script_dir);
    cm = tricks_common(); cm.configure_style();
    out_dir = cm.output_dir(script_dir);

    rho_exact = pi / 4.0;
    a = -1.0; b = 1.0;
    Ns = 2:30;
    err_io = zeros(size(Ns));
    err_bi = zeros(size(Ns));
    for k = 1:numel(Ns)
        rho_io = ioakimidis(Ns(k), a, b);
        err_io(k) = abs(rho_io - rho_exact);
        rho_bi = bisection(2*Ns(k) + 1, a, b);
        err_bi(k) = abs(rho_bi - rho_exact);
    end

    fig = figure('Units', 'inches', 'Position', [1, 1, 7.5, 4.4], 'Color', 'w');
    n_evals = 2 * Ns + 1;
    semilogy(n_evals, err_io + 1e-18, 'o-', 'Color', cm.NAVY, 'MarkerFaceColor', 'w', ...
             'MarkerSize', 6, 'LineWidth', 1.0, ...
             'DisplayName', 'Ioakimidis (non-iterative)');
    hold on;
    semilogy(n_evals, err_bi + 1e-18, 's--', 'Color', cm.CORAL, 'MarkerFaceColor', 'w', ...
             'MarkerSize', 5, 'LineWidth', 0.9, ...
             'DisplayName', 'bisection (same \# of evals)');
    yline(1e-15, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.4, 'Alpha', 0.5);
    set(gca, 'YScale', 'log');
    hold off;
    xlabel('number of $f$ evaluations  ($= 2N+1$)', 'Interpreter', 'latex');
    ylabel('$|\rho_N - \pi/4|$', 'Interpreter', 'latex');
    title('Ioakimidis: geometric, no Newton, on $f(x)=\sin(x-\pi/4)/\sqrt{1+10x^2}$', ...
          'Interpreter', 'latex');
    legend('Location', 'southwest', 'Interpreter', 'latex', 'FontSize', 10);
    grid on; box on;
    ylim([1e-16, 1]);

    exportgraphics(fig, fullfile(out_dir, 'ioakimidis_root.pdf'), 'ContentType', 'vector');
    exportgraphics(fig, fullfile(out_dir, 'ioakimidis_root.png'), 'Resolution', 300);
    close(fig);

    fprintf('[Etude 21.7]  Ioakimidis non-iterative root\n');
    fprintf('  exact rho = pi/4 = %.16f\n', rho_exact);
    for kk = [3 6 10 20 30]
        rho = ioakimidis(kk, a, b);
        fprintf('  N = %2d (%2d evals): rho = %.16f, err = %.3e\n', ...
                kk, 2*kk+1, rho, abs(rho - rho_exact));
    end
    fprintf('  figure: %s\n', fullfile(out_dir, 'ioakimidis_root.pdf'));

    if ~isempty(dump_path)
        r = struct('rho_exact', rho_exact, 'Ns', Ns, ...
                   'err_ioakimidis', err_io, 'err_bisection', err_bi);
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

function v = f_test(x)
    v = sin(x - pi/4) ./ sqrt(1.0 + 10.0 * x.^2);
end

function rho = ioakimidis(N, a, b)
    j = 0:(2*N);
    half = 0.5 * ((a + b) - (a - b) * cos(j * pi / (2 * N)));
    fj = f_test(half);
    sign_j = (-1).^j;
    weight = ones(size(j));
    weight(1) = 0.5;
    weight(end) = 0.5;
    num = sum(weight .* sign_j .* half ./ fj);
    den = sum(weight .* sign_j        ./ fj);
    rho = num / den;
end

function rho = bisection(n_evals, a, b)
    fa = f_test(a); fb = f_test(b);
    assert(fa * fb < 0);
    for k = 1:(n_evals - 2)
        c = 0.5 * (a + b);
        fc = f_test(c);
        if fc * fa < 0
            b = c; fb = fc;
        else
            a = c; fa = fc;
        end
    end
    rho = 0.5 * (a + b);
end
