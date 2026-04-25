function radiation_scattering(varargin)
%% radiation_scattering - Etude 21.5: 1-D Schroedinger scattering by a
% sech^2 potential, via Boyd's radiation-augmented basis (Eqs 19.20-19.31).
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
% Part of "Computational Etudes: A Spectral Approach"

    dump_path = parse_args(varargin{:});
    script_dir = fileparts(mfilename('fullpath'));
    addpath(script_dir);
    cm = tricks_common(); cm.configure_style();
    out_dir = cm.output_dir(script_dir);

    v = 1.0; ell = 2.0; N = 48;
    k_axis = [0.3 0.6 0.9 1.2 1.5 1.8 2.1 2.4 2.7 3.0];
    R_num = zeros(size(k_axis));
    R_exact = zeros(size(k_axis));
    drifts  = zeros(size(k_axis));
    for i = 1:numel(k_axis)
        [alpha, drift] = solve_scattering(N, ell, k_axis(i), v);
        R_num(i)  = abs(alpha)^2;
        R_exact(i) = reflection_exact(k_axis(i), v);
        drifts(i)  = abs(drift);
    end
    abs_err = abs(R_num - R_exact);

    fig = figure('Units', 'inches', 'Position', [1, 1, 11.5, 4.4], 'Color', 'w');
    tl = tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile(tl); hold on;
    semilogy(k_axis, R_exact, 'Color', cm.NAVY, 'LineWidth', 1.4, ...
             'DisplayName', '$R$ exact (Boyd Eq 19.31)');
    semilogy(k_axis, R_num, 'o', 'Color', cm.CORAL, 'MarkerFaceColor', 'w', ...
             'MarkerSize', 8, 'DisplayName', sprintf('$R$ numerical, $N=%d$, $\\ell=%d$', N, ell));
    set(gca, 'YScale', 'log');
    hold off;
    xlabel('wavenumber $k$', 'Interpreter', 'latex');
    ylabel('reflection coefficient $R$', 'Interpreter', 'latex');
    title('$\mathrm{sech}^2$ scattering: spectral $R$ matches the closed form', ...
          'Interpreter', 'latex');
    legend('Location', 'southwest', 'Interpreter', 'latex', 'FontSize', 10);
    grid on; box on;
    ylim([1e-10, 2]);

    nexttile(tl); hold on;
    semilogy(k_axis, abs_err + 1e-18, 's-', 'Color', cm.NAVY, 'MarkerFaceColor', 'w', ...
             'MarkerSize', 6, 'LineWidth', 1.0, ...
             'DisplayName', '$|R_{\mathrm{num}} - R_{\mathrm{exact}}|$');
    semilogy(k_axis, drifts + 1e-18, '^--', 'Color', cm.TEAL, 'MarkerFaceColor', 'w', ...
             'MarkerSize', 5, 'LineWidth', 0.8, ...
             'DisplayName', '$|\sum_n a_n|$ (asymptotic-constant drift)');
    set(gca, 'YScale', 'log');
    hold off;
    xlabel('wavenumber $k$', 'Interpreter', 'latex');
    ylabel('error magnitude');
    title('error stays well below the rapidly small $R$', 'Interpreter', 'latex');
    legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 10);
    grid on; box on;

    exportgraphics(fig, fullfile(out_dir, 'radiation_scattering.pdf'), 'ContentType', 'vector');
    exportgraphics(fig, fullfile(out_dir, 'radiation_scattering.png'), 'Resolution', 300);
    close(fig);

    fprintf('[Etude 21.5]  Schroedinger sech^2 scattering with radiation basis\n');
    fprintf('  v = %.1f, ell = %.1f, N = %d\n', v, ell, N);
    fprintf('   k        R_numerical       R_exact         abs error\n');
    for i = 1:numel(k_axis)
        fprintf('  %.1f  %.7e   %.7e   %.2e\n', k_axis(i), R_num(i), R_exact(i), abs_err(i));
    end
    fprintf('  figure: %s\n', fullfile(out_dir, 'radiation_scattering.pdf'));

    if ~isempty(dump_path)
        r = struct('N', N, 'ell', ell, 'v', v, 'k_axis', k_axis, ...
                   'R_numerical', R_num, 'R_exact', R_exact, ...
                   'abs_err', abs_err, 'drifts', drifts);
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

function R = reflection_exact(k, v)
    s = 2 * pi * sqrt(v + 0.25);
    R = (1 + cos(s)) / (cosh(2 * pi * k) + cos(s));
end

function [alpha, drift] = solve_scattering(N, ell, k, v)
    i = 1:N;
    t = pi * (2 * i - 1) / (2 * N);
    x = ell ./ tan(t);
    Vx = -v ./ cosh(x).^2;

    Phi = zeros(N, N);
    Phipp = zeros(N, N);
    n_idx = 0:N-3;
    for j = 1:numel(n_idx)
        nn = n_idx(j);
        Phi(:, j) = cos(nn * t).';
        s = sin(t).'; c = cos(t).';
        Phipp(:, j) = -(nn / ell^2) * s.^3 .* (nn * cos(nn * t).' .* s + 2.0 * sin(nn * t).' .* c);
    end
    H_x  = 0.5 * (1 + tanh(x));
    Hp_x = 0.5 ./ cosh(x).^2;
    Hpp_x = -tanh(x) ./ cosh(x).^2;
    Phi(:, N-1) = (H_x .* cos(k * x)).';
    Phipp(:, N-1) = (Hpp_x .* cos(k * x) - 2 * k * Hp_x .* sin(k * x) ...
                     - k^2 * H_x .* cos(k * x)).';
    Phi(:, N) = (H_x .* sin(k * x)).';
    Phipp(:, N) = (Hpp_x .* sin(k * x) + 2 * k * Hp_x .* cos(k * x) ...
                   - k^2 * H_x .* sin(k * x)).';

    M = Phipp + (k^2 - Vx.') .* Phi;
    f = (Vx .* cos(k * x)).';
    g = (Vx .* sin(k * x)).';

    a = M \ f;
    b = M \ g;

    gamma1 = a(N-1) + 1.0;
    gamma2 = a(N);
    sigma1 = b(N-1);
    sigma2 = b(N) + 1.0;

    A = [gamma1 + sigma2, sigma1 - gamma2;
         gamma2 - sigma1, gamma1 + sigma2];
    rhs = [sigma2 - gamma1; -sigma1 - gamma2];
    re_im = A \ rhs;
    alpha = complex(re_im(1), re_im(2));
    drift = sum(a(1:N-2));
end
