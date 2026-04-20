function eigen_complex_detour(varargin)
%% eigen_complex_detour - Etude 18.11: detour around the singularity.
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
% Part of "Computational Etudes: A Spectral Approach"

    dump_path = parse_args(varargin{:});
    script_dir = fileparts(mfilename('fullpath'));
    addpath(fullfile(script_dir, '..', 'ch07'));
    out_dir = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch18', 'matlab');
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end
    configure_style(); [NAVY, CORAL, TEAL, PURPLE] = colours();

    Delta_scan = [0.25, 0.50, 0.75];
    N_demo = 40; N_low = 20;

    [lam_ref, ~, ~] = solve_detoured(N_demo, 0.5);
    first6_ref = lam_ref(1:6);

    fig = figure('Units', 'inches', 'Position', [1, 1, 11.5, 4.5], 'Color', 'w');
    tl = tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile(tl); hold on;
    cols = {NAVY, CORAL, TEAL};
    for k = 1:numel(Delta_scan)
        xd = linspace(-1, 1, 401);
        yd = xd + 1i * Delta_scan(k) * (xd.^2 - 1);
        plot(real(yd), imag(yd), 'Color', cols{k}, 'LineWidth', 1.2, ...
            'DisplayName', sprintf('$\\Delta = %.2f$', Delta_scan(k)));
    end
    yline(0, 'k', 'LineWidth', 0.4, 'Alpha', 0.5, 'HandleVisibility', 'off');
    xline(0, 'k', 'LineWidth', 0.4, 'Alpha', 0.5, 'HandleVisibility', 'off');
    plot(0, 0, 'x', 'Color', PURPLE, 'MarkerSize', 10, 'LineWidth', 2, ...
         'DisplayName', 'singularity $y = 0$');
    plot([-1, 1], [0, 0], 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 8, ...
         'LineWidth', 1.5, 'HandleVisibility', 'off');
    hold off;
    xlabel('Re $y$', 'Interpreter', 'latex');
    ylabel('Im $y$', 'Interpreter', 'latex');
    title('parabolic detour $y(x) = x + i\Delta(x^2 - 1)$', 'Interpreter', 'latex');
    xlim([-1.2, 1.2]); ylim([-1.0, 0.3]);
    grid on; box on;
    legend('Location', 'southeast', 'FontSize', 9, 'Interpreter', 'latex');

    nexttile(tl); hold on;
    marks = {'o', 's', '^'};
    for k = 1:numel(Delta_scan)
        [lam, ~, ~] = solve_detoured(N_demo, Delta_scan(k));
        first_few = lam(1:min(7, numel(lam)));
        scatter(real(first_few), imag(first_few), 45, 'MarkerEdgeColor', cols{k}, ...
            'MarkerFaceColor', 'w', 'Marker', marks{k}, 'LineWidth', 1.2, ...
            'DisplayName', sprintf('$\\Delta = %.2f$, $N = %d$', Delta_scan(k), N_demo));
    end
    yline(0, 'k', 'LineWidth', 0.4, 'Alpha', 0.5, 'HandleVisibility', 'off');
    xline(0, 'k', 'LineWidth', 0.4, 'Alpha', 0.5, 'HandleVisibility', 'off');
    hold off;
    xlabel('Re $\lambda$', 'Interpreter', 'latex');
    ylabel('Im $\lambda$', 'Interpreter', 'latex');
    title('first eigenvalues in the complex plane', 'Interpreter', 'latex');
    grid on; box on;
    legend('Location', 'northeast', 'FontSize', 9, 'Interpreter', 'latex');

    exportgraphics(fig, fullfile(out_dir, 'eigen_complex_detour.pdf'), 'ContentType', 'vector');
    exportgraphics(fig, fullfile(out_dir, 'eigen_complex_detour.png'), 'Resolution', 300);
    close(fig);

    [lam_low, ~, ~] = solve_detoured(N_low, 0.5);

    fprintf('[Etude 18.11]  complex-plane detour\n');
    fprintf('  first 6 at N = %d, Delta = 0.5:\n', N_demo);
    for k = 1:6
        z = first6_ref(k);
        s = ' +'; if imag(z) < 0; s = ' -'; end
        fprintf('    %+.6f%s i %.6f\n', real(z), s, abs(imag(z)));
    end
    fprintf('  low-N check (N=%d, Delta=0.5):\n', N_low);
    for k = 1:6
        z = lam_low(k);
        s = ' +'; if imag(z) < 0; s = ' -'; end
        fprintf('    %+.6f%s i %.6f\n', real(z), s, abs(imag(z)));
    end
    fprintf('  figure: %s\n', fullfile(out_dir, 'eigen_complex_detour.pdf'));

    if ~isempty(dump_path)
        r = struct('N_demo', N_demo, 'N_low', N_low, 'Delta_scan', Delta_scan);
        r.first6_ref_Delta0p5 = [real(first6_ref), imag(first6_ref)];
        r.low_N_first6 = [real(lam_low(1:6)), imag(lam_low(1:6))];
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

function [lam, y, y_i] = solve_detoured(N, Delta)
    [D, x] = cheb_matrix(N);
    y = x + 1i * Delta * (x.^2 - 1);
    dy_dx = 1 + 2i * Delta * x;
    d2y_dx2 = 2i * Delta * ones(size(x));
    idx = 2:N;
    D_int = D(idx, idx);
    D2_int = D * D;
    D2_int = D2_int(idx, idx);
    dydx_i = dy_dx(idx);
    d2ydx2_i = d2y_dx2(idx);
    y_i = y(idx);
    L = diag(dydx_i.^(-2)) * D2_int ...
        - diag(d2ydx2_i .* dydx_i.^(-3)) * D_int ...
        + diag(1 ./ y_i);
    lam = eig(L);
    [~, order] = sort(abs(lam));
    lam = lam(order);
end

function configure_style()
    set(groot, 'defaultAxesFontName', 'CMU Serif');
    set(groot, 'defaultTextFontName', 'CMU Serif');
    set(groot, 'defaultAxesFontSize', 11);
    set(groot, 'defaultAxesLineWidth', 0.8);
    set(groot, 'defaultLineLineWidth', 1.0);
    set(groot, 'defaultAxesBox', 'on');
end

function [NAVY, CORAL, TEAL, PURPLE] = colours()
    NAVY   = [0x14, 0x2D, 0x6E] / 255;
    CORAL  = [0xE7, 0x4C, 0x3C] / 255;
    TEAL   = [0x16, 0xA0, 0x85] / 255;
    PURPLE = [0x8E, 0x44, 0xAD] / 255;
end
