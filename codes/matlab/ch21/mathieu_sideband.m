function mathieu_sideband(varargin)
%% mathieu_sideband - Etude 21.2: Mathieu's equation and sideband truncation.
% Reproduces Boyd Figs 19.1 (|a_n| for ce_15 at q = 10) and 19.2
% (eigenvalue correction delta(q) for n = 15) plus a small-n breakdown.
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
% Part of "Computational Etudes: A Spectral Approach"

    dump_path = parse_args(varargin{:});
    script_dir = fileparts(mfilename('fullpath'));
    addpath(script_dir);
    cm = tricks_common(); cm.configure_style();
    out_dir = cm.output_dir(script_dir);

    n_carrier = 15;
    N_modes = 64;
    q_demo = 10.0;
    q_axis = linspace(0.0, 25.0, 200);

    [lam_demo, vec_demo, modes] = reference_eigenpair(N_modes, q_demo, n_carrier);
    idx_carrier = find(modes == n_carrier);
    vec_demo = vec_demo * sign(vec_demo(idx_carrier));
    coeff_abs = abs(vec_demo);

    delta_full = zeros(size(q_axis));
    delta_3 = zeros(size(q_axis));
    delta_5 = zeros(size(q_axis));
    for k = 1:numel(q_axis)
        [lam, ~, ~] = reference_eigenpair(N_modes, q_axis(k), n_carrier);
        delta_full(k) = lam - n_carrier^2;
        delta_3(k) = sideband_eigenvalue(n_carrier, q_axis(k), 1) - n_carrier^2;
        delta_5(k) = sideband_eigenvalue(n_carrier, q_axis(k), 2) - n_carrier^2;
    end

    n_small = 3;
    q_small_axis = linspace(0.0, 10.0, 200);
    delta_full_small = zeros(size(q_small_axis));
    delta_3_small   = zeros(size(q_small_axis));
    delta_5_small   = zeros(size(q_small_axis));
    for k = 1:numel(q_small_axis)
        [lam, ~, ~] = reference_eigenpair(N_modes, q_small_axis(k), n_small);
        delta_full_small(k) = lam - n_small^2;
        delta_3_small(k) = sideband_eigenvalue(n_small, q_small_axis(k), 1) - n_small^2;
        delta_5_small(k) = sideband_eigenvalue(n_small, q_small_axis(k), 2) - n_small^2;
    end

    % NEW (panel b): coefficient picture for ce_3 at q = q_demo
    [~, vec_demo3, ~] = reference_eigenpair(N_modes, q_demo, n_small);
    idx3 = find(modes == n_small);
    vec_demo3 = vec_demo3 * sign(vec_demo3(idx3));
    coeff_abs3 = abs(vec_demo3);

    % NEW (panel c): cluster width vs q at three carriers
    q_cluster = linspace(0.5, 50.0, 80);
    n_cs = [3 7 15];
    cluster_widths = zeros(numel(n_cs), numel(q_cluster));
    for ci = 1:numel(n_cs)
        n_c = n_cs(ci);
        idx_c = find(modes == n_c);
        for k = 1:numel(q_cluster)
            [~, vec, ~] = reference_eigenpair(N_modes, q_cluster(k), n_c);
            ac = abs(vec(idx_c));
            if ac == 0; cluster_widths(ci, k) = 0; continue; end
            mask = abs(vec) > 1e-3 * ac;
            offsets = modes(:) - n_c;
            if any(mask)
                cluster_widths(ci, k) = floor(max(abs(offsets(mask))) / 2);
            end
        end
    end

    % NEW (panel e): convergence vs sideband size at q=q_demo, n=15
    half_widths = [1 2 3 4 5 7 10];
    [~, k_demo] = min(abs(q_axis - q_demo));
    delta_full_15 = delta_full(k_demo);
    err_vs_hw = zeros(size(half_widths));
    for k = 1:numel(half_widths)
        d_hw = sideband_eigenvalue(n_carrier, q_demo, half_widths(k)) - n_carrier^2;
        err_vs_hw(k) = abs(d_hw - delta_full_15) + 1e-18;
    end

    fig = figure('Units', 'inches', 'Position', [1, 1, 11.0, 12.0], 'Color', 'w');
    tl = tiledlayout(3, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    % (a) |a_n| for ce_15
    nexttile(tl);
    bar(modes, coeff_abs, 0.85, 'FaceColor', cm.NAVY, 'EdgeColor', cm.NAVY);
    hold on;
    xline(n_carrier, '--', 'Color', cm.CORAL, 'LineWidth', 0.8, 'Alpha', 0.7, ...
          'DisplayName', 'carrier $n=15$');
    hold off;
    xlim([0, 31]);
    xlabel('Fourier degree $n$ (cos basis, odd $n$)', 'Interpreter', 'latex');
    ylabel('$|a_n|$', 'Interpreter', 'latex');
    title(sprintf('(a) $|a_n|$ for $\\mathrm{ce}_{15}$ at $q=%.0f$', q_demo), ...
          'Interpreter', 'latex');
    legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 9);
    grid on; box on;

    % (b) cluster width vs q
    nexttile(tl); hold on;
    cluster_pal = {cm.CORAL, cm.TEAL, cm.NAVY};
    for ci = 1:numel(n_cs)
        plot(q_cluster, cluster_widths(ci, :), '-o', 'Color', cluster_pal{ci}, ...
             'MarkerFaceColor', 'w', 'MarkerSize', 3, 'LineWidth', 1.0, ...
             'DisplayName', sprintf('$n = %d$', n_cs(ci)));
    end
    yline(2, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.4, ...
          'DisplayName', '$5\times5$ box capacity (hw=2)');
    hold off;
    xlabel('$q$', 'Interpreter', 'latex');
    ylabel('cluster half-width (modes)', 'Interpreter', 'latex');
    title('(b) cluster width vs $q$ at three carriers', 'Interpreter', 'latex');
    legend('Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 9);
    grid on; box on;

    % (c) |a_n| for ce_3
    nexttile(tl);
    bar(modes, coeff_abs3, 0.85, 'FaceColor', cm.NAVY, 'EdgeColor', cm.NAVY);
    hold on;
    xline(n_small, '--', 'Color', cm.CORAL, 'LineWidth', 0.8, 'Alpha', 0.7, ...
          'DisplayName', 'carrier $n=3$');
    hold off;
    xlim([0, 31]);
    xlabel('Fourier degree $n$ (cos basis, odd $n$)', 'Interpreter', 'latex');
    ylabel('$|a_n|$', 'Interpreter', 'latex');
    title(sprintf('(c) $|a_n|$ for $\\mathrm{ce}_3$ at $q=%.0f$ (cluster spread)', q_demo), ...
          'Interpreter', 'latex');
    legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 9);
    grid on; box on;

    % (d) delta(q) for n=15
    nexttile(tl); hold on;
    plot(q_axis, delta_full, 'Color', cm.NAVY, 'LineWidth', 1.4, ...
         'DisplayName', '$\delta_{\mathrm{full}}$ (high-$N$)');
    plot(q_axis, delta_5, '--', 'Color', cm.TEAL, 'LineWidth', 1.0, ...
         'DisplayName', '$\delta_5$ (5$\times$5)');
    plot(q_axis, delta_3, ':', 'Color', cm.CORAL, 'LineWidth', 1.0, ...
         'DisplayName', '$\delta_3$ (3$\times$3)');
    hold off;
    xlabel('$q$', 'Interpreter', 'latex');
    ylabel('$\delta(q) = \lambda(q) - 15^2$', 'Interpreter', 'latex');
    title('(d) $\mathrm{ce}_{15}$: 5$\times$5 indistinguishable from full', ...
          'Interpreter', 'latex');
    legend('Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 9);
    grid on; box on;

    % (e) convergence vs sideband size
    nexttile(tl); hold on; set(gca, 'YScale', 'log');
    plot(half_widths, err_vs_hw, '-o', 'Color', cm.NAVY, ...
         'MarkerFaceColor', 'w', 'MarkerSize', 6, 'LineWidth', 1.1, ...
         'DisplayName', '$|\delta_{\mathrm{hw}} - \delta_{\mathrm{full}}|$');
    xline(2, ':', 'Color', cm.CORAL, 'LineWidth', 0.8, ...
          'DisplayName', '$5\times5$ box (hw=2)');
    hold off;
    xlabel('sideband half-width hw', 'Interpreter', 'latex');
    ylabel('absolute error in $\delta$', 'Interpreter', 'latex');
    title(sprintf('(e) convergence vs sideband size ($n=15$, $q=%.0f$)', q_demo), ...
          'Interpreter', 'latex');
    legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 9);
    grid on; box on;

    % (f) delta(q) for n=3
    nexttile(tl); hold on;
    plot(q_small_axis, delta_full_small, 'Color', cm.NAVY, 'LineWidth', 1.4, ...
         'DisplayName', '$\delta_{\mathrm{full}}$');
    plot(q_small_axis, delta_5_small, '--', 'Color', cm.TEAL, 'LineWidth', 1.0, ...
         'DisplayName', '$\delta_5$');
    plot(q_small_axis, delta_3_small, ':', 'Color', cm.CORAL, 'LineWidth', 1.0, ...
         'DisplayName', '$\delta_3$');
    hold off;
    xlabel('$q$', 'Interpreter', 'latex');
    ylabel('$\delta(q) = \lambda(q) - 3^2$', 'Interpreter', 'latex');
    title('(f) $\mathrm{ce}_3$: $q/n^2$ large, sideband breaks', ...
          'Interpreter', 'latex');
    legend('Location', 'southwest', 'Interpreter', 'latex', 'FontSize', 9);
    grid on; box on;

    exportgraphics(fig, fullfile(out_dir, 'mathieu_sideband.pdf'), 'ContentType', 'vector');
    exportgraphics(fig, fullfile(out_dir, 'mathieu_sideband.png'), 'Resolution', 300);
    close(fig);

    fprintf('[Etude 21.2]  Mathieu sideband truncation\n');
    fprintf('  ce_%d, q = %.1f\n', n_carrier, q_demo);
    [~, idx_sort] = sort(coeff_abs, 'descend');
    for j = idx_sort(1:5).'
        fprintf('    n = %2d  |a_n| = %.4e\n', modes(j), coeff_abs(j));
    end
    fprintf('  figure: %s\n', fullfile(out_dir, 'mathieu_sideband.pdf'));

    if ~isempty(dump_path)
        r = struct('n_carrier', n_carrier, 'N_modes', N_modes, 'q_demo', q_demo, ...
                   'modes_demo', modes, 'coeff_abs_demo', coeff_abs, ...
                   'q_axis', q_axis, 'delta_full', delta_full, ...
                   'delta_5', delta_5, 'delta_3', delta_3);
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

function modes = odd_modes(M)
    modes = 2 * (1:M) - 1;
end

function [M, modes] = galerkin_full(N_modes, q)
    modes = odd_modes(N_modes);
    M = diag(double(modes).^2);
    for i = 1:N_modes
        for j = 1:N_modes
            if abs(modes(i) - modes(j)) == 2
                M(i, j) = q;
            end
        end
    end
end

function [lam_pick, v_pick, modes] = reference_eigenpair(N_modes, q, n_carrier)
    [M, modes] = galerkin_full(N_modes, q);
    [V, D] = eig((M + M.') / 2);
    lam = diag(D);
    idx_carrier = find(modes == n_carrier);
    [~, dominant] = max(abs(V), [], 1);
    candidates = find(dominant == idx_carrier);
    if isempty(candidates)
        [~, k] = min(abs(lam - n_carrier^2));
    else
        [~, kk] = max(abs(V(idx_carrier, candidates)));
        k = candidates(kk);
    end
    lam_pick = lam(k);
    v_pick = V(:, k);
end

function lam = sideband_eigenvalue(n_carrier, q, half_width)
    modes = arrayfun(@(k) n_carrier + 2 * k, -half_width:half_width);
    M = diag(double(modes).^2);
    for i = 1:numel(modes)
        for j = 1:numel(modes)
            if abs(modes(i) - modes(j)) == 2
                M(i, j) = q;
            end
        end
    end
    lams = sort(eig((M + M.') / 2));
    [~, idx] = min(abs(lams - n_carrier^2));
    lam = lams(idx);
end
