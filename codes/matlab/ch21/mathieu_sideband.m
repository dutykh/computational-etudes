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

    fig = figure('Units', 'inches', 'Position', [1, 1, 13.5, 4.4], 'Color', 'w');
    tl = tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile(tl);
    bar(modes, coeff_abs, 0.85, 'FaceColor', cm.NAVY, 'EdgeColor', cm.NAVY);
    hold on;
    xline(n_carrier, '--', 'Color', cm.CORAL, 'LineWidth', 0.8, 'Alpha', 0.7, ...
          'DisplayName', 'carrier $n=15$');
    hold off;
    xlim([0, 31]);
    xlabel('Fourier degree $n$  (cos basis, odd $n$)', 'Interpreter', 'latex');
    ylabel('$|a_n|$', 'Interpreter', 'latex');
    title(sprintf('Panel A. $|a_n|$ for $\\mathrm{ce}_{15}$ at $q=%.0f$', q_demo), ...
          'Interpreter', 'latex');
    legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 9);
    grid on; box on;

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
    title('Panel B. $\mathrm{ce}_{15}$: 5$\times$5 indistinguishable from full', ...
          'Interpreter', 'latex');
    legend('Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 9);
    grid on; box on;

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
    title('Panel C. $\mathrm{ce}_3$: $q/n^2$ large, sideband breaks', ...
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
