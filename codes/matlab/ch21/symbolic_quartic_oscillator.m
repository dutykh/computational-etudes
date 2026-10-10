function symbolic_quartic_oscillator(varargin)
%% symbolic_quartic_oscillator - Etude 21.10: quartic oscillator on the line.
% Symbolic 5-term Galerkin via the rational-Chebyshev map y = ell x/sqrt(1-x^2),
% ell = 2.  Produces the secular polynomial of Boyd Eq 20.16; roots compared
% to Bender-Orszag references.
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
% Part of "Computational Etudes: A Spectral Approach"

    dump_path = parse_args(varargin{:});
    script_dir = fileparts(mfilename('fullpath'));
    addpath(script_dir);
    cm = tricks_common(); cm.configure_style();
    out_dir = cm.output_dir(script_dir);

    [coeffs_E, ~] = build_secular();
    coeffs_num = double(coeffs_E);                    % high-degree first
    fprintf('Symbolic secular determinant (high degree -> low):\n');
    for k = 1:numel(coeffs_num)
        fprintf('   E^%d: %.6e\n', numel(coeffs_num) - k, coeffs_num(k));
    end
    norm_const = coeffs_num(end);
    fprintf('\nBoyd Eq 20.16 form D(E) / D(0):\n');
    for k = 1:numel(coeffs_num)
        fprintf('   E^%d: %.6e\n', k - 1, coeffs_num(end - k + 1) / norm_const);
    end

    rts = roots(coeffs_num);
    rts_real = sort(real(rts(abs(imag(rts)) < 1e-6)));
    bender_orszag = [1.060362, 7.45570, 16.2618, 26.528, 37.92];
    fprintf('\nEigenvalues from secular determinant vs Bender-Orszag:\n');
    for k = 1:numel(rts_real)
        rel = (rts_real(k) - bender_orszag(k)) / bender_orszag(k);
        fprintf('   E_num = %11.4f    E_exact = %9.4f    rel err = %+.2e\n', ...
                rts_real(k), bender_orszag(k), rel);
    end

    fig = figure('Units', 'inches', 'Position', [1, 1, 11.5, 4.4], 'Color', 'w');
    tl = tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile(tl); hold on;
    E_lo = 0.0; E_hi = 50.0;
    Es = linspace(E_lo, E_hi, 2000);
    Ds = polyval(coeffs_num, Es);
    log_absD = log10(abs(Ds) + 1e-30);
    plot(Es, log_absD, 'Color', cm.NAVY, 'LineWidth', 1.0, ...
         'DisplayName', '$\log_{10}|D(E)|$');

    in_window = @(v) (v >= E_lo) & (v <= E_hi);
    refs_in   = bender_orszag(in_window(bender_orszag));
    roots_in  = rts_real(in_window(rts_real(:).')).';
    roots_out = rts_real(~in_window(rts_real(:).')).';

    rug_y = min(log_absD) - 0.6;
    for r = refs_in
        xline(r, ':', 'Color', cm.TEAL, 'LineWidth', 0.6, 'Alpha', 0.35, ...
              'HandleVisibility', 'off');
    end
    scatter(refs_in, rug_y * ones(size(refs_in)), 70, 'o', ...
            'MarkerEdgeColor', cm.TEAL, 'LineWidth', 1.2, ...
            'DisplayName', 'Bender-Orszag $E_n$ (reference)');
    scatter(roots_in, rug_y * ones(size(roots_in)), 80, cm.CORAL, 'x', ...
            'LineWidth', 1.5, 'DisplayName', 'numerical roots of $D$');

    if ~isempty(roots_out)
        labs = {'$E_6$', '$E_8$'};
        parts = {};
        for k = 1:min(numel(roots_out), 2)
            parts{end+1} = sprintf('%s $\\approx$ %.0f', labs{k}, roots_out(k));
        end
        annot_text = [strjoin(parts, ', '), ' (numerical mirages, see right panel)'];
        text(E_hi - 28, rug_y - 1.6, annot_text, ...
             'Color', cm.CORAL, 'FontSize', 8, 'Interpreter', 'latex');
        annotation('arrow', 'Color', cm.CORAL, 'LineWidth', 0.8);
    end

    hold off;
    xlabel('eigenvalue $E$', 'Interpreter', 'latex');
    ylabel('$\log_{10}|D(E)|$', 'Interpreter', 'latex');
    title('secular determinant zeros on the lower window $E\in[0,50]$', ...
          'Interpreter', 'latex');
    xlim([E_lo, E_hi]);
    ylim([rug_y - 2.2, max(log_absD) + 0.5]);
    legend('Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 9);
    grid on; box on;

    nexttile(tl);
    ns = [0 2 4 6 8];
    rels = abs((rts_real(:).' - bender_orszag) ./ bender_orszag);
    semilogy(ns, rels, 'o-', 'Color', cm.NAVY, 'MarkerFaceColor', 'w', ...
             'MarkerSize', 8, 'LineWidth', 1.0, ...
             'DisplayName', '$|E_n^{\mathrm{sym}} - E_n| / E_n$');
    yline(0.01, 'Color', cm.TEAL, 'LineWidth', 0.4, 'Alpha', 0.5, ...
          'DisplayName', '1% threshold');
    set(gca, 'YScale', 'log');
    xlabel('physical mode index $n$', 'Interpreter', 'latex');
    ylabel('relative error');
    title('trust the lower spectrum: bottom 3 modes good, top 2 unusable', ...
          'Interpreter', 'latex');
    legend('Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 9);
    grid on; box on;
    ylim([1e-4, 1e2]);

    exportgraphics(fig, fullfile(out_dir, 'symbolic_quartic_oscillator.pdf'), 'ContentType', 'vector');
    exportgraphics(fig, fullfile(out_dir, 'symbolic_quartic_oscillator.png'), 'Resolution', 300);
    close(fig);

    fprintf('  figure: %s\n', fullfile(out_dir, 'symbolic_quartic_oscillator.pdf'));

    if ~isempty(dump_path)
        r = struct('ell', 2.0, 'secular_coeffs', coeffs_num, ...
                   'roots', rts_real.', 'bender_orszag', bender_orszag, 'rel_err', rels);
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

function [coeffs_E, M] = build_secular()
    syms x E real
    ell = sym(2);
    a = sym('a_', [1 5]);
    p = sum(a .* x.^(2 * (0:4)));
    u = (1 - x^2) * p;
    u_x = diff(u, x);
    u_xx = diff(u, x, 2);
    R = (1 - x^2)^4 * ((1 - x^2) * u_xx - 3 * x * u_x) / ell^2 ...
        + ((1 - x^2)^2 * E - ell^4 * x^4) * u;
    eqs = sym(zeros(5, 1));
    for j = 0:4
        eqs(j + 1) = int(x^(2 * j) * R, x, -1, 1);
    end
    M = sym(zeros(5, 5));
    for i = 1:5
        for j = 1:5
            M(i, j) = simplify(diff(eqs(i), a(j)));
        end
    end
    DE = simplify(det(M));
    p_E = sym2poly(DE);
    coeffs_E = p_E;
end
