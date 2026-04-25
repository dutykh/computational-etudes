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
    Es = linspace(0, 50, 1000);
    Ds = polyval(coeffs_num, Es);
    plot(Es, Ds, 'Color', cm.NAVY, 'LineWidth', 1.2, 'DisplayName', '$D(E)$');
    yline(0, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.4, 'Alpha', 0.5);
    for r = bender_orszag
        xline(r, ':', 'Color', cm.TEAL, 'LineWidth', 0.8, 'Alpha', 0.5);
    end
    scatter(rts_real, zeros(size(rts_real)), 80, cm.CORAL, 'x', 'LineWidth', 1.5, ...
            'DisplayName', 'numerical roots of $D$');
    scatter(bender_orszag, zeros(size(bender_orszag)), 70, cm.TEAL, 'o', 'LineWidth', 1.0, ...
            'DisplayName', 'Bender-Orszag $E_n$');
    hold off;
    xlabel('eigenvalue $E$', 'Interpreter', 'latex');
    ylabel('$D(E)$', 'Interpreter', 'latex');
    title('secular determinant $D(E)$ and its roots', 'Interpreter', 'latex');
    legend('Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 9);
    grid on; box on;

    nexttile(tl);
    ns = [0 2 4 6 8];
    rels = abs((rts_real(:).' - bender_orszag) ./ bender_orszag);
    semilogy(ns, rels, 'o-', 'Color', cm.NAVY, 'MarkerFaceColor', 'w', ...
             'MarkerSize', 8, 'LineWidth', 1.0, ...
             'DisplayName', '$|E_n^{\mathrm{sym}} - E_n| / E_n$');
    yline(0.01, 'Color', cm.TEAL, 'LineWidth', 0.4, 'Alpha', 0.5, ...
          'DisplayName', '1\% threshold');
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
