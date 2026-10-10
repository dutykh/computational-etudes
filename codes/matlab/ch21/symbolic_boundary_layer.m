function symbolic_boundary_layer(varargin)
%% symbolic_boundary_layer - Etude 21.9: Carrier-Pearson boundary-layer BVP
% solved by a four-term symbolic Galerkin method (Boyd Eq 20.5-20.8).
% Reproduces Boyd Table 20.3.
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
% Part of "Computational Etudes: A Spectral Approach"

    dump_path = parse_args(varargin{:});
    script_dir = fileparts(mfilename('fullpath'));
    addpath(script_dir);
    cm = tricks_common(); cm.configure_style();
    out_dir = cm.output_dir(script_dir);

    [u4_expr, a_dict, x, eps_sym] = galerkin_4_term();
    fprintf('Symbolic four-term Galerkin solution:\n');
    fns = fieldnames(a_dict);
    for k = 1:numel(fns)
        fprintf('  %s = %s\n', fns{k}, char(simplify(a_dict.(fns{k}))));
    end

    u4_func = matlabFunction(u4_expr, 'Vars', [x, eps_sym]);
    u_exact_func = @(xx, ee) 1 - cosh(xx ./ ee) ./ cosh(1.0 ./ ee);

    eps_table = [1/20, 3/40, 1/10, 3/20, 1/5, 1/4, 3/10, 2/5, 1/2, 3/4, 1];
    x_eval = linspace(-1.0, 1.0, 4001);
    err_table = zeros(size(eps_table));
    for k = 1:numel(eps_table)
        u4_v = u4_func(x_eval, eps_table(k));
        ue_v = u_exact_func(x_eval, eps_table(k));
        err_table(k) = max(abs(u4_v - ue_v));
    end

    e_show = [0.1, 0.05];

    fig = figure('Units', 'inches', 'Position', [1, 1, 11.5, 4.4], 'Color', 'w');
    tl = tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile(tl); hold on;
    for i = 1:numel(e_show)
        e = e_show(i);
        u4_v = u4_func(x_eval, e);
        ue_v = u_exact_func(x_eval, e);
        ls = '-'; if i == 2; ls = '--'; end
        plot(x_eval, ue_v, 'Color', cm.NAVY, 'LineWidth', 1.2, 'LineStyle', ls, ...
             'DisplayName', sprintf('$u_{\\mathrm{exact}}, \\epsilon = 1/%d$', round(1/e)));
        if i == 1
            plot(x_eval(1:80:end), u4_v(1:80:end), 'o', 'Color', cm.CORAL, ...
                 'MarkerFaceColor', 'w', 'MarkerSize', 5, ...
                 'DisplayName', sprintf('$u_4, \\epsilon = 1/%d$', round(1/e)));
        else
            plot(x_eval(1:80:end), u4_v(1:80:end), 'o', 'Color', cm.CORAL, ...
                 'MarkerFaceColor', 'w', 'MarkerSize', 5, 'HandleVisibility', 'off');
        end
    end
    hold off;
    xlabel('$x$', 'Interpreter', 'latex');
    ylabel('$u(x)$', 'Interpreter', 'latex');
    title('$\epsilon^2 u'''' - u = -1$:  4-term Galerkin vs exact', 'Interpreter', 'latex');
    legend('Location', 'southwest', 'Interpreter', 'latex', 'FontSize', 9);
    grid on; box on;

    nexttile(tl);
    loglog(eps_table, err_table, 'o-', 'Color', cm.NAVY, 'MarkerFaceColor', 'w', ...
           'MarkerSize', 6, 'LineWidth', 1.0, ...
           'DisplayName', 'max $|u_4 - u_{\mathrm{exact}}|$');
    set(gca, 'XScale', 'log', 'YScale', 'log');
    xlabel('$\epsilon$', 'Interpreter', 'latex');
    ylabel('$L^\infty$ error', 'Interpreter', 'latex');
    title('Boyd Table 20.3:  4-term symbolic across $\epsilon$', 'Interpreter', 'latex');
    legend('Location', 'southwest', 'Interpreter', 'latex', 'FontSize', 10);
    grid on; box on;

    exportgraphics(fig, fullfile(out_dir, 'symbolic_boundary_layer.pdf'), 'ContentType', 'vector');
    exportgraphics(fig, fullfile(out_dir, 'symbolic_boundary_layer.png'), 'Resolution', 300);
    close(fig);

    fprintf('\n[Etude 21.9]  symbolic boundary-layer BVP, 4-term Galerkin\n');
    for k = 1:numel(eps_table)
        fprintf('   eps = %.5f   L_inf err = %.5e\n', eps_table(k), err_table(k));
    end
    fprintf('  figure: %s\n', fullfile(out_dir, 'symbolic_boundary_layer.pdf'));

    if ~isempty(dump_path)
        r = struct('epsilon', eps_table, 'Linf_error', err_table);
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

function [u4_subs, sol, x, eps_sym] = galerkin_4_term()
    syms x eps_sym real
    a0 = sym('a_0'); a2 = sym('a_2'); a4 = sym('a_4'); a6 = sym('a_6');
    p = a0 + a2 * x^2 + a4 * x^4 + a6 * x^6;
    u4 = (1 - x^2) * p;
    R4 = eps_sym^2 * diff(u4, x, 2) - u4 + 1;
    eqs = sym(zeros(4, 1));
    for j = 0:3
        eqs(j + 1) = int(x^(2 * j) * R4, x, -1, 1);
    end
    s = solve(eqs, [a0 a2 a4 a6]);
    sol = struct('a_0', s.a_0, 'a_2', s.a_2, 'a_4', s.a_4, 'a_6', s.a_6);
    u4_subs = subs(u4, [a0, a2, a4, a6], [s.a_0, s.a_2, s.a_4, s.a_6]);
end
