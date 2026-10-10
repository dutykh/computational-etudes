function tau_first_order(varargin)
%% tau_first_order - Etude 21.11: tau method on u' + u = 0, u(-1) = 1.
% Modified problem v' + v = tau T_N(x) admits an exact polynomial solution
% of degree N.  Compare with standard pseudospectral.
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
% Part of "Computational Etudes: A Spectral Approach"

    dump_path = parse_args(varargin{:});
    script_dir = fileparts(mfilename('fullpath'));
    addpath(script_dir);
    addpath(fullfile(script_dir, '..', 'ch07'));
    cm = tricks_common(); cm.configure_style();
    out_dir = cm.output_dir(script_dir);

    Ns = [4 6 8 12 16 24 32];
    tau_vals = zeros(size(Ns));
    err_tau  = zeros(size(Ns));
    err_pseu = zeros(size(Ns));
    for k = 1:numel(Ns)
        [x, v, tau] = tau_solve(Ns(k));
        tau_vals(k) = abs(tau);
        u_exact = exp(-(x + 1.0));
        err_tau(k) = max(abs(v - u_exact));
        [x2, u_ps] = standard_pseudospectral(Ns(k));
        err_pseu(k) = max(abs(u_ps - exp(-(x2 + 1.0))));
    end

    N_show = 16;
    [x_show, v_show, tau_show] = tau_solve(N_show);
    x_dense = linspace(-1, 1, 401);
    u_exact_dense = exp(-(x_dense + 1.0));

    fig = figure('Units', 'inches', 'Position', [1, 1, 11.5, 4.4], 'Color', 'w');
    tl = tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile(tl); hold on;
    semilogy(Ns, tau_vals, 'o-', 'Color', cm.NAVY, 'MarkerFaceColor', 'w', ...
             'MarkerSize', 6, 'LineWidth', 1.0, ...
             'DisplayName', '$|\tau|$ (Lanczos perturbation)');
    semilogy(Ns, err_tau, 's--', 'Color', cm.TEAL, 'MarkerFaceColor', 'w', ...
             'MarkerSize', 5, 'LineWidth', 0.9, ...
             'DisplayName', 'max $|v_N - u_\mathrm{exact}|$');
    semilogy(Ns, err_pseu, '^:', 'Color', cm.CORAL, 'MarkerFaceColor', 'w', ...
             'MarkerSize', 5, 'LineWidth', 0.9, ...
             'DisplayName', 'max $|u_N^{\mathrm{p\text{-}s}} - u_\mathrm{exact}|$');
    set(gca, 'YScale', 'log');
    hold off;
    xlabel('$N$', 'Interpreter', 'latex');
    ylabel('$|\tau|$ and pointwise error', 'Interpreter', 'latex');
    title('geometric decay of $\tau$, geometric convergence of both methods', ...
          'Interpreter', 'latex');
    legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 9);
    grid on; box on;
    ylim([1e-16, 5]);

    nexttile(tl); hold on;
    plot(x_dense, u_exact_dense, 'Color', cm.NAVY, 'LineWidth', 1.4, ...
         'DisplayName', '$u_\mathrm{exact}(x) = e^{-(x+1)}$');
    plot(x_show, v_show, 'o', 'Color', cm.TEAL, 'MarkerFaceColor', 'w', ...
         'MarkerSize', 6, 'DisplayName', sprintf('$v_{%d}(x_j)$ tau solution', N_show));
    hold off;
    xlabel('$x$', 'Interpreter', 'latex');
    ylabel('$u(x)$', 'Interpreter', 'latex');
    title(sprintf('$N=%d$: $\\tau \\approx %.3e$, max err $\\approx$%.1e', ...
                  N_show, tau_show, err_tau(Ns == N_show)), ...
          'Interpreter', 'latex');
    legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 9);
    grid on; box on;

    exportgraphics(fig, fullfile(out_dir, 'tau_first_order.pdf'), 'ContentType', 'vector');
    exportgraphics(fig, fullfile(out_dir, 'tau_first_order.png'), 'Resolution', 300);
    close(fig);

    fprintf('[Etude 21.11]  tau method on u'' + u = 0, u(-1) = 1\n');
    for k = 1:numel(Ns)
        fprintf('  N = %2d: |tau| = %.3e, err_tau = %.3e, err_p-s = %.3e\n', ...
                Ns(k), tau_vals(k), err_tau(k), err_pseu(k));
    end
    fprintf('  figure: %s\n', fullfile(out_dir, 'tau_first_order.pdf'));

    if ~isempty(dump_path)
        r = struct('Ns', Ns, 'tau', tau_vals, 'err_tau', err_tau, ...
                   'err_pseudospectral', err_pseu);
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

function [x, v, tau] = tau_solve(N)
    [D, x] = cheb_matrix(N);
    L = D + eye(N + 1);
    TN = cos(N * acos(x)).';
    A = zeros(N + 2, N + 2);
    b = zeros(N + 2, 1);
    A(1:N+1, 1:N+1) = L;
    A(1:N+1, N+2) = -TN;
    A(N+2, :) = 0;
    A(N+2, N+1) = 1;
    b(N+2) = 1;
    sol = A \ b;
    v = sol(1:N+1);                  % column, matches x's column orientation
    tau = sol(N+2);
end

function [x, u] = standard_pseudospectral(N)
    [D, x] = cheb_matrix(N);
    L = D + eye(N + 1);
    L(N+1, :) = 0;
    L(N+1, N+1) = 1;
    rhs = zeros(N + 1, 1);
    rhs(N+1) = 1;
    u = L \ rhs;                     % column
end
