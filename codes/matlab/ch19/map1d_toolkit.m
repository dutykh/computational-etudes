function map1d_toolkit()
%% map1d_toolkit - Chapter 19, Computational Etude 19.3.
%
% Demonstrates a reusable one-dimensional map abstraction through two
% validation cases:
%   * algebraic semi-infinite map  y = ell (1 + x) / (1 - x)  applied to
%     u(y) = exp(-y);
%   * tanh map  y = tanh(x)  applied to  u(y) = 1 / (1 + y^2).
%
% For each case, the first and second mapped derivatives are compared
% against known analytic values.  The abstraction ``Map1D'' stores
% forward, inverse, and first/second derivatives of f and builds the
% mapped derivative matrices from the Chebyshev D.
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)

    script_dir = fileparts(mfilename('fullpath'));
    addpath(fullfile(script_dir, '..', 'ch07'));
    out_dir = fullfile(script_dir, '..', '..', '..', ...
                       'textbook', 'figures', 'ch19', 'matlab');
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end

    Ns = [8 12 16 20 24 32 40 48 64];
    err_alg_1 = zeros(size(Ns));  err_alg_2 = zeros(size(Ns));
    err_tanh_1 = zeros(size(Ns)); err_tanh_2 = zeros(size(Ns));
    for i = 1:length(Ns)
        N = Ns(i);
        [err_alg_1(i),  err_alg_2(i)]  = convergence_alg(N, 2.0);
        [err_tanh_1(i), err_tanh_2(i)] = convergence_tanh(N);
    end

    NAVY=[20 45 110]/255; CORAL=[231 76 60]/255; TEAL=[22 160 133]/255;

    Ngrid = 24;  ELL = 2.0;  Y_LIM_ALG = 12.0;
    [~, x_demo] = cheb_matrix(Ngrid);
    x_line = linspace(-1.0, 1.0 - 1e-12, 401)';

    fig = figure('Position',[100 100 1100 760],'Color','w');

    %% (a) tanh map: y = tanh(x) with x_j ticks (bottom) + y_j ticks (right)
    subplot(2,2,1); hold on;
    y_tanh_line = tanh(x_line);
    y_tanh_grid = tanh(x_demo);
    plot(x_line, y_tanh_line, 'Color', TEAL, 'LineWidth', 1.6);
    plot(x_demo, -1.18*ones(size(x_demo)), '|', 'Color', NAVY, 'MarkerSize', 12);
    plot(1.10*ones(size(y_tanh_grid)), y_tanh_grid, '_', 'Color', TEAL, 'MarkerSize', 12);
    for j = 1:3:length(x_demo)
        plot([x_demo(j) x_demo(j) 1.10], ...
             [-1.18 y_tanh_grid(j) y_tanh_grid(j)], ...
             ':', 'Color', TEAL, 'LineWidth', 0.4);
    end
    plot([-1.20 1.20], [0 0], 'Color', [0.85 0.85 0.85], 'LineWidth', 0.4);
    plot([0 0], [-1.30 1.20], 'Color', [0.85 0.85 0.85], 'LineWidth', 0.4);
    xlim([-1.20 1.20]); ylim([-1.30 1.20]);
    xlabel('computational coordinate $x$', 'Interpreter', 'latex');
    ylabel('physical $y$', 'Interpreter', 'latex');
    title('(a) tanh map: $x_j$ (bottom) $\to y_j$ (right)', 'Interpreter', 'latex');
    legend({'$y = \tanh(x)$'}, 'Location', 'northwest', ...
        'Interpreter', 'latex', 'Box', 'off');
    box on;

    %% (b) algebraic semi-infinite map y = ell(1+x)/(1-x)
    subplot(2,2,2); hold on;
    y_alg_line = ELL*(1 + x_line)./(1 - x_line);
    visible_line = y_alg_line <= Y_LIM_ALG;
    plot(x_line(visible_line), y_alg_line(visible_line), ...
         'Color', CORAL, 'LineWidth', 1.6);
    x_finite = x_demo(1:end-1);
    y_alg_grid = ELL*(1 + x_finite)./(1 - x_finite);
    visible_grid = y_alg_grid <= Y_LIM_ALG;
    plot(x_demo, -0.7*ones(size(x_demo)), '|', 'Color', NAVY, 'MarkerSize', 12);
    plot(1.10*ones(sum(visible_grid),1), y_alg_grid(visible_grid), ...
         '_', 'Color', CORAL, 'MarkerSize', 12);
    for j = 1:length(x_finite)
        if y_alg_grid(j) <= Y_LIM_ALG && mod(j, 2) == 0
            plot([x_finite(j) x_finite(j) 1.10], ...
                 [-0.7 y_alg_grid(j) y_alg_grid(j)], ...
                 ':', 'Color', CORAL, 'LineWidth', 0.4);
        end
    end
    n_off = sum(~visible_grid) + 1;   % +1 for x=1 endpoint at infinity
    text(0.10, Y_LIM_ALG - 1.6, ...
         sprintf('$y \\to \\infty$ (%d ticks beyond view)', n_off), ...
         'Interpreter', 'latex', 'FontSize', 8, 'Color', CORAL);
    plot([-1.20 1.20], [0 0], 'Color', [0.85 0.85 0.85], 'LineWidth', 0.4);
    plot([0 0], [-1.6 Y_LIM_ALG+1.0], 'Color', [0.85 0.85 0.85], 'LineWidth', 0.4);
    xlim([-1.20 1.20]); ylim([-1.6 Y_LIM_ALG + 1.0]);
    xlabel('computational coordinate $x$', 'Interpreter', 'latex');
    ylabel('physical $y$', 'Interpreter', 'latex');
    title('(b) algebraic map ($\ell = 2$): clusters near $y = 0$', ...
          'Interpreter', 'latex');
    legend({'$y = \ell(1+x)/(1-x)$'}, 'Location', 'northwest', ...
        'Interpreter', 'latex', 'Box', 'off');
    box on;

    %% (c) First-derivative convergence
    subplot(2,2,3);
    semilogy(Ns, err_alg_1 + 1e-18, '-o', 'Color', CORAL, 'LineWidth', 1.1); hold on;
    semilogy(Ns, err_tanh_1 + 1e-18, '-s', 'Color', TEAL, 'LineWidth', 1.1);
    xlabel('$N$', 'Interpreter','latex');
    ylabel('$\|u''_N - u''\|_\infty$', 'Interpreter','latex');
    title('(c) First derivative: $D_y u$', 'Interpreter','latex');
    grid on; box on;
    legend({'algebraic', 'tanh'}, 'Location', 'best');

    %% (d) Second-derivative convergence -- the N^2-lag panel
    subplot(2,2,4);
    semilogy(Ns, err_alg_2 + 1e-18, '--o', 'Color', CORAL, 'LineWidth', 1.1); hold on;
    semilogy(Ns, err_tanh_2 + 1e-18, '--s', 'Color', TEAL, 'LineWidth', 1.1);
    xlabel('$N$', 'Interpreter','latex');
    ylabel('$\|u''''_N - u''''\|_\infty$', 'Interpreter','latex');
    title('(d) Second derivative: $D_y^2 u$ ($N^2$-lag)', 'Interpreter','latex');
    grid on; box on;
    legend({'algebraic', 'tanh'}, 'Location', 'best');

    set(fig, 'PaperPositionMode','auto');
    pos = get(fig, 'Position');
    set(fig, 'PaperUnits','points', ...
             'PaperSize',[pos(3) pos(4)], ...
             'PaperPosition',[0 0 pos(3) pos(4)]);
    print(fig, fullfile(out_dir, 'map1d_toolkit.pdf'), '-dpdf');
    print(fig, fullfile(out_dir, 'map1d_toolkit.png'), '-dpng');
    close(fig);
    fprintf('[19.3-matlab] figure saved\n');
end

function [e1, e2] = convergence_alg(N, ell)
    [Dx, x] = cheb_matrix(N);
    fp  = 2*ell ./ (1 - x).^2;
    fpp = 4*ell ./ (1 - x).^3;
    Dy  = diag(1./fp) * Dx;
    Dy2 = diag(1./fp.^2) * (Dx*Dx) - diag(fpp./fp.^3) * Dx;
    y = ell*(1+x)./(1-x);
    u = exp(-y);
    mask = (abs(y) < 1e6) & isfinite(y);
    du_ex = -exp(-y); d2u_ex = exp(-y);
    du  = Dy*u;   d2u = Dy2*u;
    e1 = max(abs(du(mask)  - du_ex(mask)));
    e2 = max(abs(d2u(mask) - d2u_ex(mask)));
end

function [e1, e2] = convergence_tanh(N)
    [Dx, x] = cheb_matrix(N);
    fp  = 1 ./ cosh(x).^2;
    fpp = -2*tanh(x) ./ cosh(x).^2;
    Dy  = diag(1./fp) * Dx;
    Dy2 = diag(1./fp.^2)*(Dx*Dx) - diag(fpp./fp.^3)*Dx;
    y = tanh(x);
    u = 1./(1 + y.^2);
    du_ex = -2*y./(1 + y.^2).^2;
    d2u_ex = (6*y.^2 - 2) ./ (1 + y.^2).^3;
    du  = Dy*u;   d2u = Dy2*u;
    e1 = max(abs(du - du_ex));
    e2 = max(abs(d2u - d2u_ex));
end
