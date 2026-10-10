function heal_branch_point()
%% heal_branch_point - Chapter 19, Computational Etude 19.5.
%
% Approximate g(X) = sqrt(1 - X^2) on [-1, 1] (with weak branch points
% at X = +/-1) directly, and its tanh-mapped analogue sech(y) on a
% wide window y in [-Y, Y].  The direct expansion converges only
% algebraically (~1/N); the tanh map recovers geometric convergence.
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)

    script_dir = fileparts(mfilename('fullpath'));
    addpath(fullfile(script_dir, '..', 'ch07'));
    out_dir = fullfile(script_dir, '..', '..', '..', ...
                       'textbook', 'figures', 'ch19', 'matlab');
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end

    Ns = [8 16 24 32 48 64 96 128];
    err_direct = zeros(size(Ns));
    err_mapped = zeros(size(Ns));
    for i = 1:length(Ns)
        err_direct(i) = chebyshev_error(Ns(i));
        err_mapped(i) = mapped_error(Ns(i), 10);
    end
    [~, a_direct] = chebyshev_error(64);
    [~, a_mapped] = mapped_error(64, 10);

    NAVY=[20 45 110]/255; CORAL=[231 76 60]/255; TEAL=[22 160 133]/255;

    fig = figure('Position',[100 100 1100 760],'Color','w');

    %% (a) Two views of the same function
    subplot(2,2,1); hold on;
    Xp = linspace(-1, 1, 401); yp = linspace(-6, 6, 401);
    h_g = plot(Xp, sqrt(max(1-Xp.^2,0)), '-', 'Color', CORAL, 'LineWidth',1.4);
    h_s = plot(yp, 1./cosh(yp), '--', 'Color', TEAL, 'LineWidth',1.4);
    xlabel('coordinate'); ylabel('value');
    title('(a) Two views of the same function');
    legend([h_g, h_s], {'$\sqrt{1 - X^2}$', '$\mathrm{sech}\,y$'}, ...
           'Interpreter','latex', 'Location','best', 'Box','off');
    grid on; box on;

    %% (b) NEW: pointwise error at fixed N = 32
    subplot(2,2,2); hold on;
    Nshow = 32;
    Xfine = linspace(-1+1e-6, 1-1e-6, 1001);
    err_dir_pw = chebyshev_pw_error(Nshow, Xfine);
    err_map_pw = mapped_pw_error(Nshow, 10, Xfine);
    h_d = semilogy(Xfine, err_dir_pw + 1e-18, '-', ...
                   'Color', CORAL, 'LineWidth', 1.4);
    h_m = semilogy(Xfine, err_map_pw + 1e-18, '-', ...
                   'Color', TEAL, 'LineWidth', 1.4);
    set(gca, 'YScale', 'log');
    xlim([-1 1]); ylim([1e-16 1e0]);
    xlabel('$X$', 'Interpreter','latex');
    ylabel('$|g(X) - g_N(X)|$', 'Interpreter','latex');
    title(sprintf('(b) Pointwise error at $N = %d$', Nshow), ...
          'Interpreter','latex');
    legend([h_d, h_m], {'direct, $\sqrt{1-X^2}$', ...
                         'tanh-mapped, $\mathrm{sech}\,y$'}, ...
           'Interpreter','latex', 'Location','south', 'Box','off');
    grid on; box on;

    %% (c) Convergence: algebraic vs geometric
    subplot(2,2,3);
    loglog(Ns, err_direct, '-o', 'Color', CORAL, 'LineWidth',1.1); hold on;
    loglog(Ns, err_mapped+1e-18, '-s', 'Color', TEAL, 'LineWidth',1.1);
    loglog(Ns, 1./Ns, ':', 'Color', NAVY, 'LineWidth',0.8);
    xlabel('$N$', 'Interpreter','latex');
    ylabel('$\|g - g_N\|_\infty$', 'Interpreter','latex');
    title('(c) Convergence: algebraic vs geometric', ...
          'Interpreter','latex');
    legend({'direct, $\sqrt{1-X^2}$', 'tanh-mapped', '$1/N$'}, ...
           'Interpreter','latex', 'Location', 'best', 'Box','off');
    grid on; box on;

    %% (d) Coefficient decay at N = 64
    subplot(2,2,4);
    semilogy(0:length(a_direct)-1, a_direct+1e-18, 'x', ...
             'Color', CORAL, 'MarkerSize', 5); hold on;
    semilogy(0:length(a_mapped)-1, a_mapped+1e-18, 'o', ...
             'Color', TEAL, 'MarkerSize', 4);
    xlabel('Chebyshev index $n$', 'Interpreter','latex');
    ylabel('$|a_n|$', 'Interpreter','latex');
    title('(d) Coefficient decay, $N = 64$', 'Interpreter','latex');
    ylim([1e-17 10]); grid on; box on;
    legend({'direct Chebyshev', 'tanh-mapped'}, 'Location', 'best');

    set(fig, 'PaperPositionMode','auto');
    pos = get(fig, 'Position');
    set(fig, 'PaperUnits','points', ...
             'PaperSize',[pos(3) pos(4)], ...
             'PaperPosition',[0 0 pos(3) pos(4)]);
    print(fig, fullfile(out_dir, 'heal_branch_point.pdf'), '-dpdf');
    print(fig, fullfile(out_dir, 'heal_branch_point.png'), '-dpng');
    close(fig);
    fprintf('[19.5-matlab] figure saved\n');
end

function a = coeffs(v)
    N = length(v) - 1;
    V = [v; v(N:-1:2)];
    A = real(fft(V))/N;
    A(1) = A(1)/2; A(N+1) = A(N+1)/2;
    a = A(1:N+1);
end

function v = cheb_eval(a, xfine, N)
    Tkm2 = ones(size(xfine)); Tkm1 = xfine;
    v = a(1)*Tkm2 + a(2)*Tkm1;
    for n = 2:N
        Tk = 2*xfine.*Tkm1 - Tkm2;
        v = v + a(n+1)*Tk;
        Tkm2 = Tkm1; Tkm1 = Tk;
    end
end

function [e, a_abs] = chebyshev_error(N)
    [~, x] = cheb_matrix(N);
    v = sqrt(max(1 - x.^2, 0));
    a = coeffs(v);
    xfine = linspace(-1+1e-10, 1-1e-10, 2001);
    val = cheb_eval(a, xfine, N);
    e = max(abs(val - sqrt(max(1-xfine.^2, 0))));
    a_abs = abs(a);
end

function [e, a_abs] = mapped_error(N, Y)
    [~, xi] = cheb_matrix(N);
    y = Y*xi;
    f = 1./cosh(y);
    a = coeffs(f);
    xfine = linspace(-1, 1, 2001);
    val = cheb_eval(a, xfine, N);
    yfine = Y*xfine;
    e = max(abs(val - 1./cosh(yfine)));
    a_abs = abs(a);
end

function err = chebyshev_pw_error(N, Xfine)
    [~, x] = cheb_matrix(N);
    a = coeffs(sqrt(max(1 - x.^2, 0)));
    val = cheb_eval(a, Xfine, N);
    err = abs(val - sqrt(max(1 - Xfine.^2, 0)));
end

function err = mapped_pw_error(N, Y, Xfine)
    [~, xi] = cheb_matrix(N);
    a = coeffs(1./cosh(Y*xi));
    Xc = max(min(Xfine, 1 - 1e-12), -1 + 1e-12);
    yfine = atanh(Xc);
    inside = abs(yfine) <= Y;
    xi_fine = zeros(size(yfine));
    xi_fine(inside) = yfine(inside) / Y;
    val = cheb_eval(a, xi_fine, N);
    err = abs(val - 1./cosh(yfine));
    err(~inside) = NaN;
end
