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
    fig = figure('Position',[100 100 1180 340],'Color','w');

    subplot(1,3,1); hold on;
    Xp = linspace(-1, 1, 401); yp = linspace(-6, 6, 401);
    plot(Xp, sqrt(max(1-Xp.^2,0)), '-', 'Color', CORAL, 'LineWidth',1.2);
    plot(yp, 1./cosh(yp), '--', 'Color', TEAL, 'LineWidth',1.2);
    legend({'sqrt(1-X^2)','sech(y)'});
    title('(a) Two views'); grid on; box on;

    subplot(1,3,2);
    loglog(Ns, err_direct, '-o', 'Color', CORAL); hold on;
    loglog(Ns, err_mapped+1e-18, '-s', 'Color', TEAL);
    loglog(Ns, 1./Ns, ':', 'Color', NAVY);
    xlabel('N'); ylabel('max error'); grid on; box on;
    title('(b) Algebraic vs geometric');
    legend({'direct','tanh-mapped','1/N'});

    subplot(1,3,3);
    semilogy(0:length(a_direct)-1, a_direct+1e-18, 'x', 'Color', CORAL); hold on;
    semilogy(0:length(a_mapped)-1, a_mapped+1e-18, 'o', 'Color', TEAL, 'MarkerSize',3);
    xlabel('Chebyshev index'); ylabel('|a_n|');
    title('(c) Coefficient decay, N=64'); ylim([1e-17 10]); grid on; box on;

    print(fig, fullfile(out_dir, 'heal_branch_point.pdf'), '-dpdf');
    print(fig, fullfile(out_dir, 'heal_branch_point.png'), '-dpng');
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
