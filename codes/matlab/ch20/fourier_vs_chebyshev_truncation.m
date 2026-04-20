function fourier_vs_chebyshev_truncation()
%% fourier_vs_chebyshev_truncation - Chapter 20, Etude 20.2.
%
% Fourier vs Chebyshev for domain truncation of sech(y).  On a truncated
% infinite interval, Fourier has a pi/2 denser interior grid than
% Chebyshev, so it wins for the main "interior" where the solution
% lives.  The Gibbs overshoot at y = +/- L is bounded by the domain-
% truncation error and is irrelevant.
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)

    script_dir = fileparts(mfilename('fullpath'));
    addpath(fullfile(script_dir, '..', 'ch07'));
    out_dir = fullfile(script_dir, '..', '..', '..', ...
                       'textbook', 'figures', 'ch20', 'matlab');
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end

    L = 10;  Ns = [8 12 16 24 32 48 64 96 128];
    err_c = arrayfun(@(N) cheb_err(N, L), Ns);
    err_f = arrayfun(@(N) fourier_err(N, L), Ns);

    NAVY=[20 45 110]/255; CORAL=[231 76 60]/255; TEAL=[22 160 133]/255;
    fig = figure('Position',[100 100 1020 340],'Color','w');

    subplot(1,2,1);
    Nshow = 32;
    [~, x_GL] = cheb_matrix(Nshow);
    y_GL = L * x_GL;
    y_F = -L + 2*L*(0:2*Nshow-1)/(2*Nshow);
    y_plot = linspace(-L, L, 401);
    plot(y_plot, 1./cosh(y_plot), '-', 'Color', NAVY, 'LineWidth',1.0); hold on;
    plot(y_F, -0.08*ones(size(y_F)), 'x', 'Color', CORAL);
    plot(y_GL, -0.18*ones(size(y_GL)), 'o', 'Color', TEAL, 'MarkerSize',4);
    xlim([-L L]); ylim([-0.25 1.1]); grid on; box on;
    xlabel('y'); title('(a) Grid densities on [-L, L]');

    subplot(1,2,2);
    semilogy(Ns, err_c, '-o', 'Color', TEAL, 'LineWidth',1.1); hold on;
    semilogy(Ns, err_f, '-s', 'Color', CORAL, 'LineWidth',1.1);
    xlabel('N'); ylabel('interior max error'); grid on; box on;
    title('(b) Interior error at fixed L=10');
    legend({'Chebyshev','Fourier'});

    print(fig, fullfile(out_dir, 'fourier_vs_chebyshev_truncation.pdf'), '-dpdf');
    print(fig, fullfile(out_dir, 'fourier_vs_chebyshev_truncation.png'), '-dpng');
    fprintf('[20.2-matlab] figure saved\n');
end

function e = cheb_err(N, L)
    [~, x] = cheb_matrix(N);
    fv = 1 ./ cosh(L * x);
    a = dct1(fv);
    y = linspace(-L+1, L-1, 4001);
    approx = cheb_eval(a, y/L, N);
    e = max(abs(approx - 1./cosh(y)));
end

function e = fourier_err(N, L)
    M = 2 * N;
    yj = -L + 2*L*(0:M-1)/M;
    fj = 1./cosh(yj);
    c = fft(fj) / M;
    k = [0:M/2-1, -M/2:-1];
    y = linspace(-L+1, L-1, 4001);
    t = pi * (y + L) / L;
    vals = zeros(size(y));
    for m = 1:M
        vals = vals + real(c(m) * exp(1i*k(m)*t));
    end
    e = max(abs(vals - 1./cosh(y)));
end

function a = dct1(v)
    N = length(v) - 1;
    V = [v(:); v(N:-1:2)'];
    A = real(fft(V))/N; A(1) = A(1)/2; A(N+1) = A(N+1)/2;
    a = A(1:N+1);
end

function v = cheb_eval(a, x, N)
    T0 = ones(size(x)); T1 = x;
    v = a(1)*T0 + a(2)*T1;
    for n = 2:N
        Tk = 2*x.*T1 - T0;
        v = v + a(n+1)*Tk;
        T0 = T1; T1 = Tk;
    end
end
