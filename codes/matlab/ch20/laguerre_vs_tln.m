function laguerre_vs_tln()
%% laguerre_vs_tln - Chapter 20, Etude 20.5.  Laguerre vs TL_n.
% Author: Dr. Denys Dutykh.

    script_dir = fileparts(mfilename('fullpath'));
    addpath(fullfile(script_dir, '..', 'ch07'));
    out_dir = fullfile(script_dir, '..', '..', '..', ...
                       'textbook', 'figures', 'ch20', 'matlab');
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end

    Ns = [8 12 16 24 32 48 64 96];
    err_lag_e = arrayfun(@(N) lag_err(@(y) exp(-y), N), Ns);
    err_lag_a = arrayfun(@(N) lag_err(@(y) 1./(1+y), N), Ns);
    err_tln_e = arrayfun(@(N) tln_err(@(y) exp(-y), N, 2.0), Ns);
    err_tln_a = arrayfun(@(N) tln_err(@(y) 1./(1+y), N, 5.0), Ns);

    NAVY=[20 45 110]/255; CORAL=[231 76 60]/255; TEAL=[22 160 133]/255;
    fig = figure('Position',[100 100 1020 340],'Color','w');

    subplot(1,2,1);
    semilogy(Ns, err_lag_e + 1e-18, '-o', 'Color', CORAL); hold on;
    semilogy(Ns, err_tln_e + 1e-18, '-s', 'Color', TEAL);
    xlabel('N'); ylabel('max error'); grid on; box on;
    title('(a) f(y) = e^{-y} (exponential)');
    legend({'Laguerre','TL_n L=2'});

    subplot(1,2,2);
    semilogy(Ns, err_lag_a + 1e-18, '-o', 'Color', CORAL); hold on;
    semilogy(Ns, err_tln_a + 1e-18, '-s', 'Color', TEAL);
    xlabel('N'); ylabel('max error'); grid on; box on;
    title('(b) f(y) = 1/(1+y) (algebraic)');
    legend({'Laguerre','TL_n L=5'});

    print(fig, fullfile(out_dir, 'laguerre_vs_tln.pdf'), '-dpdf');
    print(fig, fullfile(out_dir, 'laguerre_vs_tln.png'), '-dpng');
    fprintf('[20.5-matlab] figure saved\n');
end

function L = lag_poly(n, y)
    y = y(:)';
    L = zeros(n+1, length(y));
    L(1,:) = 1;
    if n >= 1, L(2,:) = 1 - y; end
    for k = 1:n-1
        L(k+2,:) = ((2*k+1 - y) .* L(k+1,:) - k * L(k,:)) / (k+1);
    end
end

function e = lag_err(f, N)
    [x, w] = gauss_laguerre(N + 32);
    L = lag_poly(N, x);
    fv = f(x);
    c = sum(w .* exp(x/2) .* fv .* L, 2);
    y = linspace(0.001, 60, 4001);
    Ly = lag_poly(N, y);
    approx = exp(-y/2) .* sum(c .* Ly, 1);
    e = max(abs(approx - f(y)));
end

function e = tln_err(f, N, L)
    [~, x] = cheb_matrix(N);
    y = L * (1 + x) ./ (1 - x);
    fv = zeros(size(y));
    ok = abs(x) < 1 - 1e-12;
    fv(ok) = f(y(ok));
    a = dct1(fv);
    y_fine = linspace(0.001, 60, 4001);
    xf = (y_fine - L) ./ (y_fine + L);
    approx = cheb_eval(a, xf, N);
    e = max(abs(approx - f(y_fine)));
end

function [x, w] = gauss_laguerre(K)
    n = K; i = 1:n-1;
    alpha_i = (2*(0:n-1) + 1);
    beta_i = i;
    J = diag(alpha_i) + diag(sqrt(beta_i.*beta_i), 1) + diag(sqrt(beta_i.*beta_i), -1);
    [V, D] = eig(J);
    x = diag(D)';
    w = V(1,:).^2;
end

function a = dct1(v)
    v = v(:);
    N = length(v) - 1;
    V = [v; v(N:-1:2)];
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
