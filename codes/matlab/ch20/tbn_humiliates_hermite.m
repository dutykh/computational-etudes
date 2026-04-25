function tbn_humiliates_hermite()
%% tbn_humiliates_hermite - Chapter 20, Etude 20.6.  1/(1+y^2): TB_n wins.
% Author: Dr. Denys Dutykh.

    script_dir = fileparts(mfilename('fullpath'));
    addpath(fullfile(script_dir, '..', 'ch07'));
    out_dir = fullfile(script_dir, '..', '..', '..', ...
                       'textbook', 'figures', 'ch20', 'matlab');
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end

    Ns = [8 12 16 24 32 48 64 96 128];
    err_her = arrayfun(@(N) hermite_err(N), Ns);
    err_sinc = arrayfun(@(N) sinc_err(N), Ns);
    err_tbn = arrayfun(@(N) tbn_err(N, 1.0), Ns);

    NAVY=[20 45 110]/255; CORAL=[231 76 60]/255; TEAL=[22 160 133]/255;
    fig = figure('Position',[100 100 1020 340],'Color','w');

    subplot(1,2,1);
    y = linspace(-8, 8, 401);
    plot(y, 1./(1+y.^2), '-', 'Color', NAVY, 'LineWidth',1.2);
    xlabel('y'); ylabel('f'); grid on; box on;
    title('(a) f(y) = 1/(1+y^2)');

    subplot(1,2,2);
    loglog(Ns, err_her + 1e-18, '-o', 'Color', CORAL); hold on;
    loglog(Ns, err_sinc + 1e-18, '-s', 'Color', TEAL);
    loglog(Ns, err_tbn + 1e-18, '-^', 'Color', NAVY, 'LineWidth',1.2);
    xlabel('N'); ylabel('max error'); grid on; box on;
    title('(b) Algebraic / subgeometric / geometric');
    legend({'Hermite','sinc','TB_n ell=1'});

    print(fig, fullfile(out_dir, 'tbn_humiliates_hermite.pdf'), '-dpdf');
    print(fig, fullfile(out_dir, 'tbn_humiliates_hermite.png'), '-dpng');
    fprintf('[20.6-matlab] figure saved\n');
end

function psi = hermite_psi(n, y)
    y = y(:)';
    psi = zeros(n+1, length(y));
    psi(1,:) = pi^(-0.25) * exp(-0.5*y.^2);
    if n >= 1, psi(2,:) = sqrt(2)*y.*psi(1,:); end
    for k = 1:n-1
        psi(k+2,:) = sqrt(2/(k+1))*y.*psi(k+1,:) - sqrt(k/(k+1))*psi(k,:);
    end
end

function e = hermite_err(N)
    [x, w] = gauss_hermite(N + 40);
    fv = 1 ./ (1 + x.^2);
    psi = hermite_psi(N, x);
    integ = fv .* psi .* exp(x.^2);
    c = sum(integ .* w, 2);
    y = linspace(-40, 40, 8001);
    psi_y = hermite_psi(N, y);
    approx = sum(c .* psi_y, 1);
    e = max(abs(approx - 1./(1+y.^2)));
end

function e = sinc_err(N)
    h = sqrt(pi^2 / (2*N));
    j = -floor(N/2):floor(N/2);
    yj = j * h;
    fj = 1 ./ (1 + yj.^2);
    y = linspace(-40, 40, 8001);
    z = (y(:) - yj) / h;
    S = sin(pi*z) ./ (pi*z); S(abs(z) < 1e-14) = 1;
    approx = S * fj(:);
    e = max(abs(approx(:) - 1./(1 + y(:).^2)));
end

function e = tbn_err(N, ell)
    [~, x] = cheb_matrix(N);
    y = ell * x ./ sqrt(1 - x.^2);
    fv = zeros(size(y));
    ok = abs(x) < 1 - 1e-12;
    fv(ok) = 1 ./ (1 + y(ok).^2);
    a = dct1(fv);
    y_fine = linspace(-40, 40, 8001);
    xf = y_fine ./ sqrt(ell^2 + y_fine.^2);
    approx = cheb_eval(a, xf, N);
    e = max(abs(approx - 1./(1 + y_fine.^2)));
end

function [x, w] = gauss_hermite(K)
    n = K; i = 1:n-1;
    a = sqrt(i / 2);
    J = diag(a, 1) + diag(a, -1);
    [V, D] = eig(J);
    x = diag(D)'; w = sqrt(pi) * V(1,:).^2;
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
