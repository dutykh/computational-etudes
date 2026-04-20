function hermite_width_mismatch()
%% hermite_width_mismatch - Chapter 20, Etude 20.4.  Hermite width mismatch.
% Author: Dr. Denys Dutykh.

    script_dir = fileparts(mfilename('fullpath'));
    out_dir = fullfile(script_dir, '..', '..', '..', ...
                       'textbook', 'figures', 'ch20', 'matlab');
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end

    Ns = [8 12 16 24 32 48 64];
    A_list = [0.1 0.5 2.0 8.0];
    err1 = zeros(length(A_list), length(Ns));
    err2 = zeros(length(A_list), length(Ns));
    for ai = 1:length(A_list)
        A = A_list(ai);
        f = @(y) exp(-A*y.^2);
        alpha_m = sqrt(2*A);
        for ni = 1:length(Ns)
            N = Ns(ni);
            c1 = hermite_expand(f, N, 1.0);
            err1(ai, ni) = hermite_maxerr(f, c1, 1.0);
            c2 = hermite_expand(f, N, alpha_m);
            err2(ai, ni) = hermite_maxerr(f, c2, alpha_m);
        end
    end

    NAVY=[20 45 110]/255; CORAL=[231 76 60]/255; TEAL=[22 160 133]/255;
    PURPLE=[142 68 173]/255;
    cs = {CORAL, NAVY, TEAL, PURPLE};
    fig = figure('Position',[100 100 1340 340],'Color','w');

    subplot(1,3,1);
    hold on; set(gca, 'YScale','log');
    for ai = 1:length(A_list)
        plot(Ns, err1(ai,:), '-o', 'Color', cs{ai});
    end
    grid on; box on; xlabel('N'); ylabel('max error');
    title('(a) \alpha = 1 (unscaled)');
    legend(arrayfun(@(A) sprintf('A=%g', A), A_list, 'UniformOutput', false));

    subplot(1,3,2);
    hold on; set(gca, 'YScale','log');
    for ai = 1:length(A_list)
        plot(Ns, err2(ai,:) + 1e-18, '-s', 'Color', cs{ai});
    end
    grid on; box on; xlabel('N'); ylabel('max error');
    title('(b) \alpha = sqrt(2A) (matched)');

    subplot(1,3,3);
    Ns_q = [1 2 4 8 16 32];
    f0 = @(y) pi^(-0.25) * exp(-0.5*y.^2);
    err_q = zeros(size(Ns_q));
    for ni = 1:length(Ns_q)
        c = hermite_expand(f0, Ns_q(ni), 1.0);
        err_q(ni) = hermite_maxerr(f0, c, 1.0);
    end
    semilogy(Ns_q, err_q + 1e-18, '-D', 'Color', NAVY);
    grid on; box on; xlabel('N'); ylabel('max error');
    title('(c) Quantum oscillator \psi_0');

    print(fig, fullfile(out_dir, 'hermite_width_mismatch.pdf'), '-dpdf');
    print(fig, fullfile(out_dir, 'hermite_width_mismatch.png'), '-dpng');
    fprintf('[20.4-matlab] figure saved\n');
end

function psi = hermite_psi(n, y)
    y = y(:)';
    psi = zeros(n+1, length(y));
    psi(1,:) = pi^(-0.25) * exp(-0.5*y.^2);
    if n >= 1
        psi(2,:) = sqrt(2) * y .* psi(1,:);
    end
    for k = 1:n-1
        psi(k+2,:) = sqrt(2/(k+1))*y.*psi(k+1,:) - sqrt(k/(k+1))*psi(k,:);
    end
end

function c = hermite_expand(f, N, alpha)
    [x, w] = gauss_hermite(N + 32);
    y_nodes = x / alpha;
    fv = f(y_nodes);
    psi = hermite_psi(N, x);
    integ = fv .* psi .* exp(x.^2);
    c = sum(integ .* w, 2) / sqrt(alpha);
end

function e = hermite_maxerr(f, c, alpha)
    y = linspace(-15, 15, 4001);
    psi = hermite_psi(length(c)-1, alpha*y);
    approx = sqrt(alpha) * sum(c .* psi, 1);
    e = max(abs(approx - f(y)));
end

function [x, w] = gauss_hermite(K)
    % physicists' Gauss-Hermite nodes and weights via the Golub-Welsch algorithm.
    n = K; i = 1:n-1;
    a = sqrt(i / 2);
    J = diag(a, 1) + diag(a, -1);
    [V, D] = eig(J);
    x = diag(D)';
    w = sqrt(pi) * V(1,:).^2;
end
