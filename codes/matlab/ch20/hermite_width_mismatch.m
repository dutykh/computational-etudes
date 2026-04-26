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

    % NEW: coefficient decay at A = 8
    A_pick = 8.0;
    f_pick = @(y) exp(-A_pick * y.^2);
    coeffs_unscaled = abs(hermite_expand(f_pick, 32, 1.0));
    coeffs_matched  = abs(hermite_expand(f_pick, 32, sqrt(2 * A_pick)));

    % NEW: optimal-alpha scan at fixed N = 16
    alphas = linspace(0.3, 5.0, 60);
    A_scan = [0.5 2.0 8.0];
    N_scan = 16;
    err_scan = zeros(length(A_scan), length(alphas));
    for ai = 1:length(A_scan)
        f = @(y) exp(-A_scan(ai) * y.^2);
        for k = 1:length(alphas)
            c = hermite_expand(f, N_scan, alphas(k));
            err_scan(ai, k) = hermite_maxerr(f, c, alphas(k));
        end
    end

    Ns_q = [1 2 4 8 16 32];
    f0 = @(y) pi^(-0.25) * exp(-0.5*y.^2);
    err_q = zeros(size(Ns_q));
    for ni = 1:length(Ns_q)
        c = hermite_expand(f0, Ns_q(ni), 1.0);
        err_q(ni) = hermite_maxerr(f0, c, 1.0);
    end

    fig = figure('Position',[100 100 1100 1200],'Color','w');

    % ---- (a) Width-mismatch cartoon -----------------------------------
    subplot(3,2,1);
    y_show = linspace(-6, 6, 401);
    psi0 = pi^(-0.25) .* exp(-0.5 .* y_show.^2);
    plot(y_show, psi0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.0, ...
         'DisplayName', '\psi_0 (basis envelope)'); hold on;
    for ai = 1:length(A_list)
        plot(y_show, exp(-A_list(ai) * y_show.^2), 'Color', cs{ai}, ...
             'LineWidth', 1.2, 'DisplayName', sprintf('A = %g', A_list(ai)));
    end
    grid on; box on; xlabel('y'); ylabel('amplitude');
    title('(a) target e^{-Ay^2} vs basis envelope \psi_0');
    legend('show', 'Location', 'northeast', 'FontSize', 8);

    % ---- (b) Coefficient decay at A = 8 --------------------------------
    subplot(3,2,2);
    ns = 0:length(coeffs_unscaled)-1;
    semilogy(ns, coeffs_unscaled + 1e-20, '-o', 'Color', CORAL, ...
             'MarkerFaceColor', 'w', 'LineWidth', 1.0, 'MarkerSize', 4, ...
             'DisplayName', '\alpha = 1 (slow algebraic decay)'); hold on;
    semilogy(ns, coeffs_matched + 1e-20, 'o', 'Color', NAVY, ...
             'MarkerFaceColor', NAVY, 'MarkerSize', 5, ...
             'DisplayName', '\alpha = sqrt(2A) (only a_0)');
    grid on; box on; xlabel('degree n'); ylabel('|a_n|');
    title(sprintf('(b) Coefficient decay at A = %g', A_pick));
    legend('show', 'Location', 'southwest', 'FontSize', 8);

    % ---- (c) Unscaled convergence ------------------------------------
    subplot(3,2,3);
    hold on; set(gca, 'YScale','log');
    for ai = 1:length(A_list)
        plot(Ns, err1(ai,:), '-o', 'Color', cs{ai}, ...
             'DisplayName', sprintf('A=%g', A_list(ai)));
    end
    grid on; box on; xlabel('N'); ylabel('max error');
    title('(c) \alpha = 1 (unscaled)');
    legend('show', 'Location', 'best');

    % ---- (d) Matched convergence -------------------------------------
    subplot(3,2,4);
    hold on; set(gca, 'YScale','log');
    for ai = 1:length(A_list)
        plot(Ns, err2(ai,:) + 1e-18, '-s', 'Color', cs{ai}, ...
             'MarkerFaceColor', 'w', ...
             'DisplayName', sprintf('A=%g', A_list(ai)));
    end
    grid on; box on; xlabel('N'); ylabel('max error');
    title('(d) \alpha = sqrt(2A) (matched: machine precision)');
    legend('show', 'Location', 'best');

    % ---- (e) Optimal-alpha scan -------------------------------------
    subplot(3,2,5);
    hold on; set(gca, 'YScale', 'log');
    cs_scan = {NAVY, TEAL, PURPLE};
    for ai = 1:length(A_scan)
        plot(alphas, err_scan(ai, :) + 1e-18, 'Color', cs_scan{ai}, ...
             'LineWidth', 1.0, 'DisplayName', sprintf('A = %g', A_scan(ai)));
        xline(sqrt(2 * A_scan(ai)), ':', 'Color', cs_scan{ai}, ...
              'LineWidth', 0.8, 'Alpha', 0.6, 'HandleVisibility', 'off');
    end
    grid on; box on; xlabel('basis scale \alpha'); ylabel('max error');
    title(sprintf('(e) Optimal-\\alpha scan at N = %d (dotted: sqrt(2A))', N_scan));
    legend('show', 'Location', 'northeast', 'FontSize', 8);

    % ---- (f) QHO ground state ---------------------------------------
    subplot(3,2,6);
    semilogy(Ns_q, err_q + 1e-20, '-D', 'Color', NAVY, ...
             'MarkerFaceColor', 'w', 'MarkerSize', 6, 'LineWidth', 1.1, ...
             'DisplayName', 'f = \psi_0');
    grid on; box on; xlabel('N'); ylabel('max error');
    title('(f) QHO ground state: basis matches physics');
    legend('show');

    set(fig, 'PaperPositionMode', 'auto');
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
