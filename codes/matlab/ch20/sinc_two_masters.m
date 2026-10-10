function sinc_two_masters()
%% sinc_two_masters - Chapter 20, Etude 20.3.  Sinc expansions: two masters, span and spacing.
% Author: Dr. Denys Dutykh.

    script_dir = fileparts(mfilename('fullpath'));
    out_dir = fullfile(script_dir, '..', '..', '..', ...
                       'textbook', 'figures', 'ch20', 'matlab');
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end

    N_fix = 48; hs = linspace(0.15, 1.8, 80);
    err_sweep = arrayfun(@(h) sinc_err(N_fix, h), hs);

    Ns = [8 12 16 24 32 48 64 96 128 192];
    h_opts = sqrt(pi^2 ./ (2*Ns));
    err_opt = arrayfun(@(i) sinc_err(Ns(i), h_opts(i)), 1:length(Ns));

    NAVY=[20 45 110]/255; CORAL=[231 76 60]/255; TEAL=[22 160 133]/255;
    fig = figure('Position',[100 100 1100 800],'Color','w');
    [~, jstar] = min(err_sweep);
    h_theory = sqrt(pi^2/(2*N_fix));

    % ---- (a) V-shape: empirical error vs h ---------------------------
    subplot(2,2,1);
    semilogy(hs, err_sweep, 'Color', CORAL, 'LineWidth', 1.2); hold on;
    xline(hs(jstar), ':', 'Color', NAVY, 'LineWidth', 1.0, ...
          'DisplayName', sprintf('h* = %.2f empirical', hs(jstar)));
    xline(h_theory, '--', 'Color', TEAL, 'LineWidth', 1.0, ...
          'DisplayName', sprintf('sqrt(pi^2/(2N)) = %.2f', h_theory));
    xlabel('h'); ylabel('max error'); grid on; box on;
    title(sprintf('(a) Fix N=%d, vary h', N_fix));
    legend('show', 'Location', 'northwest');

    % ---- (b) NEW: two-master decomposition ---------------------------
    subplot(2,2,2);
    hs_an = linspace(min(hs), max(hs), 400);
    E_W  = (4/pi) .* exp(-pi^2 ./ (2*hs_an));
    E_DT = 4 .* exp(-N_fix .* hs_an ./ 2);
    semilogy(hs_an, E_W,  '--', 'Color', NAVY,  'LineWidth', 1.0, ...
             'DisplayName', 'E_W (bandwidth) ~ exp(-pi^2/(2h))'); hold on;
    semilogy(hs_an, E_DT, '--', 'Color', CORAL, 'LineWidth', 1.0, ...
             'DisplayName', 'E_{DT} (grid span) ~ exp(-Nh/2)');
    semilogy(hs_an, E_W + E_DT, 'Color', TEAL, 'LineWidth', 1.4, ...
             'DisplayName', 'sum: E_W + E_{DT}');
    semilogy(hs, err_sweep, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.7, ...
             'DisplayName', 'empirical V');
    xline(hs(jstar), ':', 'Color', NAVY, 'LineWidth', 0.8, 'HandleVisibility', 'off');
    xlabel('h'); ylabel('error'); grid on; box on;
    title('(b) Two masters: bandwidth vs grid-span');
    legend('show', 'Location', 'southwest', 'FontSize', 7);
    ylim([1e-12 1e1]);

    % ---- (c) Subgeometric convergence at h_opt -----------------------
    subplot(2,2,3);
    semilogy(sqrt(Ns), err_opt, '-o', 'Color', TEAL, ...
             'MarkerFaceColor', 'w', 'LineWidth', 1.1, ...
             'DisplayName', 'sinc at h = sqrt(pi^2/(2N))'); hold on;
    semilogy(sqrt(Ns), exp(-pi*sqrt(Ns/2)), ':', 'Color', NAVY, 'LineWidth', 1.0, ...
             'DisplayName', 'exp(-pi sqrt(N/2))');
    xlabel('sqrt(N)'); ylabel('max error'); grid on; box on;
    title('(c) Subgeometric convergence at h = h*');
    legend('show', 'Location', 'northeast');

    % ---- (d) Grid cartoon -------------------------------------------
    subplot(2,2,4);
    hold on;
    for k = 1:2
        Nc = [8 32]; offs = [0.7 0.3]; cs = {CORAL, TEAL};
        hk = sqrt(pi^2/(2*Nc(k)));
        j = -floor(Nc(k)/2):floor(Nc(k)/2);
        yj = j * hk;
        scatter(yj, offs(k)*ones(size(yj)), 16, cs{k}, 'x', ...
                'DisplayName', sprintf('N=%d, h=%.2f', Nc(k), hk));
    end
    y_plot = linspace(-10, 10, 401);
    plot(y_plot, 0.05 + 0.6./cosh(y_plot), 'Color', NAVY, 'LineWidth', 1.0, ...
         'HandleVisibility', 'off');
    xlim([-10 10]); ylim([0 1.2]); box on;
    xlabel('y'); title('(d) Span and spacing grow together');
    legend('show', 'Location', 'northeast');

    set(fig, 'PaperPositionMode', 'auto');
    print(fig, fullfile(out_dir, 'sinc_two_masters.pdf'), '-dpdf');
    print(fig, fullfile(out_dir, 'sinc_two_masters.png'), '-dpng');
    fprintf('[20.3-matlab] figure saved\n');
end

function e = sinc_err(N, h)
    j = -floor(N/2):floor(N/2);
    yj = j * h;
    fj = 1 ./ cosh(yj);
    y = linspace(-20, 20, 4001);
    z = (y(:) - yj) / h;
    S = sin(pi*z) ./ (pi*z);
    S(abs(z) < 1e-14) = 1;
    approx = S * fj(:);
    e = max(abs(approx(:) - 1./cosh(y(:))));
end
