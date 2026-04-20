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
    fig = figure('Position',[100 100 1360 340],'Color','w');

    subplot(1,3,1);
    semilogy(hs, err_sweep, 'Color', CORAL, 'LineWidth',1.2); hold on;
    [~, jstar] = min(err_sweep);
    xline(hs(jstar), ':', 'Color', NAVY);
    xline(sqrt(pi^2/(2*N_fix)), '--', 'Color', TEAL);
    xlabel('h'); ylabel('max error'); grid on; box on;
    title(sprintf('(a) Fix N=%d, vary h', N_fix));

    subplot(1,3,2);
    semilogy(sqrt(Ns), err_opt, '-o', 'Color', TEAL); hold on;
    semilogy(sqrt(Ns), exp(-pi*sqrt(Ns/2)), ':', 'Color', NAVY);
    xlabel('sqrt(N)'); ylabel('max error'); grid on; box on;
    title('(b) Subgeometric convergence');

    subplot(1,3,3);
    hold on;
    for k = 1:2
        Nc = [8 32]; offs = [0.7 0.3]; cs = {CORAL, TEAL};
        hk = sqrt(pi^2/(2*Nc(k)));
        j = -floor(Nc(k)/2):floor(Nc(k)/2);
        yj = j * hk;
        scatter(yj, offs(k)*ones(size(yj)), 16, cs{k}, 'x');
    end
    y_plot = linspace(-10, 10, 401);
    plot(y_plot, 0.05 + 0.6./cosh(y_plot), 'Color', NAVY, 'LineWidth',1.0);
    xlim([-10 10]); ylim([0 1.2]); box on;
    xlabel('y'); title('(c) Span and spacing grow together');

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
