function arctan_tan_sweep()
%% arctan_tan_sweep - Chapter 19, Computational Etude 19.7.
%
% Sweep the map parameter ell of the 2-pi-periodic arctan/tan map over
% a range of values, for each N, and plot the error landscape.
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)

    script_dir = fileparts(mfilename('fullpath'));
    out_dir = fullfile(script_dir, '..', '..', '..', ...
                       'textbook', 'figures', 'ch19', 'matlab');
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end

    KAPPA = 80;
    Ns = [12 16 24 32 48 64 96];
    Ls = [0.08 0.10 0.15 0.20 0.30 0.45 0.60 0.80 1.0 1.5];
    y_eval = linspace(-pi+1e-9, pi-1e-9, 4097);
    truth  = exp(-KAPPA*(1 - cos(y_eval)));
    E = zeros(length(Ns), length(Ls));
    for i = 1:length(Ns)
        for j = 1:length(Ls)
            E(i,j) = eval_error(Ns(i), Ls(j), y_eval, truth, KAPPA);
        end
    end

    NAVY=[20 45 110]/255; CORAL=[231 76 60]/255; TEAL=[22 160 133]/255;
    ORANGE=[230 126 34]/255; PURPLE=[142 68 173]/255;
    GREY=[0.65 0.65 0.65];

    fig = figure('Position',[100 100 1100 760],'Color','w');

    %% (a) Mapped grids in physical y, N=32
    subplot(2,2,1); hold on;
    yy = linspace(-pi, pi, 401);
    h_f = plot(yy, exp(-KAPPA*(1-cos(yy))), '-', 'Color', NAVY, 'LineWidth',1.4);
    Ng = 32; x = -pi + 2*pi*(0:Ng-1)/Ng;
    legendH = h_f;
    legendL = {'$f(y)$'};
    cols = {CORAL, TEAL, ORANGE};
    offs = [-0.08, -0.18, -0.28];
    ells = [0.1, 0.3, 1.0];
    for k = 1:3
        ell = ells(k); col = cols{k}; off = offs(k);
        yc = 2*atan(ell*tan(x/2));
        h = plot(yc, off*ones(size(yc)), 'o', 'Color', col, 'MarkerSize', 4);
        legendH(end+1) = h;
        legendL{end+1} = sprintf('$\\ell = %g$', ell);
    end
    xlim([-pi pi]); ylim([-0.35 1.18]); grid on; box on;
    xlabel('physical $y$', 'Interpreter','latex');
    title('(a) Mapped grids in physical $y$, $N = 32$', ...
          'Interpreter','latex');
    legend(legendH, legendL, 'Interpreter','latex', ...
           'Location','northeast', 'Box','off');

    %% (b) NEW: pulse f(y) and broadened f(y(x)) at ell* = 0.3
    subplot(2,2,2); hold on;
    ELL_STAR = 0.3;
    x_line = linspace(-pi+1e-9, pi-1e-9, 401);
    f_y = exp(-KAPPA*(1-cos(yy)));
    f_x = exp(-KAPPA*(1-cos(2*atan(ELL_STAR*tan(x_line/2)))));
    h_orig = plot(yy, f_y, '--', 'Color', GREY, 'LineWidth', 1.0);
    h_brd  = plot(x_line, f_x, '-', 'Color', TEAL, 'LineWidth', 1.6);
    x_grid = -pi + 2*pi*(0:Ng-1)/Ng;
    h_g = plot(x_grid, -0.18*ones(size(x_grid)), 'o', ...
               'Color', TEAL, 'MarkerSize', 4);
    xlim([-pi pi]); ylim([-0.35 1.18]); grid on; box on;
    xlabel('computational coordinate $x$', 'Interpreter','latex');
    title('(b) Pulse in computational $x$ at $\ell^* = 0.3$', ...
          'Interpreter','latex');
    legend([h_orig, h_brd, h_g], ...
           {'original $f(y)$', '$\tilde f(x) = f(y(x))$, $\ell = 0.3$', ...
            'uniform $x$-grid'}, ...
           'Interpreter','latex', 'Location','northeast', 'Box','off');

    %% (c) Error landscape (N, ell)
    subplot(2,2,3);
    imagesc(log10(Ls), Ns, log10(E+1e-16));
    set(gca,'YDir','normal');
    xlabel('$\log_{10} \ell$', 'Interpreter','latex');
    ylabel('$N$', 'Interpreter','latex');
    title('(c) Error landscape', 'Interpreter','latex');
    cb = colorbar;
    cb.Label.String = 'log_{10} max error';

    %% (d) Slices at fixed N
    subplot(2,2,4);
    slice_Ns = [16 24 32 48 64];
    colours = [CORAL; ORANGE; TEAL; NAVY; PURPLE]; hold on;
    for k = 1:length(slice_Ns)
        i = find(Ns == slice_Ns(k));
        semilogy(Ls, E(i,:)+1e-16, '-o', 'Color', colours(k,:), ...
                 'LineWidth', 1.1);
    end
    set(gca,'XScale','log','YScale','log'); grid on; box on;
    xlabel('$\ell$', 'Interpreter','latex');
    ylabel('$\|f - f_N\|_\infty$', 'Interpreter','latex');
    title('(d) Slices at fixed $N$', 'Interpreter','latex');
    legend(arrayfun(@(N) sprintf('$N = %d$', N), slice_Ns, ...
                    'UniformOutput', false), ...
           'Interpreter','latex', 'Location','best');

    set(fig, 'PaperPositionMode','auto');
    pos = get(fig, 'Position');
    set(fig, 'PaperUnits','points', ...
             'PaperSize',[pos(3) pos(4)], ...
             'PaperPosition',[0 0 pos(3) pos(4)]);
    print(fig, fullfile(out_dir, 'arctan_tan_sweep.pdf'), '-dpdf');
    print(fig, fullfile(out_dir, 'arctan_tan_sweep.png'), '-dpng');
    close(fig);
    fprintf('[19.7-matlab] figure saved\n');
end

function e = eval_error(N, ell, y_eval, truth, kappa)
    x  = -pi + 2*pi*(0:N-1)/N;
    y  = 2*atan(ell*tan(x/2));
    fv = exp(-kappa*(1 - cos(y)));
    coeffs = fft(fv)/N;
    k = [0:N/2-1, -N/2:-1];
    x_eval = 2*atan(tan(y_eval/2)/ell);
    vals = zeros(size(y_eval));
    for m = 1:N
        vals = vals + real(coeffs(m) * exp(1i*k(m)*(x_eval + pi)));
    end
    e = max(abs(vals - truth));
end
