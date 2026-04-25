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
    fig = figure('Position',[100 100 1320 360],'Color','w');

    % (a) three grids
    subplot(1,3,1); hold on;
    yy = linspace(-pi, pi, 401);
    plot(yy, exp(-KAPPA*(1-cos(yy))), '-', 'Color', NAVY, 'LineWidth',1.2);
    Ng = 32; x = -pi + 2*pi*(0:Ng-1)/Ng;
    for data = {[0.1 CORAL -0.08], [0.3 TEAL -0.18], [1.0 ORANGE -0.28]}
        d = data{1};
        ell  = d(1);   col = d(2:4);   off = d(5);
        yc = 2*atan(ell*tan(x/2));
        plot(yc, off*ones(size(yc)), 'o', 'Color', col, 'MarkerSize',4);
    end
    xlim([-pi pi]); ylim([-0.35 1.1]); grid on; box on;
    title('(a) Mapped grids, N=32');

    subplot(1,3,2);
    imagesc(log10(Ls), Ns, log10(E+1e-16));
    set(gca,'YDir','normal');
    xlabel('log_{10} ell'); ylabel('N');
    title('(b) Error landscape'); colorbar;

    subplot(1,3,3);
    slice_Ns = [16 24 32 48 64];
    colours = [CORAL; ORANGE; TEAL; NAVY; PURPLE];
    hold on;
    for k = 1:length(slice_Ns)
        i = find(Ns == slice_Ns(k));
        semilogy(Ls, E(i,:)+1e-16, '-o', 'Color', colours(k,:));
    end
    set(gca,'XScale','log','YScale','log'); grid on; box on;
    xlabel('ell'); ylabel('max error');
    legend(arrayfun(@(N) sprintf('N=%d',N), slice_Ns, 'UniformOutput', false));
    title('(c) Slices at fixed N');

    print(fig, fullfile(out_dir, 'arctan_tan_sweep.pdf'), '-dpdf');
    print(fig, fullfile(out_dir, 'arctan_tan_sweep.png'), '-dpng');
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
