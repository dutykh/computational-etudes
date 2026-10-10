function periodic_pulse_two_grids()
%% periodic_pulse_two_grids - Chapter 19, Computational Etude 19.1.
%
% Periodic localised pulse f(y) = exp(-kappa (1 - cos y)) approximated on
%
%   (1) a uniform Fourier grid;
%   (2) the 2-pi-periodic arctan/tan map
%         y = 2 arctan( ell tan(x/2) ),
%       using a uniform computational grid in x.
%
% Compares max-norm interpolation error and coefficient decay.  This is
% the chapter's opening shock: the computational coordinate is part of the
% numerical method.
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
% Part of "Computational Etudes: A Spectral Approach"

    script_dir = fileparts(mfilename('fullpath'));
    out_dir = fullfile(script_dir, '..', '..', '..', ...
                       'textbook', 'figures', 'ch19', 'matlab');
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end

    KAPPA = 80; ell = 0.3;
    Ns = [8 12 16 24 32 48 64 96 128];
    y_eval = linspace(-pi + 1e-9, pi - 1e-9, 4097);
    truth  = target(y_eval, KAPPA);

    err_U = zeros(size(Ns));
    err_M = zeros(size(Ns));
    for i = 1:length(Ns)
        N = Ns(i);
        yU = -pi + 2*pi*(0:N-1)/N;
        fU = target(yU, KAPPA);
        approx_U = fourier_interp(yU, fU, y_eval);
        err_U(i) = max(abs(approx_U - truth));

        x  = -pi + 2*pi*(0:N-1)/N;
        yM = 2*atan(ell*tan(x/2));
        fM = target(yM, KAPPA);
        approx_M = mapped_interp(x, fM, y_eval, ell);
        err_M(i) = max(abs(approx_M - truth));
    end

    % coefficient decay at N=96
    N = 96;
    yU = -pi + 2*pi*(0:N-1)/N;
    cU = abs(fftshift(fft(target(yU, KAPPA))/N));
    x  = -pi + 2*pi*(0:N-1)/N;
    yM = 2*atan(ell*tan(x/2));
    cM = abs(fftshift(fft(target(yM, KAPPA))/N));
    k  = fftshift((0:N-1) - (0:N-1 >= N/2)*N);

    NAVY=[20 45 110]/255; CORAL=[231 76 60]/255; TEAL=[22 160 133]/255;
    GREY = [0.65 0.65 0.65];

    fig = figure('Position',[100 100 1100 760],'Color','w');

    %% (a) Pulse and grids in physical y
    subplot(2,2,1); hold on;
    y_line = linspace(-pi, pi, 401);
    h_f = plot(y_line, target(y_line, KAPPA), '-', 'Color', NAVY, 'LineWidth',1.4);
    Nshow = 32;
    yU32 = -pi + 2*pi*(0:Nshow-1)/Nshow;
    x32  = yU32;
    yM32 = 2*atan(ell*tan(x32/2));
    h_uU = plot(yU32, -0.08*ones(size(yU32)), 'x', 'Color', CORAL, 'MarkerSize',7);
    h_uM = plot(yM32, -0.18*ones(size(yM32)), 'o', 'Color', TEAL, 'MarkerSize',5);
    xlim([-pi pi]); ylim([-0.28 1.18]);
    xlabel('$y$', 'Interpreter','latex');
    ylabel('$f(y)$, grid points', 'Interpreter','latex');
    title('(a) Pulse and grids in physical $y$', 'Interpreter','latex');
    legend([h_f, h_uU, h_uM], {'$f(y)$', sprintf('uniform, $N=%d$', Nshow), ...
            sprintf('arctan/tan, $\\ell=%g$', ell)}, ...
            'Location','northwest', 'Interpreter','latex', 'Box','off');
    box on; grid on;

    %% (b) Pulse in computational coordinate x  -- the punchline
    subplot(2,2,2); hold on;
    x_line = linspace(-pi + 1e-9, pi - 1e-9, 1001);
    f_of_x = target(2*atan(ell*tan(x_line/2)), KAPPA);
    h_ref  = plot(y_line, target(y_line, KAPPA), '--', 'Color', GREY, 'LineWidth', 0.9);
    h_ftil = plot(x_line, f_of_x, '-', 'Color', TEAL, 'LineWidth', 1.6);
    h_xg   = plot(x32, -0.18*ones(size(x32)), 'o', 'Color', TEAL, 'MarkerSize', 5);
    xlim([-pi pi]); ylim([-0.28 1.18]);
    xlabel('$x$', 'Interpreter','latex');
    ylabel('$\tilde f(x)$, grid points', 'Interpreter','latex');
    title('(b) Pulse in computational $x$', 'Interpreter','latex');
    legend([h_ref, h_ftil, h_xg], ...
        {'original $f$ (physical $y$)', '$\tilde f(x) = f(y(x))$', ...
         sprintf('uniform $x$-grid, $N=%d$', Nshow)}, ...
        'Location','northwest', 'Interpreter','latex', 'Box','off');
    box on; grid on;

    %% (c) Convergence
    subplot(2,2,3);
    semilogy(Ns, err_U, '-o', 'Color', CORAL, 'LineWidth',1.1); hold on;
    semilogy(Ns, err_M, '-s', 'Color', TEAL,  'LineWidth',1.1);
    xlabel('$N$', 'Interpreter','latex');
    ylabel('$\|f - f_N\|_\infty$', 'Interpreter','latex');
    title('(c) Convergence', 'Interpreter','latex');
    box on; grid on;
    legend({'uniform Fourier', sprintf('arctan/tan, $\\ell=%g$', ell)}, ...
            'Location','best', 'Interpreter','latex');

    %% (d) Coefficient decay
    subplot(2,2,4);
    mask = k >= 0;
    semilogy(k(mask), cU(mask)+1e-18, 'x', 'Color', CORAL, 'MarkerSize',6); hold on;
    semilogy(k(mask), cM(mask)+1e-18, 'o', 'Color', TEAL, 'MarkerSize',4);
    ylim([1e-16 1]);
    xlabel('wave-number $k$', 'Interpreter','latex');
    ylabel('$|c_k|$', 'Interpreter','latex');
    title('(d) Coefficient decay, $N=96$', 'Interpreter','latex');
    box on; grid on;
    legend({'uniform Fourier','arctan/tan'}, 'Location','best');

    set(fig, 'PaperPositionMode','auto');
    pos = get(fig, 'Position');
    set(fig, 'PaperUnits','points', ...
             'PaperSize',[pos(3) pos(4)], ...
             'PaperPosition',[0 0 pos(3) pos(4)]);
    print(fig, fullfile(out_dir, 'periodic_pulse_two_grids.pdf'), '-dpdf');
    print(fig, fullfile(out_dir, 'periodic_pulse_two_grids.png'), '-dpng');
    close(fig);
    fprintf('[19.1-matlab] figure saved\n');
end

function v = target(y, kappa)
    v = exp(-kappa*(1 - cos(y)));
end

function vals = fourier_interp(y_nodes, f_nodes, y_eval)
    N = length(y_nodes);
    coeffs = fft(f_nodes)/N;
    k = [0:N/2-1, -N/2:-1];
    vals = zeros(size(y_eval));
    for m = 1:N
        vals = vals + real(coeffs(m) * exp(1i*k(m)*(y_eval + pi)));
    end
end

function vals = mapped_interp(x_nodes, f_nodes, y_eval, ell)
    x_eval = 2*atan(tan(y_eval/2) / ell);
    vals = fourier_interp(x_nodes, f_nodes, x_eval);
end
