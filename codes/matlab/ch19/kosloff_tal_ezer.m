function kosloff_tal_ezer()
%% kosloff_tal_ezer - Chapter 19, Computational Etude 19.8.
%
% Minimum grid spacing and spectral radius of the first-derivative
% matrix under the Kosloff-Tal-Ezer arcsine map, and Chebyshev-
% coefficient |a_N| of f(y) = y under two policies:
%   aggressive     beta = 1 - cos(1/N)   (Boyd Fig 16.3)
%   conservative   beta = 1 - cos(1/2)   (Hesthaven et al. 1999)
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)

    script_dir = fileparts(mfilename('fullpath'));
    addpath(fullfile(script_dir, '..', 'ch07'));
    out_dir = fullfile(script_dir, '..', '..', '..', ...
                       'textbook', 'figures', 'ch19', 'matlab');
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end

    Ns = 8:2:96;
    min_std = zeros(size(Ns)); min_opt = zeros(size(Ns)); min_cons = zeros(size(Ns));
    stiff_std = zeros(size(Ns)); stiff_opt = zeros(size(Ns)); stiff_cons = zeros(size(Ns));
    for i = 1:length(Ns)
        N = Ns(i);
        [y, ~, D] = kte_grid(N, 0);
        min_std(i) = min(diff(sort(y)));
        stiff_std(i) = max(abs(eig(D)));
        beta_opt = 1 - cos(1/N);
        [y, ~, ~, Dy] = kte_derivative(N, beta_opt);
        min_opt(i) = min(diff(sort(y)));
        stiff_opt(i) = max(abs(eig(Dy)));
        beta_cons = 1 - cos(0.5);
        [y, ~, ~, Dy] = kte_derivative(N, beta_cons);
        min_cons(i) = min(diff(sort(y)));
        stiff_cons(i) = max(abs(eig(Dy)));
    end

    Ns_c = 3:2:63;
    aN_agg = zeros(size(Ns_c)); aN_cons = zeros(size(Ns_c));
    for i = 1:length(Ns_c)
        N = Ns_c(i);
        [y_a, ~, ~] = kte_grid(N, 1 - cos(1/N));
        aN_agg(i) = abs(cheb_coeff_last(y_a));
        [y_c, ~, ~] = kte_grid(N, 1 - cos(0.5));
        aN_cons(i) = abs(cheb_coeff_last(y_c));
    end

    NAVY=[20 45 110]/255; CORAL=[231 76 60]/255; TEAL=[22 160 133]/255;

    fig = figure('Position',[100 100 1100 760],'Color','w');

    %% (a) NEW: KTE map shape y(xi) for the three regimes
    subplot(2,2,1); hold on;
    Ndemo = 32;
    [~, xi_demo] = cheb_matrix(Ndemo);
    x_line = linspace(-1, 1, 401);
    h_std = plot(x_line, x_line, '-', 'Color', NAVY, 'LineWidth', 1.4);
    beta_a = 1 - cos(1/Ndemo);
    y_a = asin((1 - beta_a).*x_line) / asin(1 - beta_a);
    h_agg = plot(x_line, y_a, '-', 'Color', CORAL, 'LineWidth', 1.6);
    beta_c = 1 - cos(0.5);
    y_c = asin((1 - beta_c).*x_line) / asin(1 - beta_c);
    h_con = plot(x_line, y_c, '-', 'Color', TEAL, 'LineWidth', 1.6);
    plot(xi_demo, -1.12*ones(size(xi_demo)), '|', ...
         'Color', NAVY, 'MarkerSize', 10);
    text(-1.0, -1.30, 'CGL nodes \xi_j', 'Color', NAVY, 'FontSize', 8);
    plot([-1.10 1.10], [0 0], 'Color', [0.85 0.85 0.85], 'LineWidth', 0.4);
    plot([0 0], [-1.40 1.10], 'Color', [0.85 0.85 0.85], 'LineWidth', 0.4);
    xlim([-1.10 1.10]); ylim([-1.40 1.10]);
    axis equal; box on;
    xlabel('computational coordinate $\xi$', 'Interpreter','latex');
    ylabel('physical coordinate $y$', 'Interpreter','latex');
    title('(a) KTE map  $y = \arcsin((1-\beta)\xi)/\arcsin(1-\beta)$', ...
          'Interpreter','latex');
    legend([h_std, h_agg, h_con], ...
           {'standard ($\beta = 0$): identity', ...
            sprintf('aggressive $\\beta \\sim 1/N^2$ ($N = %d$)', Ndemo), ...
            'conservative $\beta = 1 - \cos(1/2)$'}, ...
           'Interpreter','latex', 'Location','northwest', 'Box','off');

    %% (b) Minimum spacing vs N
    subplot(2,2,2);
    loglog(Ns, min_std, '-o', 'Color', NAVY); hold on;
    loglog(Ns, min_opt, '-s', 'Color', CORAL);
    loglog(Ns, min_cons, '-^', 'Color', TEAL);
    loglog(Ns, 2./Ns, ':', 'Color', [0.5 0.5 0.5]);
    xlabel('$N$', 'Interpreter','latex');
    ylabel('$\min_j |y_{j+1} - y_j|$', 'Interpreter','latex');
    grid on; box on;
    legend({'standard Chebyshev', 'KTE $\beta \sim 1/N^2$', ...
            'KTE $\beta = 1 - \cos(1/2)$', '$2/N$'}, ...
           'Interpreter','latex', 'Location','best');
    title('(b) Minimum spacing vs $N$', 'Interpreter','latex');

    %% (c) Stiffness rho(D) vs N
    subplot(2,2,3);
    loglog(Ns, stiff_std, '-o', 'Color', NAVY); hold on;
    loglog(Ns, stiff_opt, '-s', 'Color', CORAL);
    loglog(Ns, stiff_cons, '-^', 'Color', TEAL);
    xlabel('$N$', 'Interpreter','latex');
    ylabel('$\rho(D)$', 'Interpreter','latex');
    grid on; box on;
    legend({'standard Chebyshev', 'KTE $\beta \sim 1/N^2$', ...
            'KTE $\beta = 1 - \cos(1/2)$'}, ...
           'Interpreter','latex', 'Location','best');
    title('(c) Stiffness $\rho(D)$ vs $N$', 'Interpreter','latex');

    %% (d) Coefficient |a_N| of f(y) = y
    subplot(2,2,4);
    loglog(Ns_c, aN_agg+1e-18, 'o', 'Color', CORAL, 'MarkerSize', 5); hold on;
    loglog(Ns_c, aN_cons+1e-18, 's', 'Color', TEAL, 'MarkerSize', 5);
    loglog(Ns_c, 0.488./Ns_c.^2, ':', 'Color', NAVY, 'LineWidth', 0.9);
    xlabel('$N$', 'Interpreter','latex');
    ylabel('$|a_N|$ of $f(y) = y$', 'Interpreter','latex');
    grid on; box on;
    legend({'aggressive $\beta = 1 - \cos(1/N)$', ...
            'conservative $\beta = 1 - \cos(1/2)$', ...
            '$0.488/N^2$ (Boyd)'}, ...
           'Interpreter','latex', 'Location','best');
    title('(d) Accuracy destroyed by aggressive scaling', ...
          'Interpreter','latex');

    set(fig, 'PaperPositionMode','auto');
    pos = get(fig, 'Position');
    set(fig, 'PaperUnits','points', ...
             'PaperSize',[pos(3) pos(4)], ...
             'PaperPosition',[0 0 pos(3) pos(4)]);
    print(fig, fullfile(out_dir, 'kosloff_tal_ezer.pdf'), '-dpdf');
    print(fig, fullfile(out_dir, 'kosloff_tal_ezer.png'), '-dpng');
    close(fig);
    fprintf('[19.8-matlab] figure saved\n');
end

function [y, xi, D] = kte_grid(N, beta)
    [D, xi] = cheb_matrix(N);
    denom = asin(1 - beta);
    if denom == 0, y = xi; return; end
    y = asin((1 - beta).*xi) / denom;
end

function [y, xi, D, Dy] = kte_derivative(N, beta)
    [D, xi] = cheb_matrix(N);
    denom = asin(1 - beta);
    y = asin((1 - beta).*xi) / denom;
    fp = (1 - beta) ./ (sqrt(1 - (1-beta)^2 .* xi.^2) * denom);
    Dy = diag(1./fp) * D;
end

function aN = cheb_coeff_last(v)
    N = length(v) - 1;
    V = [v; v(N:-1:2)];
    A = real(fft(V))/N;
    A(1) = A(1)/2; A(N+1) = A(N+1)/2;
    aN = A(N+1);
end
