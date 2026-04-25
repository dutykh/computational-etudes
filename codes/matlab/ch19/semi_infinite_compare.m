function semi_infinite_compare()
%% semi_infinite_compare - Chapter 19, Computational Etude 19.4.
%
% Semi-infinite BVP  u'' - u = 0,  u(0)=1, u(infty)=0,  with exact
% solution u_ex(y) = exp(-y).  Two methods:
%   (I)  truncation to [0, L];
%   (II) algebraic map  y = L (1+x)/(1-x).
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)

    script_dir = fileparts(mfilename('fullpath'));
    addpath(fullfile(script_dir, '..', 'ch07'));
    out_dir = fullfile(script_dir, '..', '..', '..', ...
                       'textbook', 'figures', 'ch19', 'matlab');
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end

    Ns = [12 16 20 24 32 40 48 64];
    L_trunc = [10 20 40];     % truncation half-widths
    ell_map = [1 2 4 8];      % map parameters
    err_trunc = zeros(length(Ns), length(L_trunc));
    err_map   = zeros(length(Ns), length(ell_map));
    for i = 1:length(Ns)
        N = Ns(i);
        for j = 1:length(L_trunc)
            [y, u] = solve_truncation(N, L_trunc(j));
            err_trunc(i,j) = max(abs(u - exp(-y)));
        end
        for j = 1:length(ell_map)
            [y, u] = solve_algebraic(N, ell_map(j));
            mask = isfinite(y) & y < 1e6;
            err_map(i,j) = max(abs(u(mask) - exp(-y(mask))));
        end
    end

    NAVY=[20 45 110]/255; CORAL=[231 76 60]/255; TEAL=[22 160 133]/255;
    ORANGE=[230 126 34]/255; PURPLE=[142 68 173]/255; GOLD=[212 160 23]/255;
    fig = figure('Position',[100 100 1280 340],'Color','w');

    % (a) solution picture
    subplot(1,3,1);
    [yt, ut] = solve_truncation(24, 20);
    [ym, um] = solve_algebraic(24, 2);
    yplot = linspace(0, 12, 401);
    plot(yplot, exp(-yplot), '-', 'Color', NAVY, 'LineWidth',1.2); hold on;
    plot(yt, ut, 'o', 'Color', CORAL, 'MarkerSize',4);
    ym_c = min(ym, 12);
    plot(ym_c, um, 's', 'Color', TEAL, 'MarkerSize',4);
    xlim([0 12]); ylim([-0.05 1.1]); grid on; box on;
    xlabel('y'); ylabel('u(y)'); title('(a) Solution at N=24');
    legend({'exact','truncation L=20','algebraic ell=2'});

    subplot(1,3,2);
    colours = [CORAL; ORANGE; GOLD];
    hold on;
    for j = 1:length(L_trunc)
        semilogy(Ns, err_trunc(:,j), '-o', 'Color', colours(j,:));
    end
    set(gca, 'YScale', 'log'); grid on; box on;
    xlabel('N'); ylabel('max error'); title('(b) Truncation');
    legend(arrayfun(@(L) sprintf('L=%d', L), L_trunc, 'UniformOutput', false));

    subplot(1,3,3);
    colours = [CORAL; TEAL; PURPLE; NAVY];
    hold on;
    for j = 1:length(ell_map)
        semilogy(Ns, err_map(:,j), '-s', 'Color', colours(j,:));
    end
    set(gca, 'YScale', 'log'); grid on; box on;
    xlabel('N'); ylabel('max error'); title('(c) Algebraic map');
    legend(arrayfun(@(e) sprintf('ell=%d', e), ell_map, 'UniformOutput', false));

    print(fig, fullfile(out_dir, 'semi_infinite_compare.pdf'), '-dpdf');
    print(fig, fullfile(out_dir, 'semi_infinite_compare.png'), '-dpng');
    fprintf('[19.4-matlab] figure saved\n');
end

function [y, u] = solve_truncation(N, L)
    [D, x] = cheb_matrix(N);
    y = 0.5*L*(x+1);
    Dy = (2/L)*D;  Dy2 = Dy*Dy;
    A = Dy2(2:N, 2:N) - eye(N-1);
    rhs = -Dy2(2:N, 1)*0 - Dy2(2:N, N+1)*1;
    u_int = A \ rhs;
    u = zeros(N+1,1); u(1) = 0; u(N+1) = 1; u(2:N) = u_int;
end

function [y, u] = solve_algebraic(N, ell)
    [D, x] = cheb_matrix(N);
    fp  = 2*ell ./ (1 - x).^2;
    fpp = 4*ell ./ (1 - x).^3;
    Dy  = diag(1./fp)*D;
    Dy2 = diag(1./fp.^2)*(D*D) - diag(fpp./fp.^3)*D;
    y = ell*(1+x)./(1-x);
    A = Dy2(2:N, 2:N) - eye(N-1);
    rhs = -Dy2(2:N, 1)*0 - Dy2(2:N, N+1)*1;
    u_int = A \ rhs;
    u = zeros(N+1,1); u(1) = 0; u(N+1) = 1; u(2:N) = u_int;
end
