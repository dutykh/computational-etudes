function corner_tensor_clustering()
%% corner_tensor_clustering - Chapter 19, Computational Etude 19.6.
%
% 2D Poisson on the unit square with weak corner branch points:
% compares unmapped Chebyshev collocation against an independent
% tanh clustering in each coordinate, for a sweep of alpha.
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)

    script_dir = fileparts(mfilename('fullpath'));
    addpath(fullfile(script_dir, '..', 'ch07'));
    out_dir = fullfile(script_dir, '..', '..', '..', ...
                       'textbook', 'figures', 'ch19', 'matlab');
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end

    % reference at N=96 unmapped
    [xref, Uref] = solve_unmapped(96);
    Ns = [12 16 20 24 32 40 48 64];
    ALPHA_LIST = [1.0 2.0 3.0];
    err_unmapped = zeros(size(Ns));
    err_mapped = zeros(length(ALPHA_LIST), length(Ns));
    for i = 1:length(Ns)
        N = Ns(i);
        [xu, Uu] = solve_unmapped(N);
        err_unmapped(i) = compare(xu, Uu, xref, Uref);
        for j = 1:length(ALPHA_LIST)
            [xm, Um] = solve_mapped(N, ALPHA_LIST(j));
            err_mapped(j, i) = compare(xm, Um, xref, Uref);
        end
    end

    NAVY=[20 45 110]/255; CORAL=[231 76 60]/255; TEAL=[22 160 133]/255;
    PURPLE=[142 68 173]/255;

    % Pointwise error fields at N_err for the bottom row
    N_err = 24;
    [xu, Uu] = solve_unmapped(N_err);
    [xf, err_um] = pw_error(xu, Uu, xref, Uref);
    [xm1, Um1]  = solve_mapped(N_err, 1.0);
    [~,    err_t1] = pw_error(xm1, Um1, xref, Uref);
    [xm3, Um3]  = solve_mapped(N_err, 3.0);
    [~,    err_t3] = pw_error(xm3, Um3, xref, Uref);
    err_floor = max(1e-8, min([nonzero_min(err_um), ...
                               nonzero_min(err_t1), nonzero_min(err_t3)]));
    err_ceil  = max([err_um(:); err_t1(:); err_t3(:)]);

    fig = figure('Position',[100 100 1350 800],'Color','w');

    %% (a) Solution contour at N=32
    subplot(2,3,1);
    [xs, Us] = solve_unmapped(32);
    [~, idx] = sort(xs);
    contourf(xs(idx), xs(idx), Us(idx, idx), 15);
    xlabel('x'); ylabel('y'); axis equal tight;
    title('(a) $-\Delta u = 1$,  $u|_{\partial\Omega} = 0$', ...
          'Interpreter','latex');
    colorbar;

    %% (b) 1-D grids at N=24
    subplot(2,3,2);
    [~, xi24] = cheb_matrix(24);
    x_plain = xi24;
    x_tanh = tanh(2.0*xi24)/tanh(2.0);
    hold on;
    for k = 1:length(x_plain)
        plot([x_plain(k) x_plain(k)], [0 0.45], 'Color', NAVY);
    end
    for k = 1:length(x_tanh)
        plot([x_tanh(k) x_tanh(k)], [0.55 1.0], 'Color', CORAL);
    end
    text(-0.95, 0.22, 'standard Chebyshev', 'Color', NAVY);
    text(-0.95, 0.78, 'tanh clustered, $\alpha = 2$', ...
         'Color', CORAL, 'Interpreter','latex');
    xlim([-1.05 1.05]); ylim([0 1]); yticks([]); box on;
    xlabel('physical coordinate');
    title('(b) Grids at $N = 24$', 'Interpreter','latex');

    %% (c) Corner convergence
    subplot(2,3,3);
    cmap = [NAVY; CORAL; TEAL; PURPLE]; hold on;
    loglog(Ns, err_unmapped+1e-18, '-o', 'Color', cmap(1,:), 'LineWidth',1.2);
    for j = 1:length(ALPHA_LIST)
        loglog(Ns, err_mapped(j,:)+1e-18, '-s', 'Color', cmap(j+1,:), ...
               'MarkerFaceColor','none');
    end
    set(gca, 'XScale','log','YScale','log'); grid on; box on;
    xlabel('$N$', 'Interpreter','latex');
    ylabel('error vs reference');
    title('(c) Corner convergence', 'Interpreter','latex');
    lg = [{'unmapped'}, arrayfun(@(a) sprintf('tanh \\alpha = %.0f', a), ...
          ALPHA_LIST, 'UniformOutput', false)];
    legend(lg, 'Location', 'best');

    %% (d)/(e)/(f) Pointwise error contours, common log scale
    cmin_log = log10(err_floor);
    cmax_log = log10(err_ceil + 1e-30);

    subplot(2,3,4);
    contourf(xf, xf, log10(err_um + err_floor), 20, 'LineStyle','none');
    caxis([cmin_log cmax_log]); colormap(gca, flipud(hot(64)));
    axis equal tight; xlabel('x'); ylabel('y');
    title(sprintf('(d) $|u_N - u_{\\rm ref}|$, unmapped, $N = %d$', N_err), ...
          'Interpreter','latex');

    subplot(2,3,5);
    contourf(xf, xf, log10(err_t1 + err_floor), 20, 'LineStyle','none');
    caxis([cmin_log cmax_log]); colormap(gca, flipud(hot(64)));
    axis equal tight; xlabel('x'); ylabel('y');
    title(sprintf('(e) tanh $\\alpha = 1$, $N = %d$', N_err), ...
          'Interpreter','latex');

    subplot(2,3,6);
    contourf(xf, xf, log10(err_t3 + err_floor), 20, 'LineStyle','none');
    caxis([cmin_log cmax_log]); colormap(gca, flipud(hot(64)));
    axis equal tight; xlabel('x'); ylabel('y');
    title(sprintf('(f) tanh $\\alpha = 3$, $N = %d$', N_err), ...
          'Interpreter','latex');
    cb = colorbar;
    cb.Label.String = 'log_{10}|u_N - u_{ref}|';

    set(fig, 'PaperPositionMode','auto');
    pos = get(fig, 'Position');
    set(fig, 'PaperUnits','points', ...
             'PaperSize',[pos(3) pos(4)], ...
             'PaperPosition',[0 0 pos(3) pos(4)]);
    print(fig, fullfile(out_dir, 'corner_tensor_clustering.pdf'), '-dpdf');
    print(fig, fullfile(out_dir, 'corner_tensor_clustering.png'), '-dpng');
    close(fig);
    fprintf('[19.6-matlab] figure saved\n');
end

function [xf, err] = pw_error(x, U, xref, Uref)
    if x(1) > x(end)
        xs = flipud(x); Us = flipud(fliplr(U));
    else
        xs = x; Us = U;
    end
    if xref(1) > xref(end)
        xrs = flipud(xref); Urs = flipud(fliplr(Uref));
    else
        xrs = xref; Urs = Uref;
    end
    xf = linspace(-0.995, 0.995, 121);
    [XF, YF] = meshgrid(xf, xf);
    Vm = interp2(xs, xs, Us, XF, YF, 'linear', 0);
    Vr = interp2(xrs, xrs, Urs, XF, YF, 'linear', 0);
    err = abs(Vm - Vr);
end

function v = nonzero_min(arr)
    a = arr(arr > 0);
    if isempty(a); v = 1e-12; else; v = min(a(:)); end
end

function [x, U] = solve_unmapped(N)
    [D, x] = cheb_matrix(N);
    D2 = D*D; I = eye(N+1); nx = N+1;
    L = kron(D2, I) + kron(I, D2);
    idx = reshape(1:nx*nx, nx, nx);
    interior = false(nx, nx); interior(2:N, 2:N) = true;
    idx_int = idx(interior);
    A = L(idx_int, idx_int);
    b = -ones(length(idx_int), 1);
    u_int = A \ b;
    u = zeros(nx*nx, 1); u(idx_int) = u_int;
    U = reshape(u, nx, nx);
end

function [x_phys, U] = solve_mapped(N, alpha)
    [D, xi] = cheb_matrix(N);
    th_a = tanh(alpha);
    x_phys = tanh(alpha*xi)/th_a;
    fp  = alpha ./ (cosh(alpha*xi).^2 * th_a);
    fpp = -2*alpha^2*tanh(alpha*xi) ./ (cosh(alpha*xi).^2 * th_a);
    D1 = diag(1./fp)*D;
    D2 = diag(1./fp.^2)*(D*D) - diag(fpp./fp.^3)*D;
    I = eye(N+1); nx = N+1;
    L = kron(D2, I) + kron(I, D2);
    idx = reshape(1:nx*nx, nx, nx);
    interior = false(nx, nx); interior(2:N, 2:N) = true;
    idx_int = idx(interior);
    A = L(idx_int, idx_int);
    b = -ones(length(idx_int), 1);
    u_int = A \ b;
    u = zeros(nx*nx, 1); u(idx_int) = u_int;
    U = reshape(u, nx, nx);
end

function e = compare(x, U, xref, Uref)
    if x(1) > x(end)
        xs = flipud(x); Us = flipud(fliplr(U));
    else
        xs = x; Us = U;
    end
    if xref(1) > xref(end)
        xrs = flipud(xref); Urs = flipud(fliplr(Uref));
    else
        xrs = xref; Urs = Uref;
    end
    xf = linspace(-0.99, 0.99, 121);
    [XF, YF] = meshgrid(xf, xf);
    Vm = interp2(xs, xs, Us, XF, YF, 'linear', 0);
    Vr = interp2(xrs, xrs, Urs, XF, YF, 'linear', 0);
    e = max(abs(Vm(:) - Vr(:)));
end
