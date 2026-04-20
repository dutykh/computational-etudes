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
    fig = figure('Position',[100 100 1280 340],'Color','w');

    subplot(1,3,1);
    [xs, Us] = solve_unmapped(32);
    [~, idx] = sort(xs);
    contourf(xs(idx), xs(idx), Us(idx, idx), 15);
    xlabel('x'); ylabel('y'); axis equal tight;
    title('(a) -\Delta u = 1, Dirichlet');
    colorbar;

    subplot(1,3,2);
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
    text(-0.95, 0.22, 'standard', 'Color', NAVY);
    text(-0.95, 0.78, 'tanh \alpha=2', 'Color', CORAL);
    xlim([-1.05 1.05]); ylim([0 1]); yticks([]); box on;
    title('(b) Grids at N=24');

    subplot(1,3,3);
    cmap = [NAVY; CORAL; TEAL; PURPLE];
    hold on;
    loglog(Ns, err_unmapped+1e-18, '-o', 'Color', cmap(1,:), 'LineWidth',1.2);
    for j = 1:length(ALPHA_LIST)
        loglog(Ns, err_mapped(j,:)+1e-18, '-s', 'Color', cmap(j+1,:));
    end
    set(gca, 'XScale','log','YScale','log'); grid on; box on;
    xlabel('N'); ylabel('error vs reference');
    title('(c) Corner convergence');
    lg = [{'unmapped'}, arrayfun(@(a) sprintf('tanh \\alpha=%.1f',a), ...
          ALPHA_LIST, 'UniformOutput', false)];
    legend(lg);

    print(fig, fullfile(out_dir, 'corner_tensor_clustering.pdf'), '-dpdf');
    print(fig, fullfile(out_dir, 'corner_tensor_clustering.png'), '-dpng');
    fprintf('[19.6-matlab] figure saved\n');
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
