function map1d_toolkit()
%% map1d_toolkit - Chapter 19, Computational Etude 19.3.
%
% Demonstrates a reusable one-dimensional map abstraction through two
% validation cases:
%   * algebraic semi-infinite map  y = L (1 + x) / (1 - x)  applied to
%     u(y) = exp(-y);
%   * tanh map  y = tanh(x)  applied to  u(y) = 1 / (1 + y^2).
%
% For each case, the first and second mapped derivatives are compared
% against known analytic values.  The abstraction ``Map1D'' stores
% forward, inverse, and first/second derivatives of f and builds the
% mapped derivative matrices from the Chebyshev D.
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)

    script_dir = fileparts(mfilename('fullpath'));
    addpath(fullfile(script_dir, '..', 'ch07'));
    out_dir = fullfile(script_dir, '..', '..', '..', ...
                       'textbook', 'figures', 'ch19', 'matlab');
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end

    Ns = [8 12 16 20 24 32 40 48 64];
    err_alg_1 = zeros(size(Ns));  err_alg_2 = zeros(size(Ns));
    err_tanh_1 = zeros(size(Ns)); err_tanh_2 = zeros(size(Ns));
    for i = 1:length(Ns)
        N = Ns(i);
        [err_alg_1(i),  err_alg_2(i)]  = convergence_alg(N, 2.0);
        [err_tanh_1(i), err_tanh_2(i)] = convergence_tanh(N);
    end

    NAVY=[20 45 110]/255; CORAL=[231 76 60]/255; TEAL=[22 160 133]/255;
    fig = figure('Position',[100 100 900 340],'Color','w');

    % grid panel
    subplot(1,2,1);
    [~, x_demo] = cheb_matrix(24);
    y_alg = 2.0*(1+x_demo)./(1-x_demo);
    y_alg(y_alg > 12) = 12;
    y_tanh = tanh(x_demo);
    hold on;
    plot(x_demo, 0*x_demo, '|', 'Color', NAVY, 'MarkerSize',14);
    plot(y_alg, -0.55*ones(size(y_alg)), '|', 'Color', CORAL, 'MarkerSize',14);
    plot(y_tanh, -1.1*ones(size(y_tanh)), '|', 'Color', TEAL, 'MarkerSize',14);
    xlim([-1.2 12.2]); ylim([-1.5 0.6]); yticks([]);
    xlabel('coordinate value'); title('(a) Grid clustering, N=24');

    % error panel
    subplot(1,2,2);
    semilogy(Ns, err_alg_1 + 1e-18, '-o', 'Color', CORAL); hold on;
    semilogy(Ns, err_alg_2 + 1e-18, '--o', 'Color', CORAL);
    semilogy(Ns, err_tanh_1 + 1e-18, '-s', 'Color', TEAL);
    semilogy(Ns, err_tanh_2 + 1e-18, '--s', 'Color', TEAL);
    xlabel('N'); ylabel('mapped-deriv error');
    title('(b) Validation'); grid on; box on;
    legend({'alg D1','alg D2','tanh D1','tanh D2'});

    print(fig, fullfile(out_dir, 'map1d_toolkit.pdf'), '-dpdf');
    print(fig, fullfile(out_dir, 'map1d_toolkit.png'), '-dpng');
    fprintf('[19.3-matlab] figure saved\n');
end

function [e1, e2] = convergence_alg(N, L)
    [Dx, x] = cheb_matrix(N);
    fp  = 2*L ./ (1 - x).^2;
    fpp = 4*L ./ (1 - x).^3;
    Dy  = diag(1./fp) * Dx;
    Dy2 = diag(1./fp.^2) * (Dx*Dx) - diag(fpp./fp.^3) * Dx;
    y = L*(1+x)./(1-x);
    u = exp(-y);
    mask = (abs(y) < 1e6) & isfinite(y);
    du_ex = -exp(-y); d2u_ex = exp(-y);
    du  = Dy*u;   d2u = Dy2*u;
    e1 = max(abs(du(mask)  - du_ex(mask)));
    e2 = max(abs(d2u(mask) - d2u_ex(mask)));
end

function [e1, e2] = convergence_tanh(N)
    [Dx, x] = cheb_matrix(N);
    fp  = 1 ./ cosh(x).^2;
    fpp = -2*tanh(x) ./ cosh(x).^2;
    Dy  = diag(1./fp) * Dx;
    Dy2 = diag(1./fp.^2)*(Dx*Dx) - diag(fpp./fp.^3)*Dx;
    y = tanh(x);
    u = 1./(1 + y.^2);
    du_ex = -2*y./(1 + y.^2).^2;
    d2u_ex = (6*y.^2 - 2) ./ (1 + y.^2).^3;
    du  = Dy*u;   d2u = Dy2*u;
    e1 = max(abs(du - du_ex));
    e2 = max(abs(d2u - d2u_ex));
end
