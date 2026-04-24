function rK1_composed_map()
%% rK1_composed_map - Chapter 20, Etude 20.10.  Global expansion of r K_1(r).
% Composed map  r = arcsinh(exp y)  + TB_n basis on (-infty, +infty).
% Author: Dr. Denys Dutykh.

    script_dir = fileparts(mfilename('fullpath'));
    addpath(fullfile(script_dir, '..', 'ch07'));
    out_dir = fullfile(script_dir, '..', '..', '..', ...
                       'textbook', 'figures', 'ch20', 'matlab');
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end

    Ns = [8 12 16 20 24 32 48 64];
    L = 4;
    r_fine = logspace(-3, 1.6, 400);
    truth = r_fine .* besselk(1, r_fine);
    errs = zeros(size(Ns));
    for i = 1:length(Ns)
        a = tbn_approx(Ns(i), L);
        approx = evaluate(a, r_fine, L);
        errs(i) = max(abs(approx - truth));
    end

    NAVY=[20 45 110]/255; CORAL=[231 76 60]/255; TEAL=[22 160 133]/255;
    fig = figure('Position',[100 100 1020 340],'Color','w');

    subplot(1,2,1);
    semilogx(r_fine, truth, 'Color', NAVY, 'LineWidth',1.2); hold on;
    a16 = tbn_approx(16, L);
    semilogx(r_fine, evaluate(a16, r_fine, L), '--', 'Color', CORAL);
    xlabel('r'); ylabel('f(r)'); grid on; box on;
    title('(a) r K_1(r) (N=16 approx)');

    subplot(1,2,2);
    semilogy(Ns, errs + 1e-18, '-o', 'Color', TEAL);
    xlabel('N'); ylabel('max error'); grid on; box on;
    title(sprintf('(b) Subgeometric descent, L=%d', L));

    print(fig, fullfile(out_dir, 'rK1_composed_map.pdf'), '-dpdf');
    print(fig, fullfile(out_dir, 'rK1_composed_map.png'), '-dpng');
    fprintf('[20.10-matlab] figure saved\n');
end

function r = r_of_y(y)
    r = asinh(exp(y));
end

function y = y_of_r(r)
    y = log(sinh(r));
end

function a = tbn_approx(N, L)
    [~, x] = cheb_matrix(N);
    interior = abs(x) < 1 - 1e-12;
    y = zeros(size(x));
    y(interior) = L * x(interior) ./ sqrt(1 - x(interior).^2);
    r = r_of_y(y);
    fv = zeros(size(x));
    fv(interior) = r(interior) .* besselk(1, r(interior));
    % reorder: cheb_matrix returns x descending, so x(1) = 1 (y = +infty), x(end) = -1 (y = -infty)
    fv(end) = 1;         % r -> 0, r K_1(r) -> 1
    fv(1) = 0;           % r -> +infty, r K_1(r) -> 0
    fv = fv(:);
    V = [fv; fv(N:-1:2)];
    A = real(fft(V))/N; A(1) = A(1)/2; A(N+1) = A(N+1)/2;
    a = A(1:N+1);
end

function v = evaluate(a, r, L)
    y = y_of_r(r);
    x = y ./ sqrt(L^2 + y.^2);
    T0 = ones(size(r)); T1 = x;
    v = a(1)*T0 + a(2)*T1;
    for n = 2:length(a)-1
        Tk = 2*x.*T1 - T0;
        v = v + a(n+1)*Tk;
        T0 = T1; T1 = Tk;
    end
end
