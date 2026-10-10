function behavioral_bc_laguerre()
%% behavioral_bc_laguerre - Chapter 20, Etude 20.7.
% Boundedness is sometimes enough, sometimes not.
% Author: Dr. Denys Dutykh.

    script_dir = fileparts(mfilename('fullpath'));
    addpath(fullfile(script_dir, '..', 'ch07'));
    out_dir = fullfile(script_dir, '..', '..', '..', ...
                       'textbook', 'figures', 'ch20', 'matlab');
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end

    ell = 32;
    Ns = [10 20 30 40 60 80 120];
    good_A = zeros(size(Ns)); good_B = zeros(size(Ns));
    for i = 1:length(Ns)
        eA = sort(real(solve_A(Ns(i), ell)));
        eB = sort(real(solve_B(Ns(i), ell)));
        good_A(i) = count_good(eA);
        good_B(i) = count_good(eB);
    end

    NAVY=[20 45 110]/255; CORAL=[231 76 60]/255; TEAL=[22 160 133]/255;
    fig = figure('Position',[100 100 1020 340],'Color','w');
    subplot(1,2,1);
    plot(Ns, good_A, '-o', 'Color', CORAL); hold on;
    plot(Ns, good_B, '-s', 'Color', TEAL);
    xlabel('N'); ylabel('good eigenvalues'); grid on; box on;
    title('(a) Good-eigenvalue count'); legend({'naive','behavioural'});

    subplot(1,2,2);
    eA = solve_A(40, ell); eB = solve_B(40, ell);
    plot(real(eA(1:20)), imag(eA(1:20)), 'o', 'Color', CORAL); hold on;
    plot(real(eB(1:20)), imag(eB(1:20)), 's', 'Color', TEAL);
    yline(0, '-', 'Color', [0.5 0.5 0.5]);
    xlim([-1 20]); ylim([-3 3]); grid on; box on;
    xlabel('Re(\lambda)'); ylabel('Im(\lambda)');
    title('(b) Spectrum (N=40)');
    legend({'A naive','B behavioural'});

    print(fig, fullfile(out_dir, 'behavioral_bc_laguerre.pdf'), '-dpdf');
    print(fig, fullfile(out_dir, 'behavioral_bc_laguerre.png'), '-dpng');
    fprintf('[20.7-matlab] figure saved\n');
end

function [y, Dy, Dy2] = tln_dmatrices(N, ell)
    [Dx, x] = cheb_matrix(N);
    fp = 2*ell ./ (1-x).^2;
    fpp = 4*ell ./ (1-x).^3;
    Dy = diag(1./fp) * Dx;
    Dy2 = diag(1./fp.^2)*(Dx*Dx) - diag(fpp./fp.^3)*Dx;
    y_full = ell*(1+x)./(1-x);
    y = y_full(2:end);
    Dy = Dy(2:end, 2:end);
    Dy2 = Dy2(2:end, 2:end);
end

function e = solve_A(N, ell)
    [y, Dy, Dy2] = tln_dmatrices(N, ell);
    Y = diag(y);
    A = Y*Dy2 + (Y + eye(length(y)))*Dy;
    e = eig(-A);
end

function e = solve_B(N, ell)
    [y, Dy, Dy2] = tln_dmatrices(N, ell);
    Y = diag(y); M = diag(0.5 + 0.25*y);
    A = Y*Dy2 + Dy - M;
    e = eig(-A);
end

function n = count_good(eigs)
    eigs = eigs(eigs > -0.5);
    targets = 0:49;
    taken = false(size(targets));
    n = 0;
    for i = 1:length(eigs)
        [~, order] = sort(abs(eigs(i) - targets));
        for k = 1:length(order)
            idx = order(k);
            if taken(idx), continue; end
            rel = abs(eigs(i) - targets(idx)) / max(1, targets(idx));
            if rel < 0.05
                taken(idx) = true; n = n + 1;
            end
            break;
        end
    end
end
