function ell_diagnostic()
%% ell_diagnostic - Chapter 20, Etude 20.9.  Read the coefficients before the error.
% Boyd's Rule-of-Thumb 16: a too-small ell causes |a_n| to flatten early.
% Author: Dr. Denys Dutykh.

    script_dir = fileparts(mfilename('fullpath'));
    addpath(fullfile(script_dir, '..', 'ch07'));
    out_dir = fullfile(script_dir, '..', '..', '..', ...
                       'textbook', 'figures', 'ch20', 'matlab');
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end

    NAVY=[20 45 110]/255; CORAL=[231 76 60]/255; TEAL=[22 160 133]/255;
    ORANGE=[230 126 34]/255; PURPLE=[142 68 173]/255;

    N = 64;  L_list = [0.5 1.0 2.0 4.0 8.0 16.0];
    cs = {CORAL, ORANGE, TEAL, PURPLE, NAVY, [139 69 19]/255};
    fig = figure('Position',[100 100 1160 400],'Color','w');

    subplot(1,2,1); hold on; set(gca, 'YScale','log');
    for i = 1:length(L_list)
        a = abs(tbn_coeffs(N, L_list(i)));
        semilogy(0:N, a + 1e-18, '-o', 'MarkerSize',3, 'Color', cs{i}, ...
                 'DisplayName', sprintf('ell=%g', L_list(i)));
    end
    ylim([1e-17 10]); grid on; box on; xlabel('n'); ylabel('|a_n|');
    title('(a) Coefficient decay, N=64'); legend('Location','southwest');

    subplot(1,2,2);
    Ls = [0.3 0.5 0.7 1 1.5 2 3 4 6 8 12 16 24];
    hold on; set(gca, 'YScale','log', 'XScale','log');
    N_refs = [24 48 96]; ref_cols = {CORAL, TEAL, NAVY};
    for j = 1:3
        errs = zeros(size(Ls));
        for i = 1:length(Ls)
            a = tbn_coeffs(N_refs(j), Ls(i));
            errs(i) = sum(abs(a(floor(N_refs(j)/2)+1:end)));
        end
        loglog(Ls, errs, '-o', 'Color', ref_cols{j}, 'DisplayName', sprintf('N=%d', N_refs(j)));
    end
    grid on; box on; xlabel('ell'); ylabel('tail sum');
    title('(b) Valley of good ell broadens with N');
    legend('Location','best');

    print(fig, fullfile(out_dir, 'ell_diagnostic.pdf'), '-dpdf');
    print(fig, fullfile(out_dir, 'ell_diagnostic.png'), '-dpng');
    fprintf('[20.9-matlab] figure saved\n');
end

function a = tbn_coeffs(N, ell)
    [~, x] = cheb_matrix(N);
    y = ell * x ./ sqrt(1 - x.^2);
    fv = zeros(size(x));
    ok = abs(x) < 1 - 1e-12;
    fv(ok) = 1 ./ cosh(y(ok));
    fv = fv(:);
    V = [fv; fv(N:-1:2)];
    A = real(fft(V))/N; A(1) = A(1)/2; A(N+1) = A(N+1)/2;
    a = A(1:N+1);
end
