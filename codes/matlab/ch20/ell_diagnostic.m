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

    N = 64;  L_full = [0.5 1.0 2.0 4.0 8.0 16.0];
    cs_full = {CORAL, ORANGE, TEAL, PURPLE, NAVY, [139 69 19]/255};
    fig = figure('Position',[100 100 1100 800],'Color','w');

    % ---- (a) Three regimes -----------------------------------------
    subplot(2,2,1); hold on; set(gca, 'YScale','log');
    three_ell = [0.5 2.0 16.0]; three_cs = {CORAL, NAVY, TEAL};
    three_lab = {'small (early flatten)', 'good (clean descent)', ...
                 'large (gentle small-n slope)'};
    for k = 1:3
        a = abs(tbn_coeffs(N, three_ell(k)));
        env = envelope_max(a, 3);
        ns = 0:N;
        plot(ns, env + 1e-18, '-', 'Color', three_cs{k}, 'LineWidth', 1.2, ...
             'DisplayName', sprintf('ell=%g - %s', three_ell(k), three_lab{k}));
        plot(ns(1:4:end), env(1:4:end) + 1e-18, 'o', 'Color', three_cs{k}, ...
             'MarkerFaceColor', 'w', 'MarkerSize', 4, 'HandleVisibility', 'off');
    end
    ylim([1e-17 10]); grid on; box on; xlabel('degree n'); ylabel('envelope of |a_n|');
    title('(a) three regimes of ell (envelope view), N=64');
    legend('Location', 'southwest', 'FontSize', 8);

    % ---- (b) Full sweep, envelope only -----------------------------
    subplot(2,2,2); hold on; set(gca, 'YScale','log');
    for i = 1:length(L_full)
        a = abs(tbn_coeffs(N, L_full(i)));
        env = envelope_max(a, 3);
        plot(0:N, env + 1e-18, '-', 'Color', cs_full{i}, 'LineWidth', 0.9, ...
             'DisplayName', sprintf('ell=%g', L_full(i)));
    end
    ylim([1e-17 10]); grid on; box on; xlabel('degree n'); ylabel('envelope of |a_n|');
    title('(b) full ell sweep at N=64, envelope only');
    legend('Location', 'southwest', 'FontSize', 8, 'NumColumns', 2);

    % ---- (c) Tail sum vs ell at three resolutions ------------------
    subplot(2,2,3);
    Ls = [0.3 0.5 0.7 1 1.5 2 3 4 6 8 12 16 24];
    hold on; set(gca, 'YScale','log', 'XScale','log');
    N_refs = [24 48 96]; ref_cols = {CORAL, TEAL, NAVY};
    for j = 1:3
        errs = zeros(size(Ls));
        for i = 1:length(Ls)
            a = tbn_coeffs(N_refs(j), Ls(i));
            errs(i) = sum(abs(a(floor(N_refs(j)/2)+1:end)));
        end
        plot(Ls, errs, '-o', 'Color', ref_cols{j}, 'LineWidth', 1.0, ...
             'DisplayName', sprintf('N=%d', N_refs(j)));
    end
    grid on; box on; xlabel('ell'); ylabel('sum_{n>N/2} |a_n| (tail size)');
    title('(c) Valley of good ell broadens with N');
    legend('Location', 'best');

    % ---- (d) Valley centre and width vs N --------------------------
    subplot(2,2,4);
    Ns_extra = [24 32 48 64 96 128];
    centres = zeros(size(Ns_extra));
    widths  = zeros(size(Ns_extra));
    for k = 1:length(Ns_extra)
        N_ref = Ns_extra(k);
        errs_N = zeros(size(Ls));
        for i = 1:length(Ls)
            a = tbn_coeffs(N_ref, Ls(i));
            errs_N(i) = sum(abs(a(floor(N_ref/2)+1:end)));
        end
        emin = min(errs_N);
        good = errs_N <= 3.0 * emin;
        centres(k) = exp(mean(log(Ls(good))));
        widths(k)  = max(Ls(good)) / min(Ls(good));
    end
    hold on; set(gca, 'XScale','log', 'YScale','log');
    plot(Ns_extra, centres, '-o', 'Color', NAVY, 'MarkerFaceColor', 'w', ...
         'MarkerSize', 6, 'LineWidth', 1.1, 'DisplayName', 'valley centre');
    plot(Ns_extra, widths, '-s', 'Color', TEAL, 'MarkerFaceColor', 'w', ...
         'MarkerSize', 5, 'LineWidth', 1.0, 'DisplayName', 'valley width factor');
    grid on; box on; xlabel('N'); ylabel('value (log-log)');
    title('(d) broad-valley centre and width vs N');
    legend('Location', 'best');

    set(fig, 'PaperPositionMode', 'auto');
    print(fig, fullfile(out_dir, 'ell_diagnostic.pdf'), '-dpdf');
    print(fig, fullfile(out_dir, 'ell_diagnostic.png'), '-dpng');
    fprintf('[20.9-matlab] figure saved\n');
end

function out = envelope_max(a, win)
    n = numel(a); out = zeros(size(a));
    for i = 1:n
        lo = max(1, i - win); hi = min(n, i + win);
        out(i) = max(a(lo:hi));
    end
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
