function truncation_stalls()
%% truncation_stalls - Chapter 20, Computational Etude 20.1.
%
% Domain truncation of sech(y) on (-infty, +infty): three sweeps
% demonstrating Boyd's Rule of Thumb 14 (unbounded intervals demand TWO
% parameters, N and L).
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)

    script_dir = fileparts(mfilename('fullpath'));
    addpath(fullfile(script_dir, '..', 'ch07'));
    out_dir = fullfile(script_dir, '..', '..', '..', ...
                       'textbook', 'figures', 'ch20', 'matlab');
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end

    % (A) fix L, vary N
    L_A = 6; Ns_A = [8 12 16 24 32 48 64 96 128];
    err_A = arrayfun(@(N) cheb_trunc_err(N, L_A), Ns_A);
    % (B) fix N, vary L
    N_B = 32; Ls_B = [2 3 4 5 6 8 10 12 16 20];
    err_B = arrayfun(@(L) cheb_trunc_err(N_B, L), Ls_B);
    % (C) grow both
    pairs = [8 3; 12 4; 16 5; 24 6; 32 8; 48 10; 64 12; 96 14];
    err_C = zeros(size(pairs, 1), 1);
    for i = 1:size(pairs, 1)
        err_C(i) = cheb_trunc_err(pairs(i,1), pairs(i,2));
    end
    Ns_C = pairs(:,1);

    NAVY=[20 45 110]/255; CORAL=[231 76 60]/255; TEAL=[22 160 133]/255;

    % NEW (panel d): full 2D error surface E(N, L)
    Ns_grid = [8 12 16 24 32 48 64 96 128];
    Ls_grid = [2 3 4 5 6 8 10 12 16 20];
    err_grid = zeros(numel(Ls_grid), numel(Ns_grid));
    for i = 1:numel(Ls_grid)
        for j = 1:numel(Ns_grid)
            err_grid(i, j) = cheb_trunc_err(Ns_grid(j), Ls_grid(i));
        end
    end
    log_err = log10(err_grid + 1e-18);

    fig = figure('Position',[100 100 1100 800],'Color','w');

    subplot(2,2,1);
    semilogy(Ns_A, err_A, '-o', 'Color', CORAL, 'LineWidth',1.1, ...
             'DisplayName', sprintf('error at L=%d', L_A)); hold on;
    yline(exp(-L_A), '--', 'Color', NAVY, 'DisplayName', 'e^{-L}');
    xlabel('N'); ylabel('max error'); grid on; box on;
    title(sprintf('(a) Fix L=%d, vary N: plateau', L_A));
    legend('Location','best');

    subplot(2,2,2);
    semilogy(Ls_B, err_B, '-s', 'Color', TEAL, 'LineWidth',1.1, ...
             'MarkerFaceColor', 'w');
    xlabel('L'); ylabel('max error'); grid on; box on;
    title(sprintf('(b) Fix N=%d, vary L: sweet spot', N_B));

    subplot(2,2,3);
    semilogy(Ns_C, err_C, '-^', 'Color', NAVY, 'LineWidth',1.1);
    xlabel('N (L growing with N)'); ylabel('max error');
    grid on; box on;
    title('(c) Grow both: subgeometric descent');

    subplot(2,2,4);
    imagesc(Ns_grid, Ls_grid, log_err);
    set(gca, 'YDir', 'normal');
    colormap(gca, parula);
    cb = colorbar; cb.Label.String = 'log_{10} max error';
    hold on;
    plot(Ns_A, L_A * ones(size(Ns_A)), '-o', 'Color', CORAL, ...
         'MarkerFaceColor', 'w', 'MarkerSize', 4, 'LineWidth', 1.4, ...
         'DisplayName', '(a) L = 6, vary N');
    plot(N_B * ones(size(Ls_B)), Ls_B, '-s', 'Color', TEAL, ...
         'MarkerFaceColor', 'w', 'MarkerSize', 4, 'LineWidth', 1.4, ...
         'DisplayName', sprintf('(b) N = %d, vary L', N_B));
    Ls_C = pairs(:,2);
    plot(Ns_C, Ls_C, '-^', 'Color', 'w', 'MarkerFaceColor', 'w', ...
         'MarkerSize', 5, 'LineWidth', 1.6, ...
         'DisplayName', '(c) grow both jointly');
    hold off;
    xlabel('N'); ylabel('L');
    title('(d) error surface E(N, L) and its three slices');
    legend('Location', 'northeast', 'TextColor', 'w', 'Color', 'k', ...
           'EdgeColor', 'none', 'FontSize', 8);

    set(fig, 'PaperPositionMode', 'auto');
    print(fig, fullfile(out_dir, 'truncation_stalls.pdf'), '-dpdf');
    print(fig, fullfile(out_dir, 'truncation_stalls.png'), '-dpng');
    fprintf('[20.1-matlab] figure saved\n');
end

function e = cheb_trunc_err(N, L)
    [~, x] = cheb_matrix(N);
    y_nodes = L * x;
    fv = 1 ./ cosh(y_nodes);
    a = dct1_coeffs(fv);
    y_fine = linspace(-20, 20, 4001);
    in_window = abs(y_fine) <= L;
    xin = y_fine(in_window) / L;
    approx = zeros(size(y_fine));
    approx(in_window) = cheb_eval(a, xin, N);
    truth = 1 ./ cosh(y_fine);
    e = max(abs(approx - truth));
end

function a = dct1_coeffs(v)
    v = v(:);
    N = length(v) - 1;
    V = [v; v(N:-1:2)];
    A = real(fft(V)) / N;
    A(1) = A(1)/2; A(N+1) = A(N+1)/2;
    a = A(1:N+1);
end

function v = cheb_eval(a, xfine, N)
    T0 = ones(size(xfine)); T1 = xfine;
    v = a(1)*T0 + a(2)*T1;
    for n = 2:N
        Tk = 2*xfine.*T1 - T0;
        v = v + a(n+1)*Tk;
        T0 = T1; T1 = Tk;
    end
end
