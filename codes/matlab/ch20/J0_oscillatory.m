function J0_oscillatory()
%% J0_oscillatory - Chapter 20, Etude 20.11.  Teach the basis the oscillation first.
% Asymptotic augmentation of the TL_n basis to capture J_0's non-decaying carrier.
% Author: Dr. Denys Dutykh.

    script_dir = fileparts(mfilename('fullpath'));
    out_dir = fullfile(script_dir, '..', '..', '..', ...
                       'textbook', 'figures', 'ch20', 'matlab');
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end

    L = 4;
    y_fine = linspace(0.01, 50, 8001);
    truth = sqrt(1 + y_fine) .* besselj(0, y_fine);
    y_samp = linspace(0.01, 80, 2001)';
    f_samp = sqrt(1 + y_samp) .* besselj(0, y_samp);

    Ns = [4 6 8 10 15 20 30 40];
    err_n = zeros(size(Ns)); err_a = zeros(size(Ns));
    for i = 1:length(Ns)
        c = naive_fit(y_samp, f_samp, Ns(i), L);
        err_n(i) = max(abs(naive_eval(c, y_fine, L) - truth));
        [a, b] = aug_fit(y_samp, f_samp, Ns(i), L);
        err_a(i) = max(abs(aug_eval(a, b, y_fine, L) - truth));
    end

    NAVY=[20 45 110]/255; CORAL=[231 76 60]/255; TEAL=[22 160 133]/255;
    fig = figure('Position',[100 100 1020 340],'Color','w');

    subplot(1,2,1);
    plot(y_fine, truth, 'Color', NAVY, 'LineWidth',1.0); hold on;
    [a15, b15] = aug_fit(y_samp, f_samp, 15, L);
    M = build_TL(y_fine', 15, L);
    plot(y_fine, abs(M*a15), '--', 'Color', CORAL);
    plot(y_fine, abs(M*b15), ':',  'Color', TEAL);
    xlim([0 20]); grid on; box on;
    xlabel('y'); title('(a) sqrt(1+y) J_0(y) with amp-phase');
    legend({'target','|a(y)|','|\phi(y)|'});

    subplot(1,2,2);
    semilogy(Ns, err_n + 1e-18, '-o', 'Color', CORAL); hold on;
    semilogy(Ns, err_a + 1e-18, '-s', 'Color', TEAL);
    xlabel('N'); ylabel('max error'); grid on; box on;
    title('(b) Asymptotic augmentation wins');
    legend({'naive TL_n','augmented'});

    print(fig, fullfile(out_dir, 'J0_oscillatory.pdf'), '-dpdf');
    print(fig, fullfile(out_dir, 'J0_oscillatory.png'), '-dpng');
    fprintf('[20.11-matlab] figure saved\n');
end

function M = build_TL(y, N, L)
    x = (y - L) ./ (y + L);
    M = zeros(length(y), N + 1);
    M(:, 1) = 1;
    if N >= 1
        M(:, 2) = x;
    end
    for k = 3:N+1
        M(:, k) = 2*x.*M(:, k-1) - M(:, k-2);
    end
end

function c = naive_fit(y, f, N, L)
    M = build_TL(y, N, L);
    c = M \ f;
end

function [a, b] = aug_fit(y, f, N, L)
    M = build_TL(y, N, L);
    C = cos(y - pi/4); S = sin(y - pi/4);
    D = [M .* C, M .* S];
    ab = D \ f;
    a = ab(1:N+1); b = ab(N+2:end);
end

function v = naive_eval(c, y, L)
    M = build_TL(y', length(c)-1, L);
    v = (M * c)';
end

function v = aug_eval(a, b, y, L)
    M = build_TL(y', length(a)-1, L);
    v = ((M*a) .* cos(y' - pi/4) + (M*b) .* sin(y' - pi/4))';
end
