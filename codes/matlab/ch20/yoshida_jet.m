function yoshida_jet()
%% yoshida_jet - Chapter 20, Etude 20.8.  Pick the right half of the basis.
% Author: Dr. Denys Dutykh.

    script_dir = fileparts(mfilename('fullpath'));
    out_dir = fullfile(script_dir, '..', '..', '..', ...
                       'textbook', 'figures', 'ch20', 'matlab');
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end

    ell = 3;
    Ns = [2 3 4 5 6 8 10];
    y_fine = linspace(0.1, 20, 4001);
    [c_ref, idx_ref] = assemble(21, ell);
    v_ref = evaluate(c_ref, idx_ref, y_fine, ell);
    errs = zeros(size(Ns));
    for i = 1:length(Ns)
        [c, idx] = assemble(Ns(i), ell);
        v = evaluate(c, idx, y_fine, ell);
        errs(i) = max(abs(v - v_ref));
    end

    NAVY=[20 45 110]/255; CORAL=[231 76 60]/255; TEAL=[22 160 133]/255;
    fig = figure('Position',[100 100 1020 340],'Color','w');
    subplot(1,2,1);
    y_p = linspace(0.05, 8, 401);
    plot(y_p, evaluate(c_ref, idx_ref, y_p, ell), 'Color', NAVY, 'LineWidth',1.4); hold on;
    styles = {CORAL, '--', 2; TEAL, '-.', 3; [142 68 173]/255, ':', 4};
    for ii = 1:size(styles, 1)
        [c, idx] = assemble(styles{ii, 3}, ell);
        plot(y_p, evaluate(c, idx, y_p, ell), 'LineStyle', styles{ii, 2}, ...
             'Color', styles{ii, 1}, 'LineWidth', 1.0);
    end
    grid on; box on; xlabel('y'); ylabel('v(y)');
    title('(a) Yoshida jet'); legend({'ref (21 SB)','N=2','N=3','N=4'});

    subplot(1,2,2);
    semilogy(Ns, errs + 1e-18, '-o', 'Color', TEAL);
    xlabel('N'); ylabel('max error'); grid on; box on;
    title('(b) Odd SB basis resolves 1/y tail');

    print(fig, fullfile(out_dir, 'yoshida_jet.pdf'), '-dpdf');
    print(fig, fullfile(out_dir, 'yoshida_jet.png'), '-dpng');
    fprintf('[20.8-matlab] figure saved\n');
end

function [c, idx] = assemble(N, ell)
    idx = 2*(0:N-1) + 1;
    t = (1:N) * pi / (2*(N+1));
    y = ell ./ tan(t);
    A = zeros(N, N); rhs = y(:);
    for j = 1:N
        A(:, j) = sb_ddot(idx(j), y, ell) - (y.^2) .* sb_eval(idx(j), y, ell);
    end
    c = A \ rhs;
end

function v = evaluate(c, idx, y, ell)
    v = zeros(size(y));
    for k = 1:length(c)
        v = v + c(k) * sb_eval(idx(k), y, ell);
    end
end

function s = sb_eval(n, y, ell)
    t = pi/2 - atan(y/ell);
    s = sin((n+1) * t);
end

function s = sb_ddot(n, y, ell)
    % second derivative of SB_n(y; ell) w.r.t. y.
    t = pi/2 - atan(y/ell);
    denom = ell^2 + y.^2;
    dt = -ell ./ denom;
    dt2 = 2*ell*y ./ denom.^2;
    c = cos((n+1)*t); sv = sin((n+1)*t);
    s = -((n+1)^2) * sv .* dt.^2 + (n+1) * c .* dt2;
end
