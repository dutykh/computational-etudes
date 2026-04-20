function chebyshev_as_cosine()
%% chebyshev_as_cosine - Chapter 19, Computational Etude 19.2.
%
% Solve  u'' - q u = f  on [-1, 1] with Dirichlet BCs  in two ways:
%   (X) Chebyshev collocation in x using the dense differentiation matrix;
%   (T) Chebyshev collocation via the FFT-based cosine differentiation,
%       which is the "working in t" viewpoint of Boyd Sec 16.2.
%
% Exact solution: u_ex(x) = sin(pi x).
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)

    script_dir = fileparts(mfilename('fullpath'));
    addpath(fullfile(script_dir, '..', 'ch07'));
    out_dir = fullfile(script_dir, '..', '..', '..', ...
                       'textbook', 'figures', 'ch19', 'matlab');
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end

    q = 4;
    Ns = [8 12 16 20 24 32 40 48];
    err_x = zeros(size(Ns)); err_t = zeros(size(Ns));
    for i = 1:length(Ns)
        N = Ns(i);
        [x1, u1] = solve_x(N, q);
        err_x(i) = max(abs(u1 - u_ex(x1)));
        [x2, u2] = solve_t(N, q);
        err_t(i) = max(abs(u2 - u_ex(x2)));
    end

    NAVY=[20 45 110]/255; CORAL=[231 76 60]/255; TEAL=[22 160 133]/255;
    fig = figure('Position',[100 100 820 340],'Color','w');

    subplot(1,2,1);
    [xN, uN] = solve_x(24, q);
    xfine = linspace(-1, 1, 401);
    plot(xfine, u_ex(xfine), '-', 'Color', NAVY, 'LineWidth',1.2); hold on;
    plot(xN, uN, 'o', 'Color', CORAL, 'MarkerSize',4);
    xlabel('x'); ylabel('u(x)');
    title('(a) Manufactured BVP'); grid on; box on;
    legend({'exact','X-form N=24'});

    subplot(1,2,2);
    semilogy(Ns, err_x + 1e-18, '-o', 'Color', CORAL, 'LineWidth',1.1); hold on;
    semilogy(Ns, err_t + 1e-18, '-s', 'Color', TEAL,  'LineWidth',1.1);
    xlabel('N'); ylabel('max error');
    title('(b) Same convergence, different arithmetic'); grid on; box on;
    legend({'X-form (DM)','T-form (FFT)'});

    print(fig, fullfile(out_dir, 'chebyshev_as_cosine.pdf'), '-dpdf');
    print(fig, fullfile(out_dir, 'chebyshev_as_cosine.png'), '-dpng');
    fprintf('[19.2-matlab] figure saved\n');
end

function v = u_ex(x), v = sin(pi*x); end
function v = f_src(x), v = -pi^2*sin(pi*x) - 4*sin(pi*x); end

function [x, u] = solve_x(N, q)
    [D, x] = cheb_matrix(N);
    D2 = D*D;
    A = D2(2:N, 2:N) - q*eye(N-1);
    rhs = f_src(x(2:N));
    u_int = A \ rhs;
    u = zeros(N+1,1); u(2:N) = u_int;
end

function [x, u] = solve_t(N, q)
    x = cos(pi*(0:N)/N)';
    % build D^2 via FFT-cosine derivative applied column-wise
    D2 = zeros(N+1, N+1);
    e = zeros(N+1, 1);
    for j = 1:N+1
        e(:) = 0; e(j) = 1;
        D2(:, j) = chebfft(chebfft(e));
    end
    A = D2(2:N, 2:N) - q*eye(N-1);
    rhs = f_src(x(2:N));
    u_int = A \ rhs;
    u = zeros(N+1,1); u(2:N) = u_int;
end

function w = chebfft(v)
    N = length(v) - 1;
    x = cos(pi*(0:N)/N)';
    V = [v; v(N:-1:2)];
    U = real(fft(V));
    k = (0:N)';
    w_hat = 1i*[k(1:N); 0; k(2:N)-N] .* U;
    W = real(ifft(w_hat));
    w = zeros(N+1, 1);
    w(2:N) = -W(2:N) ./ sqrt(1 - x(2:N).^2);
    k2 = (0:N)'.^2;
    w(1) = sum(k2 .* U(1:N+1)) / N + 0.5*N*U(N+1);
    w(N+1) = sum(((-1).^((0:N)'+1)) .* k2 .* U(1:N+1))/N ...
             + 0.5*(-1)^(N+1)*N*U(N+1);
end
