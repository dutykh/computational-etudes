function hilbert_fourier(varargin)
%% hilbert_fourier - Etude 21.8: Hilbert transform via Fourier multiplier.
% Periodic test: f(x) = exp(cos x) on [-pi, pi].  Fourier series uses the
% modified Bessel coefficients I_k(1), so convergence is super-geometric.
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
% Part of "Computational Etudes: A Spectral Approach"

    dump_path = parse_args(varargin{:});
    script_dir = fileparts(mfilename('fullpath'));
    addpath(script_dir);
    cm = tricks_common(); cm.configure_style();
    out_dir = cm.output_dir(script_dir);

    Ns = [4 6 8 10 12 14 16 20 24 32];
    err = zeros(size(Ns));
    for k = 1:numel(Ns)
        N = Ns(k);
        x = -pi + 2 * pi * (0:N-1) / N;
        Hf = hilbert_via_fft(exp(cos(x)));
        Hf_ex = hilbert_exact(x);
        err(k) = max(abs(Hf - Hf_ex));
    end

    N_show = 32;
    x_show = -pi + 2 * pi * (0:N_show-1) / N_show;
    Hf_show = hilbert_via_fft(exp(cos(x_show)));
    x_dense = linspace(-pi, pi, 401);
    f_dense = exp(cos(x_dense));
    Hf_ex_dense = hilbert_exact(x_dense);

    fig = figure('Units', 'inches', 'Position', [1, 1, 11.5, 4.4], 'Color', 'w');
    tl = tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile(tl); hold on;
    plot(x_dense, f_dense, 'Color', cm.NAVY, 'LineWidth', 1.4, ...
         'DisplayName', '$f(x) = e^{\cos x}$');
    plot(x_dense, Hf_ex_dense, 'Color', cm.TEAL, 'LineWidth', 1.4, ...
         'DisplayName', '$H\{f\}(y)$ exact');
    plot(x_show, Hf_show, 'o', 'Color', cm.CORAL, 'MarkerFaceColor', 'w', ...
         'MarkerSize', 5, 'DisplayName', sprintf('$H\\{f\\}_N(y_j)$ via FFT, $N=%d$', N_show));
    yline(0, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.4, 'Alpha', 0.5);
    hold off;
    xlabel('$x$ (or $y$)', 'Interpreter', 'latex');
    ylabel('amplitude', 'Interpreter', 'latex');
    title('$f(x) = e^{\cos x}$ and its Hilbert transform on the circle', ...
          'Interpreter', 'latex');
    legend('Location', 'southwest', 'Interpreter', 'latex', 'FontSize', 9);
    grid on; box on;
    xlim([-pi - 0.1, pi + 0.1]);

    nexttile(tl); hold on;
    semilogy(Ns, err + 1e-18, 'o-', 'Color', cm.NAVY, 'MarkerFaceColor', 'w', ...
             'MarkerSize', 6, 'LineWidth', 1.0, ...
             'DisplayName', 'max $|H\{f\}_N - H\{f\}|$');
    bound = 2.0 * besseli(Ns / 2, 1.0);
    semilogy(Ns, bound, '--', 'Color', cm.CORAL, 'LineWidth', 0.8, ...
             'DisplayName', '$\sim 2 I_{N/2}(1)$ (Bessel-tail bound)');
    yline(1e-15, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.4, 'Alpha', 0.5);
    set(gca, 'YScale', 'log');
    hold off;
    xlabel('$N$ (Fourier modes per period)', 'Interpreter', 'latex');
    ylabel('max error', 'Interpreter', 'latex');
    title('super-geometric convergence on a periodic analytic test', ...
          'Interpreter', 'latex');
    legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 9);
    grid on; box on;
    ylim([1e-17, 5]);

    exportgraphics(fig, fullfile(out_dir, 'hilbert_fourier.pdf'), 'ContentType', 'vector');
    exportgraphics(fig, fullfile(out_dir, 'hilbert_fourier.png'), 'Resolution', 300);
    close(fig);

    fprintf('[Etude 21.8]  Hilbert transform via Fourier multiplier\n');
    for k = 1:numel(Ns)
        fprintf('  N = %3d: max err = %.3e\n', Ns(k), err(k));
    end
    fprintf('  figure: %s\n', fullfile(out_dir, 'hilbert_fourier.pdf'));

    if ~isempty(dump_path)
        r = struct('Ns', Ns, 'err', err);
        fid = fopen(dump_path, 'w'); fprintf(fid, '%s', jsonencode(r)); fclose(fid);
    end
end

function dump_path = parse_args(varargin)
    dump_path = ''; k = 1;
    while k <= numel(varargin)
        if strcmp(varargin{k}, '--dump'); dump_path = char(varargin{k+1}); k = k + 2;
        else; k = k + 1; end
    end
end

function Hf = hilbert_via_fft(fv)
    N = numel(fv);
    F = fft(fv);
    k = [0:N/2-1, -N/2:-1];   % integer wavenumbers in fftfreq order
    multiplier = -1i * sign(k);
    multiplier(1) = 0;
    Hf = real(ifft(multiplier .* F));
end

function v = hilbert_exact(y)
    v = zeros(size(y));
    for k = 1:40
        v = v + 2.0 * besseli(k, 1.0) * sin(k * y);
    end
end
