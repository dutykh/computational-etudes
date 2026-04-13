% fourier_cfl.m
% Chapter 17: Stability of Spectral Methods
%
% Computational Etude 17.3: Fourier CFL Scaling.
%
% Two-panel figure:
%   (a) Eigenvalues of the Fourier differentiation matrix for N=32,
%       plotted on the imaginary axis from -16i to 16i.
%   (b) Log-log plot of critical dt vs N for RK4 applied to Fourier
%       advection u_t + u_x = 0. Binary search for dt_crit.
%
% Author: Dr. Denys Dutykh
% Part of "Computational Etudes: A Spectral Approach"

clear; clc; close all;

NAVY  = [20, 45, 110]/255;
CORAL = [231, 76, 60]/255;

%% ---- Panel (a): Fourier eigenvalues for N=32 ----
N = 32;
k = [0:N/2-1, -N/2:-1];   % fftfreq equivalent
k_sorted = sort(k);
eigs_val = 1i * k_sorted;  % eigenvalues of d/dx

figure('Position', [100, 100, 1100, 450]);

subplot(1, 2, 1);
plot(real(eigs_val), imag(eigs_val), 'o', 'Color', NAVY, 'MarkerSize', 6, ...
     'MarkerFaceColor', NAVY, 'MarkerEdgeColor', 'w');
hold on;
plot([-2, 2], [0, 0], 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
plot([0, 0], [-20, 20], 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);

% Annotate k_max
text(0.5, 14, 'ik_{max} = i \cdot 16', 'Color', CORAL, 'FontSize', 10);
text(0.5, -14, '-ik_{max}', 'Color', CORAL, 'FontSize', 10);
hold off;

xlabel('Re(\lambda)'); ylabel('Im(\lambda)');
title(sprintf('(a) Fourier eigenvalues, N = %d', N));
xlim([-2, 2]); ylim([-20, 20]);
grid on; set(gca, 'GridAlpha', 0.3);

%% ---- Panel (b): dt_crit vs N ----
N_values = [16, 32, 64, 128, 256];
dt_crits = zeros(size(N_values));

for idx = 1:length(N_values)
    Nv = N_values(idx);
    dtc = find_dt_crit_rk4(Nv);
    dt_crits(idx) = dtc;
    fprintf('N = %4d, dt_crit = %.6f, dt_crit * (N/2) = %.4f\n', ...
            Nv, dtc, dtc * Nv / 2);
end

subplot(1, 2, 2);
loglog(N_values, dt_crits, 'o-', 'Color', NAVY, 'MarkerSize', 7, ...
       'LineWidth', 1.5, 'MarkerFaceColor', NAVY, ...
       'DisplayName', '\Delta t_{crit} (binary search)');
hold on;

% Reference line: 2.83 / (N/2)
N_ref = linspace(12, 300, 100);
loglog(N_ref, 2.83 ./ (N_ref/2), '--', 'Color', CORAL, 'LineWidth', 1.2, ...
       'DisplayName', '2.83 / (N/2)');

% Slope annotation
text(80, 0.06, 'slope = -1', 'FontSize', 10, 'Color', CORAL, ...
     'Rotation', -28, 'HorizontalAlignment', 'center');
hold off;

xlabel('N'); ylabel('\Delta t_{crit}');
title('(b) RK4 CFL scaling for Fourier advection');
legend('Location', 'northeast', 'FontSize', 9);
grid on; set(gca, 'GridAlpha', 0.3);

set(gcf, 'PaperPositionMode', 'auto');
print('-dpdf', 'fourier_cfl_scaling.pdf');
fprintf('\nFigure saved to fourier_cfl_scaling.pdf\n');


%% ---- Helper functions ----

function g = rk4_amplification(z)
    % RK4 amplification factor
    g = 1 + z + z.^2/2 + z.^3/6 + z.^4/24;
end

function gmax = rk4_max_amp(dt, N)
    % Maximum |g(dt * lambda)| over all Fourier eigenvalues
    k = [0:N/2-1, -N/2:-1];
    z = -1i * k * dt;   % dt * lambda_k for advection
    gmax = max(abs(rk4_amplification(z)));
end

function dtc = find_dt_crit_rk4(N)
    % Binary search for critical dt where max|g| = 1
    tol = 1e-8;
    dt_lo = 1e-10;
    dt_hi = 1.0;

    % Find unstable upper bound
    while rk4_max_amp(dt_hi, N) <= 1.0
        dt_hi = dt_hi * 2;
        if dt_hi > 100, dtc = dt_hi; return; end
    end

    for iter = 1:80
        dt_mid = 0.5 * (dt_lo + dt_hi);
        if rk4_max_amp(dt_mid, N) <= 1.0
            dt_lo = dt_mid;
        else
            dt_hi = dt_mid;
        end
        if dt_hi - dt_lo < tol * dt_lo, break; end
    end
    dtc = 0.5 * (dt_lo + dt_hi);
end
