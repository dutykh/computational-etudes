% chebyshev_stiffness.m
% Chapter 17: Stability of Spectral Methods
%
% Computational Etude 17.4: Chebyshev Stiffness and Eigenvalue Structure.
%
% Three figures:
%   (1) Eigenvalues of Chebyshev D^2 (Dirichlet BCs) for N=32.
%   (2) Physical eigenvector (mode N/4) and outlier eigenvector.
%   (3) Log-log plot of spectral radius rho(D^2) vs N, showing N^4 scaling.
%
% Author: Dr. Denys Dutykh
% Part of "Computational Etudes: A Spectral Approach"

clear; clc; close all;

NAVY  = [20, 45, 110]/255;
SKY   = [120, 150, 210]/255;
CORAL = [231, 76, 60]/255;
TEAL  = [22, 160, 133]/255;

%% ---- Figure 1: Eigenvalues of Chebyshev D^2 for N=32 ----
N = 32;
[D, x] = cheb(N);
D2 = D * D;
D2_int = D2(2:N, 2:N);

eigs_all = sort(real(eig(D2_int)));

% Physical eigenvalue predictions: -(n*pi/2)^2 for n = 1, 2, ...
n_phys = 1:(N-1);
phys_approx = -(n_phys * pi / 2).^2;

% Separate physical from outliers
threshold = -(N * pi / 2)^2 * 0.5;
phys_mask = eigs_all > threshold;
out_mask  = ~phys_mask;

figure('Position', [100, 100, 1000, 350]);
hold on;
plot(eigs_all(phys_mask), zeros(sum(phys_mask), 1), 'o', ...
     'Color', NAVY, 'MarkerSize', 7, 'MarkerFaceColor', NAVY, ...
     'MarkerEdgeColor', 'w', 'DisplayName', 'physical modes');
if any(out_mask)
    plot(eigs_all(out_mask), zeros(sum(out_mask), 1), 's', ...
         'Color', CORAL, 'MarkerSize', 7, 'MarkerFaceColor', CORAL, ...
         'MarkerEdgeColor', 'w', 'DisplayName', 'outlier modes');
end

% Mark physical approximations
for n = 1:min(7, N-1)
    plot(phys_approx(n), 0.05, 'v', 'Color', TEAL, 'MarkerSize', 5, ...
         'MarkerFaceColor', TEAL, 'HandleVisibility', 'off');
end
plot(NaN, NaN, 'v', 'Color', TEAL, 'MarkerSize', 5, ...
     'MarkerFaceColor', TEAL, 'DisplayName', '-n^2\pi^2/4 (predicted)');

plot([min(eigs_all)*1.1, max(eigs_all)*0.5], [0, 0], ...
     'Color', [0.5 0.5 0.5], 'LineWidth', 0.5, 'HandleVisibility', 'off');
hold off;

xlabel('\lambda (eigenvalue)');
title(sprintf('Eigenvalues of Chebyshev D^2 with Dirichlet BCs (N = %d)', N));
set(gca, 'YTick', []);
legend('Location', 'northwest', 'FontSize', 9);
grid on; set(gca, 'GridAlpha', 0.3);

set(gcf, 'PaperPositionMode', 'auto');
print('-dpdf', 'cheb_eigenvalues.pdf');
fprintf('Figure saved to cheb_eigenvalues.pdf\n');

%% ---- Figure 2: Physical and outlier eigenvectors ----
x_int = x(2:N);
[vecs, lam_diag] = eig(D2_int);
lam_vals = diag(lam_diag);
[lam_sorted, idx] = sort(real(lam_vals));
vecs_sorted = vecs(:, idx);

figure('Position', [100, 100, 1100, 450]);

% Physical eigenvector: mode N/4 = 8 (0-indexed: mode_idx = N/4 - 1)
mode_idx = N/4;   % 1-based in MATLAB, so mode 8 is index 8 from bottom
v_phys = vecs_sorted(:, mode_idx);
v_phys = v_phys / max(abs(v_phys));
if v_phys(floor(length(v_phys)/4)) < 0
    v_phys = -v_phys;
end
v_full = zeros(N+1, 1);
v_full(2:N) = v_phys;

subplot(1, 2, 1);
plot(x, v_full, '-', 'Color', NAVY, 'LineWidth', 1.8);
hold on;
plot(x_int, v_phys, 'o', 'Color', NAVY, 'MarkerSize', 4, ...
     'MarkerFaceColor', 'w', 'MarkerEdgeColor', NAVY);
hold off;
xlabel('x'); ylabel('v(x)');
title(sprintf('(a) Physical eigenvector (mode %d, \\lambda = %.1f)', ...
              N/4, lam_sorted(mode_idx)));
xlim([-1.05, 1.05]);
grid on; set(gca, 'GridAlpha', 0.3);

% Outlier eigenvector: most negative eigenvalue (index 1 after sorting)
v_out = real(vecs_sorted(:, 1));
v_out = v_out / max(abs(v_out));
v_full_out = zeros(N+1, 1);
v_full_out(2:N) = v_out;

subplot(1, 2, 2);
plot(x, v_full_out, '-', 'Color', CORAL, 'LineWidth', 1.8);
hold on;
plot(x_int, v_out, 'o', 'Color', CORAL, 'MarkerSize', 4, ...
     'MarkerFaceColor', 'w', 'MarkerEdgeColor', CORAL);
hold off;
xlabel('x'); ylabel('v(x)');
title(sprintf('(b) Outlier eigenvector (\\lambda = %.0f)', lam_sorted(1)));
xlim([-1.05, 1.05]);
grid on; set(gca, 'GridAlpha', 0.3);

set(gcf, 'PaperPositionMode', 'auto');
print('-dpdf', 'cheb_eigenvectors.pdf');
fprintf('Figure saved to cheb_eigenvectors.pdf\n');

%% ---- Figure 3: CFL scaling (spectral radius vs N) ----
N_values = [8, 12, 16, 20, 24, 32, 48, 64];
rho_values = zeros(size(N_values));

for idx = 1:length(N_values)
    Nv = N_values(idx);
    [Dv, ~] = cheb(Nv);
    D2v = Dv * Dv;
    D2v_int = D2v(2:Nv, 2:Nv);
    ev = eig(D2v_int);
    rho_values(idx) = max(abs(ev));
    fprintf('N = %3d, rho(D2) = %.2e, rho/N^4 = %.4f\n', ...
            Nv, rho_values(idx), rho_values(idx) / Nv^4);
end

figure('Position', [100, 100, 700, 500]);
loglog(N_values, rho_values, 'o-', 'Color', NAVY, 'MarkerSize', 7, ...
       'LineWidth', 1.5, 'MarkerFaceColor', NAVY, ...
       'DisplayName', '\rho(\tilde{D}_2)');
hold on;

% Reference line: 0.048 * N^4
N_ref = linspace(6, 80, 100);
loglog(N_ref, 0.048 * N_ref.^4, '--', 'Color', CORAL, 'LineWidth', 1.2, ...
       'DisplayName', '0.048 N^4');

% Slope annotation
text(30, 1e5, 'slope = 4', 'FontSize', 11, 'Color', CORAL, ...
     'Rotation', 40, 'HorizontalAlignment', 'center');
hold off;

xlabel('N'); ylabel('\rho(\tilde{D}_2)');
title('Spectral radius of Chebyshev D^2 (Dirichlet BCs)');
legend('Location', 'northwest', 'FontSize', 10);
grid on; set(gca, 'GridAlpha', 0.3);

set(gcf, 'PaperPositionMode', 'auto');
print('-dpdf', 'cheb_cfl_scaling.pdf');
fprintf('Figure saved to cheb_cfl_scaling.pdf\n');
