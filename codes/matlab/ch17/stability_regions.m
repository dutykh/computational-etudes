% stability_regions.m
% Chapter 17: Stability of Spectral Methods
%
% Computational Etude 17.2: Stability Regions and Spectral Eigenvalue Overlay.
%
% Two figures:
%   (1) 2x3 grid showing stability regions of six classical time-stepping
%       methods: Forward Euler, RK4, Leapfrog, Adams-Bashforth 2,
%       Crank-Nicolson, Backward Euler.
%   (2) 1x2 figure overlaying spectral eigenvalues on stability regions
%       of RK4 (Fourier) and Forward Euler (Chebyshev D^2).
%
% Author: Dr. Denys Dutykh
% Part of "Computational Etudes: A Spectral Approach"

clear; clc; close all;

NAVY       = [20, 45, 110]/255;
SKY        = [120, 150, 210]/255;
CORAL      = [231, 76, 60]/255;
TEAL       = [22, 160, 133]/255;
AMBER      = [230, 126, 34]/255;
LIGHT_NAVY = [197, 208, 232]/255;

%% ---- Figure 1: Stability regions (2x3 grid) ----
figure('Position', [50, 50, 1200, 800]);

% (a) Forward Euler
subplot(2, 3, 1);
plot_stability_contour(@(z) 1 + z, [-3, 1.5], [-2, 2], ...
    '(a) Forward Euler', NAVY, LIGHT_NAVY);

% (b) RK4
subplot(2, 3, 2);
plot_stability_contour(@(z) 1 + z + z.^2/2 + z.^3/6 + z.^4/24, ...
    [-5, 2], [-4, 4], '(b) RK4', NAVY, LIGHT_NAVY);

% (c) Leapfrog -- parametric boundary
subplot(2, 3, 3);
theta = linspace(0, 2*pi, 1000);
g = exp(1i*theta);
z_lf = (g.^2 - 1) ./ (2*g);
fill(real(z_lf), imag(z_lf), LIGHT_NAVY, 'FaceAlpha', 0.5, ...
     'EdgeColor', NAVY, 'LineWidth', 1.5);
hold on;
plot([min(real(z_lf))-0.5, max(real(z_lf))+0.5], [0, 0], ...
     'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
plot([0, 0], [-1.5, 1.5], 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
hold off;
xlim([-1.5, 1.5]); ylim([-1.5, 1.5]);
xlabel('Re(z)'); ylabel('Im(z)');
title('(c) Leapfrog');
axis equal; grid on; set(gca, 'GridAlpha', 0.2);

% (d) Adams-Bashforth 2 -- parametric boundary
subplot(2, 3, 4);
z_ab = (g.^2 - g) ./ (1.5*g - 0.5);
fill(real(z_ab), imag(z_ab), LIGHT_NAVY, 'FaceAlpha', 0.5, ...
     'EdgeColor', NAVY, 'LineWidth', 1.5);
hold on;
plot([-2, 1.5], [0, 0], 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
plot([0, 0], [-1.5, 1.5], 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
hold off;
xlim([-2, 1.5]); ylim([-1.5, 1.5]);
xlabel('Re(z)'); ylabel('Im(z)');
title('(d) Adams-Bashforth 2');
axis equal; grid on; set(gca, 'GridAlpha', 0.2);

% (e) Crank-Nicolson
subplot(2, 3, 5);
plot_stability_contour(@(z) (1 + z/2)./(1 - z/2), ...
    [-5, 1], [-4, 4], '(e) Crank-Nicolson', NAVY, LIGHT_NAVY);

% (f) Backward Euler
subplot(2, 3, 6);
plot_stability_contour(@(z) 1./(1 - z), ...
    [-1.5, 5], [-4, 4], '(f) Backward Euler', NAVY, LIGHT_NAVY);

set(gcf, 'PaperPositionMode', 'auto');
print('-dpdf', 'stability_regions.pdf');
fprintf('Figure saved to stability_regions.pdf\n');

%% ---- Figure 2: Spectra overlay ----
figure('Position', [100, 100, 1100, 500]);

% Left panel: RK4 + Fourier imaginary eigenvalues
subplot(1, 2, 1);
n_grid = 400;
xv = linspace(-5, 2, n_grid);
yv = linspace(-4, 4, n_grid);
[X, Y] = meshgrid(xv, yv);
Z = X + 1i*Y;
G = abs(1 + Z + Z.^2/2 + Z.^3/6 + Z.^4/24);
contourf(X, Y, G, [0, 1], 'FaceColor', LIGHT_NAVY, 'FaceAlpha', 0.4);
hold on;
contour(X, Y, G, [1, 1], 'LineColor', NAVY, 'LineWidth', 1.5);

% Fourier eigenvalues for advection
N_f = 32;
kf = [0:N_f/2-1, -N_f/2:-1];
lam_f = 1i * kf;

% Stable dt
dt_s = 2.8 / (N_f/2);
eigs_s = dt_s * lam_f;
plot(real(eigs_s), imag(eigs_s), 'o', 'Color', TEAL, 'MarkerSize', 4, ...
     'MarkerFaceColor', TEAL, ...
     'DisplayName', sprintf('dt = %.4f (stable)', dt_s));

% Unstable dt
dt_u = 3.5 / (N_f/2);
eigs_u = dt_u * lam_f;
plot(real(eigs_u), imag(eigs_u), 's', 'Color', CORAL, 'MarkerSize', 4, ...
     'MarkerFaceColor', CORAL, ...
     'DisplayName', sprintf('dt = %.4f (unstable)', dt_u));

plot([-5, 2], [0, 0], 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5, ...
     'HandleVisibility', 'off');
plot([0, 0], [-4, 4], 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5, ...
     'HandleVisibility', 'off');
hold off;
xlabel('Re(\Delta t \lambda)'); ylabel('Im(\Delta t \lambda)');
title(sprintf('(a) RK4 + Fourier advection (N=%d)', N_f));
axis equal; legend('Location', 'southwest', 'FontSize', 8);
grid on; set(gca, 'GridAlpha', 0.2);

% Right panel: Forward Euler + Chebyshev D^2 eigenvalues
subplot(1, 2, 2);
xv2 = linspace(-3, 1.5, n_grid);
yv2 = linspace(-2, 2, n_grid);
[X2, Y2] = meshgrid(xv2, yv2);
Z2 = X2 + 1i*Y2;
G2 = abs(1 + Z2);
contourf(X2, Y2, G2, [0, 1], 'FaceColor', LIGHT_NAVY, 'FaceAlpha', 0.4);
hold on;
contour(X2, Y2, G2, [1, 1], 'LineColor', NAVY, 'LineWidth', 1.5);

% Chebyshev D^2 eigenvalues (Dirichlet)
N_c = 16;
[D, ~] = cheb(N_c);
D2 = D * D;
D2_int = D2(2:N_c, 2:N_c);
eigs_c = eig(D2_int);

% Scale dt so a few eigenvalues fit in the disk
phys_mask = abs(real(eigs_c)) < 200;
out_mask  = ~phys_mask;
dt_c = 0.8 / max(abs(real(eigs_c(1:5))));
scaled_eigs = dt_c * eigs_c;

plot(real(scaled_eigs(phys_mask)), imag(scaled_eigs(phys_mask)), 'o', ...
     'Color', TEAL, 'MarkerSize', 5, 'MarkerFaceColor', TEAL, ...
     'DisplayName', 'physical modes');
if any(out_mask)
    plot(real(scaled_eigs(out_mask)), imag(scaled_eigs(out_mask)), 's', ...
         'Color', CORAL, 'MarkerSize', 5, 'MarkerFaceColor', CORAL, ...
         'DisplayName', 'outlier modes');
end

plot([-3, 1.5], [0, 0], 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5, ...
     'HandleVisibility', 'off');
plot([0, 0], [-2, 2], 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5, ...
     'HandleVisibility', 'off');
hold off;
xlabel('Re(\Delta t \lambda)'); ylabel('Im(\Delta t \lambda)');
title(sprintf('(b) Forward Euler + Chebyshev D^2 (N=%d)', N_c));
axis equal; legend('Location', 'northwest', 'FontSize', 8);
grid on; set(gca, 'GridAlpha', 0.2);

set(gcf, 'PaperPositionMode', 'auto');
print('-dpdf', 'spectra_overlay.pdf');
fprintf('Figure saved to spectra_overlay.pdf\n');


%% ---- Helper function: contour-based stability region ----
function plot_stability_contour(g_fun, xlims, ylims, ttl, navy, light_navy)
    n_grid = 400;
    xv = linspace(xlims(1), xlims(2), n_grid);
    yv = linspace(ylims(1), ylims(2), n_grid);
    [X, Y] = meshgrid(xv, yv);
    Z = X + 1i*Y;
    G = abs(g_fun(Z));
    contourf(X, Y, G, [0, 1], 'FaceColor', light_navy, 'FaceAlpha', 0.5);
    hold on;
    contour(X, Y, G, [1, 1], 'LineColor', navy, 'LineWidth', 1.5);
    plot([xlims(1), xlims(2)], [0, 0], 'Color', [0.5 0.5 0.5], ...
         'LineWidth', 0.5);
    plot([0, 0], [ylims(1), ylims(2)], 'Color', [0.5 0.5 0.5], ...
         'LineWidth', 0.5);
    hold off;
    xlim(xlims); ylim(ylims);
    xlabel('Re(z)'); ylabel('Im(z)');
    title(ttl);
    axis equal; grid on; set(gca, 'GridAlpha', 0.2);
end
