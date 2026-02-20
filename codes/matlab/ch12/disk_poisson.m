% disk_poisson.m - Poisson equation on the unit disk
% Chapter 12: Spectral Methods in Polar Coordinates
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi, UAE)
% Last modified: February 2026
%
% Computational Etude 12.4: Solves -Delta u = f on the disk with u = 0
% at r = 1, using an off-centre Gaussian source. Also compares radially
% symmetric vs asymmetric forcing to verify symmetry preservation.

clear; close all;

% Add Chebyshev matrix from ch07
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'ch07'));

% Output directory
fig_dir = fullfile(fileparts(mfilename('fullpath')), ...
    '..', '..', '..', 'textbook', 'figures', 'ch12', 'matlab');
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

% -----------------------------------------------------------------------
% Build the Laplacian
% -----------------------------------------------------------------------

Nr = 31;   % radial points (odd)
M  = 40;   % angular points (even, more than eigenmodes for smoother plots)
N2 = (Nr - 1) / 2;

[L, r_pos, theta] = laplacian_polar(Nr, M);

% Grid in (r, theta) space
[RR, TT] = ndgrid(r_pos, theta);
XX = RR .* cos(TT);
YY = RR .* sin(TT);

% -----------------------------------------------------------------------
% Poisson solve with off-centre Gaussian forcing
% -----------------------------------------------------------------------

% Source: localized load at (x0, y0) = (0.4, 0.2)
x0   = 0.4;
y0   = 0.2;
beta = 40.0;
f_vals = exp(-beta * ((XX - x0).^2 + (YY - y0).^2));

% Flatten and solve: L u = Delta u, and we want -Delta u = f
% so L u = -f, i.e. u = L \ (-f)
f_vec = -f_vals(:);
u_vec = L \ f_vec;
U     = reshape(u_vec, [N2, M]);

% Extend for plotting: add boundary (r=1, u=0) and wrap in theta
U_ext = zeros(N2 + 1, M + 1);
U_ext(2:end, 1:M) = U;
U_ext(:, M+1)     = U_ext(:, 1);

r_ext     = [1.0; r_pos];
theta_ext = [theta; 2*pi];
[RR_ext, TT_ext] = ndgrid(r_ext, theta_ext);
XX_ext = RR_ext .* cos(TT_ext);
YY_ext = RR_ext .* sin(TT_ext);

% -----------------------------------------------------------------------
% Figure 12.6: 3D surface plot of the Poisson solution
% -----------------------------------------------------------------------

figure('Position', [100 100 640 500]);
ax = axes;

surf(XX_ext, YY_ext, U_ext, 'EdgeColor', 'k', 'LineWidth', 0.2, ...
    'FaceAlpha', 0.9);
colormap(ax, viridis_colormap());
hold on;
% Draw disk boundary circle at z=0
th_circle = linspace(0, 2*pi, 200);
plot3(cos(th_circle), sin(th_circle), zeros(1, 200), 'k-', 'LineWidth', 1.0);
hold off;

xlabel('$x$', 'Interpreter', 'latex');
ylabel('$y$', 'Interpreter', 'latex');
zlabel('$u$', 'Interpreter', 'latex');
view(25, -55);
title('Poisson equation on the disk: off-centre Gaussian load', ...
    'FontSize', 10);

exportgraphics(gcf, fullfile(fig_dir, 'poisson_disk_surface.pdf'), ...
    'ContentType', 'vector');

% -----------------------------------------------------------------------
% Figure 12.7: Filled contour plot on the disk
% -----------------------------------------------------------------------

figure('Position', [100 100 500 460]);

contourf(XX_ext, YY_ext, U_ext, 20);
hold on;
contour(XX_ext, YY_ext, U_ext, 20, 'k', 'LineWidth', 0.3);
% Disk boundary
plot(cos(th_circle), sin(th_circle), 'k-', 'LineWidth', 1.2);
% Mark source location
plot(x0, y0, 'rp', 'MarkerSize', 10, 'MarkerFaceColor', 'r', ...
    'MarkerEdgeColor', 'k');
text(x0 + 0.08, y0 + 0.08, 'source', 'FontSize', 9, 'Color', 'r');
hold off;

colormap(viridis_colormap());
xlim([-1.15 1.15]);
ylim([-1.15 1.15]);
axis equal;
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$y$', 'Interpreter', 'latex');
title('Poisson solution: filled contour', 'FontSize', 10);
cb = colorbar;
cb.Label.String = '$u(x,y)$';
cb.Label.Interpreter = 'latex';

exportgraphics(gcf, fullfile(fig_dir, 'poisson_disk_contour.pdf'), ...
    'ContentType', 'vector');

% -----------------------------------------------------------------------
% Figure 12.8: Radial symmetry test
% -----------------------------------------------------------------------

% Case 1: radially symmetric Gaussian centred at origin
f_sym     = exp(-beta * (XX.^2 + YY.^2));
f_vec_sym = -f_sym(:);
u_vec_sym = L \ f_vec_sym;
U_sym     = reshape(u_vec_sym, [N2, M]);

% Case 2: off-centre Gaussian (already computed as U)

figure('Position', [100 100 900 380]);

% Left: radially symmetric case -- profiles should collapse
subplot(1, 2, 1);
n_angles      = 8;
angle_indices = round(linspace(1, M, n_angles + 1));
angle_indices = angle_indices(1:n_angles);
hold on;
for k = angle_indices
    plot(r_pos, U_sym(:, k), '-', 'Color', [0.12 0.47 0.71], ...
        'LineWidth', 0.8);
end
hold off;
xlabel('$r$', 'Interpreter', 'latex');
ylabel('$u(r, \theta_k)$', 'Interpreter', 'latex');
title({'Centred Gaussian: radial profiles', '(all angles coincide)'}, ...
    'FontSize', 10);
set(gca, 'XDir', 'reverse');   % r decreasing from left (r=1) to right (r~0)

% Right: off-centre case -- profiles differ
subplot(1, 2, 2);
hold on;
for ki = 1:length(angle_indices)
    k = angle_indices(ki);
    if ki <= 3
        plot(r_pos, U(:, k), 'LineWidth', 0.8, ...
            'DisplayName', sprintf('$\\theta = %.2f$', theta(k)));
    else
        plot(r_pos, U(:, k), 'LineWidth', 0.8, ...
            'HandleVisibility', 'off');
    end
end
hold off;
xlabel('$r$', 'Interpreter', 'latex');
ylabel('$u(r, \theta_k)$', 'Interpreter', 'latex');
title({'Off-centre Gaussian: radial profiles', '(angles differ)'}, ...
    'FontSize', 10);
legend('FontSize', 8, 'Location', 'northeast', 'Interpreter', 'latex');
set(gca, 'XDir', 'reverse');

exportgraphics(gcf, fullfile(fig_dir, 'radial_symmetry_test.pdf'), ...
    'ContentType', 'vector');

fprintf('Figures saved to: %s\n', fig_dir);
fprintf('  poisson_disk_surface.pdf\n');
fprintf('  poisson_disk_contour.pdf\n');
fprintf('  radial_symmetry_test.pdf\n');
fprintf('\nMax deflection: %.6f at off-centre load\n', max(U(:)));
fprintf('Max deflection: %.6f at centred load\n', max(U_sym(:)));


% =====================================================================
% Helper function
% =====================================================================

function cmap = viridis_colormap()
% VIRIDIS_COLORMAP  Approximation of matplotlib's viridis colormap.
    n = 256;
    nodes = [0.267004 0.004874 0.329415;
             0.282327 0.140926 0.457517;
             0.253935 0.265254 0.529983;
             0.206756 0.371758 0.553117;
             0.163625 0.471133 0.558148;
             0.127568 0.566949 0.550556;
             0.134692 0.658636 0.517649;
             0.266941 0.748751 0.440573;
             0.477504 0.821444 0.318195;
             0.741388 0.873449 0.149561;
             0.993248 0.906157 0.143936];
    x_nodes = linspace(0, 1, size(nodes, 1));
    xi = linspace(0, 1, n);
    cmap = interp1(x_nodes, nodes, xi);
end
