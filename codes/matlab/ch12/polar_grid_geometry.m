% polar_grid_geometry.m - Grid geometry and the cost of clustering
% Chapter 12: Spectral Methods in Polar Coordinates
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi, UAE)
% Last modified: February 2026
%
% Computational Etude 12.1: Compares naive (r in [0,1]) and doubled
% (r in [-1,1]) Chebyshev-Fourier grids on the unit disk. Visualises
% grid clustering, area distribution, and CFL time-step scaling.

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
% Figure 12.1: Naive vs doubled polar grids on the disk
% -----------------------------------------------------------------------

Nr = 25;   % must be odd for doubled grid
M  = 20;   % angular points

% Chebyshev points on [-1, 1]
[~, x_cheb] = cheb_matrix(Nr);

% --- Naive grid: map x in [-1,1] to r in [0,1] via r = (x+1)/2 ---
r_naive = (x_cheb + 1) / 2;   % includes r=0 and r=1
theta   = linspace(0, 2*pi, M+1);
theta   = theta(1:M);          % M points, no endpoint

[RR_n, TT_n] = meshgrid(r_naive, theta);
XX_n = RR_n .* cos(TT_n);
YY_n = RR_n .* sin(TT_n);

% --- Doubled grid: r in [-1,1], use only r > 0, odd Nr so no r=0 ---
N2    = (Nr - 1) / 2;
r_pos = x_cheb(2:N2+1);   % r > 0 interior points

[RR_d, TT_d] = meshgrid(r_pos, theta);
XX_d = RR_d .* cos(TT_d);
YY_d = RR_d .* sin(TT_d);

% Median radius (circle containing half the points)
r_naive_interior  = r_naive(2:end-1);
r_median_naive    = median(r_naive_interior);
r_median_doubled  = median(r_pos);

figure('Position', [100 100 900 430]);

% --- Left panel: naive grid ---
subplot(1, 2, 1);
hold on;
% Disk boundary
th_circle = linspace(0, 2*pi, 200);
plot(cos(th_circle), sin(th_circle), 'k-', 'LineWidth', 1.2);
% Grid points
plot(XX_n(:), YY_n(:), 'k.', 'MarkerSize', 2.5);
% Median circle
plot(r_median_naive*cos(th_circle), r_median_naive*sin(th_circle), ...
    '--', 'Color', [0.84 0.15 0.16], 'LineWidth', 1.0);
hold off;
area_frac_naive = r_median_naive^2 * 100;
title(sprintf('Naive grid (r \\in [0, 1])\nMedian circle: %.0f%% of area', ...
    area_frac_naive), 'FontSize', 10);
xlim([-1.15 1.15]);
ylim([-1.15 1.15]);
axis equal off;

% --- Right panel: doubled grid ---
subplot(1, 2, 2);
hold on;
plot(cos(th_circle), sin(th_circle), 'k-', 'LineWidth', 1.2);
plot(XX_d(:), YY_d(:), 'k.', 'MarkerSize', 2.5);
plot(r_median_doubled*cos(th_circle), r_median_doubled*sin(th_circle), ...
    '--', 'Color', [0.84 0.15 0.16], 'LineWidth', 1.0);
hold off;
area_frac_doubled = r_median_doubled^2 * 100;
title(sprintf('Doubled grid (r \\in [-1, 1])\nMedian circle: %.0f%% of area', ...
    area_frac_doubled), 'FontSize', 10);
xlim([-1.15 1.15]);
ylim([-1.15 1.15]);
axis equal off;

% Save figure
exportgraphics(gcf, fullfile(fig_dir, 'polar_grids.pdf'), ...
    'ContentType', 'vector');

% -----------------------------------------------------------------------
% Figure 12.2: Area distribution per ring and CFL scaling
% -----------------------------------------------------------------------

figure('Position', [100 100 900 380]);

% --- Left: area per radial ring ---
subplot(1, 2, 1);

% Naive grid
r_sorted_naive = sort(r_naive);
n_naive = length(r_sorted_naive);
areas_naive     = zeros(n_naive - 1, 1);
midpoints_naive = zeros(n_naive - 1, 1);
for i = 1:n_naive-1
    areas_naive(i)     = pi * (r_sorted_naive(i+1)^2 - r_sorted_naive(i)^2);
    midpoints_naive(i) = 0.5 * (r_sorted_naive(i) + r_sorted_naive(i+1));
end

% Doubled grid: use sorted positive-r points, plus r=0 and r=1 as boundaries
r_sorted_doubled = sort([0; r_pos; 1]);
n_doubled = length(r_sorted_doubled);
areas_doubled     = zeros(n_doubled - 1, 1);
midpoints_doubled = zeros(n_doubled - 1, 1);
for i = 1:n_doubled-1
    areas_doubled(i)     = pi * (r_sorted_doubled(i+1)^2 - r_sorted_doubled(i)^2);
    midpoints_doubled(i) = 0.5 * (r_sorted_doubled(i) + r_sorted_doubled(i+1));
end

width = 0.012;
hold on;
bar(midpoints_naive - width, areas_naive, 2*width, ...
    'FaceColor', [0.12 0.47 0.71], 'FaceAlpha', 0.7);
bar(midpoints_doubled + width, areas_doubled, 2*width, ...
    'FaceColor', [0.84 0.15 0.16], 'FaceAlpha', 0.7);
hold off;
xlabel('$r$', 'Interpreter', 'latex');
ylabel('Area of radial annulus');
legend('Naive', 'Doubled', 'FontSize', 9);
title('Area per ring', 'FontSize', 10);

% --- Right: CFL time step vs Nr ---
subplot(1, 2, 2);

Nr_values = 7:2:51;   % odd values only
dt_naive_list   = zeros(size(Nr_values));
dt_doubled_list = zeros(size(Nr_values));

for k = 1:length(Nr_values)
    N = Nr_values(k);
    [~, x] = cheb_matrix(N);

    % Naive: r = (x+1)/2 in [0,1]
    r_n      = (x + 1) / 2;
    h_min_n  = min(abs(diff(r_n)));
    dt_naive_list(k) = h_min_n^2;   % CFL ~ h_min^2

    % Doubled: use positive interior points
    n2  = (N - 1) / 2;
    r_p = x(2:n2+1);
    h_min_d = min(abs(diff(sort(r_p))));
    dt_doubled_list(k) = h_min_d^2;
end

Nr_ref = double(Nr_values);
loglog(Nr_values, dt_naive_list, 'o-', 'MarkerSize', 3.5, ...
    'Color', [0.12 0.47 0.71], 'DisplayName', 'Naive');
hold on;
loglog(Nr_values, dt_doubled_list, 's-', 'MarkerSize', 3.5, ...
    'Color', [0.84 0.15 0.16], 'DisplayName', 'Doubled');
% Reference slopes
loglog(Nr_ref, 5e2 * Nr_ref.^(-4), 'k--', 'LineWidth', 0.8, ...
    'DisplayName', '$O(N^{-4})$');
loglog(Nr_ref, 2e1 * Nr_ref.^(-2), 'k:', 'LineWidth', 0.8, ...
    'DisplayName', '$O(N^{-2})$');
hold off;
xlabel('$N_r$', 'Interpreter', 'latex');
ylabel('$\Delta t_{\max} \propto h_{\min}^2$', 'Interpreter', 'latex');
legend('Location', 'best', 'FontSize', 8, 'Interpreter', 'latex');
title('CFL time-step scaling', 'FontSize', 10);

% Save figure
exportgraphics(gcf, fullfile(fig_dir, 'polar_area_cfl.pdf'), ...
    'ContentType', 'vector');

fprintf('Figures saved to: %s\n', fig_dir);
fprintf('  polar_grids.pdf\n');
fprintf('  polar_area_cfl.pdf\n');
