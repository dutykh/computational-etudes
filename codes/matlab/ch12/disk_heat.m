% disk_heat.m - Heat equation on the unit disk
% Chapter 12: Spectral Methods in Polar Coordinates
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi, UAE)
% Last modified: February 2026
%
% Computational Etude 12.5: Solves the heat equation u_t = alpha*Delta*u + S(r,theta)
% on the unit disk with Crank-Nicolson time stepping. Models heat dissipation
% on a semiconductor wafer with a localised laser source.

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
% Parameters
% -----------------------------------------------------------------------

Nr      = 25;      % radial points (odd)
M       = 20;      % angular points (even)
alpha   = 0.1;     % thermal diffusivity
dt      = 0.01;    % time step
t_final = 20.0;
n_steps = round(t_final / dt);

N2   = (Nr - 1) / 2;
ndof = N2 * M;

% Build the Laplacian
[L, r_pos, theta] = laplacian_polar(Nr, M);

% Grid
[RR, TT] = ndgrid(r_pos, theta);
XX = RR .* cos(TT);
YY = RR .* sin(TT);

% Source: localised laser spot at (0.3, 0.0)
x_src = 0.3;
y_src = 0.0;
S_vals = 2.0 * exp(-30.0 * ((XX - x_src).^2 + (YY - y_src).^2));
S_vec  = S_vals(:);

% -----------------------------------------------------------------------
% Crank-Nicolson time stepping
% -----------------------------------------------------------------------

% (I - alpha*dt/2*L) u^{n+1} = (I + alpha*dt/2*L) u^n + dt*S
I_mat = eye(ndof);
A = I_mat - alpha * dt / 2 * L;    % LHS matrix
B = I_mat + alpha * dt / 2 * L;    % RHS matrix part

% LU factorise once (reuse at every step)
[L_lu, U_lu, P_lu] = lu(A);

% Initial condition: u = 0 (cold start)
u = zeros(ndof, 1);

% Store snapshots at selected times
snapshot_times = [0.0, 0.5, 2.0, 5.0, 10.0, 20.0];
n_snaps        = length(snapshot_times);
snapshots      = cell(n_snaps, 1);
snap_found     = false(n_snaps, 1);

% Store t=0 snapshot
snapshots{1} = u;
snap_found(1) = true;

max_temp = zeros(n_steps + 1, 1);
times    = zeros(n_steps + 1, 1);
max_temp(1) = 0.0;
times(1)    = 0.0;

t = 0.0;
for step = 1:n_steps
    rhs = B * u + dt * S_vec;
    u   = U_lu \ (L_lu \ (P_lu * rhs));
    t   = step * dt;

    max_temp(step + 1) = max(u);
    times(step + 1)    = t;

    % Check if this is a snapshot time (within tolerance)
    for s = 1:n_snaps
        if ~snap_found(s) && abs(t - snapshot_times(s)) < dt / 2
            snapshots{s} = u;
            snap_found(s) = true;
        end
    end
end

% Ensure final snapshot is stored
if ~snap_found(end)
    snapshots{end} = u;
    snap_found(end) = true;
end

% -----------------------------------------------------------------------
% Figure 12.9: Heat equation snapshots at 6 time levels
% -----------------------------------------------------------------------

% Extend grid for plotting
r_ext     = [1.0; r_pos];
theta_ext = [theta; 2*pi];
[RR_ext, TT_ext] = ndgrid(r_ext, theta_ext);
XX_ext = RR_ext .* cos(TT_ext);
YY_ext = RR_ext .* sin(TT_ext);

% Global color scale
u_max = 0;
for s = 1:n_snaps
    u_max = max(u_max, max(abs(snapshots{s})));
end
if u_max < 1e-10
    u_max = 1.0;
end

figure('Position', [100 100 1000 680]);
th_circle = linspace(0, 2*pi, 200);

for s = 1:n_snaps
    ax = subplot(2, 3, s);
    U_snap = reshape(snapshots{s}, [N2, M]);

    % Extend for plotting
    U_ext = zeros(N2 + 1, M + 1);
    U_ext(2:end, 1:M) = U_snap;
    U_ext(:, M+1)     = U_ext(:, 1);

    contourf(XX_ext, YY_ext, U_ext, 20, 'LineStyle', 'none');
    colormap(ax, hot_colormap());
    caxis([0 u_max]);
    hold on;
    plot(cos(th_circle), sin(th_circle), 'k-', 'LineWidth', 1.0);
    hold off;

    xlim([-1.15 1.15]);
    ylim([-1.15 1.15]);
    axis equal off;
    title(sprintf('$t = %.1f$', snapshot_times(s)), ...
        'Interpreter', 'latex', 'FontSize', 10);
end

% Single colorbar
h = colorbar;
h.Label.String = 'Temperature $u$';
h.Label.Interpreter = 'latex';
h.Position = [0.92 0.15 0.02 0.7];

sgtitle('Heat dissipation on a disk with localised source', 'FontSize', 11);

exportgraphics(gcf, fullfile(fig_dir, 'heat_disk_snapshots.pdf'), ...
    'ContentType', 'vector');

% -----------------------------------------------------------------------
% Figure 12.10: Temperature evolution and steady-state comparison
% -----------------------------------------------------------------------

figure('Position', [100 100 900 380]);

% Left: maximum temperature vs time
subplot(1, 2, 1);
plot(times, max_temp, '-', 'Color', [0.12 0.47 0.71], 'LineWidth', 1.0);
xlabel('Time $t$', 'Interpreter', 'latex');
ylabel('$\max_{x,y} \, u(x,y,t)$', 'Interpreter', 'latex');
title('Maximum temperature evolution', 'FontSize', 10);
grid on;
set(gca, 'GridAlpha', 0.3);

% Right: compare steady state with direct Poisson solve
% Steady state satisfies alpha*Delta*u + S = 0, i.e. Delta*u = -S/alpha
u_poisson = L \ (-S_vec / alpha);

% Plot radial profiles at theta = 0
angle_idx  = 1;
U_final    = reshape(snapshots{end}, [N2, M]);
U_poisson  = reshape(u_poisson, [N2, M]);

subplot(1, 2, 2);
plot(r_pos, U_final(:, angle_idx), '-', 'Color', [0.12 0.47 0.71], ...
    'LineWidth', 1.5, 'DisplayName', 'CN steady state');
hold on;
plot(r_pos, U_poisson(:, angle_idx), '--', 'Color', [0.84 0.15 0.16], ...
    'LineWidth', 1.5, 'DisplayName', 'Direct Poisson solve');
hold off;
xlabel('$r$', 'Interpreter', 'latex');
ylabel('$u(r, 0)$', 'Interpreter', 'latex');
title('Steady state vs Poisson solution ($\theta = 0$)', ...
    'Interpreter', 'latex', 'FontSize', 10);
legend('FontSize', 9, 'Location', 'best');
set(gca, 'XDir', 'reverse');

exportgraphics(gcf, fullfile(fig_dir, 'heat_disk_energy.pdf'), ...
    'ContentType', 'vector');

% Print summary
err_ss = max(abs(snapshots{end} - u_poisson));
fprintf('Steady-state vs Poisson max error: %.2e\n', err_ss);
fprintf('Max temperature at t=%.1f: %.6f\n', t_final, max_temp(end));
fprintf('\nFigures saved to: %s\n', fig_dir);
fprintf('  heat_disk_snapshots.pdf\n');
fprintf('  heat_disk_energy.pdf\n');


% =====================================================================
% Helper function
% =====================================================================

function cmap = hot_colormap()
% HOT_COLORMAP  Standard 'hot' colormap (black-red-yellow-white).
    cmap = hot(256);
end
