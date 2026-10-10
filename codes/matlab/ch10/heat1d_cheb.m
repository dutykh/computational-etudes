%% heat1d_cheb.m
%
% Figure 10.5: 1D heat equation on Chebyshev grid
%
% PDE: u_t = kappa * u_xx, -1 < x < 1, t > 0
% BCs: u(-1, t) = u(1, t) = 0 (Dirichlet)
% ICs: u(x, 0) = exp(-20*(x+0.5)^2) + 0.5*exp(-30*(x-0.4)^2) (two bumps)
%
% Uses Crank--Nicolson time stepping (implicit, A-stable).
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
% Part of "Computational Etudes: A Spectral Approach"
% https://github.com/dutykh/computational-etudes

clear; close all; clc;

%% Configuration
N = 128;          % Chebyshev intervals
kappa = 1.0;      % Diffusion coefficient
tmax = 0.5;       % Final time
dt = 0.005;       % Time step

% Colors
NAVY = [0.078 0.176 0.431];
SKY = [0.475 0.588 0.824];
CORAL = [0.906 0.298 0.235];
TEAL = [0.086 0.627 0.522];
PURPLE = [0.608 0.349 0.714];

% Output path
script_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch10', 'matlab');
output_file = fullfile(output_dir, 'heat1d_evolution.pdf');

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% Setup grid and operators
[D, x] = cheb_matrix(N);
D2 = D * D;

% Interior block (Dirichlet BCs)
D2_int = D2(2:N, 2:N);
I_int = eye(N-1);

% Crank--Nicolson matrices
A = I_int - 0.5 * kappa * dt * D2_int;
B = I_int + 0.5 * kappa * dt * D2_int;

% Pre-factorize A (constant across all time steps)
dA = decomposition(A, 'lu');

%% Initial condition: two bumps
u_full = exp(-20 * (x + 0.5).^2) + 0.5 * exp(-30 * (x - 0.4).^2);
u = u_full(2:N);  % Interior values

%% Time stepping
nsteps = ceil(tmax / dt);
n_save = 80;
save_every = max(1, floor(nsteps / n_save));

U_save = zeros(N+1, n_save + 1);
t_save = zeros(1, n_save + 1);
U_save(:, 1) = u_full;
t_save(1) = 0;
save_idx = 2;

t = 0;
for n = 1:nsteps
    % Solve CN system using pre-factorized A
    u = dA \ (B * u);
    t = t + dt;

    % Save snapshot
    if mod(n, save_every) == 0 && save_idx <= n_save + 1
        u_full = zeros(N+1, 1);
        u_full(2:N) = u;
        U_save(:, save_idx) = u_full;
        t_save(save_idx) = t;
        save_idx = save_idx + 1;
    end
end

%% Create figure
fig = figure('Units', 'inches', 'Position', [1, 1, 12, 6]);

% Left panel: 3D surface
subplot(1, 2, 1);
[X, T] = meshgrid(x, t_save);
surf(X, T, U_save', 'EdgeColor', 'none', 'FaceAlpha', 0.9);
colormap('hot');
xlabel('x', 'FontSize', 11);
ylabel('t', 'FontSize', 11);
zlabel('u(x, t)', 'FontSize', 11);
title('Heat Diffusion: Space-Time View', 'FontSize', 12);
view([-60, 25]);
colorbar;

% Right panel: time snapshots
subplot(1, 2, 2);
t_targets = [0.0, 0.05, 0.1, 0.2, 0.5];
colors = {NAVY, SKY, TEAL, CORAL, PURPLE};

hold on;
for i = 1:length(t_targets)
    [~, idx] = min(abs(t_save - t_targets(i)));
    plot(x, U_save(:, idx), '-', 'Color', colors{i}, 'LineWidth', 1.5, ...
         'DisplayName', sprintf('t = %.2f', t_save(idx)));
end
hold off;

xlabel('x', 'FontSize', 11);
ylabel('u(x, t)', 'FontSize', 11);
title('Heat Diffusion: Time Snapshots', 'FontSize', 12);
legend('Location', 'northeast', 'FontSize', 10);
xlim([-1, 1]);
ylim([0, 1.1]);
grid on;
box on;

%% Save figure
exportgraphics(fig, output_file, 'ContentType', 'vector');
exportgraphics(fig, strrep(output_file, '.pdf', '.png'), 'Resolution', 300);

fprintf('Figure saved to: %s\n', output_file);

close(fig);
