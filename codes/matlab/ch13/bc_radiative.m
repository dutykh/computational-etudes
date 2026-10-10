%% bc_radiative.m - Radiative cooling with nonlinear boundary condition
%
% Chapter 13: Boundary Conditions in Spectral Methods
%
% Heat equation with radiative (Stefan-Boltzmann) boundary condition:
%
%   u_t = u_xx,           x in (0, 1),  t > 0
%   u(0, t) = 1           (hot wall, Dirichlet)
%   u_x(1, t) + sigma * u(1, t)^4 = 0   (radiative cooling, nonlinear Robin)
%
% Initial condition: u(x, 0) = 1 (uniform temperature)
%
% Mapping to Chebyshev domain [-1, 1]:
%   Physical x in [0, 1] maps to xi in [-1, 1] via x = (xi + 1)/2
%   D_x = 2 * D_xi
%
% Note: x(1) = 1 in Chebyshev ordering maps to physical x = 1 (radiating end),
%       x(N+1) = -1 maps to physical x = 0 (hot wall).
%
% Time stepping: Backward Euler with Newton iteration for the nonlinear BC.
%
% Generates Figure: radiative_cooling.pdf (2 panels)
%   Left:  temperature profiles at t = 0, 0.1, 0.5, 1, 5 for sigma = 1
%   Right: u(1, t) vs t for sigma = 0.1, 1, 10
%
% Author: Dr. Denys Dutykh
%         Mathematics Department
%         Khalifa University of Science and Technology
%         Abu Dhabi, UAE
%
% Part of the book "Computational Etudes: A Spectral Approach"
% https://github.com/dutykh/computational-etudes
%
% Last modified: February 2026

clear;
close all;
clc;

% Add path to Chapter 7 functions
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'ch07'));

%% Publication-quality figure settings
set(groot, 'DefaultAxesFontSize', 10);
set(groot, 'DefaultAxesFontName', 'CMU Serif');
set(groot, 'DefaultTextFontSize', 10);
set(groot, 'DefaultTextFontName', 'CMU Serif');
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');
set(groot, 'DefaultAxesLineWidth', 0.8);
set(groot, 'DefaultLineLineWidth', 1.5);

%% Color scheme
NAVY   = [20,  45, 110] / 255;
SKY    = [120, 150, 210] / 255;
CORAL  = [231,  76,  60] / 255;
TEAL   = [22,  160, 133] / 255;
PURPLE = [142,  68, 173] / 255;
ORANGE = [230, 126,  34] / 255;

%% Output directory
script_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch13', 'matlab');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% Parameters
N       = 40;        % Chebyshev intervals
t_final = 5.0;      % Final time
dt      = 0.01;     % Time step
nsteps  = round(t_final / dt);
newton_tol = 1e-12;
newton_max = 20;

%% Setup grid and operators
[D_xi, xi] = cheb_matrix(N);

% Physical derivative: D_x = 2 * D_xi (mapping [0,1] -> [-1,1])
D  = 2 * D_xi;
D2 = D * D;

% Physical x coordinates: x = (xi + 1)/2
x_phys = (xi + 1) / 2;

%% =====================================================================
%  Part 1: Temperature profiles for sigma = 1
%  =====================================================================
sigma = 1.0;
snapshot_times = [0, 0.1, 0.5, 1.0, 5.0];
snapshots = {};

[u_history, t_history] = solve_radiative(D, D2, xi, N, dt, nsteps, sigma, ...
                                          newton_tol, newton_max);

% Extract snapshots
for k = 1:length(snapshot_times)
    [~, idx] = min(abs(t_history - snapshot_times(k)));
    snapshots{k} = u_history(:, idx); %#ok<SAGROW>
end

%% Part 2: Tip temperature for sigma = 0.1, 1, 10
sigma_values = [0.1, 1.0, 10.0];
tip_temps = cell(length(sigma_values), 1);

for j = 1:length(sigma_values)
    [u_hist, t_hist] = solve_radiative(D, D2, xi, N, dt, nsteps, ...
                                        sigma_values(j), newton_tol, newton_max);
    % Tip temperature: u at x(1) = 1 (physical x = 1, radiating end)
    tip_temps{j} = u_hist(1, :);
end

%% Create figure
fig = figure('Units', 'inches', 'Position', [1, 1, 11, 4.5]);

% ---- Left panel: temperature profiles ----
subplot(1, 2, 1);

colors_snap = {SKY, TEAL, [0.18 0.80 0.44], CORAL, NAVY};

hold on;
for k = 1:length(snapshot_times)
    plot(x_phys, snapshots{k}, '-', 'Color', colors_snap{k}, 'LineWidth', 1.5, ...
         'DisplayName', sprintf('$t = %g$', snapshot_times(k)));
end
hold off;

xlabel('$x$');
ylabel('$u(x, t)$');
title(sprintf('Temperature Profiles ($\\sigma = %g$)', sigma), 'FontSize', 11);
legend('Location', 'southwest', 'FontSize', 9);
xlim([0, 1]);
ylim([0, 1.1]);
grid on;
set(gca, 'GridAlpha', 0.3);
box on;

% ---- Right panel: tip temperature vs time ----
subplot(1, 2, 2);

colors_sigma = {TEAL, CORAL, PURPLE};

hold on;
for j = 1:length(sigma_values)
    plot(t_hist, tip_temps{j}, '-', 'Color', colors_sigma{j}, 'LineWidth', 1.5, ...
         'DisplayName', sprintf('$\\sigma = %g$', sigma_values(j)));
end
hold off;

xlabel('$t$');
ylabel('$u(1, t)$');
title('Radiating End Temperature', 'FontSize', 11);
legend('Location', 'northeast', 'FontSize', 9);
xlim([0, t_final]);
ylim([0, 1.1]);
grid on;
set(gca, 'GridAlpha', 0.3);
box on;

sgtitle('Radiative Cooling: $u_t = u_{xx}$, $u_x(1,t) + \sigma u(1,t)^4 = 0$', ...
        'FontSize', 12);

%% Save figure
exportgraphics(fig, fullfile(output_dir, 'radiative_cooling.pdf'), ...
               'ContentType', 'vector');
fprintf('Figure saved to: %s\n', fullfile(output_dir, 'radiative_cooling.pdf'));
close(fig);


%% =====================================================================
%  Local functions
%  =====================================================================

function [u_history, t_history] = solve_radiative(D, D2, xi, N, dt, nsteps, ...
                                                   sigma, newton_tol, newton_max)
% SOLVE_RADIATIVE  Solve the heat equation with radiative BC using backward Euler
%                  and Newton iteration.
%
% In Chebyshev ordering:
%   x(1) = +1   -> physical x = 1 (radiating end)
%   x(N+1) = -1 -> physical x = 0 (hot wall, Dirichlet u = 1)

    % IC: u = 1 everywhere
    u = ones(N + 1, 1);

    % Save every few steps
    n_save = min(400, nsteps);
    save_every = max(1, floor(nsteps / n_save));

    u_history = u;
    t_history = 0;

    for n = 1:nsteps
        % Backward Euler: (I - dt * D2) * u^{n+1} = u^n
        % with nonlinear BC at row 1: D(1,:)*u + sigma*u(1)^4 = 0
        % and Dirichlet BC at row N+1: u(N+1) = 1

        % Newton iteration starting from current solution
        u_new = u;

        for iter = 1:newton_max
            % Residual of the full system
            % Interior: (I - dt*D2)*u_new - u = 0
            % Row 1 (radiating): D(1,:)*u_new + sigma*u_new(1)^4 = 0
            % Row N+1 (hot wall): u_new(N+1) - 1 = 0

            I_mat = eye(N + 1);
            F = (I_mat - dt * D2) * u_new - u;

            % Overwrite boundary rows
            F(1) = D(1, :) * u_new + sigma * u_new(1)^4;
            F(N + 1) = u_new(N + 1) - 1.0;

            % Jacobian
            J = I_mat - dt * D2;

            % Row 1: d/du [D(1,:)*u + sigma*u(1)^4]
            J(1, :)  = D(1, :);
            J(1, 1)  = J(1, 1) + 4 * sigma * u_new(1)^3;

            % Row N+1: identity row
            J(N + 1, :)     = 0;
            J(N + 1, N + 1) = 1.0;

            % Newton update
            delta = J \ (-F);
            u_new = u_new + delta;

            if norm(delta, inf) < newton_tol
                break;
            end
        end

        u = u_new;

        % Save snapshot
        if mod(n, save_every) == 0
            u_history = [u_history, u]; %#ok<AGROW>
            t_history = [t_history, n * dt]; %#ok<AGROW>
        end
    end
end
