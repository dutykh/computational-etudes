%% kdv_soliton.m
%
% Chapter 11: Fourier Pseudospectral Methods for Periodic PDEs
%
% Fourier pseudospectral solver for the KdV equation with integrating factor RK4:
%   u_t + 6 u u_x + u_xxx = 0 on [0, L), periodic.
%
% Simulates the elastic collision of two solitons of different amplitudes.
% Soliton solution: u(x,t) = 2*a^2 * sech^2(a*(x - x0 - 4*a^2*t))
%
% Generates figure: kdv_snapshots.pdf (snapshots before, during, after collision)
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
% Part of "Computational Etudes: A Spectral Approach"
% https://github.com/dutykh/computational-etudes
% Date: February 2026

clear; close all; clc;

%% Configuration
N  = 512;
L  = 80.0;
a1 = 1.0;  a2 = 0.5;
x1 = 20.0; x2 = 40.0;   % Taller soliton behind
t_final = 25.0;
dt = 0.002;

NAVY  = [0.078 0.176 0.431];
SKY   = [0.471 0.588 0.824];
CORAL = [0.906 0.298 0.235];
TEAL  = [0.086 0.627 0.522];

% Output path
script_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch11', 'matlab');

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% Setup grid and wavenumbers
x = L * (0:N-1)' / N;

% Wavenumbers scaled for domain [0, L)
if mod(N, 2) == 0
    kn = [0:N/2-1, 0, -N/2+1:-1]';
else
    kn = [0:(N-1)/2, -(N-1)/2:-1]';
end
k = (2 * pi / L) * kn;

% Linear symbol: lambda_k = i * k^3
lam = 1i * k.^3;

% 2/3 dealiasing mask
cutoff = floor(N / 3);
mask = true(N, 1);
if cutoff > 0
    mask(cutoff+2:N-cutoff) = false;
end

%% Initial condition: two solitons
u0 = soliton(x, a1, x1) + soliton(x, a2, x2);
v_hat = fft(u0);

%% Time stepping: integrating factor RK4 (cumulative form)
n_steps = ceil(t_final / dt);
dt = t_final / n_steps;

n_save = 200;
save_every = max(1, floor(n_steps / n_save));
U_save = u0;
t_save = 0;

t = 0;
for n = 1:n_steps
    k1 = G_rhs(t,            v_hat,                 lam, k, mask);
    k2 = G_rhs(t + 0.5 * dt, v_hat + 0.5 * dt * k1, lam, k, mask);
    k3 = G_rhs(t + 0.5 * dt, v_hat + 0.5 * dt * k2, lam, k, mask);
    k4 = G_rhs(t + dt,       v_hat + dt * k3,         lam, k, mask);
    v_hat = v_hat + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
    t = t + dt;

    if mod(n, save_every) == 0
        u_hat = exp(lam * t) .* v_hat;
        u = real(ifft(u_hat));
        U_save = [U_save, u]; %#ok<AGROW>
        t_save = [t_save, t]; %#ok<AGROW>
    end
end

%% Create figure: 2x2 snapshots
fig = figure('Units', 'inches', 'Position', [1, 1, 12, 8]);

t_targets = [0.0, 8.0, 12.0, t_final];
titles = {'Initial condition ($t = 0$)', ...
          'Before collision', ...
          'During collision', ...
          sprintf('After collision ($t = %d$)', round(t_final))};
colors_snap = {NAVY, TEAL, CORAL, NAVY};

for i = 1:4
    subplot(2, 2, i);
    [~, idx] = min(abs(t_save - t_targets(i)));
    plot(x, U_save(:, idx), '-', 'Color', colors_snap{i}, 'LineWidth', 1.5);

    title(titles{i}, 'Interpreter', 'latex', 'FontSize', 12);
    xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 11);
    ylabel('$u(x, t)$', 'Interpreter', 'latex', 'FontSize', 11);
    xlim([0, L]);
    ylim([-0.2, 2.5]);
    grid on;
    set(gca, 'GridAlpha', 0.3);
    box on;
end

sgtitle(sprintf('KdV Two-Soliton Collision ($a_1 = %g$, $a_2 = %g$, $N = %d$)', ...
        a1, a2, N), 'Interpreter', 'latex', 'FontSize', 13);

%% Save figure
exportgraphics(fig, fullfile(output_dir, 'kdv_snapshots.pdf'), ...
               'ContentType', 'vector');
exportgraphics(fig, fullfile(output_dir, 'kdv_snapshots.png'), ...
               'Resolution', 300);
fprintf('Figure saved: %s\n', fullfile(output_dir, 'kdv_snapshots.pdf'));
close(fig);


%% -----------------------------------------------------------------------
%  Helper functions
%  -----------------------------------------------------------------------

function u = soliton(x, a, x0)
    % KdV soliton: 2*a^2 * sech^2(a*(x - x0))
    u = 2.0 * a^2 ./ cosh(a * (x - x0)).^2;
end

function g = G_rhs(t_now, v_hat_now, lam, k, mask)
    % RHS in integrating factor variables (cumulative form)
    E = exp(lam * t_now);
    u_hat = E .* v_hat_now;
    u  = real(ifft(u_hat));
    ux = real(ifft(1i * k .* u_hat));

    % Nonlinear term: -6 * u * u_x
    f_nl = -6.0 * u .* ux;
    f_nl_hat = fft(f_nl);
    f_nl_hat(~mask) = 0;

    g = exp(-lam * t_now) .* f_nl_hat;
end
