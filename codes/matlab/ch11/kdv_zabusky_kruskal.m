%% kdv_zabusky_kruskal.m
%
% Chapter 11: Fourier Pseudospectral Methods for Periodic PDEs
%
% Fourier pseudospectral solver for the KdV equation (Zabusky-Kruskal form):
%   u_t + u u_x + delta^2 u_xxx = 0 on [0, 2), periodic.
%
% with u(x, 0) = cos(pi * x) and delta^2 = 0.022^2.
%
% Replicates the Zabusky-Kruskal (1965) experiment: a cosine initial condition
% breaks into a train of solitons.
%
% Uses integrating factor RK4 (cumulative form) with 2/3 dealiasing.
%
% Generates figure: kdv_zabusky_kruskal.pdf (space-time diagram + snapshots)
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
% Part of "Computational Etudes: A Spectral Approach"
% https://github.com/dutykh/computational-etudes
% Date: February 2026

clear; close all; clc;

%% Configuration
N  = 256;
L  = 2.0;
t_final = 3.6 / pi;    % Standard ZK experiment time
dt = 5e-5;
delta2 = 0.022^2;

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

% Linear symbol: lambda_k = -i * delta^2 * k^3
lam = -1i * delta2 * k.^3;

% 2/3 dealiasing mask
cutoff = floor(N / 3);
mask = true(N, 1);
if cutoff > 0
    mask(cutoff+2:N-cutoff) = false;
end

%% Initial condition
u0 = cos(pi * x);
v_hat = fft(u0);

%% Time stepping: integrating factor RK4 (cumulative form)
n_steps = ceil(t_final / dt);
dt = t_final / n_steps;

n_save = 300;
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

%% Create figure
fig = figure('Units', 'inches', 'Position', [1, 1, 12, 5]);

% Left panel: space-time plot
subplot(1, 2, 1);
[X, T] = meshgrid(x, t_save);
pcolor(X, T, U_save');
shading interp;
colormap(subplot(1,2,1), redblue(256));
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 11);
ylabel('$t$', 'Interpreter', 'latex', 'FontSize', 11);
title('Zabusky--Kruskal Experiment', 'FontSize', 12);
cb = colorbar;
cb.Label.String = '$u(x, t)$';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 11;

% Right panel: snapshots at selected times
subplot(1, 2, 2);
hold on;

t_targets = [0.0, 0.3, 0.6, 0.9, t_final];
colors = [SKY; TEAL; 0.180 0.800 0.443; 0.953 0.612 0.071; NAVY];

for i = 1:length(t_targets)
    [~, idx] = min(abs(t_save - t_targets(i)));
    plot(x, U_save(:, idx), '-', 'Color', colors(i, :), 'LineWidth', 1.3, ...
         'DisplayName', sprintf('$t = %.2f$', t_save(idx)));
end

xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 11);
ylabel('$u(x, t)$', 'Interpreter', 'latex', 'FontSize', 11);
title('Soliton Train Formation', 'FontSize', 12);
lgd = legend('Interpreter', 'latex', 'FontSize', 9, 'Location', 'northeast');
lgd.BoxFace.ColorType = 'truecoloralpha';
lgd.BoxFace.ColorData = uint8([255; 255; 255; 230]);
xlim([0, L]);
grid on;
set(gca, 'GridAlpha', 0.3);
box on;

%% Save figure
exportgraphics(fig, fullfile(output_dir, 'kdv_zabusky_kruskal.pdf'), ...
               'ContentType', 'vector');
exportgraphics(fig, fullfile(output_dir, 'kdv_zabusky_kruskal.png'), ...
               'Resolution', 300);
fprintf('Figure saved: %s\n', fullfile(output_dir, 'kdv_zabusky_kruskal.pdf'));
close(fig);


%% -----------------------------------------------------------------------
%  Helper functions
%  -----------------------------------------------------------------------

function g = G_rhs(t_now, v_hat_now, lam, k, mask)
    % RHS in integrating factor variables (cumulative form)
    E = exp(lam * t_now);
    u_hat = E .* v_hat_now;
    u  = real(ifft(u_hat));
    ux = real(ifft(1i * k .* u_hat));

    % Nonlinear term: -u * u_x (ZK normalisation)
    f_nl = -u .* ux;
    f_nl_hat = fft(f_nl);
    f_nl_hat(~mask) = 0;

    g = exp(-lam * t_now) .* f_nl_hat;
end

function cmap = redblue(m)
    if nargin < 1, m = 256; end
    n = ceil(m / 2);
    r1 = linspace(0.2, 1, n)';
    g1 = linspace(0.3, 1, n)';
    b1 = linspace(0.8, 1, n)';
    r2 = linspace(1, 0.8, n)';
    g2 = linspace(1, 0.2, n)';
    b2 = linspace(1, 0.2, n)';
    cmap = [r1 g1 b1; r2(2:end,:) g2(2:end,:) b2(2:end,:)];
    if size(cmap, 1) > m
        cmap = cmap(1:m, :);
    end
end
