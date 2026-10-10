% catastrophe_gallery.m
% Chapter 17: Stability of Spectral Methods
%
% Computational Etude 17.1: Gallery of Deliberate Instabilities.
%
% Three-panel figure showing what happens when CFL / stability limits are
% deliberately violated in spectral discretisations:
%   (a) Fourier advection u_t + u_x = 0 with forward Euler, dt 20% above CFL.
%   (b) Chebyshev wave equation u_tt = u_xx (Dirichlet) with leapfrog,
%       dt 50% above CFL.
%   (c) Chebyshev heat equation u_t = u_xx (Dirichlet) with forward Euler,
%       dt doubled above the stability limit.
%
% Author: Dr. Denys Dutykh
% Part of "Computational Etudes: A Spectral Approach"

clear; clc; close all;

NAVY  = [20, 45, 110]/255;
SKY   = [120, 150, 210]/255;
CORAL = [231, 76, 60]/255;
TEAL  = [22, 160, 133]/255;
AMBER = [230, 126, 34]/255;

% Output path
script_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch17', 'matlab');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

%% ---- (a) Fourier advection instability ----
N = 64;
L = 2*pi;
x_f = linspace(0, L, N+1); x_f(end) = [];  % endpoint=false
k = [0:N/2-1, -N/2:-1];                      % fftfreq equivalent
k = 2*pi*k/L;

u0 = exp(-10*(x_f - pi).^2);
dt_stable = 0.5 * L / (N * 1.0);   % reference CFL ~ dx/c
dt = 1.2 * dt_stable;               % 20% above

snap_times = [0, 50, 100, 150];
ik = 1i * k;

% Stable run with RK4
snapshots_stable = cell(1, length(snap_times));
u_hat = fft(u0);
snap_idx = 1;
for n = 0:max(snap_times)
    if ismember(n, snap_times)
        snapshots_stable{snap_idx} = real(ifft(u_hat));
        snap_idx = snap_idx + 1;
    end
    % RK4 step
    k1 = -ik .* u_hat;
    k2 = -ik .* (u_hat + 0.5*dt*k1);
    k3 = -ik .* (u_hat + 0.5*dt*k2);
    k4 = -ik .* (u_hat + dt*k3);
    u_hat = u_hat + dt/6 * (k1 + 2*k2 + 2*k3 + k4);
end

% Unstable run with forward Euler
snapshots_unstable = cell(1, length(snap_times));
u_hat = fft(u0);
snap_idx = 1;
for n = 0:max(snap_times)
    if ismember(n, snap_times)
        snapshots_unstable{snap_idx} = real(ifft(u_hat));
        snap_idx = snap_idx + 1;
    end
    u_hat = u_hat + dt * (-ik .* u_hat);
end

%% ---- (b) Chebyshev wave equation instability ----
Nw = 40;
[D, x_w] = cheb(Nw);
D2 = D * D;
D2_int = D2(2:Nw, 2:Nw);   % interior points (Dirichlet BCs)
x_int_w = x_w(2:Nw);

% Critical dt for leapfrog on wave equation
eigs_w = eig(D2_int);
rho = max(abs(eigs_w));
dt_crit_w = 2.0 / sqrt(rho);
dt_w = 1.5 * dt_crit_w;     % 50% above

% Initial condition: smooth bump
u_w = exp(-20 * x_int_w.^2);
u_prev = u_w;

n_steps_w = 60;
snap_times_w = [0, 20, 40, 60];

% Store initial snapshot
snapshots_wave = cell(1, length(snap_times_w));
snap_full = zeros(Nw+1, 1);
snap_full(2:Nw) = u_w;
snap_idx = 1;
snapshots_wave{snap_idx} = snap_full; snap_idx = snap_idx + 1;

for n = 1:n_steps_w
    if n == 1
        % First step: use forward Euler for u^1
        u_new = u_w + 0.5 * dt_w^2 * (D2_int * u_w);
    else
        u_new = 2*u_w - u_prev + dt_w^2 * (D2_int * u_w);
    end
    u_prev = u_w;
    u_w = u_new;
    if ismember(n, snap_times_w)
        snap_full = zeros(Nw+1, 1);
        snap_full(2:Nw) = u_w;
        snapshots_wave{snap_idx} = snap_full;
        snap_idx = snap_idx + 1;
    end
end

%% ---- (c) Chebyshev heat equation instability ----
Nh = 20;
[D, x_h] = cheb(Nh);
D2 = D * D;
D2_int_h = D2(2:Nh, 2:Nh);
x_int_h = x_h(2:Nh);

eigs_h = eig(D2_int_h);
rho_h = max(abs(real(eigs_h)));
dt_crit_h = 2.0 / rho_h;
dt_h = 2.0 * dt_crit_h;   % doubled

u_h = sin(pi * x_int_h);
snap_times_h = [0, 1, 2, 3];
snapshots_heat = cell(1, length(snap_times_h));
snap_idx = 1;

for n = 0:max(snap_times_h)
    if ismember(n, snap_times_h)
        snap_full = zeros(Nh+1, 1);
        snap_full(2:Nh) = u_h;
        snapshots_heat{snap_idx} = snap_full;
        snap_idx = snap_idx + 1;
    end
    u_h = u_h + dt_h * (D2_int_h * u_h);
end

%% ---- Plotting ----
figure('Position', [100, 100, 1400, 420]);

% Panel (a): Fourier advection
subplot(1, 3, 1);
alphas = [1.0, 0.5, 0.7, 1.0];
hold on;
plot(x_f, snapshots_stable{1}, 'Color', NAVY, 'LineWidth', 1.5, ...
     'DisplayName', 'initial');
for i = 2:length(snap_times)
    plot(x_f, snapshots_stable{i}, '--', 'Color', [SKY, alphas(i)], ...
         'LineWidth', 1.0, 'HandleVisibility', 'off');
    u_clip = max(min(snapshots_unstable{i}, 3), -3);
    plot(x_f, u_clip, 'Color', [CORAL, alphas(i)], 'LineWidth', 1.0, ...
         'DisplayName', sprintf('FE, n=%d', snap_times(i)));
end
hold off;
xlabel('x'); ylabel('u(x, t)');
title('(a) Fourier advection, forward Euler');
ylim([-3, 3]);
legend('Location', 'northeast', 'FontSize', 8);
grid on; set(gca, 'GridAlpha', 0.3);

% Panel (b): Chebyshev wave equation
subplot(1, 3, 2);
colors_w = {NAVY, SKY, AMBER, CORAL};
t_labels = [0, snap_times_w(2:end)];
hold on;
for i = 1:length(snapshots_wave)
    snap_clip = max(min(snapshots_wave{i}, 5), -5);
    plot(x_w, snap_clip, 'Color', colors_w{i}, 'LineWidth', 1.2, ...
         'DisplayName', sprintf('n = %d', t_labels(i)));
end
hold off;
xlabel('x'); ylabel('u(x, t)');
title('(b) Chebyshev wave, leapfrog');
ylim([-5, 5]);
legend('Location', 'northwest', 'FontSize', 8);
grid on; set(gca, 'GridAlpha', 0.3);

% Panel (c): Chebyshev heat equation
subplot(1, 3, 3);
colors_h = {NAVY, SKY, AMBER, CORAL};
hold on;
for i = 1:length(snapshots_heat)
    snap_clip = max(min(snapshots_heat{i}, 5), -5);
    plot(x_h, snap_clip, 'Color', colors_h{i}, 'LineWidth', 1.2, ...
         'DisplayName', sprintf('n = %d', snap_times_h(i)));
end
hold off;
xlabel('x'); ylabel('u(x, t)');
title('(c) Chebyshev heat, forward Euler');
ylim([-5, 5]);
legend('Location', 'northeast', 'FontSize', 8);
grid on; set(gca, 'GridAlpha', 0.3);

set(gcf, 'PaperPositionMode', 'auto');
exportgraphics(gcf, fullfile(output_dir, 'catastrophe_gallery.pdf'), 'ContentType', 'vector');
exportgraphics(gcf, fullfile(output_dir, 'catastrophe_gallery.png'), 'Resolution', 300);
fprintf('Figure saved to %s\n', fullfile(output_dir, 'catastrophe_gallery.pdf'));
