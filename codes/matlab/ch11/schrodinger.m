%% schrodinger.m
%
% Time-dependent Schrodinger equation with harmonic potential (Strang splitting)
%
% PDE: i u_t = -u_xx + x^2 u,  0 < x < 2*pi,  t > 0 (periodic)
% ICs: u(x, 0) = exp(-5*(x - pi)^2) * exp(3ix)  (Gaussian wavepacket)
%
% Method: Strang splitting (second-order, unitary)
%   1. Half potential step (physical space): u -> u * exp(-i x^2 dt/2)
%   2. Full kinetic step (Fourier space):   u_hat -> u_hat * exp(-i k^2 dt)
%   3. Half potential step (physical space): u -> u * exp(-i x^2 dt/2)
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
% Part of "Computational Etudes: A Spectral Approach"
% https://github.com/dutykh/computational-etudes

clear; close all; clc;

%% Configuration
N = 256;          % Grid points
tmax = 3.0;       % Final time

% Colors
NAVY  = [0.078 0.176 0.431];
CORAL = [0.906 0.298 0.235];
TEAL  = [0.086 0.627 0.522];

% Output path
script_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch11', 'matlab');
output_file = fullfile(output_dir, 'schrodinger.pdf');

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% Setup grid
h = 2 * pi / N;
x = h * (0:N-1)';

% Wavenumbers
if mod(N, 2) == 0
    k = [0:N/2-1, 0, -N/2+1:-1]';
else
    k = [0:(N-1)/2, -(N-1)/2:-1]';
end

%% Potential: V(x) = x^2
V = x.^2;

%% Initial condition: Gaussian wavepacket with momentum
u = exp(-5.0 * (x - pi).^2) .* exp(3i * x);

%% Time stepping (Strang splitting)
dt = 0.5 * h;
nsteps = ceil(tmax / dt);
dt = tmax / nsteps;  % Adjust to hit tmax exactly

% Precompute exponentials for splitting operators
half_potential = exp(-0.5i * V * dt);    % Half potential step
full_kinetic  = exp(-1i * k.^2 * dt);   % Full kinetic step

% Storage for snapshots
n_save = 200;
save_every = max(1, floor(nsteps / n_save));
U_save = zeros(N, n_save + 1);
t_save = zeros(1, n_save + 1);
U_save(:, 1) = abs(u).^2;
t_save(1) = 0;
save_idx = 2;

% Record initial L^2 norm for unitarity check
norm0 = sqrt(h * sum(abs(u).^2));

t = 0;
for n = 1:nsteps
    % Strang splitting: potential/2 -> kinetic -> potential/2
    u = half_potential .* u;               % Half potential step
    u_hat = fft(u);
    u_hat = full_kinetic .* u_hat;         % Full kinetic step
    u = ifft(u_hat);
    u = half_potential .* u;               % Half potential step

    t = t + dt;

    % Save snapshot
    if mod(n, save_every) == 0 && save_idx <= n_save + 1
        U_save(:, save_idx) = abs(u).^2;
        t_save(save_idx) = t;
        save_idx = save_idx + 1;
    end
end

% Trim unused storage
U_save = U_save(:, 1:save_idx-1);
t_save = t_save(1:save_idx-1);

% Check unitarity preservation
norm_final = sqrt(h * sum(abs(u).^2));
fprintf('Initial L2 norm: %.10f\n', norm0);
fprintf('Final   L2 norm: %.10f\n', norm_final);
fprintf('Relative change: %.2e\n', abs(norm_final - norm0) / norm0);

%% Create figure
fig = figure('Units', 'inches', 'Position', [1, 1, 12, 6]);

% Left panel: space-time plot of |u|^2
subplot(1, 2, 1);
[X, T] = meshgrid(x, t_save);
pcolor(X, T, U_save');
shading interp;
colormap(subplot(1,2,1), 'hot');
xlabel('x', 'FontSize', 11);
ylabel('t', 'FontSize', 11);
title('Probability Density |\psi(x,t)|^2', 'FontSize', 12);
cb = colorbar;
cb.Label.String = '|\psi|^2';
cb.Label.FontSize = 11;

% Right panel: waterfall view
subplot(1, 2, 2);
hold on;

n_lines = 40;
indices = round(linspace(1, save_idx - 1, n_lines));
cmap = parula(n_lines);

for i = 1:n_lines
    idx = indices(i);
    offset = t_save(idx) * 0.3;
    plot(x, U_save(:, idx) + offset, '-', 'Color', cmap(i, :), 'LineWidth', 0.8);
end

xlabel('x', 'FontSize', 11);
ylabel('|\psi(x,t)|^2 (offset by t)', 'FontSize', 11);
title('Waterfall View: |\psi|^2 Evolution', 'FontSize', 12);
xlim([0, 2*pi]);
box on;

% Add inset for potential profile
axes('Position', [0.75, 0.75, 0.15, 0.15]);
plot(x, V, '-', 'Color', CORAL, 'LineWidth', 1.5);
xlabel('x', 'FontSize', 8);
ylabel('V(x)', 'FontSize', 8);
title('Harmonic Potential', 'FontSize', 9);
xlim([0, 2*pi]);
box on;
set(gca, 'FontSize', 7);

%% Save figure
exportgraphics(fig, output_file, 'ContentType', 'vector');
exportgraphics(fig, strrep(output_file, '.pdf', '.png'), 'Resolution', 300);

fprintf('Figure saved to: %s\n', output_file);

close(fig);
