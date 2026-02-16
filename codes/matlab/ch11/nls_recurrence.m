%% nls_recurrence.m
%
% Chapter 11: Fourier Pseudospectral Methods for Periodic PDEs
%
% Split-step Fourier solver for the focusing nonlinear Schrodinger equation:
%   i u_t + u_xx + 2 |u|^2 u = 0 on [0, 2*pi), periodic.
%
% Demonstrates modulation instability and approximate recurrence using
% Strang splitting (linear-nonlinear-linear half-steps).
%
% Initial condition: u(x, 0) = 1 + 0.1 * cos(x)
%
% Generates figure: nls_recurrence.pdf
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
% Part of "Computational Etudes: A Spectral Approach"
% https://github.com/dutykh/computational-etudes
% Date: February 2026

clear; close all; clc;

%% Configuration
N = 256;
t_final = 5.0;
dt = 5e-4;

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
x = 2 * pi * (0:N-1)' / N;

if mod(N, 2) == 0
    k = [0:N/2-1, 0, -N/2+1:-1]';
else
    k = [0:(N-1)/2, -(N-1)/2:-1]';
end

%% Initial condition: modulated plane wave
u = (1.0 + 0.1 * cos(x)) + 0i;  % ensure complex

%% Precompute linear propagator (half-step)
n_steps = ceil(t_final / dt);
dt = t_final / n_steps;
linear_half = exp(-1i * k.^2 * dt / 2.0);

%% Storage
n_save = 400;
save_every = max(1, floor(n_steps / n_save));
U_save  = abs(u);
L2_save = sqrt(mean(abs(u).^2));
t_save  = 0;

%% Time stepping: Strang splitting
t = 0;
for n = 1:n_steps
    % Half linear step
    u_hat = fft(u);
    u_hat = linear_half .* u_hat;
    u = ifft(u_hat);

    % Full nonlinear step: u -> u * exp(2i |u|^2 dt)
    u = u .* exp(2i * abs(u).^2 * dt);

    % Second half linear step
    u_hat = fft(u);
    u_hat = linear_half .* u_hat;
    u = ifft(u_hat);

    t = t + dt;

    if mod(n, save_every) == 0
        U_save  = [U_save, abs(u)];           %#ok<AGROW>
        L2_save = [L2_save, sqrt(mean(abs(u).^2))]; %#ok<AGROW>
        t_save  = [t_save, t];                 %#ok<AGROW>
    end
end

%% Create figure
fig = figure('Units', 'inches', 'Position', [1, 1, 12, 5]);

% Left panel: space-time diagram of |u|
ax1 = subplot(1, 2, 1);
[X, T] = meshgrid(x, t_save);
pcolor(X, T, U_save');
shading interp;
colormap(ax1, 'hot');
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 11);
ylabel('$t$', 'Interpreter', 'latex', 'FontSize', 11);
title('NLS: $|u(x, t)|$ (Modulation Instability)', ...
      'Interpreter', 'latex', 'FontSize', 12);
cb = colorbar;
cb.Label.String = '$|u|$';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 11;

% Right panel: L2 norm conservation
subplot(1, 2, 2);
L2_0 = L2_save(1);
semilogy(t_save, abs(L2_save - L2_0) / L2_0, '-', 'Color', NAVY, ...
         'LineWidth', 1.0);
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 11);
ylabel('$|\|u\|_2 - \|u_0\|_2| / \|u_0\|_2$', ...
       'Interpreter', 'latex', 'FontSize', 11);
title('$L^2$ Norm Conservation', 'Interpreter', 'latex', 'FontSize', 12);
grid on;
set(gca, 'GridAlpha', 0.3);
box on;

%% Save figure
exportgraphics(fig, fullfile(output_dir, 'nls_recurrence.pdf'), ...
               'ContentType', 'vector');
exportgraphics(fig, fullfile(output_dir, 'nls_recurrence.png'), ...
               'Resolution', 300);
fprintf('Figure saved: %s\n', fullfile(output_dir, 'nls_recurrence.pdf'));
close(fig);
