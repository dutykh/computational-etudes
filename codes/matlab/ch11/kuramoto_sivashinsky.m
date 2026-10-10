%% kuramoto_sivashinsky.m
%
% Chapter 11: Fourier Pseudospectral Methods for Periodic PDEs
%
% Fourier pseudospectral solver for the Kuramoto-Sivashinsky equation:
%   u_t + u u_x + u_xx + u_xxxx = 0 on [0, L), periodic.
%
% Uses per-step integrating factor RK4 (Trefethen SMIM style) to handle
% the stiff linear terms exactly. The linear symbol is L(k) = k^2 - k^4
% (energy production at low k, damping at high k).
%
% Per-step IF-RK4 approach:
%   E  = exp(L*dt/2)   (half-step propagator)
%   E2 = E.^2          (full-step propagator)
%   Na = dt * Nhat(u_hat)
%   Nb = dt * Nhat(E .* u_hat + Na/2)
%   Nc = dt * Nhat(E .* u_hat + Nb/2)
%   Nd = dt * Nhat(E2 .* u_hat + E .* Nc)
%   u_hat = E2 .* u_hat + (E2.*Na + 2*E.*(Nb+Nc) + Nd)/6
%
% Generates two figures:
%   - ks_spacetime.pdf: space-time "fingerprint" plot
%   - ks_spectrum.pdf: time-averaged energy spectrum
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
% Part of "Computational Etudes: A Spectral Approach"
% https://github.com/dutykh/computational-etudes
% Date: February 2026

clear; close all; clc;

%% Configuration
N = 256;
L = 32 * pi;
t_final = 200.0;
dt = 0.05;

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

% Linear symbol: L(k) = k^2 - k^4
lam = k.^2 - k.^4;

% 2/3 dealiasing mask
cutoff = floor(N / 3);
mask = true(N, 1);
if cutoff > 0
    mask(cutoff+2:N-cutoff) = false;
end

%% Initial condition: small random perturbation (low-pass filtered)
rng(42);
u0 = 0.01 * randn(N, 1);

% Low-pass filter
u_hat = fft(u0);
u_hat(abs(k) > 5.0 * (2 * pi / L)) = 0;
u0 = real(ifft(u_hat));

%% Per-step integrating factor (Trefethen style)
u_hat = fft(u0);
n_steps = ceil(t_final / dt);
dt = t_final / n_steps;

E  = exp(lam * dt / 2.0);   % half-step propagator
E2 = E.^2;                   % full-step propagator

%% Storage
n_save = 500;
save_every = max(1, floor(n_steps / n_save));
U_save = u0;
t_save = 0;

%% Time stepping: per-step IF-RK4
t = 0;
for n = 1:n_steps
    Na = dt * Nhat(u_hat, k, mask);
    Nb = dt * Nhat(E .* u_hat + Na / 2.0, k, mask);
    Nc = dt * Nhat(E .* u_hat + Nb / 2.0, k, mask);
    Nd = dt * Nhat(E2 .* u_hat + E .* Nc, k, mask);
    u_hat = E2 .* u_hat + (E2 .* Na + 2.0 * E .* (Nb + Nc) + Nd) / 6.0;
    t = t + dt;

    if mod(n, save_every) == 0
        u = real(ifft(u_hat));
        U_save = [U_save, u]; %#ok<AGROW>
        t_save = [t_save, t]; %#ok<AGROW>
    end
end

%% Figure 1: Space-time "fingerprint"
fig1 = figure('Units', 'inches', 'Position', [1, 1, 8, 6]);

[X, T] = meshgrid(x, t_save);
n_quarter = max(1, floor(size(U_save, 2) / 4));
vmax = prctile(abs(U_save(:, n_quarter:end)), 98, 'all');
pcolor(X, T, U_save');
shading interp;
colormap(redblue(256));
caxis([-vmax, vmax]);
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 11);
ylabel('$t$', 'Interpreter', 'latex', 'FontSize', 11);
title(sprintf('Kuramoto--Sivashinsky: Spatiotemporal Chaos ($L = %d\\pi$)', ...
      round(L / pi)), 'Interpreter', 'latex', 'FontSize', 12);
cb = colorbar;
cb.Label.String = '$u(x, t)$';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 11;

exportgraphics(fig1, fullfile(output_dir, 'ks_spacetime.pdf'), ...
               'ContentType', 'vector');
exportgraphics(fig1, fullfile(output_dir, 'ks_spacetime.png'), ...
               'Resolution', 300);
fprintf('Figure saved: %s\n', fullfile(output_dir, 'ks_spacetime.pdf'));
close(fig1);

%% Figure 2: Time-averaged energy spectrum
% Average over the second half (statistically steady state)
n_start  = floor(size(U_save, 2) / 2);
n_total  = size(U_save, 2);
spectra  = zeros(N/2, 1);

for i = n_start:n_total
    uh = fft(U_save(:, i));
    spectra = spectra + abs(uh(1:N/2)).^2 / N^2;
end
spectra = spectra / (n_total - n_start + 1);

k_pos = (2 * pi / L) * (0:N/2-1)';

fig2 = figure('Units', 'inches', 'Position', [1, 1, 6, 4.5]);
semilogy(k_pos(2:end), spectra(2:end), '-', 'Color', NAVY, 'LineWidth', 1.5);
xlabel('Wavenumber $k$', 'Interpreter', 'latex', 'FontSize', 11);
ylabel('$\langle |\hat{u}_k|^2 \rangle$', 'Interpreter', 'latex', 'FontSize', 11);
title('Time-Averaged Energy Spectrum (KS)', 'FontSize', 12);
grid on;
set(gca, 'GridAlpha', 0.3);
box on;

exportgraphics(fig2, fullfile(output_dir, 'ks_spectrum.pdf'), ...
               'ContentType', 'vector');
exportgraphics(fig2, fullfile(output_dir, 'ks_spectrum.png'), ...
               'Resolution', 300);
fprintf('Figure saved: %s\n', fullfile(output_dir, 'ks_spectrum.pdf'));
close(fig2);


%% -----------------------------------------------------------------------
%  Helper functions
%  -----------------------------------------------------------------------

function f_hat = Nhat(u_hat_now, k, mask)
    % Nonlinear term in Fourier space: -ik/2 * FFT(u^2)
    u = real(ifft(u_hat_now));
    f_hat = -0.5i * k .* fft(u.^2);
    f_hat(~mask) = 0;
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
