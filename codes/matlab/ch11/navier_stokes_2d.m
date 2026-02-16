%% navier_stokes_2d.m
%
% Chapter 11: Fourier Pseudospectral Methods for Periodic PDEs
%
% Fourier pseudospectral solver for the 2D incompressible Navier-Stokes
% equations in vorticity-streamfunction form on [0, 2*pi)^2, periodic:
%
%   omega_t + J(psi, omega) = nu * Delta omega
%   Delta psi = -omega
%
% where J(psi, omega) = psi_x omega_y - psi_y omega_x is the Jacobian.
%
% Uses RK4 with 2/3 dealiasing. The Poisson equation for the streamfunction
% is solved trivially in Fourier space: psi_hat = omega_hat / (kx^2 + ky^2).
%
% Generates two figures:
%   - ns2d_vorticity.pdf: vorticity snapshots at several times
%   - ns2d_energy.pdf: energy and enstrophy decay curves
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
% Part of "Computational Etudes: A Spectral Approach"
% https://github.com/dutykh/computational-etudes
% Date: February 2026

clear; close all; clc;

%% Configuration
N = 128;
nu = 1e-3;
t_final = 20.0;
dt = 0.005;

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
Lx = 2 * pi;
x = Lx * (0:N-1)' / N;
y = x;

% 1D wavenumbers
if mod(N, 2) == 0
    kx_1d = [0:N/2-1, 0, -N/2+1:-1]';
else
    kx_1d = [0:(N-1)/2, -(N-1)/2:-1]';
end
ky_1d = kx_1d;

% 2D wavenumber grids (using ndgrid for 'ij' indexing)
[KX, KY] = ndgrid(kx_1d, ky_1d);
K2 = KX.^2 + KY.^2;
K2(1, 1) = 1.0;  % Avoid division by zero for mean mode

% 2/3 dealiasing mask (2D)
cutoff = floor(N / 3);
mask_1d = true(N, 1);
if cutoff > 0
    mask_1d(cutoff+2:N-cutoff) = false;
end
[MX, MY] = ndgrid(mask_1d, mask_1d);
mask = MX & MY;

% Diffusive symbol
diff_symbol = -nu * K2;

%% Initial condition: random vorticity with zero mean
rng(42);
omega_hat = zeros(N, N);

for i = 1:N
    for j = 1:N
        kk = sqrt(KX(i,j)^2 + KY(i,j)^2);
        if kk > 1.0 && kk < 10.0
            amp = 1.0 / (1.0 + kk);
            phase = 2 * pi * rand();
            omega_hat(i, j) = amp * exp(1i * phase);
        end
    end
end

% Enforce real-valued omega
omega = real(ifft2(omega_hat));
omega = omega - mean(omega(:));  % Ensure zero mean
omega = omega / std(omega(:));   % Normalize to unit RMS vorticity

%% Time stepping setup
n_steps = ceil(t_final / dt);
dt = t_final / n_steps;

n_save = 200;
save_every = max(1, floor(n_steps / n_save));

% Storage
omega_save = cell(n_save + 1, 1);
omega_save{1} = omega;
t_save = 0;
save_idx = 2;

[E0, Z0] = compute_diagnostics(omega, KX, KY, K2);
energy_save = E0;
enstrophy_save = Z0;

%% RK4 time stepping
t = 0;
for n = 1:n_steps
    k1 = ns_rhs(omega,                   KX, KY, K2, diff_symbol, mask);
    k2 = ns_rhs(omega + 0.5 * dt * k1,   KX, KY, K2, diff_symbol, mask);
    k3 = ns_rhs(omega + 0.5 * dt * k2,   KX, KY, K2, diff_symbol, mask);
    k4 = ns_rhs(omega + dt * k3,          KX, KY, K2, diff_symbol, mask);
    omega = omega + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
    t = t + dt;

    if mod(n, save_every) == 0 && save_idx <= n_save + 1
        omega_save{save_idx} = omega;
        t_save = [t_save, t]; %#ok<AGROW>
        [E, Z] = compute_diagnostics(omega, KX, KY, K2);
        energy_save    = [energy_save, E];    %#ok<AGROW>
        enstrophy_save = [enstrophy_save, Z]; %#ok<AGROW>
        save_idx = save_idx + 1;
    end
end

% Trim storage
n_saved = save_idx - 1;
omega_save = omega_save(1:n_saved);

%% Figure 1: Vorticity snapshots (2x2)
fig1 = figure('Units', 'inches', 'Position', [1, 1, 10, 10]);

[XX, YY] = ndgrid(x, y);
t_targets = [0.0, 5.0, 10.0, t_final];
vmax = prctile(abs(omega_save{1}(:)), 95);

for i = 1:4
    subplot(2, 2, i);

    [~, idx] = min(abs(t_save - t_targets(i)));
    omega_snap = omega_save{idx};

    pcolor(XX, YY, omega_snap);
    shading interp;
    colormap(redblue(256));
    caxis([-vmax, vmax]);
    xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 11);
    ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 11);
    title(sprintf('$t = %g$', t_targets(i)), ...
          'Interpreter', 'latex', 'FontSize', 12);
    axis equal;
    xlim([0, Lx]);
    ylim([0, Lx]);
    colorbar;
end

sgtitle(sprintf('2D Navier--Stokes: Decaying Turbulence ($N = %d$, $\\nu = %g$)', ...
        N, nu), 'Interpreter', 'latex', 'FontSize', 13);

exportgraphics(fig1, fullfile(output_dir, 'ns2d_vorticity.pdf'), ...
               'ContentType', 'vector');
exportgraphics(fig1, fullfile(output_dir, 'ns2d_vorticity.png'), ...
               'Resolution', 300);
fprintf('Figure saved: %s\n', fullfile(output_dir, 'ns2d_vorticity.pdf'));
close(fig1);

%% Figure 2: Energy and enstrophy decay
fig2 = figure('Units', 'inches', 'Position', [1, 1, 12, 4.5]);

t_diag = t_save(1:length(energy_save));

subplot(1, 2, 1);
plot(t_diag, energy_save, '-', 'Color', NAVY, 'LineWidth', 1.5);
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 11);
ylabel('$E(t) = \frac{1}{2}\langle |\mathbf{u}|^2 \rangle$', ...
       'Interpreter', 'latex', 'FontSize', 11);
title('Kinetic Energy Decay', 'FontSize', 12);
grid on;
set(gca, 'GridAlpha', 0.3);
box on;

subplot(1, 2, 2);
plot(t_diag, enstrophy_save, '-', 'Color', CORAL, 'LineWidth', 1.5);
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 11);
ylabel('$Z(t) = \frac{1}{2}\langle \omega^2 \rangle$', ...
       'Interpreter', 'latex', 'FontSize', 11);
title('Enstrophy Decay', 'FontSize', 12);
grid on;
set(gca, 'GridAlpha', 0.3);
box on;

exportgraphics(fig2, fullfile(output_dir, 'ns2d_energy.pdf'), ...
               'ContentType', 'vector');
exportgraphics(fig2, fullfile(output_dir, 'ns2d_energy.png'), ...
               'Resolution', 300);
fprintf('Figure saved: %s\n', fullfile(output_dir, 'ns2d_energy.pdf'));
close(fig2);


%% -----------------------------------------------------------------------
%  Helper functions
%  -----------------------------------------------------------------------

function f = ns_rhs(omega_phys, KX, KY, K2, diff_symbol, mask)
    % RHS: -J(psi, omega) + nu * laplacian(omega)
    omega_hat = fft2(omega_phys);

    % Poisson solve: psi_hat = omega_hat / K2
    psi_hat = omega_hat ./ K2;
    psi_hat(1, 1) = 0;  % Zero mean streamfunction

    % Velocity: u = psi_y, v = -psi_x
    u = real(ifft2(1i * KY .* psi_hat));
    v = real(ifft2(-1i * KX .* psi_hat));

    % Vorticity gradients
    omega_x = real(ifft2(1i * KX .* omega_hat));
    omega_y = real(ifft2(1i * KY .* omega_hat));

    % Jacobian (advection) in physical space
    jac = u .* omega_x + v .* omega_y;

    % Transform, dealias, add diffusion
    jac_hat = fft2(jac);
    jac_hat(~mask) = 0;

    f_hat = -jac_hat + diff_symbol .* omega_hat;
    f = real(ifft2(f_hat));
end

function [E, Z] = compute_diagnostics(omega_phys, KX, KY, K2)
    % Compute kinetic energy and enstrophy
    omega_hat = fft2(omega_phys);
    psi_hat = omega_hat ./ K2;
    psi_hat(1, 1) = 0;

    u = real(ifft2(1i * KY .* psi_hat));
    v = real(ifft2(-1i * KX .* psi_hat));

    E = 0.5 * mean(u(:).^2 + v(:).^2);
    Z = 0.5 * mean(omega_phys(:).^2);
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
