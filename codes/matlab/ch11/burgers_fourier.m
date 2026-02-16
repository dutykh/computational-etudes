%% burgers_fourier.m
%
% Chapter 11: Fourier Pseudospectral Methods for Periodic PDEs
%
% Fourier pseudospectral solver for the viscous Burgers equation:
%   u_t + u * u_x = nu * u_xx on [0, 2*pi), periodic BCs.
%
% Initial condition: u(x, 0) = sin(x)
% Features: RK4 time stepping, 2/3 dealiasing rule.
%
% Generates two figures:
%   - burgers_evolution.pdf: solution evolution showing front steepening
%   - burgers_dealiasing.pdf: comparison of energy spectra with/without dealiasing
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
% Part of "Computational Etudes: A Spectral Approach"
% https://github.com/dutykh/computational-etudes
% Date: February 2026

clear; close all; clc;

%% Configuration
N  = 256;
nu = 0.01;
t_final = 1.5;

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

%% Figure 1: Solution evolution (with dealiasing)
[x, t_save, U_save] = burgers_solve(N, nu, t_final, 0.4, true);

fig1 = figure('Units', 'inches', 'Position', [1, 1, 12, 5]);

% Left panel: space-time plot
subplot(1, 2, 1);
[X, T] = meshgrid(x, t_save);
pcolor(X, T, U_save');
shading interp;
colormap(subplot(1,2,1), redblue(256));
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 11);
ylabel('$t$', 'Interpreter', 'latex', 'FontSize', 11);
title(sprintf('Burgers Equation ($\\nu = %g$)', nu), ...
      'Interpreter', 'latex', 'FontSize', 12);
cb = colorbar;
cb.Label.String = '$u(x, t)$';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 11;

% Right panel: snapshots
subplot(1, 2, 2);
hold on;

snapshot_times = [0.0, 0.3, 0.6, 1.0, 1.5];
colors = [SKY; TEAL; 0.180 0.800 0.443; 0.953 0.612 0.071; NAVY];

for i = 1:length(snapshot_times)
    [~, idx] = min(abs(t_save - snapshot_times(i)));
    plot(x, U_save(:, idx), '-', 'Color', colors(i, :), 'LineWidth', 1.5, ...
         'DisplayName', sprintf('$t = %.1f$', t_save(idx)));
end

xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 11);
ylabel('$u(x, t)$', 'Interpreter', 'latex', 'FontSize', 11);
title('Snapshots: Front Steepening', 'FontSize', 12);
lgd = legend('Interpreter', 'latex', 'FontSize', 9, 'Location', 'northeast');
lgd.BoxFace.ColorType = 'truecoloralpha';
lgd.BoxFace.ColorData = uint8([255; 255; 255; 230]);
xlim([0, 2*pi]);
grid on;
set(gca, 'GridAlpha', 0.3);
box on;

exportgraphics(fig1, fullfile(output_dir, 'burgers_evolution.pdf'), ...
               'ContentType', 'vector');
exportgraphics(fig1, fullfile(output_dir, 'burgers_evolution.png'), ...
               'Resolution', 300);
fprintf('Figure saved: %s\n', fullfile(output_dir, 'burgers_evolution.pdf'));
close(fig1);

%% Figure 2: Dealiasing comparison
t_compare = 1.0;
[~, t1, U1] = burgers_solve(N, nu, t_compare, 0.4, true);
[~, t2, U2] = burgers_solve(N, nu, t_compare, 0.4, false);

% Energy spectra at final time
u_hat_deal   = fft(U1(:, end));
u_hat_nodeal = fft(U2(:, end));
spectrum_deal   = abs(u_hat_deal(1:N/2)).^2 / N^2;
spectrum_nodeal = abs(u_hat_nodeal(1:N/2)).^2 / N^2;
k_pos = (0:N/2-1)';

fig2 = figure('Units', 'inches', 'Position', [1, 1, 12, 4.5]);

% Left: energy spectra
subplot(1, 2, 1);
semilogy(k_pos(2:end), spectrum_deal(2:end), '-', 'Color', NAVY, ...
         'LineWidth', 1.5, 'DisplayName', 'With 2/3 dealiasing');
hold on;
semilogy(k_pos(2:end), spectrum_nodeal(2:end), '--', 'Color', CORAL, ...
         'LineWidth', 1.5, 'DisplayName', 'Without dealiasing');
xlabel('Wavenumber $k$', 'Interpreter', 'latex', 'FontSize', 11);
ylabel('$|\hat{u}_k|^2 / N^2$', 'Interpreter', 'latex', 'FontSize', 11);
title(sprintf('Energy Spectrum at $t = %g$', t_compare), ...
      'Interpreter', 'latex', 'FontSize', 12);
legend('Interpreter', 'latex', 'FontSize', 10);
grid on;
set(gca, 'GridAlpha', 0.3);
box on;

% Right: physical space comparison
subplot(1, 2, 2);
plot(x, U1(:, end), '-', 'Color', NAVY, 'LineWidth', 1.5, ...
     'DisplayName', 'With dealiasing');
hold on;
plot(x, U2(:, end), '--', 'Color', CORAL, 'LineWidth', 1.5, ...
     'DisplayName', 'Without dealiasing');
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 11);
ylabel('$u(x, t)$', 'Interpreter', 'latex', 'FontSize', 11);
title(sprintf('Physical Space at $t = %g$', t_compare), ...
      'Interpreter', 'latex', 'FontSize', 12);
legend('Interpreter', 'latex', 'FontSize', 10);
xlim([0, 2*pi]);
grid on;
set(gca, 'GridAlpha', 0.3);
box on;

exportgraphics(fig2, fullfile(output_dir, 'burgers_dealiasing.pdf'), ...
               'ContentType', 'vector');
exportgraphics(fig2, fullfile(output_dir, 'burgers_dealiasing.png'), ...
               'Resolution', 300);
fprintf('Figure saved: %s\n', fullfile(output_dir, 'burgers_dealiasing.pdf'));
close(fig2);


%% -----------------------------------------------------------------------
%  Solver function
%  -----------------------------------------------------------------------

function [x, t_save, U_save] = burgers_solve(N, nu, t_final, CFL, dealias)
    % Fourier pseudospectral solver for viscous Burgers equation
    %   u_t + u u_x = nu u_xx on [0, 2*pi), periodic BCs.
    %
    % Returns: x (column), t_save (row), U_save (N x n_saved)

    % Grid and wavenumbers
    x = 2 * pi * (0:N-1)' / N;

    if mod(N, 2) == 0
        k = [0:N/2-1, 0, -N/2+1:-1]';
    else
        k = [0:(N-1)/2, -(N-1)/2:-1]';
    end

    % Initial condition
    u = sin(x);

    % Time step from linearised stability
    k_max = N / 2;
    dt_adv  = CFL / k_max;
    if nu > 0
        dt_diff = 0.4 / (nu * k_max^2);
    else
        dt_diff = Inf;
    end
    dt = min(dt_adv, dt_diff);
    n_steps = ceil(t_final / dt);
    dt = t_final / n_steps;

    % Diffusive symbol
    diff_symbol = -nu * k.^2;

    % 2/3 dealiasing mask
    dealias_mask = true(N, 1);
    if dealias
        cutoff = floor(N / 3);
        if cutoff > 0
            dealias_mask(cutoff+2:N-cutoff) = false;
        end
    end

    % RHS function
    rhs = @(u_phys) burgers_rhs(u_phys, k, diff_symbol, dealias_mask);

    % Storage
    n_save = min(200, n_steps);
    save_every = max(1, floor(n_steps / n_save));
    U_save = u;
    t_save = 0;

    t = 0;
    for n = 1:n_steps
        % RK4 step
        k1 = rhs(u);
        k2 = rhs(u + 0.5 * dt * k1);
        k3 = rhs(u + 0.5 * dt * k2);
        k4 = rhs(u + dt * k3);
        u = u + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
        t = t + dt;

        if mod(n, save_every) == 0
            U_save = [U_save, u]; %#ok<AGROW>
            t_save = [t_save, t]; %#ok<AGROW>
        end
    end
end

function f = burgers_rhs(u_phys, k, diff_symbol, mask)
    u_hat = fft(u_phys);
    ux = real(ifft(1i * k .* u_hat));

    % Nonlinear term in physical space, then dealias
    f_nl = -u_phys .* ux;
    f_nl_hat = fft(f_nl);
    f_nl_hat(~mask) = 0;

    % Diffusive term in Fourier space
    f_diff_hat = diff_symbol .* u_hat;

    f = real(ifft(f_nl_hat + f_diff_hat));
end


%% Red-Blue diverging colormap
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
