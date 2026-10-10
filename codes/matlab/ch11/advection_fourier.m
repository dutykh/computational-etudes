%% advection_fourier.m
%
% Chapter 11: Fourier Pseudospectral Methods for Periodic PDEs
%
% Fourier pseudospectral solver for the linear advection equation:
%   u_t + c * u_x = 0 on [0, 2*pi), periodic BCs.
%
% Compares Fourier pseudospectral + RK4 with second-order upwind FD + RK4
% to demonstrate the dispersion-free nature of spectral advection.
%
% Initial condition: u(x, 0) = exp(-10*(x - pi)^2)
%
% Generates figure: advection_comparison.pdf
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
% Part of "Computational Etudes: A Spectral Approach"
% https://github.com/dutykh/computational-etudes
% Date: February 2026

clear; close all; clc;

%% Configuration
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

%% Parameters
N = 128;
c = 1.0;
L = 2 * pi;
t_final = 2 * pi / c;  % One full period

% Grid
x  = L * (0:N-1)' / N;
dx = L / N;

% Wavenumbers (column vector, Nyquist zeroed for odd derivatives)
if mod(N, 2) == 0
    k = [0:N/2-1, 0, -N/2+1:-1]';
else
    k = [0:(N-1)/2, -(N-1)/2:-1]';
end

% Initial condition (smooth localised bump)
u0 = exp(-10.0 * (x - pi).^2);

% Time step (CFL-based)
k_max = N / 2;
dt = 0.5 / (c * k_max);  % CFL ~ 0.5
n_steps = ceil(t_final / dt);
dt = t_final / n_steps;

%% Solve with Fourier pseudospectral + RK4
u_fourier = u0;
for n = 1:n_steps
    u_fourier = rk4_step(u_fourier, dt, @(u) rhs_fourier(u, k, c));
end

%% Solve with second-order upwind FD + RK4
u_fd = u0;
for n = 1:n_steps
    u_fd = rk4_step(u_fd, dt, @(u) rhs_fd_upwind2(u, dx, c));
end

%% Exact solution (one full period = initial condition)
u_exact = u0;

%% Create figure
fig = figure('Units', 'inches', 'Position', [1, 1, 12, 4.5]);

% Left: overlay of solutions
subplot(1, 2, 1);
plot(x, u_exact, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 2.5, ...
     'DisplayName', 'Exact (initial)');
hold on;
plot(x, u_fourier, '-', 'Color', NAVY, 'LineWidth', 1.5, ...
     'DisplayName', 'Fourier PS + RK4');
plot(x, u_fd, '--', 'Color', CORAL, 'LineWidth', 1.5, ...
     'DisplayName', 'Upwind FD2 + RK4');

xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 11);
ylabel('$u(x, t)$', 'Interpreter', 'latex', 'FontSize', 11);
title(sprintf('After one period ($t = 2\\pi / c$, $N = %d$)', N), ...
      'Interpreter', 'latex', 'FontSize', 12);
legend('Interpreter', 'latex', 'FontSize', 10, 'Location', 'northeast');
xlim([0, L]);
grid on;
set(gca, 'GridAlpha', 0.3);
box on;

% Right: error profiles
subplot(1, 2, 2);
err_fourier = abs(u_fourier - u_exact);
err_fd      = abs(u_fd - u_exact);

semilogy(x, err_fd, '-', 'Color', CORAL, 'LineWidth', 1.5, ...
         'DisplayName', sprintf('Upwind FD2 (max: %.2e)', max(err_fd)));
hold on;
semilogy(x, max(err_fourier, 1e-16), '-', 'Color', NAVY, 'LineWidth', 1.5, ...
         'DisplayName', sprintf('Fourier PS (max: %.2e)', max(err_fourier)));

xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 11);
ylabel('$|u^{\mathrm{num}} - u^{\mathrm{exact}}|$', ...
       'Interpreter', 'latex', 'FontSize', 11);
title('Pointwise Error After One Period', 'FontSize', 12);
legend('Interpreter', 'latex', 'FontSize', 10, 'Location', 'south');
xlim([0, L]);
grid on;
set(gca, 'GridAlpha', 0.3);
box on;

%% Save figure
exportgraphics(fig, fullfile(output_dir, 'advection_comparison.pdf'), ...
               'ContentType', 'vector');
exportgraphics(fig, fullfile(output_dir, 'advection_comparison.png'), ...
               'Resolution', 300);
fprintf('Figure saved: %s\n', fullfile(output_dir, 'advection_comparison.pdf'));
close(fig);


%% -----------------------------------------------------------------------
%  Helper functions
%  -----------------------------------------------------------------------

function f = rhs_fourier(u, k, c)
    % Right-hand side for Fourier pseudospectral advection: -c * u_x
    u_hat = fft(u);
    ux = real(ifft(1i * k .* u_hat));
    f = -c * ux;
end

function f = rhs_fd_upwind2(u, dx, c)
    % Right-hand side for second-order upwind FD advection: -c * u_x
    N = length(u);
    ux = zeros(N, 1);
    if c > 0
        % Second-order upwind (backward biased)
        for j = 1:N
            jm1 = mod(j - 2, N) + 1;
            jm2 = mod(j - 3, N) + 1;
            ux(j) = (3.0 * u(j) - 4.0 * u(jm1) + u(jm2)) / (2.0 * dx);
        end
    else
        for j = 1:N
            jp1 = mod(j, N) + 1;
            jp2 = mod(j + 1, N) + 1;
            ux(j) = (-3.0 * u(j) + 4.0 * u(jp1) - u(jp2)) / (2.0 * dx);
        end
    end
    f = -c * ux;
end

function u_new = rk4_step(u, dt, rhs_func)
    % One step of classical RK4
    k1 = rhs_func(u);
    k2 = rhs_func(u + 0.5 * dt * k1);
    k3 = rhs_func(u + 0.5 * dt * k2);
    k4 = rhs_func(u + dt * k3);
    u_new = u + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
end
