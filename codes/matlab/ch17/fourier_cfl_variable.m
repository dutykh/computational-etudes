% fourier_cfl_variable.m
% Chapter 17: Stability of Spectral Methods
%
% Computational Etude 17.3 (variable-coefficient extension):
% Frozen-coefficient CFL prediction for variable-speed advection.
%
% Two-panel figure:
%   (a) Speed profile c(x) = 1 + 0.5*sin(x) on [0, 2*pi) with c_max
%       and c_mean annotated.
%   (b) Log-log plot of critical dt vs N for RK4 applied to Fourier
%       advection u_t + c(x)*u_x = 0. Experimental binary search
%       compared with frozen-coefficient prediction 2.83/(c_max*N/2).
%
% Author: Dr. Denys Dutykh
% Part of "Computational Etudes: A Spectral Approach"

clear; clc; close all;

NAVY  = [20, 45, 110]/255;
SKY   = [120, 150, 210]/255;
CORAL = [231, 76, 60]/255;

% Output path
script_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch17', 'matlab');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

%% ---- Panel (a): speed profile ----
c_func = @(x) 1.0 + 0.5 * sin(x);
c_max  = 1.5;
c_mean = 1.0;

x_plot = linspace(0, 2*pi, 500);
c_vals = c_func(x_plot);

figure('Position', [100, 100, 1100, 450]);

subplot(1, 2, 1);
plot(x_plot, c_vals, 'Color', NAVY, 'LineWidth', 2.0);
hold on;
yline(c_max, '--', 'Color', CORAL, 'LineWidth', 1.2);
yline(c_mean, ':', 'Color', SKY, 'LineWidth', 1.2);
hold off;

xlabel('$x$', 'Interpreter', 'latex');
ylabel('$c(x)$', 'Interpreter', 'latex');
title('(a) Variable speed profile');
xlim([0, 2*pi]); ylim([0.3, 1.7]);
legend({'$c(x) = 1 + 0.5\,\sin(x)$', '$c_{\max} = 1.5$', ...
        '$\langle c \rangle = 1.0$'}, ...
       'Interpreter', 'latex', 'Location', 'southwest', 'FontSize', 9);
grid on; set(gca, 'GridAlpha', 0.3);

%% ---- Panel (b): dt_crit vs N ----
N_values = [16, 32, 64, 128, 256];
dt_crits = zeros(size(N_values));

for idx = 1:length(N_values)
    Nv = N_values(idx);
    dtc = find_dt_crit_variable(Nv, c_func);
    dt_crits(idx) = dtc;
    fprintf('N = %4d, dt_crit = %.6f, dt_crit * c_max * (N/2) = %.4f\n', ...
            Nv, dtc, dtc * c_max * Nv / 2);
end

subplot(1, 2, 2);
loglog(N_values, dt_crits, 'o-', 'Color', NAVY, 'MarkerSize', 7, ...
       'LineWidth', 1.5, 'MarkerFaceColor', NAVY, ...
       'DisplayName', '$\Delta t_{\rm crit}$ (experiment)');
hold on;
set(gca, 'XScale', 'log', 'YScale', 'log');

% Frozen-coefficient prediction
N_ref = linspace(12, 300, 100);
loglog(N_ref, 2.83 ./ (c_max * N_ref / 2), '--', 'Color', CORAL, ...
       'LineWidth', 1.2, 'DisplayName', '$2.83\, /\, (c_{\max} \cdot N/2)$');

% Constant-coefficient reference
loglog(N_ref, 2.83 ./ (N_ref / 2), ':', 'Color', SKY, 'LineWidth', 1.2, ...
       'DisplayName', '$2.83\, /\, (N/2)$ [const.\ coeff.]');

% Slope annotation
text(80, 0.04, 'slope = -1', 'FontSize', 10, 'Color', CORAL, ...
     'Rotation', -28, 'HorizontalAlignment', 'center');
hold off;

xlabel('$N$', 'Interpreter', 'latex');
ylabel('$\Delta t_{\rm crit}$', 'Interpreter', 'latex');
title('(b) Frozen-coefficient prediction vs experiment');
legend('Interpreter', 'latex', 'Location', 'northeast', 'FontSize', 9);
grid on; set(gca, 'GridAlpha', 0.3);

% Save figure
set(gcf, 'PaperPositionMode', 'auto');
exportgraphics(gcf, fullfile(output_dir, 'fourier_cfl_variable.pdf'), ...
               'ContentType', 'vector');
exportgraphics(gcf, fullfile(output_dir, 'fourier_cfl_variable.png'), ...
               'Resolution', 300);
fprintf('\nFigure saved to %s\n', fullfile(output_dir, 'fourier_cfl_variable.pdf'));


%% ---- Helper functions ----

function du = fourier_deriv(u, k)
    % Spectral derivative via FFT
    du = real(ifft(1i * k .* fft(u)));
end

function rhs_val = variable_rhs(u, c, k)
    % RHS of u_t + c(x)*u_x = 0
    rhs_val = -c .* fourier_deriv(u, k);
end

function stable = is_stable_rk4(dt, u0, c, k, n_steps)
    % Run n_steps of RK4; return true if solution stays bounded
    u = u0;
    stable = true;
    for step = 1:n_steps
        k1 = variable_rhs(u, c, k);
        k2 = variable_rhs(u + 0.5*dt*k1, c, k);
        k3 = variable_rhs(u + 0.5*dt*k2, c, k);
        k4 = variable_rhs(u + dt*k3, c, k);
        u = u + dt/6 * (k1 + 2*k2 + 2*k3 + k4);
        if max(abs(u)) > 1e6
            stable = false;
            return;
        end
    end
end

function dtc = find_dt_crit_variable(N, c_func)
    % Binary search for critical dt of variable-coefficient advection
    tol = 1e-6;
    x = 2*pi*(0:N-1)'/N;
    c = c_func(x);
    c_max_val = max(abs(c));
    u0 = exp(cos(x));
    k = [0:N/2-1, -N/2:-1]';

    dt_lo = 0.1 / (c_max_val * N);
    dt_hi = 10.0 / N;

    % Ensure dt_hi is unstable
    while is_stable_rk4(dt_hi, u0, c, k, 500)
        dt_hi = dt_hi * 2;
        if dt_hi > 100, dtc = dt_hi; return; end
    end

    % Binary search
    for iter = 1:60
        dt_mid = 0.5 * (dt_lo + dt_hi);
        if is_stable_rk4(dt_mid, u0, c, k, 500)
            dt_lo = dt_mid;
        else
            dt_hi = dt_mid;
        end
        if dt_hi - dt_lo < tol * dt_lo, break; end
    end
    dtc = 0.5 * (dt_lo + dt_hi);
end
