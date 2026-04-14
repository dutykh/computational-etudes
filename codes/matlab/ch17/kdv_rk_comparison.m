% kdv_rk_comparison.m
% Chapter 17: Stability of Spectral Methods
%
% Computational Etude 17.7: Plain RK4 vs IF-RK4 on the KdV Equation.
%
% Two-panel figure:
%   (a) Maximum stable dt vs N (log-log) for plain RK4 and IF-RK4 on
%       KdV u_t + 6*u*u_x + u_xxx = 0.
%   (b) Solutions at t = 0.5 for N = 128 with both methods.
%
% Author: Dr. Denys Dutykh
% Part of "Computational Etudes: A Spectral Approach"

clear; clc; close all;

NAVY  = [20, 45, 110]/255;
SKY   = [120, 150, 210]/255;
CORAL = [231, 76, 60]/255;
AMBER = [230, 126, 34]/255;

% Output path
script_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch17', 'matlab');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

%% ---- Panel (a): dt_crit vs N ----
N_values = [32, 64, 128, 256];
dt_rk4   = zeros(size(N_values));
dt_ifrk4 = zeros(size(N_values));

for idx = 1:length(N_values)
    Nv = N_values(idx);

    % RK4: analytical linear stability limit
    dtc_rk4 = rk4_dt_crit_linear(Nv);
    fprintf('N = %4d: RK4 dt_crit = %.6e\n', Nv, dtc_rk4);
    dt_rk4(idx) = dtc_rk4;

    % IF-RK4: binary search
    dtc_if = find_ifrk4_dt_crit(Nv, 0.05);
    fprintf('         IF-RK4 dt_crit = %.6e\n', dtc_if);
    dt_ifrk4(idx) = dtc_if;
end

% Fit slopes
p_rk4 = polyfit(log(N_values), log(dt_rk4), 1);
slope_rk4 = p_rk4(1); c_rk4 = exp(p_rk4(2));
fprintf('\nRK4 slope: %.2f\n', slope_rk4);

p_if = polyfit(log(N_values), log(dt_ifrk4), 1);
slope_if = p_if(1); c_if = exp(p_if(2));
fprintf('IF-RK4 slope: %.2f\n', slope_if);

figure('Position', [100, 100, 1200, 500]);

subplot(1, 2, 1);
N_ref = linspace(20, 300, 100);
loglog(N_values, dt_rk4, 'o-', 'Color', NAVY, 'MarkerSize', 7, ...
       'LineWidth', 1.5, 'MarkerFaceColor', NAVY, 'DisplayName', 'Plain RK4');
hold on;
loglog(N_values, dt_ifrk4, 's-', 'Color', CORAL, 'MarkerSize', 7, ...
       'LineWidth', 1.5, 'MarkerFaceColor', CORAL, 'DisplayName', 'IF-RK4');

loglog(N_ref, c_rk4 * N_ref.^slope_rk4, '--', 'Color', SKY, ...
       'LineWidth', 1.0, 'DisplayName', sprintf('~ N^{%.1f}', slope_rk4));
loglog(N_ref, c_if * N_ref.^slope_if, ':', 'Color', AMBER, ...
       'LineWidth', 1.0, 'DisplayName', sprintf('~ N^{%.1f}', slope_if));
hold off;

xlabel('N'); ylabel('\Delta t_{crit}');
title('(a) Maximum stable \Delta t vs N');
legend('Location', 'southwest', 'FontSize', 8);
grid on; set(gca, 'GridAlpha', 0.3);

%% ---- Panel (b): Solutions at t = 0.5 for N = 128 ----
N_sol = 128;
T_sol = 0.5;

% Plain RK4 with safe dt
dt_rk4_safe = 0.8 * rk4_dt_crit_linear(N_sol);
fprintf('\nComputing RK4 solution at t=%.1f with N=%d, dt=%.6e...\n', ...
        T_sol, N_sol, dt_rk4_safe);
[u_rk4, x_rk4] = kdv_rk4_solve(N_sol, dt_rk4_safe, T_sol);

% IF-RK4 with moderate dt
dt_if_use = 1e-3;
fprintf('Computing IF-RK4 solution at t=%.1f with N=%d, dt=%.6e...\n', ...
        T_sol, N_sol, dt_if_use);
[u_if, x_if] = kdv_ifrk4_solve(N_sol, dt_if_use, T_sol);

subplot(1, 2, 2);
hold on;
if ~isempty(u_rk4)
    plot(x_rk4, u_rk4, '-', 'Color', NAVY, 'LineWidth', 1.8, ...
         'DisplayName', 'RK4');
else
    fprintf('  WARNING: RK4 solution blew up!\n');
end
if ~isempty(u_if)
    plot(x_if, u_if, '--', 'Color', CORAL, 'LineWidth', 1.5, ...
         'DisplayName', 'IF-RK4');
else
    fprintf('  WARNING: IF-RK4 solution blew up!\n');
end
hold off;

xlabel('x'); ylabel('u(x, t)');
title(sprintf('(b) KdV solution at t = %.1f, N = %d', T_sol, N_sol));
legend('Location', 'northeast', 'FontSize', 9);
grid on; set(gca, 'GridAlpha', 0.3);

set(gcf, 'PaperPositionMode', 'auto');
exportgraphics(gcf, fullfile(output_dir, 'kdv_rk_comparison.pdf'), 'ContentType', 'vector');
exportgraphics(gcf, fullfile(output_dir, 'kdv_rk_comparison.png'), 'Resolution', 300);
fprintf('\nFigure saved to %s\n', fullfile(output_dir, 'kdv_rk_comparison.pdf'));


%% ======================================================================
%  Helper functions
%  ======================================================================

function dtc = rk4_dt_crit_linear(N)
    % Analytical critical dt for RK4 on the linear part of KdV (u_xxx).
    % Eigenvalues are lambda_k = -i*k^3. RK4 stability on imaginary axis:
    % |y| <= 2*sqrt(2). Therefore dt_crit = 2*sqrt(2) / k_max^3.
    L = 2.0;
    kv = 2*pi * [0:N/2-1, -N/2:-1] / L;
    k_max = max(abs(kv));
    dtc = 2*sqrt(2) / k_max^3;
end

function [x, kv, dealias] = setup_kdv(N)
    L = 2.0;
    x = linspace(0, L, N+1); x(end) = [];
    kv = 2*pi * [0:N/2-1, -N/2:-1] / L;
    dealias = ones(1, N);
    kmax = floor(N/3);
    dealias(kmax+1 : N-kmax+1) = 0;
end

function [u, x] = kdv_ifrk4_solve(N, dt, T)
    % IF-RK4 for KdV: integrating factor removes the linear u_xxx term.
    [x, kv, dealias] = setup_kdv(N);
    u_hat = fft(cos(pi*x));

    ik  = 1i * kv;
    lam = -1i * kv.^3;

    n_steps = max(1, round(T / dt));
    dt_actual = T / n_steps;

    E  = exp(lam * dt_actual / 2);
    E2 = E.^2;

    for step = 1:n_steps
        % NL evaluations
        Na = nl_kdv(u_hat, ik, dealias);
        a  = E .* u_hat + 0.5*dt_actual * E .* Na;
        Nb = nl_kdv(a, ik, dealias);
        b  = E .* u_hat + 0.5*dt_actual * E .* Nb;
        Nc = nl_kdv(b, ik, dealias);
        c  = E2 .* u_hat + dt_actual * E .* Nc;
        Nd = nl_kdv(c, ik, dealias);
        u_hat = E2 .* u_hat + dt_actual/6 * (E2.*Na + 2*E.*(Nb+Nc) + Nd);

        if max(abs(u_hat)) > 1e8
            u = []; return;
        end
    end

    u = real(ifft(u_hat));
    if max(abs(u)) > 100
        u = [];
    end
end

function [u, x] = kdv_rk4_solve(N, dt, T)
    % Plain RK4 for KdV in Fourier space.
    [x, kv, dealias] = setup_kdv(N);
    u_hat = fft(cos(pi*x));

    ik  = 1i * kv;
    ik3 = 1i * kv.^3;

    n_steps = max(1, round(T / dt));
    dt_actual = T / n_steps;

    for step = 1:n_steps
        k1 = rhs_kdv(u_hat, ik, ik3, dealias);
        k2 = rhs_kdv(u_hat + 0.5*dt_actual*k1, ik, ik3, dealias);
        k3 = rhs_kdv(u_hat + 0.5*dt_actual*k2, ik, ik3, dealias);
        k4 = rhs_kdv(u_hat + dt_actual*k3, ik, ik3, dealias);
        u_hat = u_hat + dt_actual/6 * (k1 + 2*k2 + 2*k3 + k4);

        if max(abs(u_hat)) > 1e8
            u = []; return;
        end
    end

    u = real(ifft(u_hat));
    if max(abs(u)) > 100
        u = [];
    end
end

function f = rhs_kdv(v_hat, ik, ik3, dealias)
    % Right-hand side of KdV: -6*u*u_x - u_xxx in Fourier space
    v  = real(ifft(v_hat));
    vx = real(ifft(ik .* v_hat));
    f  = fft(-6*v.*vx) .* dealias - ik3 .* v_hat;
end

function f = nl_kdv(v_hat, ik, dealias)
    % Nonlinear part of KdV: -6*u*u_x in Fourier space
    v  = real(ifft(v_hat));
    vx = real(ifft(ik .* v_hat));
    f  = fft(-6*v.*vx) .* dealias;
end

function dtc = find_ifrk4_dt_crit(N, T)
    % Binary search for IF-RK4 critical dt
    dt_lo = 1e-6;
    dt_hi = 0.1;

    % Find unstable upper bound
    for iter = 1:20
        [u, ~] = kdv_ifrk4_solve(N, dt_hi, T);
        if isempty(u), break; end
        dt_hi = dt_hi * 2;
        if dt_hi > 10, dtc = dt_hi; return; end
    end

    % Verify dt_lo is stable
    [u, ~] = kdv_ifrk4_solve(N, dt_lo, T);
    if isempty(u), dt_lo = 1e-8; end

    for iter = 1:40
        dt_mid = sqrt(dt_lo * dt_hi);
        [u, ~] = kdv_ifrk4_solve(N, dt_mid, T);
        if ~isempty(u)
            dt_lo = dt_mid;
        else
            dt_hi = dt_mid;
        end
        if dt_hi / dt_lo < 1.02, break; end
    end
    dtc = sqrt(dt_lo * dt_hi);
end
