% fair_comparison.m
% Chapter 17: Stability of Spectral Methods
%
% Computational Etude 17.6: Fair Comparison of Time-Stepping Methods.
%
% Two-panel figure comparing Forward Euler, Backward Euler, Crank-Nicolson,
% and Integrating Factor (exact exponential) on the 1D Fourier heat equation
% u_t = nu * u_xx on [0, 2*pi) with nu = 0.1.
%
%   (a) Error vs dt (log-log).
%   (b) Error vs wall-clock time (log-log).
%
% Author: Dr. Denys Dutykh
% Part of "Computational Etudes: A Spectral Approach"

clear; clc; close all;

NAVY  = [20, 45, 110]/255;
TEAL  = [22, 160, 133]/255;
CORAL = [231, 76, 60]/255;
AMBER = [230, 126, 34]/255;

% Output path
script_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch17', 'matlab');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

N  = 64;
nu = 0.1;
T  = 1.0;

%% Exact solution in Fourier space
x = linspace(0, 2*pi, N+1); x(end) = [];
u0 = sin(x) + 0.5*cos(3*x);
u0_hat = fft(u0);
k = [0:N/2-1, -N/2:-1];
decay = exp(-nu * k.^2 * T);
u_exact = real(ifft(u0_hat .* decay));

% CFL for forward Euler: dt <= 2 / (nu * (N/2)^2)
dt_cfl = 2.0 / (nu * (N/2)^2);
fprintf('Forward Euler CFL limit: dt <= %.6e\n', dt_cfl);

%% Time step values
dt_values_all = logspace(-5, -1, 30);
dt_values_fe  = dt_values_all(dt_values_all < 0.9 * dt_cfl);

%% Run all four methods
methods = {'Forward Euler', 'Backward Euler', 'Crank-Nicolson', 'Integrating Factor'};
colors  = {NAVY, TEAL, CORAL, AMBER};
markers = {'o', 's', '^', 'd'};

results = struct();

for m = 1:4
    if m == 1
        dt_vals = dt_values_fe;
    else
        dt_vals = dt_values_all;
    end

    errs  = [];
    times = [];
    dts   = [];

    for di = 1:length(dt_vals)
        dt = dt_vals(di);
        n_steps = round(T / dt);
        dt_actual = T / n_steps;
        lam = -nu * k.^2;

        u_hat = fft(u0);
        blowup = false;

        switch m
            case 1  % Forward Euler
                tic;
                for step = 1:n_steps
                    u_hat = u_hat + dt_actual * lam .* u_hat;
                    if max(abs(u_hat)) > 1e10
                        blowup = true; break;
                    end
                end
                elapsed = toc;

            case 2  % Backward Euler
                factor = 1.0 ./ (1.0 - dt_actual * lam);
                tic;
                for step = 1:n_steps
                    u_hat = u_hat .* factor;
                end
                elapsed = toc;

            case 3  % Crank-Nicolson
                factor = (1.0 + 0.5*dt_actual*lam) ./ (1.0 - 0.5*dt_actual*lam);
                tic;
                for step = 1:n_steps
                    u_hat = u_hat .* factor;
                end
                elapsed = toc;

            case 4  % Integrating Factor (exact exponential)
                factor = exp(lam * dt_actual);
                tic;
                for step = 1:n_steps
                    u_hat = u_hat .* factor;
                end
                elapsed = toc;
        end

        if blowup, continue; end

        u_num = real(ifft(u_hat));
        err = max(abs(u_num - u_exact));
        if err < 1e-14, err = 1e-14; end

        errs  = [errs, err];     %#ok<AGROW>
        times = [times, elapsed]; %#ok<AGROW>
        dts   = [dts, dt];       %#ok<AGROW>
    end

    results(m).name   = methods{m};
    results(m).errs   = errs;
    results(m).times  = times;
    results(m).dts    = dts;
    results(m).color  = colors{m};
    results(m).marker = markers{m};
end

%% Plot
figure('Position', [100, 100, 1200, 500]);

% Panel (a): Error vs dt
subplot(1, 2, 1);
hold on;
for m = 1:4
    if ~isempty(results(m).dts)
        loglog(results(m).dts, results(m).errs, ...
               ['-', results(m).marker], 'Color', results(m).color, ...
               'MarkerSize', 5, 'LineWidth', 1.2, ...
               'MarkerFaceColor', results(m).color, ...
               'DisplayName', results(m).name);
    end
end

% Reference slopes
dt_ref = logspace(-5, -1, 50);
loglog(dt_ref, 0.5*dt_ref, '--', 'Color', [0.5 0.5 0.5], ...
       'LineWidth', 0.8, 'HandleVisibility', 'off');
text(2e-3, 2e-3, 'slope 1', 'FontSize', 8, 'Color', [0.5 0.5 0.5], ...
     'Rotation', 35);
loglog(dt_ref, 2*dt_ref.^2, ':', 'Color', [0.5 0.5 0.5], ...
       'LineWidth', 0.8, 'HandleVisibility', 'off');
text(5e-3, 8e-5, 'slope 2', 'FontSize', 8, 'Color', [0.5 0.5 0.5], ...
     'Rotation', 42);
hold off;

xlabel('\Delta t');
ylabel('||u_{num} - u_{exact}||_\infty');
title('(a) Error vs time step');
legend('Location', 'southeast', 'FontSize', 8);
ylim([1e-15, 1e1]);
grid on; set(gca, 'GridAlpha', 0.3);

% Panel (b): Error vs wall-clock time
subplot(1, 2, 2);
hold on;
for m = 1:4
    if ~isempty(results(m).times)
        loglog(results(m).times, results(m).errs, ...
               ['-', results(m).marker], 'Color', results(m).color, ...
               'MarkerSize', 5, 'LineWidth', 1.2, ...
               'MarkerFaceColor', results(m).color, ...
               'DisplayName', results(m).name);
    end
end
hold off;

xlabel('Wall-clock time (s)');
ylabel('||u_{num} - u_{exact}||_\infty');
title('(b) Error vs computational cost');
legend('Location', 'northeast', 'FontSize', 8);
ylim([1e-15, 1e1]);
grid on; set(gca, 'GridAlpha', 0.3);

set(gcf, 'PaperPositionMode', 'auto');
exportgraphics(gcf, fullfile(output_dir, 'fair_comparison.pdf'), 'ContentType', 'vector');
exportgraphics(gcf, fullfile(output_dir, 'fair_comparison.png'), 'Resolution', 300);
fprintf('\nFigure saved to %s\n', fullfile(output_dir, 'fair_comparison.pdf'));
