%% transport_variable.m
%
% Figure 11.x: Variable coefficient transport equation (periodic)
%
% PDE: u_t + c(x) * u_x = 0, 0 < x < 2*pi, t > 0 (periodic)
% Speed: c(x) = 0.3 + sin^2(x - 1)
% ICs: u(x, 0) = exp(-20*(x - pi)^2)
%
% Uses FFT for spatial differentiation and RK4 for time stepping.
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
% Part of "Computational Etudes: A Spectral Approach"
% https://github.com/dutykh/computational-etudes

clear; close all; clc;

%% Configuration
N = 256;          % Grid points
tmax = 8.0;        % Final time

% Colors
NAVY = [0.078 0.176 0.431];
CORAL = [0.906 0.298 0.235];

% Output path
script_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch11', 'matlab');
output_file = fullfile(output_dir, 'transport_variable.pdf');

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% Setup grid
h = 2 * pi / N;
x = h * (0:N-1)';

% Variable speed
c = 0.3 + sin(x - 1).^2;

%% Initial condition: Gaussian pulse
u = exp(-20 * (x - pi).^2);

%% Time stepping (RK4)
c_max = max(c);
dt = 0.4 * h / c_max;
nsteps = ceil(tmax / dt);

% RHS function
rhs = @(u) -c .* fourier_diff(u);

% Storage for waterfall plot
n_save = 120;
save_every = max(1, floor(nsteps / n_save));
U_save = zeros(N, n_save + 1);
t_save = zeros(1, n_save + 1);
U_save(:, 1) = u;
t_save(1) = 0;
save_idx = 2;

t = 0;
for n = 1:nsteps
    % RK4 time stepping
    k1 = rhs(u);
    k2 = rhs(u + 0.5 * dt * k1);
    k3 = rhs(u + 0.5 * dt * k2);
    k4 = rhs(u + dt * k3);
    u = u + (dt / 6) * (k1 + 2*k2 + 2*k3 + k4);

    t = t + dt;

    % Save snapshot
    if mod(n, save_every) == 0 && save_idx <= n_save + 1
        U_save(:, save_idx) = u;
        t_save(save_idx) = t;
        save_idx = save_idx + 1;
    end
end

%% Create figure
fig = figure('Units', 'inches', 'Position', [1, 1, 12, 6]);

% Left panel: space-time contour plot
subplot(1, 2, 1);
[X, T] = meshgrid(x, t_save(1:save_idx-1));
pcolor(X, T, U_save(:, 1:save_idx-1)');
shading interp;
colormap('jet');
xlabel('x', 'FontSize', 11);
ylabel('t', 'FontSize', 11);
title('Variable Coefficient Transport', 'FontSize', 12);
cb = colorbar;
cb.Label.String = 'u(x, t)';
cb.Label.FontSize = 11;

% Right panel: waterfall view
subplot(1, 2, 2);
hold on;

n_lines = 40;
indices = round(linspace(1, save_idx - 1, n_lines));
cmap = parula(n_lines);

for i = 1:n_lines
    idx = indices(i);
    offset = t_save(idx) * 0.05;
    plot(x, U_save(:, idx) + offset, '-', 'Color', cmap(i, :), 'LineWidth', 0.8);
end

xlabel('x', 'FontSize', 11);
ylabel('u(x, t) (offset by t)', 'FontSize', 11);
title('Waterfall View: Pulse in Variable Speed', 'FontSize', 12);
xlim([0, 2*pi]);
box on;

% Add inset for speed profile
axes('Position', [0.75, 0.75, 0.15, 0.15]);
plot(x, c, '-', 'Color', CORAL, 'LineWidth', 1.5);
xlabel('x', 'FontSize', 8);
ylabel('c(x)', 'FontSize', 8);
title('Speed Profile', 'FontSize', 9);
xlim([0, 2*pi]);
box on;
set(gca, 'FontSize', 7);

%% Save figure
exportgraphics(fig, output_file, 'ContentType', 'vector');
exportgraphics(fig, strrep(output_file, '.pdf', '.png'), 'Resolution', 300);

fprintf('Figure saved to: %s\n', output_file);

close(fig);


%% Helper function: Fourier differentiation
function w = fourier_diff(v)
    N = length(v);
    v_hat = fft(v);

    % Wavenumbers
    if mod(N, 2) == 0
        k = [0:N/2-1, 0, -N/2+1:-1]';
    else
        k = [0:(N-1)/2, -(N-1)/2:-1]';
    end

    w = real(ifft(1i * k .* v_hat));
end
