%% allen_cahn.m
%
% Allen-Cahn equation with periodic BCs using Fourier spectral IMEX scheme.
%
% PDE: u_t = u_xx + u(1 - u^2), 0 < x < 2*pi, t > 0 (periodic)
%
% The diffusion term u_xx is treated implicitly and the nonlinear reaction
% u(1 - u^2) is treated explicitly (IMEX):
%   (I - dt * D_xx) u^{n+1} = u^n + dt * (u^n - (u^n)^3)
%
% In Fourier space this becomes diagonal:
%   u_hat_k^{n+1} = (u_hat_k^n + dt * FFT(u^n - (u^n)^3)) / (1 + k^2 * dt)
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
% Part of "Computational Etudes: A Spectral Approach"
% https://github.com/dutykh/computational-etudes
% Date: February 2026

clear; close all; clc;

%% Configuration
N = 256;           % Grid points
tmax = 5.0;        % Final time
dt = 0.01;         % Time step (IMEX allows large dt)

% Colors
NAVY  = [0.078 0.176 0.431];
SKY   = [0.471 0.588 0.824];
CORAL = [0.906 0.298 0.235];
TEAL  = [0.086 0.627 0.522];

% Output path
script_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch10', 'matlab');
output_file = fullfile(output_dir, 'allen_cahn.pdf');

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

%% Initial condition: smoothed random perturbation
rng(42);  % Fixed seed for reproducibility
u = 0.1 * randn(N, 1);

% Low-pass filter to smooth the initial data
u_hat = fft(u);
filter_mask = abs(k) <= N/8;
u_hat = u_hat .* filter_mask;
u = real(ifft(u_hat));

%% IMEX multiplier: 1 / (1 + k^2 * dt)
imex_denom = 1.0 + k.^2 * dt;

%% Time stepping
nsteps = ceil(tmax / dt);

% Storage for snapshots
n_save = 200;
save_every = max(1, floor(nsteps / n_save));
U_save = zeros(N, n_save + 1);
t_save = zeros(1, n_save + 1);
U_save(:, 1) = u;
t_save(1) = 0;
save_idx = 2;

t = 0;
for n = 1:nsteps
    % Explicit reaction term: u - u^3 = u(1 - u^2)
    reaction = u - u.^3;

    % Transform to Fourier space
    u_hat = fft(u);
    reaction_hat = fft(reaction);

    % IMEX step in Fourier space:
    % u_hat^{n+1} = (u_hat^n + dt * reaction_hat) / (1 + k^2 * dt)
    u_hat = (u_hat + dt * reaction_hat) ./ imex_denom;

    % Transform back to physical space
    u = real(ifft(u_hat));

    t = t + dt;

    % Save snapshot
    if mod(n, save_every) == 0 && save_idx <= n_save + 1
        U_save(:, save_idx) = u;
        t_save(save_idx) = t;
        save_idx = save_idx + 1;
    end
end

% Trim unused storage
n_saved = save_idx - 1;
U_save = U_save(:, 1:n_saved);
t_save = t_save(1:n_saved);

%% Create figure
fig = figure('Units', 'inches', 'Position', [1, 1, 12, 5]);

% Left panel: space-time pcolormesh plot
subplot(1, 2, 1);
[X, T] = meshgrid(x, t_save);
pcolor(X, T, U_save');
shading interp;
colormap(subplot(1,2,1), redblue(256));
caxis([-1, 1]);
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 11);
ylabel('$t$', 'Interpreter', 'latex', 'FontSize', 11);
title('Allen-Cahn: Phase Separation', 'FontSize', 12);
cb = colorbar;
cb.Label.String = '$u(x, t)$';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 11;

% Right panel: snapshots at selected times
subplot(1, 2, 2);
hold on;

snapshot_times = [0.0, 0.1, 0.5, 1.0, 2.0, 5.0];
colors = [SKY; TEAL; 0.180 0.800 0.443; 0.953 0.612 0.071; CORAL; NAVY];

for i = 1:length(snapshot_times)
    [~, idx] = min(abs(t_save - snapshot_times(i)));
    t_actual = t_save(idx);
    plot(x, U_save(:, idx), '-', 'Color', colors(i, :), 'LineWidth', 1.5, ...
         'DisplayName', sprintf('$t = %.1f$', t_actual));
end

xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 11);
ylabel('$u(x, t)$', 'Interpreter', 'latex', 'FontSize', 11);
title('Snapshots at Selected Times', 'FontSize', 12);
lgd = legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 9);
lgd.BoxFace.ColorType = 'truecoloralpha';
lgd.BoxFace.ColorData = uint8([255; 255; 255; 230]);
xlim([0, 2*pi]);
ylim([-1.3, 1.3]);
yline(1.0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5, 'Alpha', 0.5);
yline(-1.0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5, 'Alpha', 0.5);
box on;

%% Save figure
exportgraphics(fig, output_file, 'ContentType', 'vector');
exportgraphics(fig, strrep(output_file, '.pdf', '.png'), 'Resolution', 300);

fprintf('Figure saved to: %s\n', output_file);

close(fig);


%% Helper function: Red-Blue diverging colormap
function cmap = redblue(m)
    % Red-Blue diverging colormap centered at white
    if nargin < 1
        m = 256;
    end
    n = ceil(m / 2);

    % Blue to white
    r1 = linspace(0.2, 1, n)';
    g1 = linspace(0.3, 1, n)';
    b1 = linspace(0.8, 1, n)';

    % White to red
    r2 = linspace(1, 0.8, n)';
    g2 = linspace(1, 0.2, n)';
    b2 = linspace(1, 0.2, n)';

    cmap = [r1 g1 b1; r2(2:end, :) g2(2:end, :) b2(2:end, :)];

    % Ensure exact size
    if size(cmap, 1) > m
        cmap = cmap(1:m, :);
    end
end
