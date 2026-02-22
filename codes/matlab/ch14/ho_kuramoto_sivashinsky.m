%% ho_kuramoto_sivashinsky.m - Kuramoto-Sivashinsky equation with ETDRK4
%
% Chapter 14: Higher-Order Boundary Value Problems
%
% Computational Etude 14.8: Kuramoto-Sivashinsky equation
%
%   u_t + u*u_x + u_xx + u_xxxx = 0
%
% on [0, 32*pi] with periodic boundary conditions. Solved using Fourier
% pseudospectral methods in space and the ETDRK4 scheme (Cox & Matthews
% 2002, with Kassam-Trefethen 2005 contour integral improvement) in time.
%
% In Fourier space, the equation becomes:
%
%   du_hat/dt = L_k * u_hat + N_hat(u)
%
% where L_k = k^2 - k^4 (energy production at low k, dissipation at
% high k) and N_hat = -(ik/2) * FFT(u^2) uses the conservation form
% u*u_x = (u^2/2)_x.
%
% The ETDRK4 scheme evaluates the phi-functions via Kassam & Trefethen's
% contour integral trick on a circle in the complex plane, avoiding
% catastrophic cancellation near z = 0.
%
% The solution exhibits spatiotemporal chaos: cell merging, splitting,
% and irregular oscillations.
%
% Generates Figure: kuramoto_sivashinsky.pdf - Space-time plot
%
% Author: Dr. Denys Dutykh
%         Mathematics Department
%         Khalifa University of Science and Technology
%         Abu Dhabi, UAE
%
% Part of the book "Computational Etudes: A Spectral Approach"
% https://github.com/dutykh/computational-etudes
%
% Last modified: February 2026

clear;
close all;
clc;

%% Publication-quality figure settings
set(groot, 'DefaultAxesFontSize', 10);
set(groot, 'DefaultAxesFontName', 'CMU Serif');
set(groot, 'DefaultTextFontSize', 10);
set(groot, 'DefaultTextFontName', 'CMU Serif');
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');
set(groot, 'DefaultAxesLineWidth', 0.8);
set(groot, 'DefaultLineLineWidth', 1.5);

%% Color scheme
NAVY   = [20,  45, 110] / 255;
SKY    = [120, 150, 210] / 255;
CORAL  = [231,  76,  60] / 255;
TEAL   = [22,  160, 133] / 255;
PURPLE = [142,  68, 173] / 255;
ORANGE = [230, 126,  34] / 255;

%% Output directory
script_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch14', 'matlab');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% Parameters
N_modes = 256;               % Number of Fourier modes
L       = 32.0 * pi;        % Domain length
t_final = 150.0;            % Final time
dt      = 0.5;              % Time step
n_snapshots = 500;          % Number of snapshots to store

fprintf('Kuramoto-Sivashinsky Equation with ETDRK4\n');
fprintf('%s\n', repmat('=', 1, 60));
fprintf('  N = %d, L = %.0f*pi, T = %.0f, dt = %.1f\n', ...
        N_modes, L / pi, t_final, dt);
fprintf('  Initial condition: u(x,0) = cos(x/16) * (1 + sin(x/16))\n');
fprintf('  Computing...\n');

%% Spatial grid
x = L * (0:N_modes - 1)' / N_modes;

%% Wavenumbers
% k = (2*pi/L) * [0, 1, 2, ..., N/2-1, -N/2, ..., -1]
k = (2.0 * pi / L) * [0:N_modes/2 - 1, -N_modes/2:-1]';

%% Linear operator: L_k = k^2 - k^4
Lop = k.^2 - k.^4;

%% 2/3 dealiasing mask
dealias = ones(N_modes, 1);
cutoff = floor(N_modes / 3);
if cutoff > 0
    dealias(cutoff + 2 : N_modes - cutoff) = 0;
end

%% Time stepping setup
n_steps = ceil(t_final / dt);
dt = t_final / n_steps;  % Adjust dt to hit t_final exactly
save_every = max(1, floor(n_steps / n_snapshots));

%% Compute ETDRK4 coefficients via Kassam-Trefethen contour integral trick
M_contour = 32;  % Number of contour quadrature points

% Full-step and half-step propagators
E  = exp(dt * Lop);
E2 = exp(dt * Lop / 2.0);

% Roots of unity for the contour integral
r = exp(2i * pi * (1:M_contour) / M_contour);  % row vector (1 x M)

% LR(i, m) = dt * L_i + r_m, evaluated on a circle of radius 1
% around each dt * L_i in the complex plane
LR = dt * Lop + r;  % broadcasting: (N x 1) + (1 x M) -> (N x M)

% Contour integral evaluations of the phi-functions
Q  = dt * real(mean((exp(LR / 2.0) - 1.0) ./ LR, 2));

f1 = dt * real(mean( ...
    (-4.0 - LR + exp(LR) .* (4.0 - 3.0 * LR + LR.^2)) ./ LR.^3, 2));

f2 = dt * real(mean( ...
    (2.0 + LR + exp(LR) .* (-2.0 + LR)) ./ LR.^3, 2));

f3 = dt * real(mean( ...
    (-4.0 - 3.0 * LR - LR.^2 + exp(LR) .* (4.0 - LR)) ./ LR.^3, 2));

%% Initial condition
u = cos(x / 16.0) .* (1.0 + sin(x / 16.0));
v = fft(u);

%% Storage for snapshots
snapshots = zeros(n_snapshots + 1, N_modes);
snap_times = zeros(n_snapshots + 1, 1);
snapshots(1, :) = u';
snap_times(1) = 0.0;
snap_count = 1;

%% ETDRK4 time integration loop
t = 0.0;
for n_step = 1:n_steps
    % Nonlinear term: N_hat = -0.5i * k .* fft(u.^2)
    Nv = nhat(v, k, dealias);
    a  = E2 .* v + Q .* Nv;
    Na = nhat(a, k, dealias);
    b  = E2 .* v + Q .* Na;
    Nb = nhat(b, k, dealias);
    c  = E2 .* a + Q .* (2.0 * Nb - Nv);
    Nc = nhat(c, k, dealias);
    v  = E .* v + f1 .* Nv + 2.0 * f2 .* (Na + Nb) + f3 .* Nc;

    t = t + dt;

    if mod(n_step, save_every) == 0 && snap_count < n_snapshots + 1
        snap_count = snap_count + 1;
        u = real(ifft(v));
        snapshots(snap_count, :) = u';
        snap_times(snap_count) = t;
    end
end

% Trim storage
snapshots = snapshots(1:snap_count, :);
snap_times = snap_times(1:snap_count);

fprintf('  %d snapshots saved, t in [%.1f, %.1f]\n', ...
        snap_count, snap_times(1), snap_times(end));
fprintf('  max |u| = %.4f\n', max(abs(snapshots(:))));

%% Create figure: Space-time pcolormesh plot
fig = figure('Units', 'inches', 'Position', [1, 1, 8, 6]);

[X, T] = meshgrid(x, snap_times);

% Use a symmetric color range based on the 98th percentile
% (after the initial transient has passed)
n_transient = floor(snap_count / 5);
vmax = prctile(abs(snapshots(n_transient:end, :)), 98, 'all');

% Custom RdBu_r colormap (diverging blue-white-red)
colormap(rdbu_colormap(256));

pcolor(X / pi, T, snapshots);
shading interp;
caxis([-vmax, vmax]);

xlabel('$x / \pi$');
ylabel('$t$');
title(sprintf(['Kuramoto--Sivashinsky Equation: ETDRK4 ', ...
       '($L = 32\\pi$, $N = %d$)'], N_modes), 'FontSize', 11);

cb = colorbar;
cb.Label.String = '$u(x, t)$';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 10;

xlim([0, L / pi]);
ylim([0, t_final]);
box on;

sgtitle('Computational \''Etude 14.8: Kuramoto--Sivashinsky Equation', ...
        'FontSize', 13, 'Interpreter', 'latex');

%% Save figure
exportgraphics(fig, fullfile(output_dir, 'kuramoto_sivashinsky.pdf'), ...
               'ContentType', 'vector');
fprintf('\nFigure saved to: %s\n', fullfile(output_dir, 'kuramoto_sivashinsky.pdf'));

close(fig);

fprintf('\n%s\n', repmat('=', 1, 60));


%% =====================================================================
%  Local functions
%  =====================================================================

function nl = nhat(v_hat, k, dealias)
% NHAT  Evaluate dealiased nonlinear term in Fourier space.
%   N_hat = -0.5i * k .* fft(u.^2) with 2/3 dealiasing.
    u_phys = real(ifft(v_hat));
    nl = -0.5i * k .* fft(u_phys.^2);
    nl = nl .* dealias;
end


function cmap = rdbu_colormap(m)
% RDBU_COLORMAP  Red-Blue diverging colormap similar to matplotlib's RdBu_r.
    if nargin < 1, m = 256; end
    nodes = [0.0196 0.1882 0.3804;   % dark blue
             0.1294 0.4000 0.6745;   % blue
             0.2627 0.5765 0.7647;   % light blue
             0.5725 0.7725 0.8706;   % very light blue
             0.8196 0.8980 0.9412;   % pale blue
             0.9686 0.9686 0.9686;   % near white
             0.9922 0.8588 0.7804;   % pale red
             0.9569 0.6471 0.5098;   % light red
             0.8392 0.3765 0.3020;   % red
             0.6980 0.0941 0.1686;   % dark red
             0.4039 0.0000 0.1216];  % very dark red
    x_nodes = linspace(0, 1, size(nodes, 1));
    xi = linspace(0, 1, m);
    cmap = interp1(x_nodes, nodes, xi);
end
