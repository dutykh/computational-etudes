%% bc_allen_cahn.m - Allen-Cahn equation with fixed and driven boundaries
%
% Chapter 13: Boundary Conditions in Spectral Methods
%
% Allen-Cahn equation with non-periodic boundary conditions:
%
%   u_t = eps * u_xx + u - u^3,  x in (-1, 1),  t in (0, 60)
%
% Parameters: eps = 0.01, N = 80
% IC: u(x, 0) = 0.53*x + 0.47*sin(-1.5*pi*x)
%
% Part A: Fixed boundaries u(-1) = -1, u(1) = 1
%   - Lifting: w(x) = x, v = u - w with v(+-1) = 0
%   - IMEX: implicit diffusion, explicit nonlinearity
%
% Part B: Driven boundary u(-1) = -1, u(1) = 1 + sin^2(t/5)
%   - Direct boundary overwrite at each time step
%
% Generates:
%   Figure 1: allen_cahn_fixed.pdf  - space-time pcolor plot (Part A)
%   Figure 2: allen_cahn_driven.pdf - space-time pcolor plot (Part B)
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

% Add path to Chapter 7 functions
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'ch07'));

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
output_dir = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch13', 'matlab');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% Parameters
eps_ac = 0.01;      % Diffusion coefficient
N      = 80;        % Number of Chebyshev intervals
t_final = 60.0;     % Final time
dt     = 0.05;      % Time step
nsteps = round(t_final / dt);

%% Setup grid and operators
[D, x] = cheb_matrix(N);
D2 = D * D;

% Interior block (rows/cols 2..N for homogeneous Dirichlet)
D2_int = D2(2:N, 2:N);
I_int  = eye(N - 1);

% IMEX implicit matrix: (I - eps*dt*D2_int)
A = I_int - eps_ac * dt * D2_int;

% Pre-factorize for repeated solves
[L_lu, U_lu, P_lu] = lu(A);

%% Initial condition
u0 = 0.53 * x + 0.47 * sin(-1.5 * pi * x);

%% Snapshot storage
n_save = 400;
save_every = max(1, floor(nsteps / n_save));

%% =====================================================================
%  Part A: Fixed boundaries u(-1) = -1, u(1) = 1
%  =====================================================================
fprintf('Part A: Fixed boundaries...\n');

% Lifting: w(x) = x, so v = u - x with v(+-1) = 0
% u = v + x, u^3 = (v + x)^3
% u_t = v_t, u_xx = v_xx (since w_xx = 0)
% v_t = eps * v_xx + (v + x) - (v + x)^3

u = u0;
U_save_A = zeros(N + 1, n_save + 1);
t_save_A = zeros(1, n_save + 1);
U_save_A(:, 1) = u;
t_save_A(1) = 0;
save_idx = 2;

v = u - x;  % Transform to lifted variable

for n = 1:nsteps
    % Reconstruct u for nonlinear term
    u_full = v + x;

    % Explicit nonlinear term: u - u^3
    nl = u_full - u_full.^3;

    % IMEX step on interior: (I - eps*dt*D2) v^{n+1}_int = v^n_int + dt * nl_int
    rhs = v(2:N) + dt * nl(2:N);
    v_int = U_lu \ (L_lu \ (P_lu * rhs));

    % Update interior, boundaries stay zero
    v(2:N) = v_int;
    v(1) = 0;      % v(1) = 0 at x = 1
    v(N+1) = 0;    % v(N+1) = 0 at x = -1

    % Save snapshot
    t = n * dt;
    if mod(n, save_every) == 0 && save_idx <= n_save + 1
        U_save_A(:, save_idx) = v + x;
        t_save_A(save_idx) = t;
        save_idx = save_idx + 1;
    end
end

n_saved_A = save_idx - 1;
U_save_A = U_save_A(:, 1:n_saved_A);
t_save_A = t_save_A(1:n_saved_A);

%% Figure 1: Allen-Cahn with fixed boundaries
fig1 = figure('Units', 'inches', 'Position', [1, 1, 7, 5]);

[X, T] = meshgrid(x, t_save_A);
pcolor(X, T, U_save_A');
shading interp;
colormap(redblue(256));
caxis([-1, 1]);
xlabel('$x$');
ylabel('$t$');
title('Allen--Cahn: Fixed Boundaries $u(\pm 1) = \pm 1$', 'FontSize', 12);
cb = colorbar;
cb.Label.String = '$u(x, t)$';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 11;

exportgraphics(fig1, fullfile(output_dir, 'allen_cahn_fixed.pdf'), ...
               'ContentType', 'vector');
fprintf('Figure 1 saved to: %s\n', fullfile(output_dir, 'allen_cahn_fixed.pdf'));
close(fig1);

%% =====================================================================
%  Part B: Driven boundary u(-1) = -1, u(1) = 1 + sin^2(t/5)
%  =====================================================================
fprintf('Part B: Driven boundaries...\n');

u = u0;
U_save_B = zeros(N + 1, n_save + 1);
t_save_B = zeros(1, n_save + 1);
U_save_B(:, 1) = u;
t_save_B(1) = 0;
save_idx = 2;

for n = 1:nsteps
    t = n * dt;

    % Explicit nonlinear term: u - u^3
    nl = u - u.^3;

    % IMEX step on interior: (I - eps*dt*D2) u^{n+1}_int = u^n_int + dt * nl_int
    rhs = u(2:N) + dt * nl(2:N);

    % Add implicit contribution from boundary rows
    % D2_int * u_int has coupling to boundary values through D2
    % We need: (I - eps*dt*D2_int) * u_int = u_int + dt*nl_int
    %   + eps*dt * (D2(2:N, 1)*u(1) + D2(2:N, N+1)*u(N+1))
    % But boundary values change, so include their contribution
    bc_left  = -1.0;                      % u(-1) = -1
    bc_right = 1.0 + sin(t / 5.0)^2;     % u(1) = 1 + sin^2(t/5)

    rhs = rhs + eps_ac * dt * (D2(2:N, 1) * bc_right + D2(2:N, N+1) * bc_left);

    u_int = U_lu \ (L_lu \ (P_lu * rhs));

    % Update with boundary overwrite
    u(2:N) = u_int;
    u(1)   = bc_right;   % x(1) = 1
    u(N+1) = bc_left;    % x(N+1) = -1

    % Save snapshot
    if mod(n, save_every) == 0 && save_idx <= n_save + 1
        U_save_B(:, save_idx) = u;
        t_save_B(save_idx) = t;
        save_idx = save_idx + 1;
    end
end

n_saved_B = save_idx - 1;
U_save_B = U_save_B(:, 1:n_saved_B);
t_save_B = t_save_B(1:n_saved_B);

%% Figure 2: Allen-Cahn with driven boundaries
fig2 = figure('Units', 'inches', 'Position', [1, 1, 7, 5]);

[X, T] = meshgrid(x, t_save_B);
pcolor(X, T, U_save_B');
shading interp;
colormap(redblue(256));
caxis([-1, 2]);
xlabel('$x$');
ylabel('$t$');
title('Allen--Cahn: Driven Boundary $u(1) = 1 + \sin^2(t/5)$', 'FontSize', 12);
cb = colorbar;
cb.Label.String = '$u(x, t)$';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 11;

exportgraphics(fig2, fullfile(output_dir, 'allen_cahn_driven.pdf'), ...
               'ContentType', 'vector');
fprintf('Figure 2 saved to: %s\n', fullfile(output_dir, 'allen_cahn_driven.pdf'));
close(fig2);


%% =====================================================================
%  Helper function: Red-Blue diverging colormap
%  =====================================================================
function cmap = redblue(m)
    if nargin < 1, m = 256; end
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

    if size(cmap, 1) > m
        cmap = cmap(1:m, :);
    end
end
