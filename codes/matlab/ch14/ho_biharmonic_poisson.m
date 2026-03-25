%% ho_biharmonic_poisson.m - Biharmonic Poisson BVP with manufactured solution
%
% Chapter 14: Higher-Order Spectral Methods
%
% Biharmonic Poisson BVP on a clamped square plate:
%
%   Delta^2 u = f(x,y),   (x,y) in (-1,1)^2
%   u = 0, du/dn = 0      on the boundary
%
% Manufactured exact solution:
%   u(x,y) = sin^2(pi*x) * sin^2(pi*y)
%
% Right-hand side:
%   f(x,y) = 4*pi^4 * [4*cos(2*pi*x)*cos(2*pi*y)
%                       - cos(2*pi*x) - cos(2*pi*y)]
%
% Generates Figure: biharmonic_poisson.pdf
%   Left panel:  contourf of numerical solution at N = 24
%   Right panel: convergence (max error vs N, semilogy)
%
% Author: Dr. Denys Dutykh
%         Mathematics Department
%         Khalifa University of Science and Technology
%         Abu Dhabi, UAE
%
% Part of the book "Computational Etudes: A Spectral Approach"
% https://github.com/dutykh/computational-etudes

clear; close all; clc;

% Add path to Chebyshev matrix from Chapter 7
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'ch07'));

% Publication-quality settings
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultAxesTickLabelInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0, 'DefaultAxesFontSize', 10);
set(0, 'DefaultTextFontSize', 11);

% Book color scheme
NAVY   = [20, 45, 110] / 255;
SKY    = [120, 150, 210] / 255;
CORAL  = [231, 76, 60] / 255;
TEAL   = [22, 160, 133] / 255;
PURPLE = [142, 68, 173] / 255;
ORANGE = [230, 126, 34] / 255;

% Output directory
script_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch14', 'matlab');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

%% Convergence study
N_values = 6:2:24;
errors = zeros(size(N_values));

fprintf('Biharmonic Poisson BVP with Manufactured Solution\n');
fprintf('%s\n', repmat('=', 1, 55));
fprintf('  u(x,y) = sin^2(pi*x) * sin^2(pi*y)\n');
fprintf('  f(x,y) = 4*pi^4*(4*cos(2*pi*x)*cos(2*pi*y)\n');
fprintf('                    - cos(2*pi*x) - cos(2*pi*y))\n\n');
fprintf('  Convergence study:\n');
fprintf('  %s\n', repmat('-', 1, 40));
fprintf('  %4s   %12s   %12s\n', 'N', 'Interior pts', 'Max error');
fprintf('  %s\n', repmat('-', 1, 40));

U_plot = [];
x_plot = [];

for idx = 1:length(N_values)
    N = N_values(idx);
    [U, x, err] = solve_biharmonic_poisson(N);
    errors(idx) = err;
    fprintf('  %4d   %12d   %12.4e\n', N, (N - 1)^2, err);

    if N == N_values(end)
        U_plot = U;
        x_plot = x;
    end
end

fprintf('  %s\n\n', repmat('-', 1, 40));

%% Create figure: two panels
fig = figure('Units', 'inches', 'Position', [1, 1, 11, 4.5]);

% Left panel: contourf of solution
subplot(1, 2, 1);
[XX, YY] = meshgrid(x_plot, x_plot);
max_val = max(abs(U_plot(:)));
levels = linspace(-max_val, max_val, 21);
contourf(XX, YY, U_plot, levels, 'LineColor', 'none');
colormap(gca, rdbu_colormap(256));
caxis([-max_val, max_val]);
hold on;
contour(XX, YY, U_plot, 8, 'Color', NAVY, 'LineWidth', 0.6);
plot([-1, 1, 1, -1, -1], [-1, -1, 1, 1, -1], 'k-', 'LineWidth', 1.0);
hold off;
colorbar;
xlim([-1.05, 1.05]);
ylim([-1.05, 1.05]);
axis equal;
set(gca, 'XTick', [-1, 0, 1], 'YTick', [-1, 0, 1]);
xlabel('$x$');
ylabel('$y$');
title(sprintf('Numerical solution ($N = %d$)', N_values(end)));

% Right panel: convergence
subplot(1, 2, 2);
semilogy(N_values, errors, 'o-', 'Color', CORAL, 'MarkerSize', 5, ...
         'LineWidth', 1.5, 'MarkerFaceColor', CORAL);
hold on;
yline(2.2e-16, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);
hold off;
xlabel('$N$');
ylabel('Max error');
title('Spectral convergence');
legend('$\|u - u_{\mathrm{exact}}\|_\infty$', ...
       'Machine $\varepsilon$', ...
       'Location', 'northeast');
grid on;
set(gca, 'XTick', N_values);

sgtitle(['Biharmonic Poisson BVP: ', ...
         '$\Delta^2 u = f$, ', ...
         '$u = \partial_n u = 0$ on $\partial\Omega$'], ...
        'FontSize', 13);

%% Save figure
exportgraphics(fig, fullfile(output_dir, 'biharmonic_poisson.pdf'), ...
               'ContentType', 'vector');
fprintf('  Figure saved to: %s\n', fullfile(output_dir, 'biharmonic_poisson.pdf'));
fprintf('%s\n', repmat('=', 1, 55));

close(fig);


%% =====================================================================
%  Local functions
%  =====================================================================

function [U, x, err] = solve_biharmonic_poisson(N)
% SOLVE_BIHARMONIC_POISSON  Solve the biharmonic Poisson BVP on [-1,1]^2.

    [D4, D2int, x] = build_clamped_biharmonic(N);

    M = N - 1;
    I_M = eye(M);

    % 2D biharmonic via Kronecker products
    L = kron(I_M, D4) + kron(D4, I_M) + 2 * kron(D2int, D2int);

    % Interior grid
    xi = x(2:N);
    [XX, YY] = meshgrid(xi, xi);

    % RHS and solve
    f_vec = f_rhs(XX(:), YY(:));
    u_vec = L \ f_vec;

    % Error
    u_ex_vec = u_exact(XX(:), YY(:));
    err = max(abs(u_vec - u_ex_vec));

    % Reshape to full grid
    U = zeros(N + 1, N + 1);
    U(2:N, 2:N) = reshape(u_vec, [M, M]);
end

function [D4, D2int, x] = build_clamped_biharmonic(N)
% BUILD_CLAMPED_BIHARMONIC  1D clamped biharmonic operator via polynomial trick.

    [D, x] = cheb_matrix(N);

    D2 = D * D;
    D3 = D2 * D;
    D4full = D3 * D;

    x2 = x.^2;
    s = zeros(N + 1, 1);
    s(2:N) = 1.0 ./ (1.0 - x2(2:N));
    S = diag(s);

    D4mat = (diag(1.0 - x2) * D4full - 8.0 * diag(x) * D3 - 12.0 * D2) * S;

    D4 = D4mat(2:N, 2:N);
    D2int = D2(2:N, 2:N);
end

function val = u_exact(x, y)
% U_EXACT  Manufactured exact solution.
    val = sin(pi * x).^2 .* sin(pi * y).^2;
end

function val = f_rhs(x, y)
% F_RHS  Right-hand side: Delta^2 u_exact.
    c2x = cos(2 * pi * x);
    c2y = cos(2 * pi * y);
    val = 4 * pi^4 * (4 * c2x .* c2y - c2x - c2y);
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
