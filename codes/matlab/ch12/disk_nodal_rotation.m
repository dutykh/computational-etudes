% disk_nodal_rotation.m - Eigenvalue degeneracy and nodal rotations
% Chapter 12: Spectral Methods in Polar Coordinates
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi, UAE)
% Last modified: February 2026
%
% Computational Etude 12.3: For a doubly-degenerate eigenvalue (m >= 1),
% the two eigenvectors span a 2D eigenspace. Forming u_phi = cos(phi)*u1 +
% sin(phi)*u2 produces continuously rotating nodal line patterns.

clear; close all;

% Add Chebyshev matrix from ch07
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'ch07'));

% Output directory
fig_dir = fullfile(fileparts(mfilename('fullpath')), ...
    '..', '..', '..', 'textbook', 'figures', 'ch12', 'matlab');
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

% -----------------------------------------------------------------------
% Build Laplacian and compute the first degenerate pair
% -----------------------------------------------------------------------

Nr = 25;
M  = 20;
N2 = (Nr - 1) / 2;

[L, r_pos, theta] = laplacian_polar(Nr, M);

[V, D_eig] = eig(-L);
eigenvalues = real(diag(D_eig));

[eigenvalues, idx] = sort(eigenvalues);
V = V(:, idx);

% The first degenerate pair is modes 2 and 3 (eigenvalue j_{1,1}^2 ~ 14.682)
u1 = real(V(:, 2));
u2 = real(V(:, 3));

% Build the coarse grid for nodal line plotting
U1 = reshape(u1, [N2, M]);
U2 = reshape(u2, [N2, M]);

% Extend: add boundary (r=1: u=0) row at top and periodic wrap in theta
U1_ext = zeros(N2 + 1, M + 1);
U2_ext = zeros(N2 + 1, M + 1);
U1_ext(2:end, 1:M) = U1;
U2_ext(2:end, 1:M) = U2;
U1_ext(:, M+1) = U1_ext(:, 1);
U2_ext(:, M+1) = U2_ext(:, 1);

r_ext     = [1.0; r_pos];
theta_ext = [theta; 2*pi];
[RR, TT]  = ndgrid(r_ext, theta_ext);
XX = RR .* cos(TT);
YY = RR .* sin(TT);

% -----------------------------------------------------------------------
% Figure 12.5: Nodal line rotation for 6 values of phi
% -----------------------------------------------------------------------

phi_values = linspace(0, 5*pi/6, 6);   % 0, pi/6, pi/3, pi/2, 2pi/3, 5pi/6

figure('Position', [100 100 900 630]);
th_circle = linspace(0, 2*pi, 200);

for k = 1:6
    phi = phi_values(k);
    subplot(2, 3, k);

    % Linear combination
    U_phi = cos(phi) * U1_ext + sin(phi) * U2_ext;

    hold on;
    % Disk boundary
    plot(cos(th_circle), sin(th_circle), 'k-', 'LineWidth', 1.0);

    % Light filled contour for context
    contourf(XX, YY, U_phi, 20, 'LineStyle', 'none');
    colormap(gca, rdbu_colormap());
    alpha(0.3);

    % Nodal lines (zero contour)
    contour(XX, YY, U_phi, [0 0], 'Color', [0.12 0.47 0.71], ...
        'LineWidth', 1.5);
    hold off;

    xlim([-1.1 1.1]);
    ylim([-1.1 1.1]);
    axis equal off;

    if k == 1
        title('$\varphi = 0$', 'Interpreter', 'latex', 'FontSize', 10);
    else
        title(sprintf('$\\varphi = %d\\pi/6$', k-1), ...
            'Interpreter', 'latex', 'FontSize', 10);
    end
end

sgtitle(sprintf('Nodal line rotation in a degenerate eigenspace ($\\lambda \\approx %.4f$)', ...
    eigenvalues(2)), 'Interpreter', 'latex', 'FontSize', 11);

exportgraphics(gcf, fullfile(fig_dir, 'nodal_rotation.pdf'), ...
    'ContentType', 'vector');

fprintf('Figure saved to: %s/nodal_rotation.pdf\n', fig_dir);
fprintf('Degenerate eigenvalue: lambda = %.10f (modes 2 and 3)\n', eigenvalues(2));


% =====================================================================
% Helper function
% =====================================================================

function cmap = rdbu_colormap()
% RDBU_COLORMAP  Red-Blue diverging colormap similar to matplotlib's RdBu_r.
    n = 256;
    nodes = [0.0196 0.1882 0.3804;
             0.1294 0.4000 0.6745;
             0.2627 0.5765 0.7647;
             0.5725 0.7725 0.8706;
             0.8196 0.8980 0.9412;
             0.9686 0.9686 0.9686;
             0.9922 0.8588 0.7804;
             0.9569 0.6471 0.5098;
             0.8392 0.3765 0.3020;
             0.6980 0.0941 0.1686;
             0.4039 0.0000 0.1216];
    x_nodes = linspace(0, 1, size(nodes, 1));
    xi = linspace(0, 1, n);
    cmap = interp1(x_nodes, nodes, xi);
end
