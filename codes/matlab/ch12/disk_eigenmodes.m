% disk_eigenmodes.m - Eigenmodes of the Laplacian on the unit disk
% Chapter 12: Spectral Methods in Polar Coordinates
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi, UAE)
% Last modified: February 2026
%
% Computational Etude 12.2: Computes eigenmodes of -Delta u = lambda u
% on the disk with u = 0 on the boundary. Compares eigenvalues with
% Bessel function zeros to demonstrate spectral accuracy.

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
% Build the Laplacian and compute eigenmodes
% -----------------------------------------------------------------------

Nr = 25;   % radial points (odd)
M  = 20;   % angular points (even)
N2 = (Nr - 1) / 2;

[L, r_pos, theta] = laplacian_polar(Nr, M);

% Solve the eigenvalue problem -L v = lambda v
[V, D_eig] = eig(-L);
eigenvalues = real(diag(D_eig));

% Sort by eigenvalue
[eigenvalues, idx] = sort(eigenvalues);
V = V(:, idx);

% -----------------------------------------------------------------------
% Figure 12.3: Six representative eigenmodes as surface plots
% -----------------------------------------------------------------------

% Mode indices (1-based): modes 1, 2, 4, 6, 7, 9
mode_indices = [1, 2, 4, 6, 7, 9];

% Extend grid for plotting: add boundary (r=1, u=0) and wrap in theta
r_ext     = [1.0; r_pos];
theta_ext = [theta; 2*pi];
[RR, TT]  = ndgrid(r_ext, theta_ext);
XX = RR .* cos(TT);
YY = RR .* sin(TT);

figure('Position', [100 100 1100 720]);

for panel = 1:6
    mi = mode_indices(panel);
    ax = subplot(2, 3, panel);

    % Reshape eigenvector: ordered as
    % [u(r_1,th_1), u(r_1,th_2), ..., u(r_1,th_M), u(r_2,th_1), ...]
    v = real(V(:, mi));
    U = reshape(v, [N2, M]);

    % Normalise
    U = U / max(abs(U(:)));

    % Extend: boundary row of zeros at r=1 (top) and periodic wrap in theta
    U_ext = zeros(N2 + 1, M + 1);
    U_ext(2:end, 1:M) = U;
    U_ext(:, M+1)     = U_ext(:, 1);   % periodic wrap

    surf(XX, YY, U_ext, 'EdgeColor', 'k', 'LineWidth', 0.3, ...
        'FaceAlpha', 0.85);
    colormap(ax, rdbu_colormap());
    caxis([-1 1]);
    zlim([-1.1 1.1]);
    xlim([-1.1 1.1]);
    ylim([-1.1 1.1]);
    view(25, -60);
    axis off;

    title(sprintf('Mode %d\n\\lambda = %.10f', mi, eigenvalues(mi)), 'FontSize', 9);
end

sgtitle('Eigenmodes of the Laplacian on the unit disk', 'FontSize', 12);

exportgraphics(gcf, fullfile(fig_dir, 'disk_eigenmodes.pdf'), ...
    'ContentType', 'vector');

% -----------------------------------------------------------------------
% Figure 12.4: Eigenvalue convergence vs Bessel zeros
% -----------------------------------------------------------------------

figure('Position', [100 100 900 380]);

% --- Left: computed vs exact eigenvalues ---
subplot(1, 2, 1);

n_show = 25;
lam_computed = eigenvalues(1:n_show);

% Collect exact Bessel eigenvalues and sort
% besselj_zeros: compute zeros of J_m using fzero
bessel_lam = [];
for m = 0:9
    z = bessel_zeros(m, 5);
    if m == 0
        mult = 1;
    else
        mult = 2;
    end
    for k = 1:length(z)
        bessel_lam = [bessel_lam; repmat(z(k)^2, mult, 1)]; %#ok<AGROW>
    end
end
bessel_lam = sort(bessel_lam);
bessel_lam = bessel_lam(1:n_show);

plot(1:n_show, lam_computed, 'o', 'MarkerSize', 6, ...
    'MarkerFaceColor', 'none', 'MarkerEdgeColor', [0.12 0.47 0.71]);
hold on;
plot(1:n_show, bessel_lam, 'x', 'MarkerSize', 5, ...
    'Color', [0.84 0.15 0.16]);
hold off;
xlabel('Mode index');
ylabel('Eigenvalue $\lambda$', 'Interpreter', 'latex');
legend('Computed', 'Bessel zeros $j_{m,n}^2$', ...
    'FontSize', 9, 'Interpreter', 'latex');
title('Eigenvalue comparison', 'FontSize', 10);

% --- Right: convergence of first eigenvalue with Nr ---
subplot(1, 2, 2);

Nr_values = [9, 13, 17, 21, 25, 29, 33];
M_fixed   = 20;
lam_exact_1 = bessel_zeros(0, 1)^2;   % j_{0,1}^2

errors = zeros(size(Nr_values));
for k = 1:length(Nr_values)
    N = Nr_values(k);
    [L_test, ~, ~] = laplacian_polar(N, M_fixed);
    eigs_test = sort(real(eig(-L_test)));
    errors(k) = abs(eigs_test(1) - lam_exact_1);
end

semilogy(Nr_values, errors, 'o-', 'MarkerSize', 5, ...
    'Color', [0.12 0.47 0.71]);
xlabel('$N_r$ (radial points)', 'Interpreter', 'latex');
ylabel('$|\lambda_1 - j_{0,1}^2|$', 'Interpreter', 'latex');
title('Exponential convergence', 'FontSize', 10);
grid on;
set(gca, 'GridAlpha', 0.3);

exportgraphics(gcf, fullfile(fig_dir, 'eigenvalue_convergence.pdf'), ...
    'ContentType', 'vector');

% Print eigenvalue comparison table
fprintf('%5s  %16s  %16s  %12s\n', 'Mode', 'Computed', 'Exact (Bessel)', 'Error');
fprintf('%s\n', repmat('-', 1, 58));
for i = 1:min(15, n_show)
    err = abs(lam_computed(i) - bessel_lam(i));
    fprintf('%5d  %16.10f  %16.10f  %12.2e\n', i, lam_computed(i), bessel_lam(i), err);
end

fprintf('\nFigures saved to: %s\n', fig_dir);
fprintf('  disk_eigenmodes.pdf\n');
fprintf('  eigenvalue_convergence.pdf\n');


% =====================================================================
% Helper functions
% =====================================================================

function z = bessel_zeros(m, n)
% BESSEL_ZEROS  Compute the first n positive zeros of J_m(x).
    z = zeros(n, 1);
    % Initial guesses using McMahon's asymptotic expansion
    for k = 1:n
        if m == 0
            x0 = (k - 0.25) * pi;
        else
            x0 = (k + m/2 - 0.25) * pi;
        end
        % Refine with fzero
        z(k) = fzero(@(x) besselj(m, x), x0);
    end
end

function cmap = rdbu_colormap()
% RDBU_COLORMAP  Red-Blue diverging colormap similar to matplotlib's RdBu_r.
    n = 256;
    % Key color nodes for RdBu_r (blue to white to red)
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
    xi = linspace(0, 1, n);
    cmap = interp1(x_nodes, nodes, xi);
end
