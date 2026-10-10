%% higher_order_derivatives.m
%
% Demonstrates spectral differentiation for higher-order derivatives (up to 4th
% order) and compares D² construction methods: matrix squaring vs direct formula.
%
% Test function: u(x) = exp(-sin(2x)) on [0, 2π)
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
% Part of "Computational Etudes: A Spectral Approach"
% https://github.com/dutykh/computational-etudes

clear; close all; clc;

%% Configuration
N_DEMO = 32;   % Grid points for demonstration
N_FINE = 500;  % Fine grid for exact solutions

% Book color scheme
NAVY = [0.078 0.176 0.431];
SKY = [0.471 0.588 0.824];
CORAL = [0.906 0.298 0.235];
TEAL = [0.086 0.627 0.522];
PURPLE = [0.557 0.267 0.678];
ORANGE = [0.902 0.494 0.133];

% Output path
script_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch05', 'matlab');

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% Test function and exact derivatives
% u(x) = exp(-sin(2x))
u_func = @(x) exp(-sin(2*x));
u1_exact = @(x) -2*cos(2*x) .* exp(-sin(2*x));
u2_exact = @(x) 4*(sin(2*x) + cos(2*x).^2) .* exp(-sin(2*x));
u3_exact = @(x) 8*cos(2*x) .* (sin(2*x).^2 - 3*sin(2*x)) .* exp(-sin(2*x));
u4_exact = @(x) 16*(-sin(2*x).^3 + 3*sin(2*x).^2 + 5*sin(2*x).*cos(2*x).^2 ...
                    - 3*cos(2*x).^2 - sin(2*x).^2.*cos(2*x).^2) .* exp(-sin(2*x));

%% Construct differentiation matrix
[D, x] = spectral_matrix_periodic(N_DEMO);

% Sample test function
u = u_func(x);

% Compute numerical derivatives via matrix powers
u1_num = D * u;
u2_num = D * D * u;
u3_num = D * D * D * u;
u4_num = D * D * D * D * u;

% Fine grid for exact solutions
x_fine = linspace(0, 2*pi, N_FINE)';

% Compute errors
err1 = max(abs(u1_num - u1_exact(x)));
err2 = max(abs(u2_num - u2_exact(x)));
err3 = max(abs(u3_num - u3_exact(x)));
err4 = max(abs(u4_num - u4_exact(x)));

%% Create main figure (2x2 panels)
fig1 = figure('Units', 'inches', 'Position', [1, 1, 10, 8]);

% Panel 1: First derivative
subplot(2, 2, 1);
plot(x_fine, u1_exact(x_fine), '-', 'Color', NAVY, 'LineWidth', 1.5); hold on;
plot(x, u1_num, 'o', 'Color', CORAL, 'MarkerSize', 5);
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$u''(x)$', 'Interpreter', 'latex');
title('First derivative', 'FontSize', 11);
legend({'Exact', 'Spectral'}, 'Location', 'best');
set(gca, 'XTick', [0, pi/2, pi, 3*pi/2, 2*pi]);
set(gca, 'XTickLabel', {'$0$', '$\pi/2$', '$\pi$', '$3\pi/2$', '$2\pi$'}, ...
    'TickLabelInterpreter', 'latex');
xlim([0, 2*pi]);
text(0.95, 0.05, sprintf('Max error: %.2e', err1), 'Units', 'normalized', ...
    'HorizontalAlignment', 'right', 'FontSize', 9, 'BackgroundColor', 'white');
grid on; box off;

% Panel 2: Second derivative
subplot(2, 2, 2);
plot(x_fine, u2_exact(x_fine), '-', 'Color', NAVY, 'LineWidth', 1.5); hold on;
plot(x, u2_num, 'o', 'Color', TEAL, 'MarkerSize', 5);
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$u''''(x)$', 'Interpreter', 'latex');
title('Second derivative', 'FontSize', 11);
legend({'Exact', 'Spectral'}, 'Location', 'best');
set(gca, 'XTick', [0, pi/2, pi, 3*pi/2, 2*pi]);
set(gca, 'XTickLabel', {'$0$', '$\pi/2$', '$\pi$', '$3\pi/2$', '$2\pi$'}, ...
    'TickLabelInterpreter', 'latex');
xlim([0, 2*pi]);
text(0.95, 0.05, sprintf('Max error: %.2e', err2), 'Units', 'normalized', ...
    'HorizontalAlignment', 'right', 'FontSize', 9, 'BackgroundColor', 'white');
grid on; box off;

% Panel 3: Third derivative
subplot(2, 2, 3);
plot(x_fine, u3_exact(x_fine), '-', 'Color', NAVY, 'LineWidth', 1.5); hold on;
plot(x, u3_num, 'o', 'Color', PURPLE, 'MarkerSize', 5);
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$u''''''(x)$', 'Interpreter', 'latex');
title('Third derivative', 'FontSize', 11);
legend({'Exact', 'Spectral'}, 'Location', 'best');
set(gca, 'XTick', [0, pi/2, pi, 3*pi/2, 2*pi]);
set(gca, 'XTickLabel', {'$0$', '$\pi/2$', '$\pi$', '$3\pi/2$', '$2\pi$'}, ...
    'TickLabelInterpreter', 'latex');
xlim([0, 2*pi]);
text(0.95, 0.05, sprintf('Max error: %.2e', err3), 'Units', 'normalized', ...
    'HorizontalAlignment', 'right', 'FontSize', 9, 'BackgroundColor', 'white');
grid on; box off;

% Panel 4: Fourth derivative
subplot(2, 2, 4);
plot(x_fine, u4_exact(x_fine), '-', 'Color', NAVY, 'LineWidth', 1.5); hold on;
plot(x, u4_num, 'o', 'Color', ORANGE, 'MarkerSize', 5);
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$u''''''''(x)$', 'Interpreter', 'latex');
title('Fourth derivative', 'FontSize', 11);
legend({'Exact', 'Spectral'}, 'Location', 'best');
set(gca, 'XTick', [0, pi/2, pi, 3*pi/2, 2*pi]);
set(gca, 'XTickLabel', {'$0$', '$\pi/2$', '$\pi$', '$3\pi/2$', '$2\pi$'}, ...
    'TickLabelInterpreter', 'latex');
xlim([0, 2*pi]);
text(0.95, 0.05, sprintf('Max error: %.2e', err4), 'Units', 'normalized', ...
    'HorizontalAlignment', 'right', 'FontSize', 9, 'BackgroundColor', 'white');
grid on; box off;

sgtitle(sprintf('Higher-Order Spectral Derivatives of $u(x) = e^{-\\sin(2x)}$ ($N = %d$)', N_DEMO), ...
    'Interpreter', 'latex', 'FontSize', 12);

% Save figure 1
output_file1 = fullfile(output_dir, 'higher_order_derivatives.pdf');
exportgraphics(fig1, output_file1, 'ContentType', 'vector');
exportgraphics(fig1, strrep(output_file1, '.pdf', '.png'), 'Resolution', 300);
fprintf('Figure saved to: %s\n', output_file1);

%% Create D² matrix squaring demonstration
N = 16;  % Grid size for visualization
[D, x_d2] = spectral_matrix_periodic(N);
D2 = D * D;

% Test second derivative accuracy
u_d2 = u_func(x_d2);
u2_num_d2 = D2 * u_d2;
u2_exact_d2 = u2_exact(x_d2);
d2_error = max(abs(u2_num_d2 - u2_exact_d2));

% Fine grid for plotting
x_fine_d2 = linspace(0, 2*pi, 500)';

fig2 = figure('Units', 'inches', 'Position', [1, 1, 12, 4]);

vmax = max(abs(D2(:)));

% Panel 1: D² matrix structure
subplot(1, 3, 1);
imagesc(D2, [-vmax, vmax]);
colormap(gca, redblue(256));
colorbar;
title('$D^2 = D \cdot D$ (matrix structure)', 'Interpreter', 'latex', 'FontSize', 11);
xlabel('Column $j$', 'Interpreter', 'latex');
ylabel('Row $i$', 'Interpreter', 'latex');
axis square;

% Panel 2: Eigenvalue spectrum
subplot(1, 3, 2);
eigvals = eig(D2);
eigvals_sorted = sort(real(eigvals));
k_vals = (-N/2 + 1):(N/2);
theoretical_eigvals = -k_vals.^2;
plot(eigvals_sorted, 'o', 'Color', NAVY, 'MarkerSize', 6); hold on;
plot(sort(theoretical_eigvals), 's', 'Color', CORAL, 'MarkerSize', 4, ...
    'MarkerFaceColor', 'none');
xlabel('Index', 'FontSize', 11);
ylabel('Eigenvalue', 'FontSize', 11);
title('Eigenvalues of $D^2$', 'Interpreter', 'latex', 'FontSize', 11);
legend({'Numerical', 'Theoretical $-k^2$'}, 'Location', 'best');
grid on; box off;

% Panel 3: Second derivative accuracy
subplot(1, 3, 3);
plot(x_fine_d2, u2_exact(x_fine_d2), '-', 'Color', NAVY, 'LineWidth', 1.5); hold on;
plot(x_d2, u2_num_d2, 'o', 'Color', TEAL, 'MarkerSize', 6);
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$u''''(x)$', 'Interpreter', 'latex');
title(sprintf('Second derivative (error = %.2e)', d2_error), 'FontSize', 11);
set(gca, 'XTick', [0, pi/2, pi, 3*pi/2, 2*pi]);
set(gca, 'XTickLabel', {'$0$', '$\pi/2$', '$\pi$', '$3\pi/2$', '$2\pi$'}, ...
    'TickLabelInterpreter', 'latex');
xlim([0, 2*pi]);
legend({'Exact', '$D^2 u$'}, 'Location', 'best', 'Interpreter', 'latex');
grid on; box off;

sgtitle(sprintf('Matrix Squaring for Second Derivatives ($N = %d$)', N), ...
    'Interpreter', 'latex', 'FontSize', 12);

% Save figure 2
output_file2 = fullfile(output_dir, 'd2_comparison.pdf');
exportgraphics(fig2, output_file2, 'ContentType', 'vector');
exportgraphics(fig2, strrep(output_file2, '.pdf', '.png'), 'Resolution', 300);
fprintf('Figure saved to: %s\n', output_file2);

%% Generate convergence table
fprintf('\n');
fprintf('======================================================================\n');
fprintf('Convergence Table: Higher-Order Derivatives\n');
fprintf('Test function: u(x) = exp(-sin(2x))\n');
fprintf('======================================================================\n');
fprintf('%6s %14s %14s %14s %14s\n', 'N', 'Error u''', 'Error u''''', ...
    'Error u''''''', 'Error u''''''''');
fprintf('----------------------------------------------------------------------\n');

N_values = [8, 16, 32, 64];
for N = N_values
    [D, x] = spectral_matrix_periodic(N);
    u = u_func(x);

    u1_num = D * u;
    u2_num = D * D * u;
    u3_num = D * D * D * u;
    u4_num = D * D * D * D * u;

    e1 = max(abs(u1_num - u1_exact(x)));
    e2 = max(abs(u2_num - u2_exact(x)));
    e3 = max(abs(u3_num - u3_exact(x)));
    e4 = max(abs(u4_num - u4_exact(x)));

    fprintf('%6d %14.2e %14.2e %14.2e %14.2e\n', N, e1, e2, e3, e4);
end
fprintf('======================================================================\n');

close all;

%% Helper function: Red-blue colormap
function cmap = redblue(m)
    if nargin < 1, m = 256; end
    r = [linspace(0, 1, m/2), ones(1, m/2)];
    g = [linspace(0, 1, m/2), linspace(1, 0, m/2)];
    b = [ones(1, m/2), linspace(1, 0, m/2)];
    cmap = [r', g', b'];
end
