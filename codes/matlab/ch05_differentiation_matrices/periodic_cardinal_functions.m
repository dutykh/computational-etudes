%% periodic_cardinal_functions.m
%
% Visualizes periodic cardinal functions (discrete Dirichlet kernel) for
% equispaced nodes on a periodic domain.
%
% The periodic cardinal function is:
%     phi_j(x) = sin(N(x - x_j)/2) / (N sin((x - x_j)/2))
%
% These functions satisfy phi_j(x_k) = delta_{jk} (Kronecker delta), making
% them the "cardinal functions" for trigonometric interpolation on periodic
% domains.
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
% Part of "Computational Etudes: A Spectral Approach"
% https://github.com/dutykh/computational-etudes

clear; close all; clc;

%% Configuration
N = 16;           % Number of grid points (should be even)
N_FINE = 1000;    % Number of points for smooth plotting

% Book color scheme
NAVY = [0.078 0.176 0.431];    % #142D6E
SKY = [0.471 0.588 0.824];     % #7896D2
CORAL = [0.906 0.298 0.235];   % #E74C3C
TEAL = [0.086 0.627 0.522];    % #16A085
PURPLE = [0.557 0.267 0.678];  % #8E44AD
ORANGE = [0.902 0.494 0.133];  % #E67E22

% Output path
script_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch05', 'matlab');
output_file = fullfile(output_dir, 'periodic_cardinal_functions.pdf');

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% Grid points
h = 2 * pi / N;
x_nodes = h * (0:N-1)';

% Fine grid for plotting
x_fine = linspace(0, 2*pi, N_FINE)';

%% Create figure
fig = figure('Units', 'inches', 'Position', [1, 1, 9, 4.5]);

% Indices to plot (spread across domain)
indices_to_plot = [0, 4, 8, 12];  % 0-indexed for labels
colors = {CORAL, TEAL, SKY, PURPLE};

% Plot cardinal functions
hold on;
for ii = 1:length(indices_to_plot)
    j = indices_to_plot(ii);
    phi_j = periodic_cardinal(x_fine, x_nodes(j+1), N);
    plot(x_fine, phi_j, 'Color', colors{ii}, 'LineWidth', 1.5, ...
         'DisplayName', sprintf('$\\phi_{%d}(x)$', j));
end

% Mark all nodes on x-axis
plot(x_nodes, zeros(size(x_nodes)), 'o', 'Color', NAVY, 'MarkerSize', 5, ...
     'MarkerFaceColor', NAVY, 'DisplayName', 'Nodes');

% Highlight peak values
for ii = 1:length(indices_to_plot)
    j = indices_to_plot(ii);
    plot(x_nodes(j+1), 1.0, 'o', 'Color', colors{ii}, 'MarkerSize', 8, ...
         'MarkerFaceColor', colors{ii}, 'MarkerEdgeColor', 'white', ...
         'LineWidth', 1.5, 'HandleVisibility', 'off');
end

% Reference lines
yline(0, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
yline(1, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5, 'Alpha', 0.5);

% Axis labels and title
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 11);
ylabel('$\phi_j(x)$', 'Interpreter', 'latex', 'FontSize', 11);
title(sprintf('Periodic Cardinal Functions ($N = %d$)', N), ...
      'Interpreter', 'latex', 'FontSize', 12);

% Set x-axis ticks at multiples of pi/2
set(gca, 'XTick', [0, pi/2, pi, 3*pi/2, 2*pi]);
set(gca, 'XTickLabel', {'$0$', '$\pi/2$', '$\pi$', '$3\pi/2$', '$2\pi$'}, ...
    'TickLabelInterpreter', 'latex');

xlim([0, 2*pi]);
ylim([-0.25, 1.15]);

% Legend
leg = legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 9);
set(leg, 'NumColumns', 2);

% Add annotation about cardinal property
text(x_nodes(5) + 0.8, 0.85, '$\phi_j(x_k) = \delta_{jk}$', ...
     'Interpreter', 'latex', 'FontSize', 10, 'Color', [0.5 0.5 0.5]);

% Clean styling
box off;
grid on;
set(gca, 'GridLineStyle', '-', 'GridAlpha', 0.2, 'GridColor', [0.5 0.5 0.5]);
set(gca, 'FontSize', 10);

hold off;

%% Save figure
exportgraphics(fig, output_file, 'ContentType', 'vector');
exportgraphics(fig, strrep(output_file, '.pdf', '.png'), 'Resolution', 300);

fprintf('Figure saved to: %s\n', output_file);

%% Verify cardinal property
fprintf('\nVerification of cardinal property for N = %d:\n', N);
fprintf('--------------------------------------------------\n');
fprintf('%4s %12s %25s\n', 'j', 'phi_j(x_j)', 'max |phi_j(x_k)|, k!=j');
fprintf('--------------------------------------------------\n');
for j = 0:N-1
    phi_at_nodes = periodic_cardinal(x_nodes, x_nodes(j+1), N);
    phi_at_own = phi_at_nodes(j+1);
    phi_at_others = phi_at_nodes([1:j, j+2:N]);
    max_at_others = max(abs(phi_at_others));
    fprintf('%4d %12.6f %25.2e\n', j, phi_at_own, max_at_others);
end
fprintf('--------------------------------------------------\n');

close(fig);

%% Helper function: Periodic cardinal function
function phi = periodic_cardinal(x, x_j, N)
    % Compute the periodic cardinal function phi_j(x) centered at x_j.
    %
    % phi_j(x) = sin(N(x - x_j)/2) / (N sin((x - x_j)/2))
    %
    % Inputs:
    %   x   - Points at which to evaluate (column vector)
    %   x_j - Center point (node location)
    %   N   - Number of grid points
    %
    % Output:
    %   phi - Values of phi_j at x

    theta = (x - x_j) / 2.0;

    % Handle singularity at theta = 0 (x = x_j)
    phi = zeros(size(x));

    small = abs(sin(theta)) < 1e-14;

    % Non-singular points
    phi(~small) = sin(N * theta(~small)) ./ (N * sin(theta(~small)));

    % Singular points (at x_j and periodic copies)
    phi(small) = 1.0;
end
