%% convergence_comparison.m
%
% Compares the convergence of polynomial interpolation on equispaced vs
% Chebyshev nodes for the Runge function.
%
% For smooth functions with singularities at distance d from [-1,1]:
%   - Chebyshev: error ~ O(rho^{-N}) where rho > 1
%   - Equispaced: error can GROW with N
%
% The Runge function has poles at ±0.2i, giving a critical rho ≈ 1.04.
%
% Author: Dr. Denys Dutykh
%         Mathematics Department
%         Khalifa University of Science and Technology
%         Abu Dhabi, UAE
%
% Part of the book "Computational Études: A Spectral Approach"

clear; close all; clc;

%% Configuration
N_FINE = 2000;

% Book color scheme
NAVY = [20, 45, 110] / 255;
SKY = [120, 150, 210] / 255;
CORAL = [231, 76, 60] / 255;
TEAL = [22, 160, 133] / 255;

% Output path
script_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch04', 'matlab');
output_file = fullfile(output_dir, 'convergence_comparison.pdf');

%% The Runge function
runge = @(x) 1 ./ (1 + 25*x.^2);

%% Evaluation grid
x_fine = linspace(-1, 1, N_FINE)';

%% Range of polynomial degrees
N_values = 2:2:50;

errors_equi = zeros(size(N_values));
errors_cheb = zeros(size(N_values));

for i = 1:length(N_values)
    N = N_values(i);

    % Equispaced
    x_equi = linspace(-1, 1, N+1)';
    f_equi = runge(x_equi);
    p_equi = lagrange_interp(x_equi, f_equi, x_fine);
    errors_equi(i) = max(abs(runge(x_fine) - p_equi));

    % Chebyshev
    j = 0:N;
    x_cheb = cos(j * pi / N)';
    f_cheb = runge(x_cheb);
    p_cheb = lagrange_interp(x_cheb, f_cheb, x_fine);
    errors_cheb(i) = max(abs(runge(x_fine) - p_cheb));
end

%% Theoretical convergence rate
% Runge function poles at ±0.2i
rho = abs(0.2i + sqrt(-1.04));

%% Create figure
fig = figure('Units', 'inches', 'Position', [1, 1, 8, 5.5]);
hold on; box on;

% Plot errors
semilogy(N_values, errors_equi, 'o-', 'Color', CORAL, 'LineWidth', 1.5, ...
         'MarkerSize', 5, 'DisplayName', 'Equispaced nodes');
semilogy(N_values, errors_cheb, 's-', 'Color', TEAL, 'LineWidth', 1.5, ...
         'MarkerSize', 5, 'DisplayName', 'Chebyshev nodes');

% Theoretical rate for Chebyshev
idx_fit = N_values >= 20;
C_fit = mean(errors_cheb(idx_fit) .* rho.^N_values(idx_fit));
N_theory = linspace(10, 50, 100);
semilogy(N_theory, C_fit * rho.^(-N_theory), '--', 'Color', TEAL, 'LineWidth', 1.5, ...
         'DisplayName', sprintf('Theory: $O(\\rho^{-N})$, $\\rho \\approx %.2f$', rho));

% Reference line
yline(1, ':', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5, ...
      'DisplayName', 'Error = 1 (useless)');

xlabel('Polynomial degree $N$', 'Interpreter', 'latex');
ylabel('Maximum interpolation error $\|f - p_N\|_\infty$', 'Interpreter', 'latex');
title('Convergence: Equispaced vs Chebyshev Interpolation', 'FontSize', 11);
xlim([0, 52]);
ylim([1e-15, 1e4]);
legend('Location', 'northeast', 'Interpreter', 'latex');

% Annotations
text(30, 1e3, 'Divergence!', 'FontSize', 10, 'Color', CORAL);
text(35, 1e-8, 'Convergence', 'FontSize', 10, 'Color', TEAL);

% Formula box
text(0.05, 0.05, '$f(x) = \frac{1}{1 + 25x^2}$', 'Units', 'normalized', ...
     'Interpreter', 'latex', 'FontSize', 12, ...
     'BackgroundColor', 'white', 'EdgeColor', NAVY, 'LineWidth', 0.5);

set(gca, 'FontSize', 10, 'TickLabelInterpreter', 'latex');
grid on;
set(gca, 'GridAlpha', 0.2);

%% Save figure
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

exportgraphics(fig, output_file, 'ContentType', 'vector');
fprintf('Figure saved to: %s\n', output_file);

png_file = strrep(output_file, '.pdf', '.png');
exportgraphics(fig, png_file, 'Resolution', 300);

%% Print summary
fprintf('\nConvergence Comparison (Runge function):\n');
fprintf('%s\n', repmat('-', 1, 55));
fprintf('%4s %15s %15s %12s\n', 'N', 'Equispaced', 'Chebyshev', 'Ratio');
fprintf('%s\n', repmat('-', 1, 55));
for i = 1:5:length(N_values)
    N = N_values(i);
    err_e = errors_equi(i);
    err_c = errors_cheb(i);
    ratio = err_e / err_c;
    fprintf('%4d %15.2e %15.2e %12.1e\n', N, err_e, err_c, ratio);
end
fprintf('%s\n', repmat('-', 1, 55));

fprintf('\nAnalysis:\n');
fprintf('  Runge function poles at z = +/- 0.2i\n');
fprintf('  Chebyshev convergence parameter: rho = %.4f\n', rho);
fprintf('  Expected rate: error ~ rho^{-N} = %.4f^N\n', 1/rho);

% Estimate rho from data
log_errors = log(errors_cheb(10:end));
N_fit = N_values(10:end);
p = polyfit(N_fit, log_errors, 1);
rho_estimated = exp(-p(1));
fprintf('  Estimated from data: rho ~ %.4f\n', rho_estimated);

%% Lagrange interpolation function
function p = lagrange_interp(x_nodes, f_nodes, x_eval)
    n = length(x_nodes);
    p = zeros(size(x_eval));

    for k = 1:n
        L_k = ones(size(x_eval));
        for j = 1:n
            if j ~= k
                L_k = L_k .* (x_eval - x_nodes(j)) / (x_nodes(k) - x_nodes(j));
            end
        end
        p = p + f_nodes(k) * L_k;
    end
end
