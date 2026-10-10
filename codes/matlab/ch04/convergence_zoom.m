%% convergence_zoom.m
%
% Zoomed view of Chebyshev interpolation convergence for the Runge function,
% without the divergent equispaced curve.
%
% This allows better visualization of the geometric convergence rate O(rho^{-N})
% where rho ~ 1.22 for the Runge function.
%
% The Runge function f(x) = 1/(1 + 25x^2) has poles at +/-0.2i, giving:
%     rho = |0.2i + sqrt((0.2i)^2 - 1)| ~ 1.22
%
% Author: Dr. Denys Dutykh
%         Mathematics Department
%         Khalifa University of Science and Technology
%         Abu Dhabi, UAE
%
% Part of the book "Computational Etudes: A Spectral Approach"

clear; close all; clc;

%% Configuration
N_FINE = 2000;

% Book color scheme
NAVY = [20, 45, 110] / 255;
TEAL = [22, 160, 133] / 255;

% Output path
script_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch04', 'matlab');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

%% The Runge function
runge = @(x) 1 ./ (1 + 25*x.^2);

%% Evaluation grid
x_fine = linspace(-1, 1, N_FINE)';
f_fine = runge(x_fine);

%% Compute errors for Chebyshev interpolation
N_values = 2:2:70;
errors_cheb = zeros(size(N_values));

for idx = 1:length(N_values)
    N = N_values(idx);
    % Chebyshev-Gauss-Lobatto nodes
    j = (0:N)';
    x_cheb = cos(j * pi / N);
    f_cheb = runge(x_cheb);
    p_cheb = lagrange_interpolate(x_cheb, f_cheb, x_fine);
    errors_cheb(idx) = max(abs(f_fine - p_cheb));
end

%% Theoretical convergence parameter
rho = abs(0.2i + sqrt(-1.04 + 0i));

%% Fit constant using data from N >= 20
idx_fit = N_values >= 20;
C_fit = mean(errors_cheb(idx_fit) .* rho.^N_values(idx_fit));

%% Create figure
fig = figure('Position', [100, 100, 700, 500]);

semilogy(N_values, errors_cheb, 's-', 'Color', TEAL, 'LineWidth', 1.5, ...
         'MarkerSize', 5, 'MarkerFaceColor', TEAL, ...
         'DisplayName', 'Chebyshev interpolation');
hold on;

% Theoretical convergence rate
N_theory = linspace(10, 70, 100);
semilogy(N_theory, C_fit * rho.^(-N_theory), '--', 'Color', NAVY, ...
         'LineWidth', 2, ...
         'DisplayName', sprintf('Theory: C \\cdot \\rho^{-N}, \\rho \\approx %.2f', rho));

% Machine precision
yline(1e-15, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8, ...
      'DisplayName', 'Machine precision');
hold off;

xlabel('Polynomial degree $N$', 'Interpreter', 'latex');
ylabel('Maximum interpolation error $\|f - p_N\|_\infty$', 'Interpreter', 'latex');
title('Chebyshev Interpolation Convergence (Runge Function)', 'FontSize', 11);
xlim([0, 72]);
ylim([1e-16, 1e1]);
legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 10);
grid on; set(gca, 'GridAlpha', 0.2);

% Save figure
exportgraphics(fig, fullfile(output_dir, 'convergence_zoom.pdf'), 'ContentType', 'vector');
exportgraphics(fig, fullfile(output_dir, 'convergence_zoom.png'), 'Resolution', 300);
fprintf('Figure saved to %s\n', fullfile(output_dir, 'convergence_zoom.pdf'));

%% Print summary
fprintf('\nChebyshev Convergence for Runge Function:\n');
fprintf('Theoretical rho = %.6f\n', rho);
log_errors = log(errors_cheb(11:end));
N_fit_vals = N_values(11:end);
p = polyfit(N_fit_vals, log_errors, 1);
rho_est = exp(-p(1));
fprintf('Estimated from data: rho ~ %.6f\n', rho_est);


%% ---- Helper function ----
function p = lagrange_interpolate(x_nodes, f_nodes, x_eval)
    % Lagrange interpolation using barycentric formula (stable version)
    n = length(x_nodes);
    m = length(x_eval);
    p = zeros(m, 1);

    for k = 1:n
        L_k = ones(m, 1);
        for j = 1:n
            if j ~= k
                L_k = L_k .* (x_eval - x_nodes(j)) / (x_nodes(k) - x_nodes(j));
            end
        end
        p = p + f_nodes(k) * L_k;
    end
end
