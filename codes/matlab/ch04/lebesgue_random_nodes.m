%% lebesgue_random_nodes.m
%
% Computational experiment: Lebesgue constant for random nodes on [-1, 1].
%
% This script investigates the growth of the Lebesgue constant when
% interpolation nodes are chosen randomly (uniformly distributed) on [-1, 1].
% Through Monte Carlo simulation, we discover the statistical behavior and
% derive an empirical asymptotic formula.
%
% Why the median: two independently drawn nodes can fall arbitrarily close
% to each other, and a single such draw makes the Lebesgue constant
% astronomically large. The M values therefore have a very heavy right
% tail, and their sample mean is dominated by the worst draws: it jumps by
% orders of magnitude from one N to the next. The median and the quartiles
% describe the typical draw and do not feel the tail, so the empirical law
% Lambda_N ~ a * exp(b * N) is fitted to the median, by linear least
% squares on its logarithm.
%
% Key finding: over the degrees computed, the median Lebesgue constant of
% random nodes grows exponentially and faster than that of equispaced
% nodes, so that the special structure of Chebyshev points (endpoint
% clustering) is essential for stable interpolation.
%
% Author: Dr. Denys Dutykh
%         Mathematics Department
%         Khalifa University of Science and Technology
%         Abu Dhabi, UAE
%
% Part of the book "Computational Etudes: A Spectral Approach"
%
% Created: January 2026

clear; close all; clc;

%% Configuration
N_FINE = 2000;          % Number of points for Lebesgue function evaluation
M_TRIALS = 200;         % Number of Monte Carlo trials per N
N_MAX = 30;             % Maximum polynomial degree
SEED = 42;              % Random seed for reproducibility

% Book color scheme
NAVY = [20, 45, 110] / 255;
SKY = [120, 150, 210] / 255;
CORAL = [231, 76, 60] / 255;
TEAL = [22, 160, 133] / 255;
GOLD = [243, 156, 18] / 255;    % New color for random nodes

% Output path
script_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch04', 'matlab');
output_file = fullfile(output_dir, 'lebesgue_random_nodes.pdf');

%% Set random seed
rng(SEED);

%% Fine grid for Lebesgue function evaluation
x_fine = linspace(-1, 1, N_FINE)';

%% Initialize arrays
N_values = 2:N_MAX;
n_N = length(N_values);

Lambda_mean = zeros(n_N, 1);
Lambda_min = zeros(n_N, 1);
Lambda_max = zeros(n_N, 1);
Lambda_median = zeros(n_N, 1);
Lambda_p25 = zeros(n_N, 1);
Lambda_p75 = zeros(n_N, 1);
Lambda_equi = zeros(n_N, 1);
Lambda_cheb = zeros(n_N, 1);

%% Run Monte Carlo simulation
fprintf('=======================================================================\n');
fprintf('Computational Experiment: Lebesgue Constant for Random Nodes\n');
fprintf('=======================================================================\n\n');
fprintf('Running Monte Carlo simulation with M = %d trials per N...\n', M_TRIALS);
fprintf('N ranges from 2 to %d\n\n', N_MAX);

for i = 1:n_N
    N = N_values(i);

    % Monte Carlo samples
    samples = zeros(M_TRIALS, 1);
    for m = 1:M_TRIALS
        x_rand = sort(2 * rand(N+1, 1) - 1);  % Uniform on [-1, 1]
        samples(m) = max(lebesgue_function(x_rand, x_fine));
    end

    % Statistics
    Lambda_mean(i) = mean(samples);
    Lambda_min(i) = min(samples);
    Lambda_max(i) = max(samples);
    Lambda_median(i) = median(samples);
    Lambda_p25(i) = prctile(samples, 25);
    Lambda_p75(i) = prctile(samples, 75);

    % Deterministic distributions
    x_equi = linspace(-1, 1, N+1)';
    Lambda_equi(i) = max(lebesgue_function(x_equi, x_fine));

    j = 0:N;
    x_cheb = cos(j * pi / N)';
    Lambda_cheb(i) = max(lebesgue_function(x_cheb, x_fine));

    % Progress report
    if mod(i, 5) == 0 || N == N_MAX
        fprintf('  N = %2d: Lambda_median = %10.4e, IQR = [%.2e, %.2e], range = [%.2e, %.2e]\n', ...
                N, Lambda_median(i), Lambda_p25(i), Lambda_p75(i), Lambda_min(i), Lambda_max(i));
    end
end

%% Fit the empirical exponential law to the medians
fprintf('\n-----------------------------------------------------------------------\n');
fprintf('Fitting the empirical law Lambda_N ~ a * exp(b * N) to the medians...\n');

% Use data for N >= 5 to avoid small-N effects
mask = N_values >= 5;
N_fit_data = N_values(mask)';

% Least-squares straight line through the points (N, ln Lambda_N): every
% decade of Lambda_N weighs the same and the fitted curve is positive
[a_fit, b_fit, r2_fit] = fit_exponential(N_fit_data, Lambda_median(mask));

% The same straight-line fit for equispaced nodes, for comparison
[a_equi, b_equi, r2_equi] = fit_exponential(N_fit_data, Lambda_equi(mask));

fprintf('\nRandom nodes (median): Lambda_N ~ %.4f * exp(%.4f * N)  (R^2 = %.4f in the (N, ln Lambda_N) plane)\n', ...
        a_fit, b_fit, r2_fit);
fprintf('Equispaced nodes:      Lambda_N ~ %.4f * exp(%.4f * N)  (R^2 = %.4f in the (N, ln Lambda_N) plane)\n', ...
        a_equi, b_equi, r2_equi);

%% Create figure
fig = figure('Units', 'inches', 'Position', [1, 1, 11, 4.5]);

% -------------------------------------------------------------------------
% Left panel: Growth of Lebesgue constant with confidence bands
% -------------------------------------------------------------------------
subplot(1, 2, 1);
hold on; box on;
set(gca, 'YScale', 'log');  % Force logarithmic scale before fill()

% Shaded region for random nodes (min to max)
fill([N_values, fliplr(N_values)], [Lambda_min', fliplr(Lambda_max')], ...
     GOLD, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', 'Random (min-max)');

% Shaded region for interquartile range
fill([N_values, fliplr(N_values)], [Lambda_p25', fliplr(Lambda_p75')], ...
     GOLD, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Random (IQR)');

% Median line for random nodes
h1 = semilogy(N_values, Lambda_median, 'o-', 'Color', GOLD, 'LineWidth', 1.5, ...
              'MarkerSize', 4, 'DisplayName', 'Random (median)');

% Equispaced for comparison
h2 = semilogy(N_values, Lambda_equi, 's--', 'Color', CORAL, 'LineWidth', 1.2, ...
              'MarkerSize', 3, 'DisplayName', 'Equispaced');
h2.Color(4) = 0.7;

% Chebyshev for reference
h3 = semilogy(N_values, Lambda_cheb, '^--', 'Color', SKY, 'LineWidth', 1.2, ...
              'MarkerSize', 3, 'DisplayName', 'Chebyshev');
h3.Color(4) = 0.7;

% Exponential law fitted to the median
N_curve = linspace(5, N_MAX, 100);
Lambda_curve = a_fit * exp(b_fit * N_curve);
fit_label = sprintf('$\\Lambda_N \\approx %.3g\\,e^{%.3fN}$ (fit to median)', a_fit, b_fit);
h4 = semilogy(N_curve, Lambda_curve, '-', 'Color', NAVY, 'LineWidth', 1.5, ...
              'DisplayName', fit_label);
h4.Color(4) = 0.8;

xlabel('$N$ (polynomial degree)', 'Interpreter', 'latex');
ylabel('Lebesgue constant $\Lambda_N$', 'Interpreter', 'latex');
title('Growth of Lebesgue Constant', 'FontSize', 11);
xlim([0, N_MAX + 2]);
legend([h1, h2, h3, h4], {'Random (median)', 'Equispaced', 'Chebyshev', fit_label}, ...
       'Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 8);

% Growth law of the Chebyshev constant
text(N_MAX - 5, Lambda_cheb(end) * 1.5, '$O(\ln N)$', ...
     'Interpreter', 'latex', 'FontSize', 10, 'Color', SKY);

set(gca, 'FontSize', 10, 'TickLabelInterpreter', 'latex');
grid on;
set(gca, 'GridAlpha', 0.2);

% -------------------------------------------------------------------------
% Right panel: Distribution at fixed N (log-scale histogram)
% -------------------------------------------------------------------------
subplot(1, 2, 2);
hold on; box on;

N_hist = 15;
M_hist = 500;
samples_hist = zeros(M_hist, 1);
for m = 1:M_hist
    x_rand = sort(2 * rand(N_hist+1, 1) - 1);
    samples_hist(m) = max(lebesgue_function(x_rand, x_fine));
end

% Compute statistics
mean_val = mean(samples_hist);
median_val = median(samples_hist);
equi_val = max(lebesgue_function(linspace(-1, 1, N_hist+1)', x_fine));
j = 0:N_hist;
cheb_val = max(lebesgue_function(cos(j * pi / N_hist)', x_fine));

% Use log-scale histogram to show the spread across orders of magnitude
log_samples = log10(samples_hist);

% Create histogram with log-spaced bins
bins = linspace(floor(min(log_samples)), ceil(max(log_samples)), 35);
h_hist = histogram(log_samples, bins, 'FaceColor', GOLD, 'FaceAlpha', 0.7, 'EdgeColor', 'white');

% The bars fill the lower 62 % of the panel and the vertical lines stop at
% 72 %, which leaves the top of the panel to the legend
y_top = max(h_hist.Values) / 0.62;
ylim([0, y_top]);
xlim([0, bins(end) + 0.5]);
y_line = 0.72 * y_top;

% Add vertical lines for statistics (on log scale)
l1 = plot(log10(median_val) * [1, 1], [0, y_line], '--', 'Color', TEAL, 'LineWidth', 2, ...
          'DisplayName', sprintf('Median = $10^{%.1f}$', log10(median_val)));
l2 = plot(log10(mean_val) * [1, 1], [0, y_line], '-', 'Color', NAVY, 'LineWidth', 2, ...
          'DisplayName', sprintf('Mean = $10^{%.1f}$', log10(mean_val)));
l3 = plot(log10(equi_val) * [1, 1], [0, y_line], ':', 'Color', CORAL, 'LineWidth', 2, ...
          'DisplayName', sprintf('Equispaced = $10^{%.1f}$', log10(equi_val)));
l4 = plot(log10(cheb_val) * [1, 1], [0, y_line], '-.', 'Color', SKY, 'LineWidth', 2, ...
          'DisplayName', sprintf('Chebyshev = %.1f', cheb_val));

xlabel('$\log_{10}(\Lambda_{15})$', 'Interpreter', 'latex');
ylabel('Count', 'Interpreter', 'latex');
title(sprintf('Distribution of $\\Lambda_N$ for $N = %d$ (%d trials)', N_hist, M_hist), ...
      'FontSize', 11, 'Interpreter', 'latex');

% Legend, titled with the spread of the samples
spread = log10(max(samples_hist)) - log10(min(samples_hist));
lgd = legend([l1, l2, l3, l4], 'Location', 'north', 'Interpreter', 'latex', ...
             'FontSize', 8, 'NumColumns', 2);
title(lgd, sprintf('Spread of the samples: %.1f orders of magnitude', spread), ...
      'Interpreter', 'latex', 'FontWeight', 'normal');

set(gca, 'FontSize', 10, 'TickLabelInterpreter', 'latex');
grid on;
set(gca, 'GridAlpha', 0.2);

%% Save figure
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

exportgraphics(fig, output_file, 'ContentType', 'vector');
fprintf('\nFigure saved to: %s\n', output_file);

png_file = strrep(output_file, '.pdf', '.png');
exportgraphics(fig, png_file, 'Resolution', 300);

%% Print summary table
fprintf('\n==============================================================================\n');
fprintf('Summary Statistics (medians over the draws; the mean is shown for contrast)\n');
fprintf('==============================================================================\n');
fprintf('%4s %12s %12s %12s %12s %12s\n', 'N', 'Median', '(mean)', 'Min', 'Max', 'Equispaced');
fprintf('------------------------------------------------------------------------------\n');
for i = 1:3:n_N
    N = N_values(i);
    fprintf('%4d %12.4e %12.4e %12.4e %12.4e %12.4e\n', ...
            N, Lambda_median(i), Lambda_mean(i), Lambda_min(i), Lambda_max(i), Lambda_equi(i));
end
fprintf('------------------------------------------------------------------------------\n');

fprintf('\n==============================================================================\n');
fprintf('EMPIRICAL GROWTH LAW OF THE MEDIAN\n');
fprintf('==============================================================================\n');
fprintf('\nFor random nodes uniformly distributed on [-1, 1], N from %d to %d:\n', ...
        N_fit_data(1), N_MAX);
fprintf('\n    Lambda_N ~ %.4f * exp(%.4f * N) = %.4f * %.4f^N\n', a_fit, b_fit, a_fit, exp(b_fit));
fprintf('\nFor equispaced nodes, the same fit gives Lambda_N ~ %.4f * %.4f^N.\n', a_equi, exp(b_equi));
fprintf('For Chebyshev nodes, Lambda_%d = %.4f.\n', N_MAX, Lambda_cheb(end));
fprintf('==============================================================================\n');

%% Helper function: exponential law fitted in the (N, ln Lambda_N) plane
function [a, b, r_squared] = fit_exponential(N_values, Lambda_values)
    log_Lambda = log(Lambda_values(:));
    p = polyfit(N_values(:), log_Lambda, 1);
    b = p(1);
    a = exp(p(2));
    residual = log_Lambda - polyval(p, N_values(:));
    r_squared = 1 - sum(residual.^2) / sum((log_Lambda - mean(log_Lambda)).^2);
end

%% Helper function: Lebesgue function
function Lambda = lebesgue_function(x_nodes, x_eval)
    N = length(x_nodes) - 1;
    Lambda = zeros(size(x_eval));

    for k = 1:N+1
        L_k = ones(size(x_eval));
        for j = 1:N+1
            if j ~= k
                L_k = L_k .* (x_eval - x_nodes(j)) / (x_nodes(k) - x_nodes(j));
            end
        end
        Lambda = Lambda + abs(L_k);
    end
end
