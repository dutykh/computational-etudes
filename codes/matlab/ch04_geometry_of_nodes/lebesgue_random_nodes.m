%% lebesgue_random_nodes.m
%
% Computational experiment: Lebesgue constant for random nodes on [-1, 1].
%
% This script investigates the growth of the Lebesgue constant when
% interpolation nodes are chosen randomly (uniformly distributed) on [-1, 1].
% Through Monte Carlo simulation, we discover the statistical behavior and
% derive an empirical asymptotic formula.
%
% Key finding: Random nodes exhibit exponential growth similar to equispaced
% nodes, demonstrating that the special structure of Chebyshev points
% (endpoint clustering) is essential for stable interpolation.
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
Lambda_std = zeros(n_N, 1);
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
    Lambda_std(i) = std(samples);
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
        fprintf('  N = %2d: Lambda_mean = %12.4f, std = %10.4f, range = [%.2f, %.2f]\n', ...
                N, Lambda_mean(i), Lambda_std(i), Lambda_min(i), Lambda_max(i));
    end
end

%% Fit empirical formula
fprintf('\n-----------------------------------------------------------------------\n');
fprintf('Fitting empirical asymptotic formula...\n');

% Use data for N >= 5 to avoid small-N effects
mask = N_values >= 5;
N_fit_data = N_values(mask)';
Lambda_fit_data = Lambda_mean(mask);

% Simple exponential fit in log space: log(Lambda) = a + b*N
log_Lambda = log(Lambda_fit_data);
p = polyfit(N_fit_data, log_Lambda, 1);
b_fit = p(1);
a_fit = exp(p(2));

fprintf('\nEmpirical formula: Lambda_N ~ %.4f * exp(%.4f * N)\n', a_fit, b_fit);
fprintf('                 or Lambda_N ~ %.4f * %.4f^N\n', a_fit, exp(b_fit));

% Compare with equispaced
log_equi = log(Lambda_equi(mask));
p_equi = polyfit(N_fit_data, log_equi, 1);
b_equi = p_equi(1);

fprintf('\nFor comparison:\n');
fprintf('  Random nodes growth rate:    exp(%.4f * N) ~ %.4f^N\n', b_fit, exp(b_fit));
fprintf('  Equispaced nodes growth rate: exp(%.4f * N) ~ %.4f^N\n', b_equi, exp(b_equi));

%% Create figure
fig = figure('Units', 'inches', 'Position', [1, 1, 11, 4.5]);

% -------------------------------------------------------------------------
% Left panel: Growth of Lebesgue constant with confidence bands
% -------------------------------------------------------------------------
subplot(1, 2, 1);
hold on; box on;

% Shaded region for random nodes (min to max)
fill([N_values, fliplr(N_values)], [Lambda_min', fliplr(Lambda_max')], ...
     GOLD, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', 'Random (min-max)');

% Shaded region for interquartile range
fill([N_values, fliplr(N_values)], [Lambda_p25', fliplr(Lambda_p75')], ...
     GOLD, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Random (IQR)');

% Mean line for random nodes
h1 = semilogy(N_values, Lambda_mean, 'o-', 'Color', GOLD, 'LineWidth', 1.5, ...
              'MarkerSize', 4, 'DisplayName', 'Random (mean)');

% Equispaced for comparison
h2 = semilogy(N_values, Lambda_equi, 's--', 'Color', CORAL, 'LineWidth', 1.2, ...
              'MarkerSize', 3, 'DisplayName', 'Equispaced');
h2.Color(4) = 0.7;

% Chebyshev for reference
h3 = semilogy(N_values, Lambda_cheb, '^--', 'Color', SKY, 'LineWidth', 1.2, ...
              'MarkerSize', 3, 'DisplayName', 'Chebyshev');
h3.Color(4) = 0.7;

% Fitted curve
N_curve = linspace(5, N_MAX, 100);
Lambda_curve = a_fit * exp(b_fit * N_curve);
h4 = semilogy(N_curve, Lambda_curve, '-', 'Color', NAVY, 'LineWidth', 1.5, ...
              'DisplayName', sprintf('Fit: $\\Lambda_N \\sim %.2f e^{%.3fN}$', a_fit, b_fit));
h4.Color(4) = 0.8;

xlabel('$N$ (polynomial degree)', 'Interpreter', 'latex');
ylabel('Lebesgue constant $\Lambda_N$', 'Interpreter', 'latex');
title('Growth of Lebesgue Constant', 'FontSize', 11);
xlim([0, N_MAX + 2]);
legend([h1, h2, h3, h4], {'Random (mean)', 'Equispaced', 'Chebyshev', ...
        sprintf('Fit: $\\Lambda_N \\sim %.2f e^{%.3fN}$', a_fit, b_fit)}, ...
       'Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 8);

% Add text annotations
text(N_MAX - 5, Lambda_mean(end) * 2, '$O(e^{bN})$', ...
     'Interpreter', 'latex', 'FontSize', 10, 'Color', GOLD);
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
histogram(log_samples, bins, 'FaceColor', GOLD, 'FaceAlpha', 0.7, 'EdgeColor', 'white');

% Add vertical lines for statistics (on log scale)
xline(log10(median_val), '--', 'Color', TEAL, 'LineWidth', 2, ...
      'DisplayName', sprintf('Median = $10^{%.1f}$', log10(median_val)));
xline(log10(mean_val), '-', 'Color', NAVY, 'LineWidth', 2, ...
      'DisplayName', sprintf('Mean = $10^{%.1f}$', log10(mean_val)));
xline(log10(equi_val), ':', 'Color', CORAL, 'LineWidth', 2, ...
      'DisplayName', sprintf('Equispaced = $10^{%.1f}$', log10(equi_val)));
xline(log10(cheb_val), '-.', 'Color', SKY, 'LineWidth', 2, ...
      'DisplayName', sprintf('Chebyshev = %.1f', cheb_val));

xlabel('$\log_{10}(\Lambda_{15})$', 'Interpreter', 'latex');
ylabel('Count', 'Interpreter', 'latex');
title(sprintf('Distribution of $\\Lambda_N$ for $N = %d$ (%d trials)', N_hist, M_hist), ...
      'FontSize', 11, 'Interpreter', 'latex');
legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 8);

% Add annotation about the spread
spread = log10(max(samples_hist)) - log10(min(samples_hist));
text(0.05, 0.95, sprintf('Spread: %.1f orders\nof magnitude', spread), ...
     'Units', 'normalized', 'HorizontalAlignment', 'left', ...
     'VerticalAlignment', 'top', 'FontSize', 9, 'Interpreter', 'latex', ...
     'BackgroundColor', 'white');

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
fprintf('\n=======================================================================\n');
fprintf('Summary Statistics\n');
fprintf('=======================================================================\n');
fprintf('%4s %12s %12s %12s %12s %12s\n', 'N', 'Mean', 'Std', 'Min', 'Max', 'Equispaced');
fprintf('-----------------------------------------------------------------------\n');
for i = 1:3:n_N
    N = N_values(i);
    fprintf('%4d %12.4f %12.4f %12.4f %12.4f %12.4f\n', ...
            N, Lambda_mean(i), Lambda_std(i), Lambda_min(i), Lambda_max(i), Lambda_equi(i));
end
fprintf('-----------------------------------------------------------------------\n');

fprintf('\n=======================================================================\n');
fprintf('EMPIRICAL ASYMPTOTIC FORMULA\n');
fprintf('=======================================================================\n');
fprintf('\nFor random nodes uniformly distributed on [-1, 1]:\n');
fprintf('\n    Lambda_N ~ %.4f * exp(%.4f * N)\n', a_fit, b_fit);
fprintf('\nor approximately:\n');
fprintf('\n    Lambda_N ~ %.4f * %.4f^N\n', a_fit, exp(b_fit));
fprintf('\nThis confirms exponential growth similar to equispaced nodes,\n');
fprintf('demonstrating that random node placement provides no advantage\n');
fprintf('over structured distributions without the clustering property\n');
fprintf('of Chebyshev points near the interval endpoints.\n');
fprintf('=======================================================================\n');

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
