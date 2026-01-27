%% cheb_convergence.m - Spectral convergence rates demonstration
%
% Demonstrates spectral convergence rates for four functions of increasing
% smoothness. Shows the fundamental relationship between function regularity
% and convergence rate of spectral differentiation.
%
% Test functions:
% 1. |x|^(5/2) - Third derivative of bounded variation (algebraic: ~N^{-2.5})
% 2. exp(-1/(1-x^2)) - Smooth bump function, C^inf but not analytic (superalgebraic)
% 3. tanh(5x) - Analytic in [-1,1], poles at +/-i*pi/10 (exponential)
% 4. x^8 - Polynomial of degree 8 (exact for N >= 8)
%
% This script generates Figure 6.5 for Chapter 6.
%
% Author: Dr. Denys Dutykh
%         Mathematics Department
%         Khalifa University of Science and Technology
%         Abu Dhabi, UAE
%
% Part of the book "Computational Etudes: A Spectral Approach"
% https://github.com/dutykh/computational-etudes
%
% Last modified: January 2026

clear;
close all;
clc;

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
NAVY = [20, 45, 110] / 255;
SKY = [120, 150, 210] / 255;
CORAL = [231, 76, 60] / 255;
TEAL = [22, 160, 133] / 255;
PURPLE = [142, 68, 173] / 255;
ORANGE = [230, 126, 34] / 255;

%% Output directory
script_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch06', 'matlab');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% Test functions and their derivatives
% Function 1: |x|^(5/2)
f1 = @(x) abs(x).^2.5;
f1_prime = @(x) 2.5 * abs(x).^1.5 .* sign(x);

% Function 2: exp(-1/(1-x^2)) - bump function
f2 = @(x) bump_func(x);
f2_prime = @(x) bump_deriv(x);

% Function 3: tanh(5x)
f3 = @(x) tanh(5 * x);
f3_prime = @(x) 5 ./ cosh(5 * x).^2;

% Function 4: x^8
f4 = @(x) x.^8;
f4_prime = @(x) 8 * x.^7;

%% Store test functions
test_funcs = {f1, f2, f3, f4};
test_derivs = {f1_prime, f2_prime, f3_prime, f4_prime};
names = {'$|x|^{5/2}$', '$e^{-1/(1-x^2)}$', '$\tanh(5x)$', '$x^8$'};
smoothness = {'$C^2$, 3rd deriv in BV', '$C^\infty$, not analytic', 'Analytic', 'Polynomial'};
rates = {'$O(N^{-2.5})$', 'faster than any $N^{-k}$', '$O(\rho^{-N})$', 'Exact for $N \geq 8$'};
colors = {CORAL, TEAL, PURPLE, ORANGE};

%% Range of N values
N_values = [4, 6, 8, 10, 12, 16, 20, 24, 32, 40, 48, 64, 80, 96, 128];

%% Create 2x2 figure
fig = figure('Position', [100, 100, 1000, 800]);

error_data = cell(4, 1);

for idx = 1:4
    subplot(2, 2, idx);

    errors = zeros(size(N_values));

    for k = 1:length(N_values)
        N = N_values(k);
        [D, x] = cheb_matrix(N);
        v = test_funcs{idx}(x);
        w = D * v;
        w_exact = test_derivs{idx}(x);

        % Compute max error (avoiding boundary for bump function)
        if idx == 2  % Bump function
            interior = (x > -0.99) & (x < 0.99);
            if any(interior)
                error_val = max(abs(w(interior) - w_exact(interior)));
            else
                error_val = max(abs(w - w_exact));
            end
        else
            error_val = max(abs(w - w_exact));
        end

        errors(k) = max(error_val, 1e-16);
        error_data{idx}(k, :) = [N, errors(k)];
    end

    % Plot convergence
    semilogy(N_values, errors, 'o-', 'Color', colors{idx}, 'LineWidth', 1.5, ...
        'MarkerSize', 5, 'MarkerFaceColor', colors{idx}, 'MarkerEdgeColor', 'white');

    % Add reference lines
    if idx == 1  % |x|^(5/2) - algebraic O(N^-2.5)
        N_ref = N_values(N_values >= 16);
        err_ref = errors(N_values >= 16);
        if ~isempty(err_ref)
            C = err_ref(1) * N_ref(1)^2.5;
            hold on;
            semilogy(N_ref, C * N_ref.^(-2.5), '--', 'Color', [0.5, 0.5, 0.5], ...
                'LineWidth', 1);
            hold off;
        end
    elseif idx == 4  % x^8 - exact for N >= 8
        hold on;
        yline(1e-14, '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1);
        hold off;
    end

    xlabel('$N$');
    ylabel('Max Error');
    title(sprintf('%s: %s', names{idx}, smoothness{idx}));
    xlim([0, 135]);
    ylim([1e-16, 1e2]);
    grid on;
    box off;

    % Add convergence rate annotation
    text(0.95, 0.95, rates{idx}, 'Units', 'normalized', ...
        'FontSize', 9, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', ...
        'BackgroundColor', 'white', 'EdgeColor', [0.8, 0.8, 0.8]);
end

% Main title
sgtitle('Spectral Convergence: Four Functions of Increasing Smoothness', ...
    'FontSize', 13);

%% Save figure
output_file = fullfile(output_dir, 'convergence_waterfall');
exportgraphics(fig, [output_file, '.pdf'], 'ContentType', 'vector');
exportgraphics(fig, [output_file, '.png'], 'Resolution', 300);
fprintf('Figure saved to: %s.pdf\n', output_file);

%% Print convergence table
fprintf('\n%s\n', repmat('=', 1, 80));
fprintf('Convergence Table: Chebyshev Spectral Differentiation\n');
fprintf('%s\n', repmat('=', 1, 80));

fprintf('%6s', 'N');
for idx = 1:4
    fprintf('  %14s', strrep(strrep(names{idx}, '$', ''), '\', ''));
end
fprintf('\n');
fprintf('%s\n', repmat('-', 1, 80));

for k = 1:length(N_values)
    fprintf('%6d', N_values(k));
    for idx = 1:4
        err = error_data{idx}(k, 2);
        if err < 1e-14
            fprintf('  %14s', '< 1e-14');
        else
            fprintf('  %14.2e', err);
        end
    end
    fprintf('\n');
end
fprintf('%s\n', repmat('=', 1, 80));

%% Local functions
function y = bump_func(x)
    y = zeros(size(x));
    mask = abs(x) < 1;
    y(mask) = exp(-1 ./ (1 - x(mask).^2));
end

function y = bump_deriv(x)
    y = zeros(size(x));
    mask = abs(x) < 1;
    inner = 1 - x(mask).^2;
    y(mask) = exp(-1 ./ inner) .* (-2 * x(mask)) ./ inner.^2;
end
