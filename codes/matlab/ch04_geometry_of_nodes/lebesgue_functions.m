%% lebesgue_functions.m
%
% Visualizes and compares Lebesgue functions for different node distributions.
%
% The Lebesgue function is Lambda_N(x) = sum |L_k(x)|, and its maximum
% is the Lebesgue constant Lambda_N.
%
% Asymptotic growth rates:
%   - Chebyshev: Lambda_N ~ (2/pi) ln(N)
%   - Legendre:  Lambda_N ~ O(sqrt(N))
%   - Equispaced: Lambda_N ~ 2^N / (N ln N)
%
% Author: Dr. Denys Dutykh
%         Mathematics Department
%         Khalifa University of Science and Technology
%         Abu Dhabi, UAE
%
% Part of the book "Computational Études: A Spectral Approach"

clear; close all; clc;

%% Configuration
N_FINE = 1000;

% Book color scheme
NAVY = [20, 45, 110] / 255;
SKY = [120, 150, 210] / 255;
CORAL = [231, 76, 60] / 255;
TEAL = [22, 160, 133] / 255;

% Output path
script_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch04', 'matlab');
output_file = fullfile(output_dir, 'lebesgue_functions.pdf');

%% Fine grid
x_fine = linspace(-1, 1, N_FINE)';

%% Create figure
fig = figure('Units', 'inches', 'Position', [1, 1, 11, 4.5]);

% Left panel: Lebesgue functions for N = 10
subplot(1, 2, 1);
hold on; box on;

N = 10;

% Equispaced
x_equi = linspace(-1, 1, N+1)';
Lambda_equi = lebesgue_function(x_equi, x_fine);

% Chebyshev
j = 0:N;
x_cheb = cos(j * pi / N)';
Lambda_cheb = lebesgue_function(x_cheb, x_fine);

% Legendre-Gauss-Lobatto (approximation using Legendre polynomial roots)
x_leg = legendre_nodes(N);
Lambda_leg = lebesgue_function(x_leg, x_fine);

plot(x_fine, Lambda_equi, 'Color', CORAL, 'LineWidth', 1.5, ...
     'DisplayName', sprintf('Equispaced ($\\Lambda_{%d}$ = %.1f)', N, max(Lambda_equi)));
plot(x_fine, Lambda_leg, 'Color', TEAL, 'LineWidth', 1.5, ...
     'DisplayName', sprintf('Legendre ($\\Lambda_{%d}$ = %.2f)', N, max(Lambda_leg)));
plot(x_fine, Lambda_cheb, 'Color', SKY, 'LineWidth', 1.5, ...
     'DisplayName', sprintf('Chebyshev ($\\Lambda_{%d}$ = %.2f)', N, max(Lambda_cheb)));

yline(1, '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5, 'HandleVisibility', 'off');

xlabel('$x$', 'Interpreter', 'latex');
ylabel('$\Lambda_N(x)$', 'Interpreter', 'latex');
title(sprintf('Lebesgue Functions ($N = %d$)', N), 'FontSize', 11);
xlim([-1, 1]);
ylim([0, 35]);
legend('Location', 'north', 'Interpreter', 'latex');

text(-0.6, 25, 'Peaks at boundaries', 'FontSize', 9, 'Color', [0.5, 0.5, 0.5]);

set(gca, 'FontSize', 10, 'TickLabelInterpreter', 'latex');
grid on;
set(gca, 'GridAlpha', 0.2);

% Right panel: Growth of Lebesgue constants
subplot(1, 2, 2);
hold on; box on;

N_values = 2:30;

Lambda_equi_vals = zeros(size(N_values));
Lambda_cheb_vals = zeros(size(N_values));
Lambda_leg_vals = zeros(size(N_values));

for i = 1:length(N_values)
    N = N_values(i);

    % Equispaced
    x_equi = linspace(-1, 1, N+1)';
    Lambda_equi_vals(i) = max(lebesgue_function(x_equi, x_fine));

    % Chebyshev
    j = 0:N;
    x_cheb = cos(j * pi / N)';
    Lambda_cheb_vals(i) = max(lebesgue_function(x_cheb, x_fine));

    % Legendre
    x_leg = legendre_nodes(N);
    Lambda_leg_vals(i) = max(lebesgue_function(x_leg, x_fine));
end

semilogy(N_values, Lambda_equi_vals, 'o-', 'Color', CORAL, 'LineWidth', 1.5, ...
         'MarkerSize', 4, 'DisplayName', 'Equispaced');
semilogy(N_values, Lambda_leg_vals, 's-', 'Color', TEAL, 'LineWidth', 1.5, ...
         'MarkerSize', 4, 'DisplayName', 'Legendre');
semilogy(N_values, Lambda_cheb_vals, '^-', 'Color', SKY, 'LineWidth', 1.5, ...
         'MarkerSize', 4, 'DisplayName', 'Chebyshev');

% Asymptotic curve for Chebyshev
N_asy = linspace(5, 30, 100);
Lambda_cheb_asy = (2/pi) * log(N_asy) + 0.6;
semilogy(N_asy, Lambda_cheb_asy, '--', 'Color', SKY, 'LineWidth', 1, 'Alpha', 0.5, ...
         'DisplayName', '$(2/\pi)\ln N$');

xlabel('$N$ (polynomial degree)', 'Interpreter', 'latex');
ylabel('Lebesgue constant $\Lambda_N$', 'Interpreter', 'latex');
title('Growth of Lebesgue Constants', 'FontSize', 11);
xlim([0, 32]);
legend('Location', 'northwest', 'Interpreter', 'latex');

% Add growth rate labels
text(25, 1e6, '$O(2^N)$', 'Interpreter', 'latex', 'FontSize', 10, 'Color', CORAL);
text(25, 4, '$O(\sqrt{N})$', 'Interpreter', 'latex', 'FontSize', 10, 'Color', TEAL);
text(25, 2.5, '$O(\ln N)$', 'Interpreter', 'latex', 'FontSize', 10, 'Color', SKY);

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
fprintf('\nLebesgue Constants Comparison:\n');
fprintf('%s\n', repmat('-', 1, 65));
fprintf('%4s %15s %15s %15s\n', 'N', 'Equispaced', 'Legendre', 'Chebyshev');
fprintf('%s\n', repmat('-', 1, 65));
for i = 1:5:length(N_values)
    N = N_values(i);
    fprintf('%4d %15.2f %15.4f %15.4f\n', N, ...
            Lambda_equi_vals(i), Lambda_leg_vals(i), Lambda_cheb_vals(i));
end
fprintf('%s\n', repmat('-', 1, 65));

%% Helper functions
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

function x = legendre_nodes(N)
    % Legendre-Gauss-Lobatto nodes
    if N == 0
        x = 0;
        return;
    end
    if N == 1
        x = [-1; 1];
        return;
    end

    % Interior nodes are zeros of P'_N
    % Use Chebyshev nodes as initial guess and Newton iteration
    j = 0:N;
    x = cos(j * pi / N)';

    % Newton iteration for LGL nodes
    for iter = 1:10
        [P, dP] = legendre_poly(N, x);
        % LGL nodes satisfy (1-x^2)P'_N(x) = 0, i.e., x = ±1 or P'_N(x) = 0
        % For interior nodes, we use P'_N(x) = 0
        % Update formula based on P_N and P'_N
        L = P ./ dP;
        x(2:N) = x(2:N) - L(2:N);
    end
    x = sort(x);
end

function [P, dP] = legendre_poly(N, x)
    % Evaluate Legendre polynomial P_N and its derivative at x
    P_prev = ones(size(x));
    P_curr = x;

    if N == 0
        P = P_prev;
        dP = zeros(size(x));
        return;
    end
    if N == 1
        P = P_curr;
        dP = ones(size(x));
        return;
    end

    for n = 1:N-1
        P_next = ((2*n+1)*x.*P_curr - n*P_prev) / (n+1);
        P_prev = P_curr;
        P_curr = P_next;
    end
    P = P_curr;

    % Derivative: P'_N(x) = N*(x*P_N - P_{N-1}) / (x^2 - 1)
    dP = N * (x.*P - P_prev) ./ (x.^2 - 1 + eps);
    % Handle endpoints
    dP(abs(x - 1) < 1e-10) = N*(N+1)/2;
    dP(abs(x + 1) < 1e-10) = (-1)^(N+1) * N*(N+1)/2;
end
