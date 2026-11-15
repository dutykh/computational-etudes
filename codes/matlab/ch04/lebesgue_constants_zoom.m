%% lebesgue_constants_zoom.m
%
% Zoomed comparison of Lebesgue constant growth for Legendre-Gauss-Lobatto and
% Chebyshev-Gauss-Lobatto nodes.
%
% This figure omits the exponentially growing equispaced curve, so that the two
% logarithmic growths can be told apart:
%     Chebyshev-Gauss-Lobatto: Lambda_N = (2/pi)(ln N + gamma + ln(8/pi)) + o(1)
%     Legendre-Gauss-Lobatto:  Lambda_N = O(ln N), close to (2/pi) ln(N+1) + 0.685,
%                              the estimate of Hesthaven drawn from numerical experiments
% The Legendre-Gauss-Lobatto nodes are the Fekete points of [-1, 1] (Fejer),
% and their constant stays about 0.28 below the Chebyshev one. Neither family
% is optimal: every node set has Lambda_N > (2/pi) ln(N+1) + 0.5212 (Brutman,
% 1978), and both curves stay a little above this bound.
%
% Author: Dr. Denys Dutykh
%         Mathematics Department
%         Khalifa University of Science and Technology
%         Abu Dhabi, UAE
%
% Part of the book "Computational Etudes: A Spectral Approach"

clear; close all; clc;

%% Configuration
NAVY = [20, 45, 110] / 255;
SKY  = [120, 150, 210] / 255;
TEAL = [22, 160, 133] / 255;

% Output path
script_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch04', 'matlab');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

%% Compute Lebesgue constants
N_values = 2:50;
Lambda_cheb = zeros(size(N_values));
Lambda_leg  = zeros(size(N_values));

for idx = 1:length(N_values)
    N = N_values(idx);

    % Chebyshev-Gauss-Lobatto nodes
    j = (0:N)';
    x_cheb = cos(j * pi / N);
    Lambda_cheb(idx) = lebesgue_constant(x_cheb);

    % Legendre-Gauss-Lobatto nodes
    x_leg = legendre_gauss_lobatto(N);
    Lambda_leg(idx) = lebesgue_constant(x_leg);
end

%% Create figure
fig = figure('Position', [100, 100, 700, 500]);

plot(N_values, Lambda_leg, 's-', 'Color', TEAL, 'LineWidth', 1.5, ...
     'MarkerSize', 5, 'MarkerFaceColor', TEAL, 'DisplayName', 'Legendre--Gauss--Lobatto');
hold on;
plot(N_values, Lambda_cheb, '^-', 'Color', SKY, 'LineWidth', 1.5, ...
     'MarkerSize', 5, 'MarkerFaceColor', SKY, 'DisplayName', 'Chebyshev--Gauss--Lobatto');

% Logarithmic asymptotics of the two families
N_asy = linspace(2, 50, 200);
euler_gamma = 0.57721566490153286;

% Chebyshev-Gauss-Lobatto: (2/pi)(ln N + gamma + ln(8/pi)) + o(1)
Lambda_cheb_asy = (2/pi) * (log(N_asy) + euler_gamma + log(8/pi));
plot(N_asy, Lambda_cheb_asy, '--', 'Color', SKY, 'LineWidth', 1.5, ...
     'DisplayName', '$(2/\pi)(\ln N + \gamma + \ln(8/\pi))$');

% Legendre-Gauss-Lobatto: the estimate (2/pi) ln(N+1) + 0.685 of Hesthaven
Lambda_leg_asy = (2/pi) * log(N_asy + 1) + 0.685;
plot(N_asy, Lambda_leg_asy, '--', 'Color', TEAL, 'LineWidth', 1.5, ...
     'DisplayName', '$(2/\pi)\ln(N+1) + 0.685$');

% Lower bound for every node set, hence for the optimal nodes (Brutman)
Lambda_low = (2/pi) * log(N_asy + 1) + 0.5212;
plot(N_asy, Lambda_low, ':', 'Color', [0.33, 0.33, 0.33], 'LineWidth', 1.5, ...
     'DisplayName', '$(2/\pi)\ln(N+1) + 0.5212$, lower bound for all nodes');
hold off;

xlabel('$N$ (polynomial degree)', 'Interpreter', 'latex');
ylabel('Lebesgue constant $\Lambda_N$', 'Interpreter', 'latex');
title('Lebesgue Constants: Legendre vs Chebyshev', 'FontSize', 11);
xlim([0, 52]);
ylim([1.0, 4.4]);
legend('Interpreter', 'latex', 'Location', 'northwest', 'FontSize', 10);
grid on; set(gca, 'GridAlpha', 0.2);

% Save figure
exportgraphics(fig, fullfile(output_dir, 'lebesgue_constants_zoom.pdf'), 'ContentType', 'vector');
exportgraphics(fig, fullfile(output_dir, 'lebesgue_constants_zoom.png'), 'Resolution', 300);
fprintf('Figure saved to %s\n', fullfile(output_dir, 'lebesgue_constants_zoom.pdf'));

%% Print comparison table
fprintf('\nLebesgue Constants Comparison:\n');
fprintf('%4s %12s %12s %12s %13s\n', 'N', 'Legendre', 'Chebyshev', 'Difference', 'Lower bound');
for idx = 1:8:length(N_values)
    N = N_values(idx);
    low = (2/pi) * log(N + 1) + 0.5212;
    fprintf('%4d %12.4f %12.4f %12.4f %13.4f\n', N, Lambda_leg(idx), Lambda_cheb(idx), ...
            Lambda_cheb(idx) - Lambda_leg(idx), low);
end


%% ---- Helper functions ----

function Lambda = lebesgue_constant(x_nodes)
    % Compute the Lebesgue constant for given nodes
    x_eval = linspace(-1, 1, 2000)';
    N = length(x_nodes) - 1;
    Lambda_func = zeros(size(x_eval));

    for k = 1:N+1
        L_k = ones(size(x_eval));
        for j = 1:N+1
            if j ~= k
                L_k = L_k .* (x_eval - x_nodes(j)) / (x_nodes(k) - x_nodes(j));
            end
        end
        Lambda_func = Lambda_func + abs(L_k);
    end

    Lambda = max(Lambda_func);
end

function x = legendre_gauss_lobatto(N)
    % Legendre-Gauss-Lobatto nodes on [-1, 1]
    % LGL nodes are: x_0 = -1, x_N = 1, interior = zeros of P'_N(x)
    if N == 0
        x = 0;
        return;
    end
    if N == 1
        x = [-1; 1];
        return;
    end

    % Start with Chebyshev-Gauss-Lobatto as initial guess
    j = 0:N;
    x = -cos(j * pi / N)';
    x(1) = -1;
    x(N+1) = 1;

    % Newton iteration for interior nodes (zeros of P'_N)
    for iter = 1:20
        [~, dP, d2P] = legendre_poly_d2(N, x(2:N));
        dx = dP ./ d2P;
        x(2:N) = x(2:N) - dx;
        if max(abs(dx)) < 1e-15, break; end
    end
    x = sort(x);
end

function [P, dP, d2P] = legendre_poly_d2(N, x)
    % Evaluate P_N, P'_N, P''_N via three-term recurrence
    P_prev = ones(size(x));   dP_prev = zeros(size(x));  d2P_prev = zeros(size(x));
    P_curr = x;               dP_curr = ones(size(x));   d2P_curr = zeros(size(x));

    if N == 0, P = P_prev; dP = dP_prev; d2P = d2P_prev; return; end
    if N == 1, P = P_curr; dP = dP_curr; d2P = d2P_curr; return; end

    for n = 1:N-1
        P_next  = ((2*n+1)*x.*P_curr  - n*P_prev)  / (n+1);
        dP_next = ((2*n+1)*(P_curr + x.*dP_curr) - n*dP_prev) / (n+1);
        d2P_next = ((2*n+1)*(2*dP_curr + x.*d2P_curr) - n*d2P_prev) / (n+1);

        P_prev = P_curr;    P_curr = P_next;
        dP_prev = dP_curr;  dP_curr = dP_next;
        d2P_prev = d2P_curr; d2P_curr = d2P_next;
    end
    P = P_curr; dP = dP_curr; d2P = d2P_curr;
end
