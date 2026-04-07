% trap_algebraic_decay.m
% Chapter 16: Integration of Periodic Functions
%
% Computational Etude 16.3: Algebraic Convergence on |sin(x/2)|^k.
%
% Two periodic functions with different regularity at the endpoints
% of [0, 2*pi]:
%   f_2(x) = |sin(x/2)|       --> O(N^{-2}) convergence
%   f_3(x) = |sin(x/2)|^3     --> O(N^{-4}) convergence
%
% Generates Figure 16.3: algebraic_decay.pdf
%
% Author: Dr. Denys Dutykh
%         Mathematics Department
%         Khalifa University of Science and Technology
%         Abu Dhabi, UAE
%
% Part of the book "Computational Etudes: A Spectral Approach"
% https://github.com/dutykh/computational-etudes

clear; close all; clc;

set(0, 'DefaultAxesFontSize', 11);
set(0, 'DefaultTextFontSize', 11);
set(0, 'DefaultLineLineWidth', 1.5);
set(0, 'DefaultAxesLineWidth', 0.8);

NAVY  = [20, 45, 110] / 255;
CORAL = [231, 76, 60] / 255;

output_dir = fullfile(fileparts(mfilename('fullpath')), ...
    '..', '..', '..', 'textbook', 'figures', 'ch16', 'matlab');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

% -------------------------------------------------------------------------
% The two test functions
% -------------------------------------------------------------------------
f2 = @(x) abs(sin(x / 2));
f3 = @(x) abs(sin(x / 2)).^3;
I2_exact = 4;
I3_exact = 8 / 3;

N_values = 2.^(3:11);
err2 = zeros(size(N_values));
err3 = zeros(size(N_values));
for k = 1:length(N_values)
    N = N_values(k);
    theta = 2*pi*(0:N-1)/N;
    err2(k) = abs((2*pi/N) * sum(f2(theta)) - I2_exact);
    err3(k) = abs((2*pi/N) * sum(f3(theta)) - I3_exact);
end

% Theoretical asymptotic constants from Weideman 2002
theory_f2 = (pi^2 / 3) * double(N_values).^(-2);
theory_f3 = (pi^4 / 30) * double(N_values).^(-4);

% Linear-fit slopes (skip first two points)
p2 = polyfit(log10(double(N_values(3:end))), log10(err2(3:end)), 1);
p3 = polyfit(log10(double(N_values(3:end))), log10(err3(3:end)), 1);
fprintf('Linear-fit slope for f_2: %.3f (theory -2)\n', p2(1));
fprintf('Linear-fit slope for f_3: %.3f (theory -4)\n', p3(1));

% -------------------------------------------------------------------------
% Plot
% -------------------------------------------------------------------------
fig = figure('Position', [100, 100, 700, 500]);
loglog(N_values, err2, 'o-', 'Color', CORAL, 'MarkerSize', 5, ...
    'LineWidth', 1.2, 'MarkerFaceColor', CORAL, ...
    'DisplayName', '$f_2 = |\sin(x/2)|$');
hold on;
loglog(N_values, theory_f2, '--', 'Color', CORAL, 'LineWidth', 0.8, ...
    'DisplayName', '$\pi^2/(3 N^2)$');
loglog(N_values, err3, 's-', 'Color', NAVY, 'MarkerSize', 5, ...
    'LineWidth', 1.2, 'MarkerFaceColor', NAVY, ...
    'DisplayName', '$f_3 = |\sin(x/2)|^3$');
loglog(N_values, theory_f3, '--', 'Color', NAVY, 'LineWidth', 0.8, ...
    'DisplayName', '$\pi^4/(30 N^4)$');
hold off;
xlabel('Number of nodes $N$', 'Interpreter', 'latex');
ylabel('Absolute error $|I_N - I|$', 'Interpreter', 'latex');
title('Algebraic convergence: regularity of the periodic extension', ...
    'Interpreter', 'latex');
legend('Location', 'southwest', 'Interpreter', 'latex', ...
    'EdgeColor', [0.5 0.5 0.5]);
grid on; set(gca, 'GridAlpha', 0.3);
xlim([6, 3000]);

exportgraphics(gcf, fullfile(output_dir, 'algebraic_decay.pdf'), 'ContentType', 'vector');
exportgraphics(gcf, fullfile(output_dir, 'algebraic_decay.png'), 'Resolution', 300);
fprintf('\nFigure saved to %s\n', fullfile(output_dir, 'algebraic_decay.pdf'));
