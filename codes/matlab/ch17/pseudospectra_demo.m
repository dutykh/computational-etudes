% pseudospectra_demo.m
% Chapter 17: Stability of Spectral Methods
%
% Computational Etude 17.5: Pseudospectra of Normal vs Nonnormal Matrices.
%
% Two-panel figure comparing pseudospectra:
%   (a) Normal diagonal matrix with eigenvalues -1, -2, ..., -16.
%   (b) Nonnormal upper triangular matrix with same diagonal but A(i,i+1) = 10.
%
% Pseudospectra computed via resolvent norm:
%   ||(zI - A)^{-1}|| = 1 / sigma_min(zI - A)
%
% Author: Dr. Denys Dutykh
% Part of "Computational Etudes: A Spectral Approach"

clear; clc; close all;

NAVY  = [20, 45, 110]/255;
SKY   = [120, 150, 210]/255;
CORAL = [231, 76, 60]/255;
DARK  = [10, 26, 64]/255;

N = 16;

% Normal matrix: diagonal with eigenvalues -1, -2, ..., -N
A_normal = diag(-(1:N));

% Nonnormal matrix: same diagonal, upper triangular with A(i,i+1) = 10
A_nonnormal = double(A_normal);
for i = 1:N-1
    A_nonnormal(i, i+1) = 10.0;
end

eigenvalues = -(1:N);

%% Compute pseudospectra of normal matrix
fprintf('Computing pseudospectra of normal matrix...\n');
x_range_n = [-18, 4]; y_range_n = [-8, 8];
n_grid = 250;
[X_n, Y_n, log_res_n] = compute_pseudospectra(A_normal, ...
    x_range_n, y_range_n, n_grid);

%% Compute pseudospectra of nonnormal matrix
fprintf('Computing pseudospectra of nonnormal matrix...\n');
x_range_nn = [-18, 10]; y_range_nn = [-12, 12];
[X_nn, Y_nn, log_res_nn] = compute_pseudospectra(A_nonnormal, ...
    x_range_nn, y_range_nn, n_grid);

%% Plot
levels = [1, 2, 3];
colors = {SKY, NAVY, DARK};

figure('Position', [50, 50, 1300, 500]);

% Panel (a): Normal matrix
subplot(1, 2, 1);
hold on;
for li = 1:length(levels)
    contour(X_n, Y_n, log_res_n, [levels(li), levels(li)], ...
            'LineColor', colors{li}, 'LineWidth', 0.8 + 0.2*li);
end
plot(eigenvalues, zeros(1, N), 'o', 'Color', CORAL, 'MarkerSize', 5, ...
     'MarkerFaceColor', CORAL, 'MarkerEdgeColor', 'w', ...
     'DisplayName', 'eigenvalues');
plot([x_range_n(1), x_range_n(2)], [0, 0], 'Color', [0.5 0.5 0.5], ...
     'LineWidth', 0.5, 'HandleVisibility', 'off');
plot([0, 0], [y_range_n(1), y_range_n(2)], 'Color', [0.5 0.5 0.5], ...
     'LineWidth', 0.5, 'HandleVisibility', 'off');
hold off;
xlabel('Re(z)'); ylabel('Im(z)');
title('(a) Normal matrix (A = diag(-1, ..., -16))');
legend('Location', 'northwest', 'FontSize', 8);
axis equal; grid on; set(gca, 'GridAlpha', 0.15);

% Panel (b): Nonnormal matrix
subplot(1, 2, 2);
hold on;
for li = 1:length(levels)
    contour(X_nn, Y_nn, log_res_nn, [levels(li), levels(li)], ...
            'LineColor', colors{li}, 'LineWidth', 0.8 + 0.2*li);
end
plot(eigenvalues, zeros(1, N), 'o', 'Color', CORAL, 'MarkerSize', 5, ...
     'MarkerFaceColor', CORAL, 'MarkerEdgeColor', 'w', ...
     'DisplayName', 'eigenvalues');
plot([x_range_nn(1), x_range_nn(2)], [0, 0], 'Color', [0.5 0.5 0.5], ...
     'LineWidth', 0.5, 'HandleVisibility', 'off');
plot([0, 0], [y_range_nn(1), y_range_nn(2)], 'Color', [0.5 0.5 0.5], ...
     'LineWidth', 0.5, 'HandleVisibility', 'off');

% Annotate the bulge
text(2, 3, {'pseudospectral', 'bulge into RHP'}, ...
     'FontSize', 9, 'Color', CORAL, 'FontAngle', 'italic', ...
     'HorizontalAlignment', 'center');
hold off;

xlabel('Re(z)'); ylabel('Im(z)');
title('(b) Nonnormal (A_{i,i+1} = 10, same diagonal)');
legend('Location', 'northwest', 'FontSize', 8);
axis equal; grid on; set(gca, 'GridAlpha', 0.15);

set(gcf, 'PaperPositionMode', 'auto');
print('-dpdf', 'pseudospectra_comparison.pdf');
fprintf('Figure saved to pseudospectra_comparison.pdf\n');


%% ---- Helper function: compute pseudospectra ----
function [X, Y, log_resolvent] = compute_pseudospectra(A, x_range, y_range, n_grid)
    xv = linspace(x_range(1), x_range(2), n_grid);
    yv = linspace(y_range(1), y_range(2), n_grid);
    [X, Y] = meshgrid(xv, yv);

    n = size(A, 1);
    log_resolvent = zeros(size(X));

    for i = 1:n_grid
        for j = 1:n_grid
            z = X(i,j) + 1i*Y(i,j);
            M = z * eye(n) - A;
            s = svd(M);
            smin = s(end);
            if smin > 0
                log_resolvent(i,j) = log10(1.0 / smin);
            else
                log_resolvent(i,j) = 15.0;  % effectively infinity
            end
        end
    end
end
