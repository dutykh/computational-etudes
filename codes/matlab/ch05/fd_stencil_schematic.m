%% fd_stencil_schematic.m
%
% Creates a schematic illustration of finite difference stencils in 1D,
% showing how the stencil width increases with the order of accuracy,
% culminating in the spectral method that uses all nodes.
%
% Illustrates:
%     - 3-point stencil (2nd order): uses x_{i-1}, x_i, x_{i+1}
%     - 5-point stencil (4th order): uses x_{i-2}, ..., x_{i+2}
%     - 7-point stencil (6th order): uses x_{i-3}, ..., x_{i+3}
%     - Spectral stencil: uses ALL nodes
%
% Author: Dr. Denys Dutykh
%         Mathematics Department
%         Khalifa University of Science and Technology
%         Abu Dhabi, UAE
%
% Part of the book "Computational Études: A Spectral Approach"
% https://github.com/dutykh/computational-etudes

clear; close all; clc;

%% Configuration
% Book color scheme
NAVY = [20, 45, 110] / 255;
SKY = [120, 150, 210] / 255;
CORAL = [231, 76, 60] / 255;
TEAL = [22, 160, 133] / 255;
LIGHT_GRAY = [0.73, 0.73, 0.73];
VERY_LIGHT = [0.87, 0.87, 0.87];

% Output path
script_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch05', 'matlab');
output_file = fullfile(output_dir, 'fd_stencil_schematic.pdf');

%% Grid setup
n_points = 11;      % Total points (odd for symmetry)
center_idx = 6;     % Index of central point (1-indexed)

% Stencil configurations: {half_width, label, is_spectral}
stencils = {
    {1, '2nd order FD (3 points)', false};
    {2, '4th order FD (5 points)', false};
    {3, '6th order FD (7 points)', false};
    {-1, 'Spectral (all points)', true};
};

%% Create figure
fig = figure('Units', 'inches', 'Position', [1, 1, 11, 7.5]);

% Point positions
x_left = 0.08;
x_right = 0.72;
x_positions = linspace(x_left, x_right, n_points);
h = x_positions(2) - x_positions(1);

for row = 1:4
    half_width = stencils{row}{1};
    label = stencils{row}{2};
    is_spectral = stencils{row}{3};

    subplot(4, 1, row);
    hold on;

    y_base = 0.0;

    xlim([0, 1]);
    ylim([-1, 1]);
    axis off;

    % Draw horizontal grid line
    plot([x_positions(1) - 0.02, x_positions(end) + 0.02], [y_base, y_base], ...
         '--', 'Color', VERY_LIGHT, 'LineWidth', 1.2);

    % Determine active indices
    if is_spectral
        active_indices = 1:n_points;
    else
        active_indices = (center_idx - half_width):(center_idx + half_width);
    end

    % Draw shaded region
    if is_spectral
        rect_left = x_positions(1) - 0.012;
        rect_right = x_positions(end) + 0.012;
    else
        rect_left = x_positions(center_idx - half_width) - 0.012;
        rect_right = x_positions(center_idx + half_width) + 0.012;
    end

    rectangle('Position', [rect_left, y_base - 0.45, rect_right - rect_left, 0.9], ...
              'FaceColor', [SKY, 0.15], 'EdgeColor', SKY, 'LineWidth', 1.2);

    % Draw grid points
    for i = 1:n_points
        x = x_positions(i);

        if i == center_idx
            color = CORAL;
            sz = 180;
        elseif ismember(i, active_indices)
            color = NAVY;
            sz = 130;
        else
            color = LIGHT_GRAY;
            sz = 60;
        end

        scatter(x, y_base, sz, color, 'filled', 'MarkerEdgeColor', 'white', 'LineWidth', 2);
    end

    % Continuation dots
    text(x_positions(1) - 0.025, y_base, '$\cdots$', 'HorizontalAlignment', 'right', ...
         'VerticalAlignment', 'middle', 'FontSize', 16, 'Color', LIGHT_GRAY, ...
         'Interpreter', 'latex');
    text(x_positions(end) + 0.025, y_base, '$\cdots$', 'HorizontalAlignment', 'left', ...
         'VerticalAlignment', 'middle', 'FontSize', 16, 'Color', LIGHT_GRAY, ...
         'Interpreter', 'latex');

    % Labels below points
    label_y = y_base - 0.7;

    if is_spectral
        text(x_positions(1), label_y, '$x_0$', 'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'top', 'FontSize', 12, 'Color', NAVY, 'Interpreter', 'latex');
        text(x_positions(center_idx), label_y, '$x_i$', 'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'top', 'FontSize', 12, 'Color', CORAL, 'Interpreter', 'latex');
        text(x_positions(end), label_y, '$x_N$', 'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'top', 'FontSize', 12, 'Color', NAVY, 'Interpreter', 'latex');
    else
        left_idx = center_idx - half_width;
        right_idx = center_idx + half_width;

        if half_width == 1
            left_lbl = '$x_{i-1}$'; right_lbl = '$x_{i+1}$';
        elseif half_width == 2
            left_lbl = '$x_{i-2}$'; right_lbl = '$x_{i+2}$';
        else
            left_lbl = '$x_{i-3}$'; right_lbl = '$x_{i+3}$';
        end

        text(x_positions(left_idx), label_y, left_lbl, 'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'top', 'FontSize', 12, 'Color', NAVY, 'Interpreter', 'latex');
        text(x_positions(center_idx), label_y, '$x_i$', 'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'top', 'FontSize', 12, 'Color', CORAL, 'Interpreter', 'latex');
        text(x_positions(right_idx), label_y, right_lbl, 'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'top', 'FontSize', 12, 'Color', NAVY, 'Interpreter', 'latex');
    end

    % Method label on the right
    text(0.88, y_base, label, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
         'FontSize', 12, 'Color', NAVY, 'FontWeight', 'bold');

    % Bracket above stencil
    bracket_y = y_base + 0.55;
    tick_len = 0.12;
    plot([rect_left, rect_left], [bracket_y - tick_len, bracket_y], '-', 'Color', NAVY, 'LineWidth', 1.5);
    plot([rect_right, rect_right], [bracket_y - tick_len, bracket_y], '-', 'Color', NAVY, 'LineWidth', 1.5);
    plot([rect_left, rect_right], [bracket_y, bracket_y], '-', 'Color', NAVY, 'LineWidth', 1.5);

    % Grid spacing (first row only)
    if row == 1
        hx1 = x_positions(center_idx);
        hx2 = x_positions(center_idx + 1);
        hy = y_base + 0.85;

        annotation('doublearrow', ...
                   [0.08 + hx1 * 0.64, 0.08 + hx2 * 0.64], ...
                   [0.91, 0.91], ...
                   'Color', TEAL, 'LineWidth', 2, ...
                   'Head1Style', 'vback2', 'Head2Style', 'vback2', ...
                   'Head1Length', 8, 'Head2Length', 8);

        text((hx1 + hx2) / 2, hy + 0.15, '$h$', 'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'bottom', 'FontSize', 14, 'Color', TEAL, ...
             'FontWeight', 'bold', 'Interpreter', 'latex');
    end
end

% Main title
sgtitle('Finite Difference Stencils: From Local to Global', ...
        'FontSize', 15, 'FontWeight', 'bold', 'Color', NAVY);

% Bottom annotation
annotation('textbox', [0.1, 0.01, 0.8, 0.04], ...
           'String', 'As stencil width increases, accuracy improves. The spectral method is the limit: all nodes contribute to $u''(x_i)$.', ...
           'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
           'FontSize', 11, 'FontAngle', 'italic', 'Color', [0.5, 0.5, 0.5], ...
           'EdgeColor', 'none', 'Interpreter', 'latex');

%% Save figure
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

exportgraphics(fig, output_file, 'ContentType', 'vector');
fprintf('Figure saved to: %s\n', output_file);

png_file = strrep(output_file, '.pdf', '.png');
exportgraphics(fig, png_file, 'Resolution', 300);
fprintf('PNG saved to: %s\n', png_file);
