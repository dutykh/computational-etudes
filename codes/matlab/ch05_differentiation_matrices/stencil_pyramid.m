%% stencil_pyramid.m
%
% Visualizes Fornberg's recursive algorithm for computing finite difference
% weights as a "stencil pyramid" showing how weights are computed progressively.
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
% Part of "Computational Études: A Spectral Approach"
% https://github.com/dutykh/computational-etudes

clear; close all; clc;

%% Configuration
n_levels = 5;  % Number of levels (m = 0, 1, 2, 3, 4)
box_width = 0.8;
box_height = 0.5;
h_spacing = 1.0;
v_spacing = 1.2;

% Colors
NAVY = [0.078 0.176 0.431];
SKY = [0.471 0.588 0.824];
CORAL = [0.906 0.298 0.235];
TEAL = [0.086 0.627 0.522];
LIGHT_GRAY = [0.91 0.91 0.91];

% Output path
script_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(script_dir, '..', '..', '..', 'textbook', 'figures', 'ch05', 'matlab');
output_file = fullfile(output_dir, 'stencil_pyramid.pdf');

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% Create figure
fig = figure('Units', 'inches', 'Position', [1, 1, 9, 6]);
ax = axes;
hold on;

boxes = struct();

for m = 0:n_levels-1
    n_boxes = m + 1;
    y = (n_levels - 1 - m) * v_spacing;

    % Center the row
    total_width = n_boxes * box_width + (n_boxes - 1) * (h_spacing - box_width);
    x_start = -total_width / 2 + box_width / 2;

    for k = 0:n_boxes-1
        x = x_start + k * h_spacing;

        % Choose color based on position
        if m == n_levels - 1
            color = NAVY;
            text_color = 'white';
        elseif m == n_levels - 2
            color = SKY;
            text_color = 'white';
        else
            color = LIGHT_GRAY;
            text_color = NAVY;
        end

        % Draw box
        rectangle('Position', [x - box_width/2, y - box_height/2, box_width, box_height], ...
                  'Curvature', 0.2, 'FaceColor', color, 'EdgeColor', NAVY, 'LineWidth', 1.5);

        % Add text label
        if m <= 2
            label = sprintf('\\delta_%d^{(%d,%d)}', m, k, m);
        else
            label = sprintf('\\delta^{(%d,%d)}', k, m);
        end
        text(x, y, label, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
             'FontSize', 9, 'Color', text_color, 'FontWeight', 'bold', 'Interpreter', 'tex');

        % Store position
        boxes(m+1, k+1).x = x;
        boxes(m+1, k+1).y = y;
    end
end

% Draw arrows
for m = 1:n_levels-1
    n_boxes = m + 1;
    for k = 0:n_boxes-1
        x_to = boxes(m+1, k+1).x;
        y_to = boxes(m+1, k+1).y;

        if k < m
            % Arrow from previous level (existing node)
            x_from = boxes(m, k+1).x;
            y_from = boxes(m, k+1).y;
            annotation('arrow', ...
                       'X', [0.5 + x_from/12, 0.5 + x_to/12], ...
                       'Y', [0.5 + (y_from - box_height/2 - 0.05)/8, 0.5 + (y_to + box_height/2 + 0.05)/8], ...
                       'Color', TEAL, 'LineWidth', 1.2);
        end

        if k == m && m > 0
            % Arrow for new node
            x_from = boxes(m, m).x;
            y_from = boxes(m, m).y;
            annotation('arrow', ...
                       'X', [0.5 + (x_from + box_width/4)/12, 0.5 + (x_to - box_width/4)/12], ...
                       'Y', [0.5 + (y_from - box_height/2 - 0.05)/8, 0.5 + (y_to + box_height/2 + 0.05)/8], ...
                       'Color', CORAL, 'LineWidth', 1.5);
        end
    end
end

% Add level labels
for m = 0:n_levels-1
    y = (n_levels - 1 - m) * v_spacing;
    text(-4.5, y, sprintf('m = %d', m), 'FontSize', 11, 'Color', NAVY);
    text(-3.2, y, sprintf('(%d node%s)', m+1, char('s' * (m > 0))), ...
         'FontSize', 9, 'Color', [0.5 0.5 0.5]);
end

% Title
title("Fornberg's Recursive Algorithm: The Stencil Pyramid", ...
      'FontSize', 13, 'FontWeight', 'bold');

% Explanation
text(0, -0.8, {'Each box represents a weight \delta_j^{(k,m)} for node k in a stencil of m+1 nodes.', ...
               'Arrows show dependencies: weights at level m depend on level m-1.'}, ...
     'HorizontalAlignment', 'center', 'FontSize', 10, 'Color', [0.5 0.5 0.5], ...
     'FontStyle', 'italic');

% Axis settings
xlim([-5.5, 5.5]);
ylim([-1.5, (n_levels - 1) * v_spacing + 1]);
axis equal off;

%% Save figure
exportgraphics(fig, output_file, 'ContentType', 'vector');
exportgraphics(fig, strrep(output_file, '.pdf', '.png'), 'Resolution', 300);

fprintf('Figure saved to: %s\n', output_file);

close(fig);
