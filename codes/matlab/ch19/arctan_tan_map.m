function arctan_tan_map()
%% arctan_tan_map - Chapter 19, illustration of the arctan/tan map (eq. 19.2).
%
% Two-panel figure:
%  (a) the map y = 2*atan( ell * tan(x/2) ) for several values of ell,
%      with the identity (ell = 1) drawn dashed in grey;
%  (b) the images of a uniform N = 24 x-grid under the map, plotted as
%      tick marks on five stacked horizontal bands -- visualising that
%      ell < 1 clusters the grid near y = 0 and ell > 1 clusters it
%      near y = +/- pi.
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
% Part of "Computational Etudes: A Spectral Approach"

    script_dir = fileparts(mfilename('fullpath'));
    out_dir = fullfile(script_dir, '..', '..', '..', ...
                       'textbook', 'figures', 'ch19', 'matlab');
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end

    NAVY   = [ 20  45 110]/255;
    CORAL  = [231  76  60]/255;
    TEAL   = [ 22 160 133]/255;
    ORANGE = [230 126  34]/255;
    PURPLE = [142  68 173]/255;
    GREY   = [0.45 0.45 0.45];

    EPS_END = 1e-12;
    ells   = [0.1, 0.3, 1.0, 3.0, 10.0];
    cols   = {CORAL; TEAL; GREY; ORANGE; PURPLE};
    styles = {'-' , '-' , '--', '-' , '-' };

    fig = figure('Position', [100 100 1100 440], 'Color', 'w');

    %% --- Panel (a): the map curves y(x) ---------------------------------
    subplot(1,2,1); hold on;
    x_line = linspace(-pi + EPS_END, pi - EPS_END, 1001);
    legendH = gobjects(1, length(ells));
    for k = 1:length(ells)
        ell = ells(k);
        y_line = 2*atan(ell * tan(x_line/2));
        legendH(k) = plot(x_line, y_line, styles{k}, ...
            'Color', cols{k}, 'LineWidth', 1.6);
    end
    plot([-pi pi], [0 0], 'Color', [0.7 0.7 0.7], 'LineWidth', 0.4);
    plot([0 0], [-pi pi], 'Color', [0.7 0.7 0.7], 'LineWidth', 0.4);
    grid on; box on; axis equal;
    xlim([-pi pi]); ylim([-pi pi]);
    set(gca, 'XTick', [-pi -pi/2 0 pi/2 pi]);
    set(gca, 'XTickLabel', {'$-\pi$','$-\pi/2$','$0$','$\pi/2$','$\pi$'});
    set(gca, 'YTick', [-pi -pi/2 0 pi/2 pi]);
    set(gca, 'YTickLabel', {'$-\pi$','$-\pi/2$','$0$','$\pi/2$','$\pi$'});
    set(gca, 'TickLabelInterpreter', 'latex');
    xlabel('computational coordinate $x$', 'Interpreter', 'latex');
    ylabel('physical coordinate $y$', 'Interpreter', 'latex');
    title('(a) the map $y = 2\,\mathrm{atan}(\ell\,\tan(x/2))$', ...
          'Interpreter', 'latex');
    legend(legendH, ...
        arrayfun(@(e) sprintf('$\\ell = %g$', e), ells, 'UniformOutput', false), ...
        'Location', 'northwest', 'Box', 'off', 'Interpreter', 'latex');

    %% --- Panel (b): clustering of a uniform N = 24 x-grid ---------------
    subplot(1,2,2); hold on;
    N = 24;
    x_uniform = -pi + 2*pi*((0:N-1) + 0.5)/N;
    bands_y = length(ells):-1:1;     % stack top-down
    for k = 1:length(ells)
        ell = ells(k);
        band = bands_y(k);
        col = cols{k};
        y_mapped = 2*atan(ell * tan(x_uniform/2));
        plot([-pi pi], [band band], 'Color', [0.85 0.85 0.85], 'LineWidth', 0.6);
        for j = 1:N
            plot([y_mapped(j) y_mapped(j)], [band-0.18 band+0.18], ...
                'Color', col, 'LineWidth', 1.5);
        end
        text(-pi-0.18, band, sprintf('$\\ell = %g$', ell), ...
             'HorizontalAlignment','right', 'VerticalAlignment','middle', ...
             'Color', NAVY, 'FontSize', 10, 'Interpreter', 'latex');
    end
    xlim([-pi pi]); ylim([0.4 length(ells)+0.6]);
    set(gca, 'XTick', [-pi -pi/2 0 pi/2 pi]);
    set(gca, 'XTickLabel', {'$-\pi$','$-\pi/2$','$0$','$\pi/2$','$\pi$'});
    set(gca, 'TickLabelInterpreter', 'latex');
    set(gca, 'YTick', []);
    box on;
    xlabel('physical coordinate $y$', 'Interpreter', 'latex');
    title(sprintf('(b) image of a uniform $N = %d$ $x$-grid', N), ...
          'Interpreter', 'latex');

    set(fig, 'PaperPositionMode', 'auto');
    pos = get(fig, 'Position');
    set(fig, 'PaperUnits', 'points', ...
             'PaperSize', [pos(3) pos(4)], ...
             'PaperPosition', [0 0 pos(3) pos(4)]);
    print(fig, fullfile(out_dir, 'arctan_tan_map.pdf'), '-dpdf');
    print(fig, fullfile(out_dir, 'arctan_tan_map.png'), '-dpng');
    close(fig);

    fprintf('[19.0-matlab] saved figure to %s\n', ...
            fullfile(out_dir, 'arctan_tan_map.pdf'));
    for k = 1:length(ells)
        ell = ells(k);
        y_mapped = sort(2*atan(ell * tan(x_uniform/2)));
        gaps = diff(y_mapped);
        fprintf('  ell=%5.2f  min gap=%.4f  max gap=%.4f  ratio=%.2f\n', ...
                ell, min(gaps), max(gaps), max(gaps)/min(gaps));
    end
end
