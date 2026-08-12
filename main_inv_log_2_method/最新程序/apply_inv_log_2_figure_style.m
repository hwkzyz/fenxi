function apply_inv_log_2_figure_style(ax, widthMm, heightMm)
%apply_inv_log_2_figure_style  Consistent paper-style figure formatting.

if nargin < 2 || isempty(widthMm)
    widthMm = 170;
end
if nargin < 3 || isempty(heightMm)
    heightMm = 85;
end

fig = ancestor(ax, 'figure');
set(fig, 'Color', 'w', 'Units', 'centimeters');
set(fig, 'Position', [2, 2, widthMm / 10, heightMm / 10]);
set(fig, 'PaperUnits', 'centimeters', ...
    'PaperPosition', [0, 0, widthMm / 10, heightMm / 10], ...
    'PaperSize', [widthMm / 10, heightMm / 10]);
palette = inv_log_2_nature_palette();
natureOrder = [ ...
    palette.proposed; ...
    palette.wrong; ...
    palette.gray; ...
    palette.lightBlue; ...
    palette.purple; ...
    palette.orange; ...
    palette.green];

allAxes = findall(fig, 'Type', 'axes');
for iax = 1:numel(allAxes)
    set(allAxes(iax), ...
        'FontName', 'Times New Roman', ...
        'FontSize', 9, ...
        'LineWidth', 0.8, ...
        'ColorOrder', natureOrder, ...
        'Box', 'on', ...
        'TickDir', 'in', ...
        'TickLength', [0.015, 0.015], ...
        'Layer', 'top', ...
        'XColor', palette.axis, ...
        'YColor', palette.axis);
    allAxes(iax).XLabel.FontName = 'Times New Roman';
    allAxes(iax).YLabel.FontName = 'Times New Roman';
    allAxes(iax).Title.FontName = 'Times New Roman';
    allAxes(iax).XLabel.FontSize = 9;
    allAxes(iax).YLabel.FontSize = 9;
    allAxes(iax).Title.FontSize = 9;
    allAxes(iax).Title.FontWeight = 'normal';
end

allText = findall(fig, 'Type', 'text');
for it = 1:numel(allText)
    if allText(it).FontSize < 9
        set(allText(it), 'FontName', 'Times New Roman', 'FontSize', 9);
    else
        set(allText(it), 'FontName', 'Times New Roman');
    end
end

allLegends = findall(fig, 'Type', 'legend');
for il = 1:numel(allLegends)
    set(allLegends(il), 'FontName', 'Times New Roman', 'FontSize', 9, ...
        'Interpreter', 'tex', 'Box', 'off');
end

allColorbars = findall(fig, 'Type', 'colorbar');
for ic = 1:numel(allColorbars)
    set(allColorbars(ic), 'FontName', 'Times New Roman', 'FontSize', 9);
end
end
