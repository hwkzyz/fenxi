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

allAxes = findall(fig, 'Type', 'axes');
for iax = 1:numel(allAxes)
    set(allAxes(iax), ...
        'FontName', 'Times New Roman', ...
        'FontSize', 8.5, ...
        'LineWidth', 0.8, ...
        'Box', 'on', ...
        'TickDir', 'in', ...
        'TickLength', [0.015, 0.015], ...
        'Layer', 'top');
end

allText = findall(fig, 'Type', 'text');
for it = 1:numel(allText)
    set(allText(it), 'FontName', 'Times New Roman', 'FontSize', 8.5);
end

allLegends = findall(fig, 'Type', 'legend');
for il = 1:numel(allLegends)
    set(allLegends(il), 'FontName', 'Times New Roman', 'FontSize', 8, ...
        'Interpreter', 'none');
end

allColorbars = findall(fig, 'Type', 'colorbar');
for ic = 1:numel(allColorbars)
    set(allColorbars(ic), 'FontName', 'Times New Roman', 'FontSize', 8);
end
end
