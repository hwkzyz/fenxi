function h = plot_inv_log_2_bar(ax, cats, values, colors, width)
%plot_inv_log_2_bar  Draw a styled bar chart with per-bar colors.

if nargin < 5 || isempty(width)
    width = 0.72;
end
if nargin < 4 || isempty(colors)
    colors = inv_log_2_nature_palette("proposed");
end

h = bar(ax, cats, values, width, ...
    'FaceColor', 'flat', ...
    'EdgeColor', [0.35, 0.35, 0.35], ...
    'LineWidth', 0.45);
if size(colors, 1) == 1
    h.CData = repmat(colors, numel(values), 1);
else
    h.CData = colors;
end
end
