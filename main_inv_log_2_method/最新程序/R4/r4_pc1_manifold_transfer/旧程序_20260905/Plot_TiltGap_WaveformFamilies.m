function Output = Plot_TiltGap_WaveformFamilies()
%PLOT_TILTGAP_WAVEFORMFAMILIES Plot raw COMSOL waveform families by state.
% The left panel changes clearance at a fixed FE tilt. The right panel
% changes FE tilt at a fixed clearance. No interpolation or normalization
% is applied: every trace is a directly exported COMSOL waveform.

thisDir = fileparts(mfilename('fullpath'));
mainDir = fileparts(thisDir);
addpath(mainDir, '-begin');
addpath(fullfile(mainDir, 'local_func'), '-begin');

gapList = (0.5:0.2:1.5).';
angleList = [0.5 1 1.5 2 2.5 3 3.5];
fixedTiltDeg = 0.5;
fixedGapMm = 0.5;
outDir = fullfile(thisDir, 'output', 'waveform_families');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

[~, xByGap, yByGap] = load_stacked_gap_curves(tilt_file(thisDir, fixedTiltDeg), gapList);
gapIndex = find(abs(gapList - fixedGapMm) < 1e-12, 1, 'first');
assert(~isempty(gapIndex), 'The selected fixed clearance is not in gapList.');

xByTilt = cell(numel(angleList), 1);
yByTilt = cell(numel(angleList), 1);
for i = 1:numel(angleList)
    [~, xCell, yCell] = load_stacked_gap_curves(tilt_file(thisDir, angleList(i)), gapList);
    xByTilt{i} = xCell{gapIndex};
    yByTilt{i} = yCell{gapIndex};
end

fig = figure('Color', 'w', 'Units', 'centimeters', ...
    'Position', [2 2 17.0 7.6], 'Name', 'COMSOL waveform families', ...
    'NumberTitle', 'off');
tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

axGap = nexttile;
plot_family(axGap, xByGap, yByGap, gapList, ...
    'Clearance-dependent simulated waveforms', ...
    sprintf('FE tilt, \\alpha = %.1f^\\circ', fixedTiltDeg), ...
    'Clearance, g (mm)');
panel_label(axGap, '(a)');

axTilt = nexttile;
plot_family(axTilt, xByTilt, yByTilt, angleList, ...
    'Tilt-dependent simulated waveforms', ...
    sprintf('Clearance, g = %.1f mm', fixedGapMm), ...
    'FE tilt, \\alpha (deg)');
panel_label(axTilt, '(b)');

drawnow;
pngPath = fullfile(outDir, 'FigS_TiltGap_WaveformFamilies.png');
pdfPath = fullfile(outDir, 'FigS_TiltGap_WaveformFamilies.pdf');
exportgraphics(fig, pngPath, 'Resolution', 600);
exportgraphics(fig, pdfPath, 'ContentType', 'vector');

Output = struct('figure', fig, 'png_path', string(pngPath), ...
    'pdf_path', string(pdfPath), 'fixed_tilt_deg', fixedTiltDeg, ...
    'fixed_gap_mm', fixedGapMm, 'gap_list_mm', gapList, ...
    'tilt_list_deg', angleList);
fprintf('Saved raw COMSOL waveform families to:\n%s\n%s\n', pngPath, pdfPath);
end

function plot_family(ax, xCell, yCell, values, heading, conditionText, colorLabel)
colors = parula(numel(values));
hold(ax, 'on');
for i = 1:numel(values)
    plot(ax, xCell{i}, yCell{i}, 'Color', colors(i, :), 'LineWidth', 1.15);
end
hold(ax, 'off');

set(ax, 'FontName', 'Times New Roman', 'FontSize', 8.5, ...
    'LineWidth', 0.65, 'TickDir', 'out', 'Box', 'on', 'Layer', 'top');
xlabel(ax, 'Circumferential coordinate, x (mm)', 'FontName', 'Times New Roman', 'FontSize', 9);
ylabel(ax, 'Simulated sensor response (V)', 'FontName', 'Times New Roman', 'FontSize', 9);
title(ax, {heading, conditionText}, 'FontName', 'Times New Roman', ...
    'FontWeight', 'normal', 'FontSize', 9);
xlim(ax, common_x_limits(xCell));
grid(ax, 'off');

colormap(ax, colors);
clim(ax, [min(values), max(values)]);
cb = colorbar(ax);
cb.FontName = 'Times New Roman';
cb.FontSize = 8;
cb.LineWidth = 0.55;
cb.Label.String = colorLabel;
cb.Label.FontName = 'Times New Roman';
cb.Label.FontSize = 8.5;
cb.Ticks = values;
cb.TickLabels = compose('%.1f', values);
end

function limits = common_x_limits(xCell)
limits = [max(cellfun(@min, xCell)), min(cellfun(@max, xCell))];
assert(limits(2) > limits(1), 'The selected COMSOL curves have no common x-domain.');
end

function panel_label(ax, labelText)
text(ax, -0.17, 1.10, labelText, 'Units', 'normalized', ...
    'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', ...
    'Clipping', 'off');
end

function filePath = tilt_file(thisDir, angleDeg)
fileName = sprintf('直叶片2mm_不同间隙0.5_0.2_1.5_倾斜角%g.txt', angleDeg);
filePath = fullfile(thisDir, 'input_fe_tilt', fileName);
assert(exist(filePath, 'file') == 2, 'Missing COMSOL file: %s', filePath);
end
