function Output = Plot_AllTiltGap_WaveformFamilies()
%PLOT_ALLTILTGAP_WAVEFORMFAMILIES Plot every clearance waveform at each FE tilt.
% Each trace is a raw COMSOL export. The common axes make the clearance
% response family directly comparable across the seven independent tilt states.

thisDir = fileparts(mfilename('fullpath'));
mainDir = fileparts(thisDir);
addpath(mainDir, '-begin');
addpath(fullfile(mainDir, 'local_func'), '-begin');

gapList = (0.5:0.2:1.5).';
angleList = [0.5 1 1.5 2 2.5 3 3.5];
nAngle = numel(angleList);
outDir = fullfile(thisDir, 'output', 'waveform_families');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

xData = cell(nAngle, numel(gapList));
yData = cell(nAngle, numel(gapList));
peakX = nan(nAngle, numel(gapList));
peakY = nan(nAngle, numel(gapList));
for iAngle = 1:nAngle
    [~, xCell, yCell] = load_stacked_gap_curves(tilt_file(thisDir, angleList(iAngle)), gapList);
    xData(iAngle, :) = reshape(xCell, 1, []);
    yData(iAngle, :) = reshape(yCell, 1, []);
    for iGap = 1:numel(gapList)
        [peakX(iAngle, iGap), peakY(iAngle, iGap)] = ...
            sg_peak(xData{iAngle, iGap}, yData{iAngle, iGap});
    end
end
xLimits = enclosing_limits(xData);
yLimits = padded_limits(yData, 0.035);
colors = parula(numel(gapList));

fig = figure('Color', 'w', 'Units', 'centimeters', ...
    'Position', [2 1 17.0 13.4], 'Name', 'All tilt-gap waveform families', ...
    'NumberTitle', 'off');
layout = tiledlayout(fig, 3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
axesList = gobjects(nAngle, 1);

for iAngle = 1:nAngle
    ax = nexttile(layout, iAngle);
    axesList(iAngle) = ax;
    hold(ax, 'on');
    for iGap = 1:numel(gapList)
        plot(ax, xData{iAngle, iGap}, yData{iAngle, iGap}, ...
            'Color', colors(iGap, :), 'LineWidth', 0.95);
    end
    hold(ax, 'off');
    style_waveform_axes(ax, xLimits, yLimits, ...
        sprintf('FE tilt, \\alpha = %.1f^\\circ', angleList(iAngle)));

    if mod(iAngle-1, 3) == 0
        ylabel(ax, 'Sensor response (V)', 'FontName', 'Times New Roman', 'FontSize', 9);
    else
        ax.YTickLabel = [];
    end
    if iAngle > 3
        xlabel(ax, 'Circumferential coordinate, x (mm)', ...
            'FontName', 'Times New Roman', 'FontSize', 9);
    else
        ax.XTickLabel = [];
    end
end

keyAxes = nexttile(layout, 8);
imagesc(keyAxes, [gapList.'; gapList.']);
axis(keyAxes, 'off');
colormap(keyAxes, colors);
clim(keyAxes, [min(gapList), max(gapList)]);
cb = colorbar(keyAxes, 'Location', 'eastoutside');
cb.FontName = 'Times New Roman';
cb.FontSize = 8;
cb.LineWidth = 0.55;
cb.Ticks = gapList;
cb.TickLabels = compose('%.1f', gapList);
cb.Label.String = 'Clearance, g (mm)';
cb.Label.FontName = 'Times New Roman';
cb.Label.FontSize = 9;
text(keyAxes, 0.5, -0.35, 'Line colour', 'Units', 'normalized', ...
    'HorizontalAlignment', 'center', 'FontName', 'Times New Roman', 'FontSize', 8.5);

emptyAxes = nexttile(layout, 9);
axis(emptyAxes, 'off');
text(emptyAxes, 0.5, 0.58, {'Raw COMSOL waveforms', 'No interpolation or normalization'}, ...
    'Units', 'normalized', 'HorizontalAlignment', 'center', ...
    'FontName', 'Times New Roman', 'FontSize', 8.5);

drawnow;
pngPath = fullfile(outDir, 'FigS_AllTilt_ClearanceWaveformFamilies.png');
pdfPath = fullfile(outDir, 'FigS_AllTilt_ClearanceWaveformFamilies.pdf');
emfPath = fullfile(outDir, 'FigS_AllTilt_ClearanceWaveformFamilies.emf');
exportgraphics(fig, pngPath, 'Resolution', 600);
exportgraphics(fig, pdfPath, 'ContentType', 'vector');
print(fig, emfPath, '-dmeta');

[angleGrid, gapGrid] = ndgrid(angleList, gapList);
PeakTable = table(angleGrid(:), gapGrid(:), peakX(:), peakY(:), ...
    'VariableNames', {'tilt_deg', 'gap_mm', 'sg_peak_x_mm', 'sg_peak_response_V'});
peakCsvPath = fullfile(outDir, 'SG_peak_locations_all_tilt_gap.csv');
writetable(PeakTable, peakCsvPath);

Output = struct('figure', fig, 'png_path', string(pngPath), ...
    'pdf_path', string(pdfPath), 'emf_path', string(emfPath), ...
    'peak_table', PeakTable, 'peak_csv_path', string(peakCsvPath), ...
    'gap_list_mm', gapList, 'tilt_list_deg', angleList, ...
    'x_limits_mm', xLimits, 'y_limits_V', yLimits);
fprintf('Saved all tilt-gap COMSOL waveform families to:\n%s\n%s\n%s\n', ...
    pngPath, pdfPath, emfPath);
end

function [xPeak, yPeak] = sg_peak(x, y)
% Smooth with SG, then refine the local maximum by quadratic fitting.
x = x(:);
y = y(:);
window = min(11, numel(y));
if mod(window, 2) == 0
    window = window - 1;
end
ySmooth = smoothdata(y, 'sgolay', window);
[~, iMax] = max(ySmooth);
idx = max(1, iMax-2):min(numel(x), iMax+2);
if numel(idx) >= 3
    p = polyfit(x(idx), ySmooth(idx), 2);
    xCandidate = -p(2)/(2*p(1));
    if p(1) < 0 && xCandidate >= x(idx(1)) && xCandidate <= x(idx(end))
        xPeak = xCandidate;
        yPeak = polyval(p, xPeak);
        return;
    end
end
xPeak = x(iMax);
yPeak = ySmooth(iMax);
end

function style_waveform_axes(ax, xLimits, yLimits, heading)
set(ax, 'FontName', 'Times New Roman', 'FontSize', 8, ...
    'LineWidth', 0.65, 'TickDir', 'in', 'Box', 'on', 'Layer', 'top');
xlim(ax, xLimits);
ylim(ax, yLimits);
title(ax, heading, 'FontName', 'Times New Roman', 'FontWeight', 'normal', 'FontSize', 9);
grid(ax, 'off');
end

function limits = enclosing_limits(dataCell)
mins = cellfun(@(v) min(v(:)), dataCell);
maxs = cellfun(@(v) max(v(:)), dataCell);
limits = [min(mins(:)), max(maxs(:))];
end

function limits = padded_limits(dataCell, fraction)
mins = cellfun(@(v) min(v(:)), dataCell);
maxs = cellfun(@(v) max(v(:)), dataCell);
dataLimits = [min(mins(:)), max(maxs(:))];
padding = fraction * diff(dataLimits);
limits = dataLimits + [-padding padding];
end

function filePath = tilt_file(thisDir, angleDeg)
fileName = sprintf('直叶片2mm_不同间隙0.5_0.2_1.5_倾斜角%g.txt', angleDeg);
filePath = fullfile(thisDir, 'input_fe_tilt', fileName);
assert(exist(filePath, 'file') == 2, 'Missing COMSOL file: %s', filePath);
end
