%% Export Figure 03 plot-only assets for an all-editable-text PowerPoint
% The PNG files contain only data graphics, region fills and reference lines.
% Titles, tick labels, axis labels, legends and annotations are added as
% editable PowerPoint text boxes by Build_Figure03_RefinedEditablePPT.py.

clear; clc; close all;

scriptDir = fileparts(mfilename('fullpath'));
packageRoot = fileparts(scriptDir);
addpath(packageRoot);
cfg = Setup_Paths_20251222();

sourceFile = fullfile(cfg.paths.preparedFoundation, ...
    'step04_low_speed_template', 'LowSpeed_Template_SourceData_20251222.mat');
templateFile = cfg.files.lowSpeedTemplate;
rawFile = fullfile(cfg.paths.rawDatasetRoot, cfg.case.lowSpeedCase, '4-1-1000.mat');
outDir = fullfile(packageRoot, 'results', 'pptx', ...
    'Figure03_Composite_refined_editable_output', ...
    'split_png_elements', 'image01');
if exist(outDir, 'dir') ~= 7
    mkdir(outDir);
end

Ssrc = load(sourceFile, 'point_cloud');
Stpl = load(templateFile, 'Template');
pcIndex = find([Ssrc.point_cloud.sensor_id] == 1 & ...
    [Ssrc.point_cloud.blade_id] == 1, 1, 'first');
tplIndex = find([Stpl.Template.SensorBlade.sensor_id] == 1 & ...
    [Stpl.Template.SensorBlade.blade_id] == 1, 1, 'first');
P = Ssrc.point_cloud(pcIndex);
Tpl = Stpl.Template.SensorBlade(tplIndex);
thresholdV = Tpl.threshold;

pulseIds = unique(P.pulse_index(:), 'stable');
nPulse = numel(pulseIds);
peakV = nan(nPulse, 1);
widthMm = nan(nPulse, 1);
baselineDeltaV = nan(nPulse, 1);
for ip = 1:nPulse
    mask = P.pulse_index == pulseIds(ip);
    x = P.x_mm(mask);
    v = P.v(mask);
    valid = isfinite(x) & isfinite(v);
    x = x(valid);
    v = v(valid);
    peakV(ip) = max(v);
    high = v > thresholdV;
    if nnz(high) >= 2
        widthMm(ip) = max(x(high)) - min(x(high));
    end
    left = x < prctile(x, 15);
    right = x > prctile(x, 85);
    baselineDeltaV(ip) = abs(median(v(left), 'omitnan') - ...
        median(v(right), 'omitnan'));
end
score = robust_square_score_local(peakV) + ...
    robust_square_score_local(widthMm) + ...
    robust_square_score_local(baselineDeltaV);
[~, representativeIndex] = min(score);
representativePulseId = pulseIds(representativeIndex);

pulseMask = P.pulse_index == representativePulseId;
xPulse = P.x_mm(pulseMask);
vPulse = P.v(pulseMask);
tPulse = P.t_s(pulseMask);
[xPulse, order] = sort(xPulse(:));
vPulse = vPulse(order);
tPulse = tPulse(order);

xOuter = [-7.5, 7.5];
insideWindow = xPulse >= xOuter(1) & xPulse <= xOuter(2) & ...
    isfinite(vPulse) & isfinite(tPulse);
xPulse = xPulse(insideWindow);
vPulse = vPulse(insideWindow);
tPulse = tPulse(insideWindow);

calMask = vPulse > thresholdV;
xLeft = min(xPulse(calMask));
xRight = max(xPulse(calMask));
bufferWidthMm = 0.90;
leftBaseEnd = xLeft - bufferWidthMm;
rightBaseStart = xRight + bufferWidthMm;
baseMask = xPulse <= leftBaseEnd | xPulse >= rightBaseStart;
baselineV = median(vPulse(baseMask), 'omitnan');

xCal = xPulse(calMask);
vCal = vPulse(calMask);
[peakObserved, peakIndex] = max(vCal);
initial = [log(max(peakObserved - baselineV, 0.1)), ...
    log(max((xRight - xLeft) / 2, 0.3)), log(2.5 - 1), xCal(peakIndex)];
objective = @(q) mean((vCal - super_gaussian_local(xCal, q, baselineV)).^2);
options = optimset('Display', 'off', 'MaxIter', 1200, ...
    'MaxFunEvals', 5000, 'TolX', 1e-9, 'TolFun', 1e-10);
qFit = fminsearch(objective, initial, options);
xFit = linspace(xLeft, xRight, 900).';
vFit = super_gaussian_local(xFit, qFit, baselineV);

Sraw = load(rawFile, 'jilu01');
rawIndex = Sraw.jilu01(:, 1);
rawVoltage = Sraw.jilu01(:, 2);
rawTimeS = (rawIndex - rawIndex(1)) / cfg.machine.sampleRateHz;
selectedCenterTimeS = median(tPulse);
topHalfSpanS = 0.036;
topMask = rawTimeS >= selectedCenterTimeS - topHalfSpanS & ...
    rawTimeS <= selectedCenterTimeS + topHalfSpanS;
topTimeMs = 1000 * (rawTimeS(topMask) - selectedCenterTimeS);
topVoltage = rawVoltage(topMask);
topPlotStep = max(1, floor(numel(topTimeMs) / 12000));
topTimeMs = topTimeMs(1:topPlotStep:end);
topVoltage = topVoltage(1:topPlotStep:end);
selectedTimeWindowMs = 1000 * ([min(tPulse), max(tPulse)] - selectedCenterTimeS);

blue = [25, 103, 184] / 255;
lightBlue = [229, 241, 252] / 255;
orange = [232, 104, 24] / 255;
lightOrange = [255, 243, 224] / 255;
lightGray = [241, 243, 246] / 255;
dark = [0.12, 0.14, 0.18];

%% Top plot: no raster text
figTop = figure('Color', 'none', 'Units', 'pixels', ...
    'Position', [100, 100, 1800, 390], 'Visible', 'off', ...
    'Renderer', 'painters');
axTop = axes(figTop, 'Position', [0.015, 0.025, 0.97, 0.95]);
hold(axTop, 'on');
ylTop = [-0.55, 4.20];
patch(axTop, [selectedTimeWindowMs(1), selectedTimeWindowMs(2), ...
    selectedTimeWindowMs(2), selectedTimeWindowMs(1)], ...
    [ylTop(1), ylTop(1), ylTop(2), ylTop(2)], lightBlue, ...
    'EdgeColor', blue, 'LineStyle', '--', 'FaceAlpha', 0.52);
plot(axTop, topTimeMs, topVoltage, '-', 'Color', dark, 'LineWidth', 0.75);
yline(axTop, thresholdV, '--', 'Color', orange, 'LineWidth', 1.0);
xlim(axTop, [-36, 36]);
ylim(axTop, ylTop);
set(axTop, 'XTick', -30:10:30, 'XTickLabel', [], ...
    'YTick', 0:4, 'YTickLabel', [], 'FontSize', 1, ...
    'LineWidth', 0.85, 'TickDir', 'in', 'Box', 'on', ...
    'Layer', 'top', 'Color', 'none');
exportgraphics(figTop, fullfile(outDir, ...
    'image01_chart_top_data_only.png'), ...
    'Resolution', 220, 'BackgroundColor', 'none');
close(figTop);

%% Lower plot: no raster text
figLower = figure('Color', 'none', 'Units', 'pixels', ...
    'Position', [100, 100, 1800, 455], 'Visible', 'off', ...
    'Renderer', 'painters');
axPulse = axes(figLower, 'Position', [0.015, 0.025, 0.97, 0.95]);
hold(axPulse, 'on');
yl = [-0.55, 4.15];
patch_region_local(axPulse, xOuter(1), leftBaseEnd, yl, lightGray);
patch_region_local(axPulse, leftBaseEnd, xLeft, yl, lightOrange);
patch_region_local(axPulse, xLeft, xRight, yl, lightBlue);
patch_region_local(axPulse, xRight, rightBaseStart, yl, lightOrange);
patch_region_local(axPulse, rightBaseStart, xOuter(2), yl, lightGray);
displayStep = max(1, floor(numel(xPulse) / 900));
displayIndex = 1:displayStep:numel(xPulse);
plot(axPulse, xPulse(displayIndex), vPulse(displayIndex), 'o', ...
    'Color', [0.40, 0.40, 0.40], 'MarkerFaceColor', 'w', ...
    'MarkerSize', 2.6, 'LineWidth', 0.55);
calDisplayIndex = find(calMask);
calDisplayIndex = calDisplayIndex(1:max(1, ...
    floor(numel(calDisplayIndex) / 180)):end);
plot(axPulse, xPulse(calDisplayIndex), vPulse(calDisplayIndex), 'o', ...
    'Color', orange, 'MarkerFaceColor', orange, ...
    'MarkerSize', 3.0, 'LineWidth', 0.5);
plot(axPulse, xFit, vFit, '-', 'Color', blue, 'LineWidth', 1.8);
yline(axPulse, thresholdV, '--', 'Color', orange, 'LineWidth', 0.9);
xline(axPulse, xLeft, '--', 'Color', orange, 'LineWidth', 0.9);
xline(axPulse, xRight, '--', 'Color', orange, 'LineWidth', 0.9);
xline(axPulse, leftBaseEnd, ':', 'Color', [0.55, 0.55, 0.55], 'LineWidth', 0.9);
xline(axPulse, rightBaseStart, ':', 'Color', [0.55, 0.55, 0.55], 'LineWidth', 0.9);
xlim(axPulse, xOuter);
ylim(axPulse, yl);
set(axPulse, 'XTick', -6:2:6, 'XTickLabel', [], ...
    'YTick', 0:4, 'YTickLabel', [], 'FontSize', 1, ...
    'LineWidth', 0.85, 'TickDir', 'in', 'Box', 'on', ...
    'Layer', 'top', 'Color', 'none');
exportgraphics(figLower, fullfile(outDir, ...
    'image01_chart_lower_data_only.png'), ...
    'Resolution', 220, 'BackgroundColor', 'none');
close(figLower);

save(fullfile(outDir, 'Figure03_PlotAsset_Metadata.mat'), ...
    'thresholdV', 'xLeft', 'xRight', 'leftBaseEnd', 'rightBaseStart', ...
    'xOuter', 'selectedTimeWindowMs', 'representativePulseId');

fprintf('Saved plot-only assets in: %s\n', outDir);

function score = robust_square_score_local(values)
values = values(:);
center = median(values, 'omitnan');
scale = 1.4826 * median(abs(values - center), 'omitnan');
if ~isfinite(scale) || scale <= eps
    scale = max(std(values, 'omitnan'), eps);
end
score = ((values - center) ./ scale) .^ 2;
score(~isfinite(score)) = inf;
end

function voltage = super_gaussian_local(x, q, baseline)
u0 = exp(q(1));
w = exp(q(2));
beta = 1 + exp(q(3));
xc = q(4);
voltage = baseline + u0 .* exp(-abs((x - xc) ./ w) .^ beta);
end

function patch_region_local(ax, x1, x2, yLimits, color)
patch(ax, [x1, x2, x2, x1], ...
    [yLimits(1), yLimits(1), yLimits(2), yLimits(2)], color, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.72);
end
