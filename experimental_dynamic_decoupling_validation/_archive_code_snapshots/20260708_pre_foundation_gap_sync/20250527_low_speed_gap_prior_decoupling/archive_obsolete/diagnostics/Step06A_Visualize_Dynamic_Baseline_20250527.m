%% Step06A: Visualize dynamic continuous-waveform baseline estimation
% Diagnostic only.  This script visualizes the same dynamic baseline logic
% used by Step06: for each sensor, use the first two revolutions of the
% continuous waveform, reject pulse regions, trim each no-pulse segment
% inward, and compute the baseline as the mean of the remaining background.

clear; clc; close all;

%% 1. Paths and parameters
thisDir = fileparts(mfilename('fullpath'));
outDir = fullfile(thisDir, 'outputs');
diagDir = fullfile(outDir, 'Step06A_dynamic_baseline_diagnostics');
if ~exist(diagDir, 'dir')
    mkdir(diagDir);
end

highMapFile = fullfile(outDir, 'Step06_HighMap_20250527_B1_S136.mat');
if ~isfile(highMapFile)
    error('Run Step06 first. Missing file: %s', highMapFile);
end

selectedSensors = [1, 3, 6];
edgeFraction = 0.10;
voltageAutoScaleThreshold = 20;
dynamicBaselineRevCount = 2;
dynamicBaselinePulsePadS = 8e-4;
dynamicBaselinePulseThresholdFraction = 0.12;
dynamicBaselineSegmentTrimFraction = 0.12;
dynamicBaselineMinSegmentS = 3e-4;
maxBackgroundPlotPoints = 12000;

H = load(highMapFile, 'highMap');
highMap = H.highMap;
experiment_case = highMap.experiment_case;
rotFreqHz = highMap.rotFreqHz;
fs = experiment_case.Config.pinlv;

%% 2. Continuous-waveform baseline per sensor and window edge check per lap
baselineRows = {};
edgeRows = {};
detail = struct();
for is = 1:numel(selectedSensors)
    sid = selectedSensors(is);
    if sid > numel(experiment_case.Raw_Stream) || isempty(experiment_case.Raw_Stream(sid).T)
        continue;
    end
    tRaw = experiment_case.Raw_Stream(sid).T(:);
    vRaw = experiment_case.Raw_Stream(sid).V(:);
    [baselineInfo, backgroundMask, windowMask, pulseMask] = estimate_dynamic_sensor_baseline_local( ...
        tRaw, vRaw, fs, rotFreqHz, dynamicBaselineRevCount, ...
        dynamicBaselinePulsePadS, dynamicBaselinePulseThresholdFraction, ...
        dynamicBaselineSegmentTrimFraction, dynamicBaselineMinSegmentS, ...
        voltageAutoScaleThreshold);

    baselineRows{end+1, 1} = table(sid, baselineInfo.baselineOriginalUnit, ...
        baselineInfo.baselineMv, baselineInfo.backgroundMedianMv, ...
        baselineInfo.backgroundFraction, baselineInfo.windowStartS, ...
        baselineInfo.windowEndS, baselineInfo.windowPeakCount, baselineInfo.voltageScale, ...
        'VariableNames', {'sensorId', 'baselineOriginalUnit', 'baselineMv', ...
        'backgroundMedianMv', 'backgroundFraction', 'baselineWindowStartS', ...
        'baselineWindowEndS', 'baselineWindowPeakCount', 'voltageScale'});

    detail.(sprintf('CH%d', sid)) = struct('sensorId', sid, 'tRaw', tRaw, ...
        'vRawMv', baselineInfo.voltageScale * vRaw, ...
        'backgroundMask', backgroundMask, 'windowMask', windowMask, ...
        'pulseMask', pulseMask, 'baselineInfo', baselineInfo);

    dataIdx = find([experiment_case.Extracted_Data.sensor_id] == sid, 1, 'first');
    if isempty(dataIdx)
        continue;
    end
    laps = experiment_case.Extracted_Data(dataIdx).Laps;
    for ilap = 1:numel(laps)
        D = laps(ilap);
        v = D.v_points(:);
        e = max(3, round(edgeFraction * numel(v)));
        firstMedianMv = baselineInfo.voltageScale * median(v(1:e), 'omitnan');
        lastMedianMv = baselineInfo.voltageScale * median(v(end-e+1:end), 'omitnan');
        edgeRows{end+1, 1} = table(sid, ilap, min(D.t_points), max(D.t_points), ...
            firstMedianMv, lastMedianMv, abs(firstMedianMv - lastMedianMv), ...
            baselineInfo.baselineMv, ...
            'VariableNames', {'sensorId', 'lapId', 'timeStart', 'timeEnd', ...
            'firstEdgeMedianMv', 'lastEdgeMedianMv', 'edgeMedianDifferenceMv', ...
            'continuousBaselineMv'});
    end
end

baselineTable = vertcat(baselineRows{:});
edgeCheckTable = vertcat(edgeRows{:});
baselineCsv = fullfile(diagDir, 'Step06A_Dynamic_Continuous_Baseline_Summary.csv');
edgeCsv = fullfile(diagDir, 'Step06A_Dynamic_Window_Edge_Check.csv');
writetable(baselineTable, baselineCsv);
writetable(edgeCheckTable, edgeCsv);

%% 3. Detailed continuous-baseline figure
fig = figure('Name', 'Step06A dynamic continuous baseline', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 28, 18], 'NumberTitle', 'off');
tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for is = 1:numel(selectedSensors)
    sid = selectedSensors(is);
    key = sprintf('CH%d', sid);
    if ~isfield(detail, key)
        continue;
    end
    D = detail.(key);
    t = D.tRaw;
    y = D.vRawMv;
    win = D.windowMask;
    bgIdx = find(D.backgroundMask);
    bgIdx = thin_indices_local(bgIdx, maxBackgroundPlotPoints);

    nexttile; hold on; box on;
    plot(t(win), y(win), 'k-', 'LineWidth', 0.8, 'DisplayName', 'raw continuous');
    scatter(t(bgIdx), y(bgIdx), 8, [0.10, 0.45, 0.80], 'filled', ...
        'DisplayName', 'trimmed no-pulse');
    yline(D.baselineInfo.baselineMv, '-', ...
        sprintf('mean %.1f mV', D.baselineInfo.baselineMv), ...
        'Color', [0.10, 0.35, 0.75], 'LineWidth', 1.1);
    yline(D.baselineInfo.backgroundMedianMv, '--', 'median', ...
        'Color', [0.55, 0.35, 0.20], 'LineWidth', 1.0);
    xlabel('Time (s)');
    ylabel('Raw voltage (mV)');
    title(sprintf('CH%d first %d rev baseline window', sid, dynamicBaselineRevCount));
    legend('Location', 'best');

    nexttile; hold on; box on;
    plot(t(win), y(win) - D.baselineInfo.baselineMv, ...
        'Color', [0.10, 0.45, 0.80], 'LineWidth', 0.8);
    yline(0, '--k', 'corrected baseline');
    xlabel('Time (s)');
    ylabel('Corrected voltage (mV)');
    title(sprintf('CH%d after continuous-baseline subtraction', sid));
end

set(findall(fig, '-property', 'FontName'), 'FontName', 'Times New Roman');
set(findall(fig, '-property', 'FontSize'), 'FontSize', 9);
detailFile = fullfile(diagDir, 'Step06A_Dynamic_Continuous_Baseline_Detailed.png');
exportgraphics(fig, detailFile, 'Resolution', 300);

%% 4. Baseline and edge-check trends
fig2 = figure('Name', 'Step06A dynamic baseline trends', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [3, 3, 20, 13], 'NumberTitle', 'off');
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile; hold on; box on;
bar(categorical("CH" + string(baselineTable.sensorId)), baselineTable.baselineMv);
ylabel('Baseline (mV)');
title('Continuous no-pulse mean baseline');

nexttile; hold on; box on;
for sid = selectedSensors
    mask = edgeCheckTable.sensorId == sid;
    plot(edgeCheckTable.lapId(mask), edgeCheckTable.edgeMedianDifferenceMv(mask), 'o-', ...
        'LineWidth', 1.0, 'DisplayName', sprintf('CH%d', sid));
end
yline(50, '--k', '50 mV check line');
xlabel('Lap index');
ylabel('|first edge - last edge| (mV)');
title('Pulse-window edge consistency check only');
legend('Location', 'best');

set(findall(fig2, '-property', 'FontName'), 'FontName', 'Times New Roman');
set(findall(fig2, '-property', 'FontSize'), 'FontSize', 9);
trendFile = fullfile(diagDir, 'Step06A_Dynamic_Continuous_Baseline_Trends.png');
exportgraphics(fig2, trendFile, 'Resolution', 300);

fprintf('\nStep06A complete.\n');
fprintf('Detailed continuous-baseline figure:\n  %s\n', detailFile);
fprintf('Trend figure:\n  %s\n', trendFile);
fprintf('Baseline table:\n  %s\n', baselineCsv);
fprintf('Edge-check table:\n  %s\n', edgeCsv);
disp(baselineTable);

%% Local functions
function [baselineInfo, backgroundMask, windowMask, pulseVicinity] = estimate_dynamic_sensor_baseline_local(tRaw, vRaw, fs, rotFreqHz, revCount, pulsePadS, pulseThresholdFraction, segmentTrimFraction, minSegmentS, voltageAutoScaleThreshold)
tRaw = tRaw(:);
vRaw = vRaw(:);
windowStart = tRaw(1);
windowEnd = min(tRaw(end), windowStart + revCount / rotFreqHz);
windowMask = tRaw >= windowStart & tRaw <= windowEnd & isfinite(vRaw);
vWin = vRaw(windowMask);
fullMedian = median(vWin, 'omitnan');
upper = prctile(vWin, 99.5);
threshold = fullMedian + pulseThresholdFraction * max(upper - fullMedian, eps);
pulseCore = vRaw > threshold & windowMask;
padN = max(1, round(pulsePadS * fs));
pulseVicinity = conv(double(pulseCore), ones(2 * padN + 1, 1), 'same') > 0;
backgroundMask0 = windowMask & ~pulseVicinity & isfinite(vRaw);
backgroundMask = trim_background_segments_local(backgroundMask0, fs, segmentTrimFraction, minSegmentS);
if nnz(backgroundMask) < 0.10 * nnz(backgroundMask0)
    backgroundMask = backgroundMask0;
end
bg = vRaw(backgroundMask);
baselineOriginal = mean(bg, 'omitnan');
baselineMedianOriginal = median(bg, 'omitnan');
voltageScale = 1;
if max(abs(vWin - baselineOriginal), [], 'omitnan') < voltageAutoScaleThreshold
    voltageScale = 1000;
end
[~, locs] = findpeaks(vWin, 'MinPeakHeight', threshold, ...
    'MinPeakDistance', max(1, round(0.5 / (rotFreqHz * 6) * fs)));
baselineInfo = struct();
baselineInfo.baselineOriginalUnit = baselineOriginal;
baselineInfo.baselineMv = voltageScale * baselineOriginal;
baselineInfo.backgroundMedianMv = voltageScale * baselineMedianOriginal;
baselineInfo.backgroundFraction = nnz(backgroundMask) / nnz(windowMask);
baselineInfo.windowStartS = windowStart;
baselineInfo.windowEndS = windowEnd;
baselineInfo.windowPeakCount = numel(locs);
baselineInfo.voltageScale = voltageScale;
end

function trimmedMask = trim_background_segments_local(mask, fs, trimFraction, minSegmentS)
mask = mask(:);
trimmedMask = false(size(mask));
idx = find(mask);
if isempty(idx)
    return;
end
breaks = [1; find(diff(idx) > 1) + 1; numel(idx) + 1];
minSegmentN = max(3, round(minSegmentS * fs));
for iseg = 1:numel(breaks)-1
    segIdx = idx(breaks(iseg):breaks(iseg+1)-1);
    n = numel(segIdx);
    if n < minSegmentN
        continue;
    end
    trimN = floor(trimFraction * n);
    if 2 * trimN >= n - 2
        trimN = max(0, floor((n - 2) / 2));
    end
    trimmedMask(segIdx((1 + trimN):(n - trimN))) = true;
end
end

function idxThin = thin_indices_local(idx, maxPoints)
idx = idx(:);
if numel(idx) <= maxPoints
    idxThin = idx;
else
    pick = unique(round(linspace(1, numel(idx), maxPoints)));
    idxThin = idx(pick);
end
end
