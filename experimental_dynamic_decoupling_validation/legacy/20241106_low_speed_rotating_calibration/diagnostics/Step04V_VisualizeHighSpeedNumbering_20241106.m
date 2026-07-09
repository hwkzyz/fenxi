%% Step04V_VisualizeHighSpeedNumbering_20241106
% Visualize the saved high-speed numbering artifact. This is separate from
% Step04 because Step04 reads raw high-speed data and is comparatively slow.

clear; close all; clc;

thisFileDir = fileparts(mfilename('fullpath'));
routeDir = fileparts(thisFileDir);
addpath(routeDir);

%% Visualization options
showIdentificationWindowWaveforms = true;
windowIdsToShow = [];      % [] => auto-pick first/middle/last planned identification windows
maxWindowCountToShow = 3;
showPulseLabels = true;
showCrossSensorLapOverlay = true;
showLowVsHighFingerprintComparison = true;
showPeakFitAudit = true;

P = NewFlow_Config_20241106();
if exist(P.files.highSpeedNumbering, 'file') ~= 2
    error('Missing Step04 artifact:\n  %s', P.files.highSpeedNumbering);
end
if exist(P.files.lowSpeedFingerprint, 'file') ~= 2
    error('Missing Step01 artifact:\n  %s', P.files.lowSpeedFingerprint);
end

loaded = load(P.files.highSpeedNumbering, 'HighSpeedNumbering');
HighSpeedNumbering = loaded.HighSpeedNumbering;
loadedRef = load(P.files.lowSpeedFingerprint, 'LowSpeedReference');
LowSpeedReference = loadedRef.LowSpeedReference;

figDir = fullfile(P.view.figureDir, '04_high_speed_numbering');
if P.view.saveFigures && exist(figDir, 'dir') ~= 7
    mkdir(figDir);
end

sensorIds = [HighSpeedNumbering.sensor.sensor_id];
bladeIds = 1:HighSpeedNumbering.blade_count;
physicalToLocal = nan(numel(sensorIds), numel(bladeIds));
scoreMat = nan(numel(sensorIds), numel(bladeIds));
for is = 1:numel(sensorIds)
    physicalToLocal(is, :) = HighSpeedNumbering.sensor(is).physical_to_local_blade_ids;
    if isfield(HighSpeedNumbering.sensor(is), 'revolution_shift_scores') && ...
            ~isempty(HighSpeedNumbering.sensor(is).revolution_shift_scores)
        scoreMat(is, :) = mean(HighSpeedNumbering.sensor(is).revolution_shift_scores, 1, 'omitnan');
    end
end

fig = figure('Name', 'Step04 high-speed numbering audit', 'Color', 'w');
tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
imagesc(bladeIds, sensorIds, physicalToLocal);
set(gca, 'YDir', 'normal');
xlabel('Physical blade ID');
ylabel('Sensor ID');
title('Physical blade -> local pulse slot');
colorbar;
box on;

nexttile;
plot(0:(HighSpeedNumbering.blade_count - 1), scoreMat.', '-o', 'LineWidth', 1.2);
xlabel('Cyclic shift');
ylabel('Fingerprint correlation');
title('Seed-lap shift scores');
legend(compose('CH%d', sensorIds), 'Location', 'best');
grid on; box on;

nexttile;
bar(categorical(compose('CH%d', sensorIds)), [HighSpeedNumbering.sensor.best_score]);
ylabel('Best score');
title('Numbering confidence');
grid on; box on;

nexttile;
bar(categorical(compose('CH%d', sensorIds)), [HighSpeedNumbering.sensor.score_margin]);
ylabel('Best - second-best score');
title('Score margin');
grid on; box on;

if P.view.saveFigures
    exportgraphics(fig, fullfile(figDir, 'Step04V_HighSpeedNumbering_20241106.png'), 'Resolution', 300);
end

figRows = figure('Name', 'Step04 selected row propagation', 'Color', 'w');
tiledlayout(figRows, numel(sensorIds), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for is = 1:numel(sensorIds)
    nexttile;
    rows = HighSpeedNumbering.sensor(is).selected_rows_by_physical;
    plot(1:size(rows, 1), rows, '-o', 'LineWidth', 1.0);
    xlabel('Propagated lap index');
    ylabel('Pulse row');
    title(sprintf('CH%d row propagation by physical blade', sensorIds(is)));
    legend(compose('B%d', bladeIds), 'Location', 'eastoutside');
    grid on; box on;
end
if P.view.saveFigures
    exportgraphics(figRows, fullfile(figDir, 'Step04V_SelectedRowsByPhysicalBlade_20241106.png'), 'Resolution', 300);
end

if isfield(HighSpeedNumbering.sensor, 'revolution_shift_scores')
    figSeed = figure('Name', 'Step04 candidate-lap locking audit', 'Color', 'w');
    tiledlayout(figSeed, numel(sensorIds), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    for is = 1:numel(sensorIds)
        nexttile;
        scoreByLap = HighSpeedNumbering.sensor(is).revolution_shift_scores;
        imagesc(0:(HighSpeedNumbering.blade_count - 1), 1:size(scoreByLap, 1), scoreByLap);
        xlabel('Cyclic shift');
        ylabel('Selected revolution');
        title(sprintf('CH%d per-revolution shift scores', sensorIds(is)));
        colorbar;
        box on;
    end
    if P.view.saveFigures
        exportgraphics(figSeed, fullfile(figDir, 'Step04V_CandidateLapShiftScores_20241106.png'), 'Resolution', 300);
    end
end

if showIdentificationWindowWaveforms
    plot_identification_window_waveforms_local( ...
        HighSpeedNumbering, P, figDir, windowIdsToShow, maxWindowCountToShow, ...
        showPulseLabels, showCrossSensorLapOverlay, ...
        showLowVsHighFingerprintComparison, showPeakFitAudit, LowSpeedReference);
end

fprintf('\n=== Step04V: visualize high-speed numbering ===\n');
fprintf('Source: %s\n', P.files.highSpeedNumbering);
fprintf('Figures: %s\n', figDir);

function plot_identification_window_waveforms_local( ...
        HighSpeedNumbering, P, figDir, windowIdsToShow, maxWindowCountToShow, ...
        showPulseLabels, showCrossSensorLapOverlay, ...
        showLowVsHighFingerprintComparison, showPeakFitAudit, LowSpeedReference)
sensorIds = [HighSpeedNumbering.sensor.sensor_id];
bladeIds = 1:HighSpeedNumbering.blade_count;
lapCount = numel(HighSpeedNumbering.selected_revolution_ids);
if lapCount < P.waveform.windowLaps
    warning('Step04V: selected laps (%d) are fewer than Step05 windowLaps (%d).', ...
        lapCount, P.waveform.windowLaps);
    return;
end

numWindows = floor((lapCount - P.waveform.windowLaps) / P.waveform.slidingStepLaps) + 1;
windowPlan = build_identification_window_plan_local(lapCount, P.waveform.windowLaps, ...
    P.waveform.slidingStepLaps, numWindows);
windowShowIds = choose_window_ids_local(numWindows, windowIdsToShow, maxWindowCountToShow);
windowPlan = windowPlan(windowShowIds);

probeCache = load_probe_cache_local(P, sensorIds);
timeWindow = build_window_plot_time_window_local(HighSpeedNumbering, probeCache, windowPlan, P);
fileRanges = build_case_file_ranges_local(P.data.highSpeedDir, P.machine.oprChannel, P.machine.sampleRateHz);
raw = load_raw_subset_local(P.data.highSpeedDir, fileRanges, sensorIds, P.machine.sampleRateHz, timeWindow);
bladeColors = lines(HighSpeedNumbering.blade_count);

fig = figure('Name', 'Step04 identification-window waveform numbering audit', 'Color', 'w');
tiledlayout(fig, numel(sensorIds), numel(windowPlan), 'TileSpacing', 'compact', 'Padding', 'compact');

for is = 1:numel(sensorIds)
    sid = sensorIds(is);
    sensorNumbering = HighSpeedNumbering.sensor(is);
    probe = probeCache(is);
    for iw = 1:numel(windowPlan)
        W = windowPlan(iw);
        nexttile;
        render_window_waveform_local(raw(sid), probe, sensorNumbering, bladeIds, bladeColors, W, P, showPulseLabels, sid);
    end
end

if P.view.saveFigures
    exportgraphics(fig, fullfile(figDir, 'Step04V_IdentificationWindowWaveforms_20241106.png'), 'Resolution', 300);
end

if showCrossSensorLapOverlay
    plot_cross_sensor_lap_overlay_local( ...
        HighSpeedNumbering, P, figDir, windowPlan, probeCache, raw, bladeColors, showPulseLabels);
end
if showLowVsHighFingerprintComparison
    plot_low_vs_high_fingerprint_local( ...
        HighSpeedNumbering, LowSpeedReference, P, figDir, windowPlan, probeCache, raw, bladeColors);
end
if showPeakFitAudit
    plot_peak_fit_audit_local( ...
        HighSpeedNumbering, P, figDir, windowPlan, probeCache, raw, bladeColors);
end
end

function windowPlan = build_identification_window_plan_local(lapCount, windowLaps, slidingStepLaps, numWindows)
windowPlan = repmat(struct('window_id', NaN, 'lap_range', [], 'lap_ids', []), numWindows, 1);
for w = 1:numWindows
    lapStart = 1 + (w - 1) * slidingStepLaps;
    lapRange = lapStart:(lapStart + windowLaps - 1);
    windowPlan(w).window_id = w;
    windowPlan(w).lap_range = lapRange;
    windowPlan(w).lap_ids = lapRange;
end
if isempty(windowPlan) && lapCount >= windowLaps
    error('Failed to build identification window plan.');
end
end

function windowIds = choose_window_ids_local(numWindows, explicitIds, maxWindowCountToShow)
if ~isempty(explicitIds)
    windowIds = unique(explicitIds(:).', 'stable');
    windowIds = windowIds(windowIds >= 1 & windowIds <= numWindows);
    return;
end
if numWindows <= maxWindowCountToShow
    windowIds = 1:numWindows;
    return;
end
autoIds = unique(round(linspace(1, numWindows, maxWindowCountToShow)));
windowIds = autoIds(:).';
end

function probeCache = load_probe_cache_local(P, sensorIds)
probeCache = repmat(struct('sensor_id', NaN, 'jilublade', []), numel(sensorIds), 1);
for i = 1:numel(sensorIds)
    sid = sensorIds(i);
    probeFile = fullfile(P.data.highSpeedPulseDir, sprintf('jilublade_probe%d.mat', sid));
    if exist(probeFile, 'file') ~= 2
        error('Missing high-speed CH%d blade pulse file:\n  %s', sid, probeFile);
    end
    loaded = load(probeFile, 'jilublade');
    probeCache(i).sensor_id = sid;
    probeCache(i).jilublade = loaded.jilublade;
end
end

function timeWindow = build_window_plot_time_window_local(HighSpeedNumbering, probeCache, windowPlan, P)
timeWindow = [inf, -inf];
for is = 1:numel(HighSpeedNumbering.sensor)
    rowsByPhysical = HighSpeedNumbering.sensor(is).selected_rows_by_physical;
    jilublade = probeCache(is).jilublade;
    for iw = 1:numel(windowPlan)
        lapRange = windowPlan(iw).lap_range;
        rows = rowsByPhysical(lapRange, :);
        rows = rows(isfinite(rows));
        if isempty(rows)
            continue;
        end
        t0 = min(jilublade(rows, 1)) - P.waveform.pulsePadSec;
        t1 = max(jilublade(rows, 2)) + P.waveform.pulsePadSec;
        timeWindow(1) = min(timeWindow(1), t0);
        timeWindow(2) = max(timeWindow(2), t1);
    end
end
if ~all(isfinite(timeWindow))
    error('Could not build identification-window waveform plot time window.');
end
end

function render_window_waveform_local(rawSensor, probe, sensorNumbering, bladeIds, bladeColors, W, P, showPulseLabels, sid)
lapRange = W.lap_range;
rowsByPhysical = sensorNumbering.selected_rows_by_physical(lapRange, :);
rows = rowsByPhysical(isfinite(rowsByPhysical));
if isempty(rows)
    title(sprintf('CH%d W%02d: no rows', sid, W.window_id));
    axis off;
    return;
end

t0 = min(probe.jilublade(rows, 1)) - P.waveform.pulsePadSec;
t1 = max(probe.jilublade(rows, 2)) + P.waveform.pulsePadSec;
keep = rawSensor.T >= t0 & rawSensor.T <= t1;
tPlot = rawSensor.T(keep);
vPlot = rawSensor.V(keep);
plot((tPlot - t0) * 1000, vPlot, 'k-', 'LineWidth', 0.9);
hold on; box on; grid on;

lapStartRelMs = nan(numel(lapRange), 1);
lapEndRelMs = nan(numel(lapRange), 1);
for iLap = 1:numel(lapRange)
    lapRows = rowsByPhysical(iLap, :);
    lapRows = lapRows(isfinite(lapRows));
    if isempty(lapRows)
        continue;
    end
    lapStartRelMs(iLap) = (min(probe.jilublade(lapRows, 1)) - t0) * 1000;
    lapEndRelMs(iLap) = (max(probe.jilublade(lapRows, 2)) - t0) * 1000;
    if iLap > 1
        xline(lapStartRelMs(iLap), ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);
    end
end

yTop = max(vPlot, [], 'omitnan');
yBottom = min(vPlot, [], 'omitnan');
ySpan = max(yTop - yBottom, eps);
textY = yTop + 0.08 * ySpan;

for iLap = 1:numel(lapRange)
    if isfinite(lapStartRelMs(iLap))
        text(lapStartRelMs(iLap), yTop + 0.18 * ySpan, sprintf('Lap %d', lapRange(iLap)), ...
            'FontSize', 8, 'Color', [0.35 0.35 0.35], 'VerticalAlignment', 'bottom');
    end
    for bladeId = bladeIds
        rowId = rowsByPhysical(iLap, bladeId);
        if ~isfinite(rowId)
            continue;
        end
        tPeakMs = (probe.jilublade(rowId, 3) - t0) * 1000;
        vPeak = interpolate_peak_local(tPlot, vPlot, probe.jilublade(rowId, 3));
        plot(tPeakMs, vPeak, 'o', 'MarkerSize', 5.5, ...
            'MarkerFaceColor', bladeColors(bladeId, :), ...
            'MarkerEdgeColor', bladeColors(bladeId, :));
        if showPulseLabels
            text(tPeakMs, textY, sprintf('B%d', bladeId), ...
                'HorizontalAlignment', 'center', 'FontSize', 8, ...
                'Color', bladeColors(bladeId, :), 'FontWeight', 'bold');
        end
    end
end

xlabel('Time in window (ms)');
ylabel('Voltage');
title(sprintf('CH%d W%02d laps %d-%d', sid, W.window_id, lapRange(1), lapRange(end)));
ylim([yBottom - 0.05 * ySpan, yTop + 0.30 * ySpan]);
end

function plot_cross_sensor_lap_overlay_local( ...
        HighSpeedNumbering, P, figDir, windowPlan, probeCache, raw, bladeColors, showPulseLabels)
sensorIds = [HighSpeedNumbering.sensor.sensor_id];
bladeIds = 1:HighSpeedNumbering.blade_count;
lapPairs = [];
for iw = 1:numel(windowPlan)
    for il = 1:numel(windowPlan(iw).lap_range)
        lapPairs = [lapPairs; windowPlan(iw).window_id, windowPlan(iw).lap_range(il)]; %#ok<AGROW>
    end
end
if isempty(lapPairs)
    return;
end

for iPair = 1:size(lapPairs, 1)
    windowId = lapPairs(iPair, 1);
    lapId = lapPairs(iPair, 2);
    fig = figure('Name', sprintf('Step04 cross-sensor lap overlay W%02d Lap %d', windowId, lapId), ...
        'Color', 'w', 'Position', [100, 100, 1500, 420]);
    render_cross_sensor_single_lap_local( ...
        HighSpeedNumbering, probeCache, raw, sensorIds, bladeIds, bladeColors, ...
        lapId, windowId, P, showPulseLabels);
    if P.view.saveFigures
        exportgraphics(fig, fullfile(figDir, ...
            sprintf('Step04V_CrossSensorLapOverlay_W%02d_Lap%02d_20241106.png', windowId, lapId)), ...
            'Resolution', 300);
    end
end
end

function plot_low_vs_high_fingerprint_local( ...
        HighSpeedNumbering, LowSpeedReference, P, figDir, windowPlan, probeCache, raw, bladeColors)
sensorIds = [HighSpeedNumbering.sensor.sensor_id];
bladeIds = 1:HighSpeedNumbering.blade_count;
lapPairs = [];
for iw = 1:numel(windowPlan)
    for il = 1:numel(windowPlan(iw).lap_range)
        lapPairs = [lapPairs; windowPlan(iw).window_id, windowPlan(iw).lap_range(il)]; %#ok<AGROW>
    end
end
if isempty(lapPairs)
    return;
end

for iPair = 1:size(lapPairs, 1)
    windowId = lapPairs(iPair, 1);
    lapId = lapPairs(iPair, 2);
    fig = figure('Name', sprintf('Step04 low-vs-high fingerprint W%02d Lap %d', windowId, lapId), ...
        'Color', 'w');
    tiledlayout(fig, 2, ceil(numel(sensorIds) / 2), 'TileSpacing', 'compact', 'Padding', 'compact');
    for is = 1:numel(sensorIds)
        sid = sensorIds(is);
        nexttile;
        lowVec = read_low_speed_fingerprint_local(LowSpeedReference, sid, HighSpeedNumbering.blade_count);
        highVec = build_high_speed_physical_peak_vector_local( ...
            raw(sid), probeCache(is).jilublade, HighSpeedNumbering.sensor(is), lapId, P);
        render_low_vs_high_fingerprint_local(lowVec, highVec, bladeIds, bladeColors, sid, windowId, lapId);
    end
    if P.view.saveFigures
        exportgraphics(fig, fullfile(figDir, ...
            sprintf('Step04V_LowVsHighFingerprint_W%02d_Lap%02d_20241106.png', windowId, lapId)), ...
            'Resolution', 300);
    end
end
end

function plot_peak_fit_audit_local( ...
        HighSpeedNumbering, P, figDir, windowPlan, probeCache, raw, bladeColors)
sensorIds = [HighSpeedNumbering.sensor.sensor_id];
bladeIds = 1:HighSpeedNumbering.blade_count;
lapPairs = [];
for iw = 1:numel(windowPlan)
    for il = 1:numel(windowPlan(iw).lap_range)
        lapPairs = [lapPairs; windowPlan(iw).window_id, windowPlan(iw).lap_range(il)]; %#ok<AGROW>
    end
end
if isempty(lapPairs)
    return;
end

for iPair = 1:size(lapPairs, 1)
    windowId = lapPairs(iPair, 1);
    lapId = lapPairs(iPair, 2);
    fig = figure('Name', sprintf('Step04 peak-fit audit W%02d Lap %d', windowId, lapId), ...
        'Color', 'w', 'Position', [80, 60, 1800, 900]);
    tiledlayout(fig, numel(sensorIds), numel(bladeIds), 'TileSpacing', 'compact', 'Padding', 'compact');
    for is = 1:numel(sensorIds)
        sid = sensorIds(is);
        for bladeId = bladeIds
            nexttile;
            render_peak_fit_single_pulse_local( ...
                raw(sid), probeCache(is).jilublade, HighSpeedNumbering.sensor(is), ...
                lapId, bladeId, sid, windowId, P, bladeColors(bladeId, :));
        end
    end
    if P.view.saveFigures
        exportgraphics(fig, fullfile(figDir, ...
            sprintf('Step04V_PeakFitAudit_W%02d_Lap%02d_20241106.png', windowId, lapId)), ...
            'Resolution', 300);
    end
end
end

function render_cross_sensor_single_lap_local( ...
        HighSpeedNumbering, probeCache, raw, sensorIds, bladeIds, bladeColors, ...
        lapId, windowId, P, showPulseLabels)
sensorCount = numel(sensorIds);
rowsBySensor = nan(sensorCount, numel(bladeIds));
tStart = inf;
tEnd = -inf;

for is = 1:sensorCount
    rowsThis = HighSpeedNumbering.sensor(is).selected_rows_by_physical(lapId, :);
    rowsBySensor(is, :) = rowsThis;
    rowsFinite = rowsThis(isfinite(rowsThis));
    if isempty(rowsFinite)
        continue;
    end
    jilublade = probeCache(is).jilublade;
    tStart = min(tStart, min(jilublade(rowsFinite, 1)) - P.waveform.pulsePadSec);
    tEnd = max(tEnd, max(jilublade(rowsFinite, 2)) + P.waveform.pulsePadSec);
end

if ~isfinite(tStart) || ~isfinite(tEnd) || tEnd <= tStart
    title(sprintf('W%02d Lap %d: no valid rows', windowId, lapId));
    axis off;
    return;
end

hold on; box on; grid on;
sensorOffsets = sensorCount:-1:1;
markerXY = nan(sensorCount, numel(bladeIds), 2);

for is = 1:sensorCount
    sid = sensorIds(is);
    keep = raw(sid).T >= tStart & raw(sid).T <= tEnd;
    tPlot = raw(sid).T(keep);
    vPlot = raw(sid).V(keep);
    if isempty(tPlot)
        continue;
    end
    vNorm = normalize_trace_local(vPlot);
    yBase = sensorOffsets(is);
    plot((tPlot - tStart) * 1000, yBase + 0.72 * vNorm, 'k-', 'LineWidth', 0.9);

    for bladeId = bladeIds
        rowId = rowsBySensor(is, bladeId);
        if ~isfinite(rowId)
            continue;
        end
        tPeak = probeCache(is).jilublade(rowId, 3);
        tPeakMs = (tPeak - tStart) * 1000;
        vPeakNorm = interpolate_peak_local(tPlot, vNorm, tPeak);
        yPeak = yBase + 0.72 * vPeakNorm;
        markerXY(is, bladeId, :) = [tPeakMs, yPeak];
        plot(tPeakMs, yPeak, 'o', 'MarkerSize', 5.5, ...
            'MarkerFaceColor', bladeColors(bladeId, :), ...
            'MarkerEdgeColor', bladeColors(bladeId, :));
        if showPulseLabels
            text(tPeakMs, yPeak + 0.10, sprintf('B%d', bladeId), ...
                'HorizontalAlignment', 'center', 'FontSize', 8, ...
                'Color', bladeColors(bladeId, :), 'FontWeight', 'bold');
        end
    end
end

for bladeId = bladeIds
    xy = squeeze(markerXY(:, bladeId, :));
    valid = all(isfinite(xy), 2);
    if nnz(valid) >= 2
        plot(xy(valid, 1), xy(valid, 2), '-', 'Color', bladeColors(bladeId, :), ...
            'LineWidth', 0.9);
    end
end

yticks(1:sensorCount);
yticklabels(compose('CH%d', sensorIds(end:-1:1)));
xlabel('Time in lap window (ms)');
ylabel('Sensor (offset)');
title(sprintf('W%02d Lap %d: cross-sensor numbered waveform overlay', windowId, lapId));
ylim([0.5, sensorCount + 1.15]);
xlim([0, (tEnd - tStart) * 1000]);
end

function render_low_vs_high_fingerprint_local(lowVec, highVec, bladeIds, bladeColors, sid, windowId, lapId)
lowNorm = normalize_vector_max_local(lowVec);
highNorm = normalize_vector_max_local(highVec);
plot(bladeIds, lowNorm, 'k--o', 'LineWidth', 1.2, 'MarkerSize', 5);
hold on; box on; grid on;
plot(bladeIds, highNorm, '-', 'Color', [0.15 0.45 0.85], 'LineWidth', 1.4);
for bladeId = bladeIds
    plot(bladeId, highNorm(bladeId), 'o', 'MarkerSize', 6, ...
        'MarkerFaceColor', bladeColors(bladeId, :), ...
        'MarkerEdgeColor', bladeColors(bladeId, :));
    text(bladeId, highNorm(bladeId) + 0.06, sprintf('B%d', bladeId), ...
        'HorizontalAlignment', 'center', 'FontSize', 8, ...
        'Color', bladeColors(bladeId, :), 'FontWeight', 'bold');
end
corrValue = corr_scalar_local(lowNorm, highNorm);
xlabel('Physical blade ID');
ylabel('Normalized peak');
title(sprintf('CH%d W%02d Lap %d corr=%.3f', sid, windowId, lapId, corrValue));
legend({'Low-speed reference', 'High-speed current lap'}, 'Location', 'best');
ylim([min([lowNorm(:); highNorm(:); 0]) - 0.05, max([lowNorm(:); highNorm(:); 1]) + 0.18]);
end

function fp = read_low_speed_fingerprint_local(LowSpeedReference, sid, bladeCount)
if ~isKey(LowSpeedReference.fingerprints, sid)
    error('Low-speed fingerprint missing CH%d.', sid);
end
raw = LowSpeedReference.fingerprints(sid);
if numel(raw) < bladeCount
    error('Low-speed fingerprint for CH%d has only %d values.', sid, numel(raw));
end
fp = raw(1:bladeCount);
end

function highVec = build_high_speed_physical_peak_vector_local(rawSensor, jilublade, sensorNumbering, lapId, P)
rowsByPhysical = sensorNumbering.selected_rows_by_physical;
if lapId > size(rowsByPhysical, 1)
    highVec = nan(1, P.machine.bladeCount);
    return;
end
rows = rowsByPhysical(lapId, :);
highVec = nan(1, P.machine.bladeCount);
for bladeId = 1:P.machine.bladeCount
    rowId = rows(bladeId);
    if ~isfinite(rowId)
        continue;
    end
    [t0, t1] = build_dynamic_segment_window_local(jilublade, rowId, P.waveform.pulsePadSec);
    keep = rawSensor.T >= t0 & rawSensor.T <= t1;
    if nnz(keep) < max(8, P.numbering.matchPolyDegree + 2)
        continue;
    end
    fitInfo = fit_polynomial_peak_trace_local( ...
        rawSensor.T(keep), rawSensor.V(keep), P.numbering.matchPolyDegree);
    highVec(bladeId) = fitInfo.peakValue;
end
end

function [t0, t1] = build_dynamic_segment_window_local(jilublade, rowId, pulsePadSec)
t0 = jilublade(rowId, 1) - pulsePadSec;
t1 = jilublade(rowId, 2) + pulsePadSec;
end

function fitInfo = fit_polynomial_peak_trace_local(t, v, degree)
fitInfo = struct( ...
    'peakValue', NaN, ...
    'peakTime', NaN, ...
    'rawPeakValue', NaN, ...
    'rawPeakTime', NaN, ...
    't', [], ...
    'v', [], ...
    'tFine', [], ...
    'vFine', [], ...
    'fitSucceeded', false, ...
    'degreeUsed', NaN);
t = t(:);
v = v(:);
valid = isfinite(t) & isfinite(v);
t = t(valid);
v = v(valid);
if numel(t) < max(3, degree + 1) || range(t) <= eps
    return;
end
[t, order] = sort(t);
v = v(order);
[rawPeakValue, rawIdx] = max(v);
fitInfo.rawPeakValue = rawPeakValue;
fitInfo.rawPeakTime = t(rawIdx);
fitInfo.t = t;
fitInfo.v = v;
realDegree = min(degree, numel(t) - 1);
fitInfo.degreeUsed = realDegree;
[coef, ~, mu] = polyfit(t, v, realDegree);
tFine = linspace(min(t), max(t), 200).';
vFine = polyval(coef, tFine, [], mu);
[peakValue, peakIdx] = max(vFine);
fitInfo.peakValue = peakValue;
fitInfo.peakTime = tFine(peakIdx);
fitInfo.tFine = tFine;
fitInfo.vFine = vFine;
fitInfo.fitSucceeded = true;
end

function render_peak_fit_single_pulse_local( ...
        rawSensor, jilublade, sensorNumbering, lapId, bladeId, sid, windowId, P, bladeColor)
rowsByPhysical = sensorNumbering.selected_rows_by_physical;
if lapId > size(rowsByPhysical, 1)
    axis off;
    title(sprintf('CH%d B%d no lap', sid, bladeId));
    return;
end
rowId = rowsByPhysical(lapId, bladeId);
if ~isfinite(rowId)
    axis off;
    title(sprintf('CH%d B%d missing', sid, bladeId));
    return;
end
[t0, t1] = build_dynamic_segment_window_local(jilublade, rowId, P.waveform.pulsePadSec);
keep = rawSensor.T >= t0 & rawSensor.T <= t1;
if nnz(keep) < max(8, P.numbering.matchPolyDegree + 2)
    axis off;
    title(sprintf('CH%d B%d too few pts', sid, bladeId));
    return;
end
fitInfo = fit_polynomial_peak_trace_local(rawSensor.T(keep), rawSensor.V(keep), P.numbering.matchPolyDegree);
tMs = (fitInfo.t - t0) * 1000;
plot(tMs, fitInfo.v, 'k-', 'LineWidth', 0.9);
hold on; box on; grid on;
if fitInfo.fitSucceeded
    plot((fitInfo.tFine - t0) * 1000, fitInfo.vFine, '-', 'Color', [0.10 0.45 0.90], 'LineWidth', 1.1);
    plot((fitInfo.peakTime - t0) * 1000, fitInfo.peakValue, 'o', ...
        'MarkerSize', 5.5, 'MarkerFaceColor', bladeColor, 'MarkerEdgeColor', bladeColor);
end
plot((fitInfo.rawPeakTime - t0) * 1000, fitInfo.rawPeakValue, 'x', ...
    'MarkerSize', 5.5, 'LineWidth', 1.0, 'Color', [0.55 0.55 0.55]);
xline((jilublade(rowId, 1) - t0) * 1000, ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.7);
xline((jilublade(rowId, 2) - t0) * 1000, ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.7);
title(sprintf('CH%d B%d fit %.3f raw %.3f', sid, bladeId, fitInfo.peakValue, fitInfo.rawPeakValue), ...
    'FontSize', 8);
xlabel('Pulse window time (ms)');
ylabel('Voltage');
yl = ylim;
text(0.02 * max(tMs), yl(2) - 0.08 * range(yl), sprintf('W%02d L%d deg%d', windowId, lapId, fitInfo.degreeUsed), ...
    'FontSize', 7, 'Color', [0.25 0.25 0.25], 'VerticalAlignment', 'top');
end

function y = normalize_vector_max_local(x)
x = x(:);
mx = max(x);
if ~isfinite(mx) || mx <= eps
    y = nan(size(x));
else
    y = x / mx;
end
end

function value = corr_scalar_local(a, b)
if numel(a) < 3 || numel(b) < 3 || any(~isfinite(a)) || any(~isfinite(b)) || ...
        std(a) <= eps || std(b) <= eps
    value = NaN;
    return;
end
C = corrcoef(a(:), b(:));
value = C(1, 2);
end

function y = normalize_trace_local(v)
v = v(:);
if isempty(v) || ~any(isfinite(v))
    y = zeros(size(v));
    return;
end
vMin = min(v, [], 'omitnan');
vMax = max(v, [], 'omitnan');
if ~isfinite(vMin) || ~isfinite(vMax) || abs(vMax - vMin) <= eps
    y = zeros(size(v));
else
    y = (v - vMin) / (vMax - vMin);
end
end

function vPeak = interpolate_peak_local(tPlot, vPlot, tPeak)
if isempty(tPlot) || isempty(vPlot)
    vPeak = NaN;
    return;
end
[~, idx] = min(abs(tPlot - tPeak));
vPeak = vPlot(idx);
end

function fileRanges = build_case_file_ranges_local(caseDir, oprChannel, sampleRateHz)
files = dir(fullfile(caseDir, sprintf('Probe%d_Data_*.dat', oprChannel)));
if isempty(files)
    files = dir(fullfile(caseDir, sprintf('Probe%d_Data_*.txt', oprChannel)));
end
if isempty(files)
    files = dir(fullfile(caseDir, sprintf('4-%d-*.mat', oprChannel)));
end
if isempty(files)
    error('No OPR raw files found in %s.', caseDir);
end
fileIds = nan(numel(files), 1);
for i = 1:numel(files)
    token = regexp(files(i).name, 'Data_(\d+)', 'tokens', 'once');
    if isempty(token)
        token = regexp(files(i).name, '4-\d+-(\d+)\.mat', 'tokens', 'once');
    end
    fileIds(i) = str2double(token{1});
end
[fileIds, order] = sort(fileIds);
files = files(order);
fileRanges = repmat(struct('file_id', NaN, 'offset', NaN, 't_start', NaN, 't_end', NaN), numel(files), 1);
lastEnd = 0;
for i = 1:numel(files)
    [tOpr, ~] = load_raw_channel_local(caseDir, oprChannel, fileIds(i), sampleRateHz);
    if i == 1
        offset = 0;
    else
        offset = lastEnd + 1 / sampleRateHz - tOpr(1);
    end
    fileRanges(i).file_id = fileIds(i);
    fileRanges(i).offset = offset;
    fileRanges(i).t_start = tOpr(1) + offset;
    fileRanges(i).t_end = tOpr(end) + offset;
    lastEnd = fileRanges(i).t_end;
end
end

function raw = load_raw_subset_local(caseDir, fileRanges, sensorIds, sampleRateHz, timeWindow)
raw(max(sensorIds)) = struct('T', [], 'V', []);
useFiles = fileRanges([fileRanges.t_end] >= timeWindow(1) & [fileRanges.t_start] <= timeWindow(2));
for k = 1:numel(useFiles)
    fileId = useFiles(k).file_id;
    offset = useFiles(k).offset;
    for sid = sensorIds
        [tLocal, vLocal] = load_raw_channel_local(caseDir, sid, fileId, sampleRateHz);
        tGlobal = tLocal(:) + offset;
        keep = tGlobal >= timeWindow(1) & tGlobal <= timeWindow(2);
        raw(sid).T = [raw(sid).T; tGlobal(keep)]; %#ok<AGROW>
        raw(sid).V = [raw(sid).V; vLocal(keep)]; %#ok<AGROW>
    end
end
end

function [tSec, v] = load_raw_channel_local(caseDir, channelId, fileId, sampleRateHz)
patterns = { ...
    sprintf('Probe%d_Data_%d.dat', channelId, fileId), ...
    sprintf('Probe%d_Data_%d.txt', channelId, fileId), ...
    sprintf('4-%d-%d.mat', channelId, fileId)};
filePath = '';
for i = 1:numel(patterns)
    candidate = fullfile(caseDir, patterns{i});
    if exist(candidate, 'file') == 2
        filePath = candidate;
        break;
    end
end
if isempty(filePath)
    error('Missing raw file for CH%d file %d in %s.', channelId, fileId, caseDir);
end
varName = sprintf('jilu%02d', channelId);
vars = whos('-file', filePath);
if ~ismember(varName, {vars.name})
    error('Raw MAT file does not contain %s: %s', varName, filePath);
end
loaded = load(filePath, varName);
raw = loaded.(varName);
raw(raw(:, 1) == 0, :) = [];
tSec = raw(:, 1) / sampleRateHz;
v = raw(:, 2);
end

