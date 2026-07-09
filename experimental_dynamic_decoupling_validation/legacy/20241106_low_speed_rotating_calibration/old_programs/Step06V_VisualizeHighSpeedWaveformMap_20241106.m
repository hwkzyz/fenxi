%% Step06V_VisualizeHighSpeedWaveformMap_20241106
% Diagnostic visualization for Step06 waveform extraction and filtering.
% Focus:
%   1) whether raw high-speed waveforms land near the low-speed template domain,
%   2) whether filtered points cover the useful template slope region,
%   3) which sensors are effectively contributing to Step07.
%
% This visualizer intentionally uses raw voltage overlays as the main
% judgment basis. Shape-only normalized comparison is not shown here
% because it can hide amplitude and baseline mismatches.

clear; close all; clc;

P = NewFlow_Config_20241106();
targetBladeSensorPairs = [5 7; 6 5];
systematicSensorId = 7;
maxTargetLapsToPlot = 20;

if exist(P.files.waveformLibrary, 'file') ~= 2
    error('Missing Step06 artifact:\n  %s', P.files.waveformLibrary);
end
if exist(P.files.filteredWaveformLibrary, 'file') ~= 2
    error(['Missing Step06 filtered artifact:\n  %s\n\n' ...
        'Run Step06_BuildHighSpeedWaveformMap_20241106 first.'], ...
        P.files.filteredWaveformLibrary);
end

loaded = load(P.files.waveformLibrary, 'WaveformLibrary');
WaveformLibrary = loaded.WaveformLibrary;
loadedFiltered = load(P.files.filteredWaveformLibrary, 'FilteredWaveformLibrary');
FilteredWaveformLibrary = loadedFiltered.FilteredWaveformLibrary;

figDir = fullfile(P.view.figureDir, '05_waveform_library');
if P.view.saveFigures && exist(figDir, 'dir') ~= 7
    mkdir(figDir);
end

bladeIds = [WaveformLibrary.Blade.blade_id];
sensorIds = [WaveformLibrary.Blade(1).Sensor.sensor_id];
nBlade = numel(bladeIds);
nSensor = numel(sensorIds);

rawMedianMin = nan(nSensor, nBlade);
rawMedianMax = nan(nSensor, nBlade);
rawMedianCenter = nan(nSensor, nBlade);
templateDomainLeft = nan(nSensor, nBlade);
templateDomainRight = nan(nSensor, nBlade);
templateCenter = nan(nSensor, nBlade);
rawOverlapRatio = nan(nSensor, nBlade);
filteredPointMedian = nan(nSensor, nBlade);
filteredCoverageRatio = nan(nSensor, nBlade);
windowPointRatio = nan(nSensor, nBlade);

for ib = 1:nBlade
    for is = 1:nSensor
        S = WaveformLibrary.Blade(ib).Sensor(is);
        xmins = arrayfun(@(x) min_or_nan_local(x.x_rel), S.Lap);
        xmaxs = arrayfun(@(x) max_or_nan_local(x.x_rel), S.Lap);
        rawMedianMin(is, ib) = median(xmins, 'omitnan');
        rawMedianMax(is, ib) = median(xmaxs, 'omitnan');
        rawMedianCenter(is, ib) = median(0.5 * (xmins + xmaxs), 'omitnan');
        templateDomainLeft(is, ib) = S.template_domain(1);
        templateDomainRight(is, ib) = S.template_domain(2);
        templateCenter(is, ib) = S.template_xc;
        rawOverlapRatio(is, ib) = overlap_ratio_local( ...
            [rawMedianMin(is, ib), rawMedianMax(is, ib)], ...
            [templateDomainLeft(is, ib), templateDomainRight(is, ib)]);

        if numel(FilteredWaveformLibrary.Blade(ib).Window) >= 1
            Bdl = FilteredWaveformLibrary.Blade(ib).Window(1).Bundle;
            meta = FilteredWaveformLibrary.Blade(ib).Window(1).SensorMeta(is);
            keep = Bdl.S == sensorIds(is);
            filteredPointMedian(is, ib) = nnz(keep);
            filteredCoverageRatio(is, ib) = overlap_ratio_local( ...
                [min_or_nan_local(Bdl.X(keep)), max_or_nan_local(Bdl.X(keep))], ...
                [templateDomainLeft(is, ib), templateDomainRight(is, ib)]);
            if meta.raw_points > 0
                windowPointRatio(is, ib) = meta.valid_points / meta.raw_points;
            end
        end
    end
end

figSummary = figure('Name', 'Step06V waveform extraction summary', 'Color', 'w', ...
    'Position', [80, 80, 1450, 860], 'NumberTitle', 'off');
tiledlayout(figSummary, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
imagesc(bladeIds, sensorIds, rawMedianCenter);
set(gca, 'YDir', 'normal');
xlabel('Physical blade ID');
ylabel('Sensor ID');
title('Raw high-speed median x-window center (mm)');
colorbar; box on;

nexttile;
imagesc(bladeIds, sensorIds, rawOverlapRatio);
set(gca, 'YDir', 'normal');
xlabel('Physical blade ID');
ylabel('Sensor ID');
title('Raw window / template-domain overlap ratio');
colorbar; box on;

nexttile;
imagesc(bladeIds, sensorIds, windowPointRatio);
set(gca, 'YDir', 'normal');
xlabel('Physical blade ID');
ylabel('Sensor ID');
title('Filtered-point / raw-point ratio');
colorbar; box on;

nexttile;
imagesc(bladeIds, sensorIds, filteredPointMedian);
set(gca, 'YDir', 'normal');
xlabel('Physical blade ID');
ylabel('Sensor ID');
title('Filtered points in window 1');
colorbar; box on;

if P.view.saveFigures
    exportgraphics(figSummary, fullfile(figDir, 'Step06V_WaveformExtractionSummary_20241106.png'), 'Resolution', 300);
end

figDiag = figure('Name', 'Step06V raw waveform and filtered-point diagnosis', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [1, 1, 30, 22], 'NumberTitle', 'off');
tiledlayout(figDiag, nSensor, nBlade, 'TileSpacing', 'compact', 'Padding', 'compact');
maxLapsToPlot = 3;
for is = 1:nSensor
    for ib = 1:nBlade
        nexttile;
        S = WaveformLibrary.Blade(ib).Sensor(is);
        plot(S.template_x, S.template_v, 'k-', 'LineWidth', 1.3, 'DisplayName', 'low-speed template'); hold on;
        cmap = lines(maxLapsToPlot);
        for k = 1:min(maxLapsToPlot, numel(S.Lap))
            [xRel, order] = sort(S.Lap(k).x_rel(:));
            v = S.Lap(k).V(order);
            plot(xRel, v, '-', 'Color', cmap(k, :), 'LineWidth', 0.9, ...
                'DisplayName', sprintf('lap %d', k));
        end
        xline(S.template_domain(1), ':', 'Color', [0.35 0.35 0.35], 'HandleVisibility', 'off');
        xline(S.template_domain(2), ':', 'Color', [0.35 0.35 0.35], 'HandleVisibility', 'off');
        xline(0, '--', 'Color', [0.65 0.1 0.1], 'LineWidth', 0.8, 'HandleVisibility', 'off');
        if numel(FilteredWaveformLibrary.Blade(ib).Window) >= 1
            Bdl = FilteredWaveformLibrary.Blade(ib).Window(1).Bundle;
            keep = Bdl.S == sensorIds(is);
            if any(keep)
                scatter(Bdl.X(keep), Bdl.V(keep), 5, Bdl.T_rel(keep), 'filled', ...
                    'DisplayName', 'filtered W1');
            end
        end
        title(sprintf('B%d CH%d', bladeIds(ib), sensorIds(is)), 'FontWeight', 'normal');
        if ib == 1
            ylabel('Voltage (V)');
        end
        if is == nSensor
            xlabel('x relative to template center (mm)');
        end
        if is == 1 && ib == 1
            legend('Location', 'best', 'Box', 'off');
        end
        grid on; box on;
    end
end
if P.view.saveFigures
    exportgraphics(figDiag, fullfile(figDir, 'Step06V_RawTemplateFilteredDiagnosis_20241106.png'), 'Resolution', 300);
end

diagnosisRows = repmat(struct( ...
    'BladeID', NaN, ...
    'SensorID', NaN, ...
    'RawMedianXMinMM', NaN, ...
    'RawMedianXMaxMM', NaN, ...
    'RawMedianCenterMM', NaN, ...
    'TemplateDomainLeftMM', NaN, ...
    'TemplateDomainRightMM', NaN, ...
    'TemplateXCMM', NaN, ...
    'RawTemplateOverlapRatio', NaN, ...
    'FilteredWindow1PointCount', NaN, ...
    'FilteredRawRatio', NaN), nBlade * nSensor, 1);

rowIdx = 0;
for ib = 1:nBlade
    for is = 1:nSensor
        rowIdx = rowIdx + 1;
        diagnosisRows(rowIdx).BladeID = bladeIds(ib);
        diagnosisRows(rowIdx).SensorID = sensorIds(is);
        diagnosisRows(rowIdx).RawMedianXMinMM = rawMedianMin(is, ib);
        diagnosisRows(rowIdx).RawMedianXMaxMM = rawMedianMax(is, ib);
        diagnosisRows(rowIdx).RawMedianCenterMM = rawMedianCenter(is, ib);
        diagnosisRows(rowIdx).TemplateDomainLeftMM = templateDomainLeft(is, ib);
        diagnosisRows(rowIdx).TemplateDomainRightMM = templateDomainRight(is, ib);
        diagnosisRows(rowIdx).TemplateXCMM = templateCenter(is, ib);
        diagnosisRows(rowIdx).RawTemplateOverlapRatio = rawOverlapRatio(is, ib);
        diagnosisRows(rowIdx).FilteredWindow1PointCount = filteredPointMedian(is, ib);
        diagnosisRows(rowIdx).FilteredRawRatio = windowPointRatio(is, ib);
    end
end
DiagnosisTable = struct2table(diagnosisRows);
writetable(DiagnosisTable, fullfile(figDir, 'Step06V_WaveformDiagnosis_20241106.csv'));

[TargetLapTable, TargetSummaryTable] = build_target_pair_offset_tables_local( ...
    WaveformLibrary, targetBladeSensorPairs);
writetable(TargetLapTable, fullfile(figDir, 'Step06V_TargetPairLapOffsets_20241106.csv'));
writetable(TargetSummaryTable, fullfile(figDir, 'Step06V_TargetPairSummary_20241106.csv'));
plot_target_pair_offset_diagnosis_local( ...
    WaveformLibrary, FilteredWaveformLibrary, targetBladeSensorPairs, ...
    maxTargetLapsToPlot, figDir, P.view.saveFigures);

[SystematicOffsetTable, SystematicSensorSummary] = build_systematic_offset_tables_local(WaveformLibrary);
writetable(SystematicOffsetTable, fullfile(figDir, 'Step06V_SystematicOffsetByBladeSensor_20241106.csv'));
writetable(SystematicSensorSummary, fullfile(figDir, 'Step06V_SystematicOffsetBySensor_20241106.csv'));
plot_systematic_sensor_offset_local( ...
    SystematicOffsetTable, SystematicSensorSummary, systematicSensorId, ...
    figDir, P.view.saveFigures);

fprintf('\n=== Step06V: visualize waveform library ===\n');
fprintf('Source: %s\n', P.files.waveformLibrary);
fprintf('Filtered source: %s\n', P.files.filteredWaveformLibrary);
fprintf('Diagnosis CSV: %s\n', fullfile(figDir, 'Step06V_WaveformDiagnosis_20241106.csv'));
fprintf('Target pair lap offsets: %s\n', fullfile(figDir, 'Step06V_TargetPairLapOffsets_20241106.csv'));
fprintf('Systematic sensor offsets: %s\n', fullfile(figDir, 'Step06V_SystematicOffsetBySensor_20241106.csv'));
fprintf('Figures: %s\n', figDir);

function [LapTable, SummaryTable] = build_target_pair_offset_tables_local(WaveformLibrary, targetPairs)
rows = [];
summaryRows = [];
for i = 1:size(targetPairs, 1)
    bladeId = targetPairs(i, 1);
    sensorId = targetPairs(i, 2);
    S = get_blade_sensor_local(WaveformLibrary, bladeId, sensorId);
    lapCenters = nan(numel(S.Lap), 1);
    for k = 1:numel(S.Lap)
        x = S.Lap(k).x_rel(:);
        v = S.Lap(k).V(:);
        xmin = min_or_nan_local(x);
        xmax = max_or_nan_local(x);
        center = 0.5 * (xmin + xmax);
        lapCenters(k) = center;
        rows = [rows; bladeId, sensorId, k, S.Lap(k).row_id, ...
            xmin, xmax, center, max_or_nan_local(v), min_or_nan_local(v), numel(x)]; %#ok<AGROW>
    end
    summaryRows = [summaryRows; bladeId, sensorId, numel(S.Lap), ...
        median(lapCenters, 'omitnan'), mean(lapCenters, 'omitnan'), ...
        std(lapCenters, 'omitnan'), min_or_nan_local(lapCenters), max_or_nan_local(lapCenters), ...
        S.template_domain(1), S.template_domain(2)]; %#ok<AGROW>
end
LapTable = array2table(rows, 'VariableNames', { ...
    'BladeID', 'SensorID', 'LapID', 'RowID', 'RawXMinMM', 'RawXMaxMM', ...
    'RawWindowCenterRelativeToTemplateMM', 'RawVMax', 'RawVMin', 'RawPointCount'});
SummaryTable = array2table(summaryRows, 'VariableNames', { ...
    'BladeID', 'SensorID', 'LapCount', 'MedianCenterRelativeToTemplateMM', ...
    'MeanCenterRelativeToTemplateMM', 'StdCenterRelativeToTemplateMM', ...
    'MinCenterRelativeToTemplateMM', 'MaxCenterRelativeToTemplateMM', ...
    'TemplateDomainLeftMM', 'TemplateDomainRightMM'});
end

function [OffsetTable, SensorSummary] = build_systematic_offset_tables_local(WaveformLibrary)
bladeIds = [WaveformLibrary.Blade.blade_id];
sensorIds = [WaveformLibrary.Blade(1).Sensor.sensor_id];
rows = [];
for ib = 1:numel(bladeIds)
    for is = 1:numel(sensorIds)
        S = WaveformLibrary.Blade(ib).Sensor(is);
        centers = nan(numel(S.Lap), 1);
        for k = 1:numel(S.Lap)
            x = S.Lap(k).x_rel(:);
            centers(k) = 0.5 * (min_or_nan_local(x) + max_or_nan_local(x));
        end
        rows = [rows; bladeIds(ib), sensorIds(is), numel(S.Lap), ...
            median(centers, 'omitnan'), mean(centers, 'omitnan'), ...
            std(centers, 'omitnan'), min_or_nan_local(centers), max_or_nan_local(centers)]; %#ok<AGROW>
    end
end
OffsetTable = array2table(rows, 'VariableNames', { ...
    'BladeID', 'SensorID', 'LapCount', 'MedianCenterRelativeToTemplateMM', ...
    'MeanCenterRelativeToTemplateMM', 'StdCenterRelativeToTemplateMM', ...
    'MinCenterRelativeToTemplateMM', 'MaxCenterRelativeToTemplateMM'});

sensorIds = unique(OffsetTable.SensorID).';
summaryRows = [];
for sensorId = sensorIds
    keep = OffsetTable.SensorID == sensorId;
    v = OffsetTable.MedianCenterRelativeToTemplateMM(keep);
    summaryRows = [summaryRows; sensorId, nnz(keep), ...
        mean(v, 'omitnan'), median(v, 'omitnan'), std(v, 'omitnan'), ...
        mean(abs(v), 'omitnan'), max(abs(v), [], 'omitnan')]; %#ok<AGROW>
end
SensorSummary = array2table(summaryRows, 'VariableNames', { ...
    'SensorID', 'BladeCount', 'MeanBladeMedianCenterMM', 'MedianBladeMedianCenterMM', ...
    'StdBladeMedianCenterMM', 'MeanAbsBladeMedianCenterMM', 'MaxAbsBladeMedianCenterMM'});
end

function plot_target_pair_offset_diagnosis_local(WaveformLibrary, FilteredWaveformLibrary, ...
    targetPairs, maxLapsToPlot, figDir, saveFigures)
if isempty(targetPairs)
    return;
end
fig = figure('Name', 'Step06V target high-speed offset diagnosis', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [1, 1, 30, 9 * size(targetPairs, 1)], ...
    'NumberTitle', 'off');
tiledlayout(fig, size(targetPairs, 1), 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:size(targetPairs, 1)
    bladeId = targetPairs(i, 1);
    sensorId = targetPairs(i, 2);
    S = get_blade_sensor_local(WaveformLibrary, bladeId, sensorId);
    nexttile;
    plot(S.template_x, S.template_v, 'k-', 'LineWidth', 1.5, 'DisplayName', 'low-speed template'); hold on;
    nPlot = min(maxLapsToPlot, numel(S.Lap));
    for k = 1:nPlot
        [xRel, order] = sort(S.Lap(k).x_rel(:));
        plot(xRel, S.Lap(k).V(order), '-', 'Color', [0.70 0.78 0.86], ...
            'LineWidth', 0.8, 'HandleVisibility', 'off');
    end
    Bdl = get_filtered_bundle_local(FilteredWaveformLibrary, bladeId);
    if ~isempty(Bdl)
        keep = Bdl.S == sensorId;
        scatter(Bdl.X(keep), Bdl.V(keep), 8, [0.85 0.20 0.10], 'filled', ...
            'DisplayName', 'filtered W1');
    end
    xline(S.template_domain(1), ':', 'Color', [0.30 0.30 0.30], 'HandleVisibility', 'off');
    xline(S.template_domain(2), ':', 'Color', [0.30 0.30 0.30], 'HandleVisibility', 'off');
    xline(0, '--', 'Color', [0.65 0.10 0.10], 'LineWidth', 0.9, 'HandleVisibility', 'off');
    grid on; box on;
    xlabel('x relative to low-speed template center (mm)');
    ylabel('Raw voltage (V)');
    title(sprintf('B%d CH%d raw waveform overlay', bladeId, sensorId), 'FontWeight', 'normal');
    legend('Location', 'best', 'Box', 'off');

    nexttile;
    centers = nan(numel(S.Lap), 1);
    for k = 1:numel(S.Lap)
        x = S.Lap(k).x_rel(:);
        centers(k) = 0.5 * (min_or_nan_local(x) + max_or_nan_local(x));
    end
    plot(1:numel(centers), centers, 'o-', 'Color', [0.10 0.35 0.70], ...
        'MarkerFaceColor', [0.10 0.35 0.70], 'LineWidth', 1.2);
    yline(0, '--', 'Color', [0.65 0.10 0.10], 'LineWidth', 0.9);
    yline(median(centers, 'omitnan'), '-', 'Color', [0.10 0.10 0.10], 'LineWidth', 1.0);
    grid on; box on;
    xlabel('Lap ID');
    ylabel('Window center offset (mm)');
    title(sprintf('median = %.3f mm, std = %.3f mm', ...
        median(centers, 'omitnan'), std(centers, 'omitnan')), 'FontWeight', 'normal');
end
if saveFigures
    exportgraphics(fig, fullfile(figDir, 'Step06V_TargetPairHighSpeedOffset_20241106.png'), 'Resolution', 300);
end
end

function plot_systematic_sensor_offset_local(OffsetTable, SensorSummary, targetSensorId, figDir, saveFigures)
fig = figure('Name', 'Step06V systematic sensor offset diagnosis', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [1, 1, 28, 14], 'NumberTitle', 'off');
tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
sensorIds = unique(OffsetTable.SensorID).';
hold on;
for i = 1:numel(sensorIds)
    sensorId = sensorIds(i);
    keep = OffsetTable.SensorID == sensorId;
    x = i + linspace(-0.12, 0.12, nnz(keep)).';
    if sensorId == targetSensorId
        c = [0.85 0.20 0.10];
    else
        c = [0.20 0.45 0.70];
    end
    scatter(x, OffsetTable.MedianCenterRelativeToTemplateMM(keep), 34, c, 'filled');
end
yline(0, '--', 'Color', [0.35 0.35 0.35]);
set(gca, 'XTick', 1:numel(sensorIds), 'XTickLabel', compose('CH%d', sensorIds));
xlabel('Sensor');
ylabel('Blade median window-center offset (mm)');
title('Per-blade high-speed offset relative to low-speed template', 'FontWeight', 'normal');
grid on; box on;

nexttile;
[~, order] = sort(SensorSummary.SensorID);
barValues = SensorSummary.MeanAbsBladeMedianCenterMM(order);
barColors = repmat([0.25 0.45 0.65], numel(order), 1);
targetIdx = find(SensorSummary.SensorID(order) == targetSensorId, 1);
if ~isempty(targetIdx)
    barColors(targetIdx, :) = [0.85 0.20 0.10];
end
b = bar(barValues, 'FaceColor', 'flat');
b.CData = barColors;
set(gca, 'XTick', 1:numel(order), 'XTickLabel', compose('CH%d', SensorSummary.SensorID(order)));
xlabel('Sensor');
ylabel('Mean abs blade median offset (mm)');
title('Sensor-level systematic offset score', 'FontWeight', 'normal');
grid on; box on;

if saveFigures
    exportgraphics(fig, fullfile(figDir, 'Step06V_SystematicSensorOffset_20241106.png'), 'Resolution', 300);
end
end

function S = get_blade_sensor_local(WaveformLibrary, bladeId, sensorId)
bladeIds = [WaveformLibrary.Blade.blade_id];
ib = find(bladeIds == bladeId, 1);
if isempty(ib)
    error('WaveformLibrary does not contain B%d.', bladeId);
end
sensorIds = [WaveformLibrary.Blade(ib).Sensor.sensor_id];
is = find(sensorIds == sensorId, 1);
if isempty(is)
    error('WaveformLibrary B%d does not contain CH%d.', bladeId, sensorId);
end
S = WaveformLibrary.Blade(ib).Sensor(is);
end

function Bdl = get_filtered_bundle_local(FilteredWaveformLibrary, bladeId)
Bdl = [];
bladeIds = [FilteredWaveformLibrary.Blade.blade_id];
ib = find(bladeIds == bladeId, 1);
if isempty(ib) || isempty(FilteredWaveformLibrary.Blade(ib).Window)
    return;
end
Bdl = FilteredWaveformLibrary.Blade(ib).Window(1).Bundle;
end

function ratio = overlap_ratio_local(a, b)
if any(~isfinite([a(:); b(:)])) || a(2) <= a(1) || b(2) <= b(1)
    ratio = NaN;
    return;
end
left = max(a(1), b(1));
right = min(a(2), b(2));
overlap = max(0, right - left);
ratio = overlap / max(b(2) - b(1), eps);
end

function value = min_or_nan_local(x)
if isempty(x)
    value = NaN;
else
    value = min(x);
end
end

function value = max_or_nan_local(x)
if isempty(x)
    value = NaN;
else
    value = max(x);
end
end
