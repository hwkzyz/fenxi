%% Step03_SelectHighSpeedRegion_20241106
% Select one reusable high-speed region. This step only writes region
% metadata; it does not number blades or extract waveforms.

clear; close all; clc;

%% Parameters to tune
analysisSensors = [2 3 5 7];

regionMode = 'manual_start';     % 'manual_start' | 'explicit_range' | 'region_plan'
startTimeSec = 75.0;
timeRangeSec = [74.998 75.369];
targetLaps = 20;
regionId = 2;
startMode = 'peak';              % only for region_plan: 'peak' | 'start'
regionPlanFile = fullfile(fileparts(mfilename('fullpath')), 'outputs', 'Step05_BTT_WindowPlan_20241106.csv');

viewEnable = true;
saveFigures = true;

%% Paths
routeDir = fileparts(mfilename('fullpath'));
timeLabel = time_label_local(startTimeSec);
outDir = fullfile(routeDir, 'output', 'new_flow', '03_region_selection');
outFile = fullfile(outDir, sprintf('HighSpeedRegion_%s_20241106.mat', timeLabel));
figureDir = fullfile(routeDir, 'output', 'new_flow', 'figures', '03_region_selection');

if exist(outDir, 'dir') ~= 7
    mkdir(outDir);
end

RegionSelection = build_region_selection_local( ...
    analysisSensors, regionMode, startTimeSec, timeRangeSec, targetLaps, ...
    regionId, startMode, regionPlanFile);

save(outFile, 'RegionSelection');
writetable(struct2table(RegionSelection), strrep(outFile, '.mat', '.csv'));

fprintf('\n=== Step03: select high-speed region ===\n');
fprintf('Mode: %s\n', RegionSelection.selection_source);
fprintf('Sensors: %s\n', mat2str(RegionSelection.analysis_sensors));
fprintf('Start time: %.6f s\n', RegionSelection.start_time_sec);
fprintf('Region range: [%.6f, %.6f] s\n', ...
    RegionSelection.region_start_sec, RegionSelection.end_time_sec);
fprintf('Saved: %s\n', outFile);

if viewEnable
    visualize_region_selection_local(RegionSelection, targetLaps, saveFigures, figureDir);
end

function S = build_region_selection_local( ...
        analysisSensors, regionMode, startTimeSec, timeRangeSec, targetLaps, ...
        regionId, startMode, regionPlanFile)
modeText = lower(strtrim(regionMode));
S = struct( ...
    'selection_source', modeText, ...
    'region_id', NaN, ...
    'start_mode', 'manual', ...
    'analysis_sensors', analysisSensors, ...
    'target_laps', targetLaps, ...
    'start_time_sec', startTimeSec, ...
    'region_start_sec', startTimeSec, ...
    'end_time_sec', NaN, ...
    'peak_time_sec', NaN, ...
    'dominant_order', NaN, ...
    'dominant_freq_hz', NaN, ...
    'region_plan_file', regionPlanFile);

switch modeText
    case 'manual_start'
        S.start_mode = 'manual';
        S.start_time_sec = startTimeSec;
        S.region_start_sec = startTimeSec;

    case 'explicit_range'
        if numel(timeRangeSec) ~= 2 || any(~isfinite(timeRangeSec))
            error('timeRangeSec must be a finite [t0, t1] pair.');
        end
        S.start_mode = 'start';
        S.start_time_sec = timeRangeSec(1);
        S.region_start_sec = timeRangeSec(1);
        S.end_time_sec = timeRangeSec(2);

    case 'region_plan'
        if exist(regionPlanFile, 'file') ~= 2
            error('Region plan file not found: %s', regionPlanFile);
        end
        T = readtable(regionPlanFile);
        row = T(T.regionId == regionId, :);
        if height(row) ~= 1
            error('Region plan must contain exactly one row for regionId=%d.', regionId);
        end

        S.region_id = regionId;
        S.start_mode = lower(strtrim(startMode));
        S.region_start_sec = row.bttStartSec(1);
        S.end_time_sec = row.bttEndSec(1);
        S.peak_time_sec = row.peakBttTimeSec(1);
        S.dominant_order = row.dominantOrder(1);
        S.dominant_freq_hz = row.dominantFreqHz(1);

        if strcmpi(S.start_mode, 'peak')
            S.start_time_sec = S.peak_time_sec;
        elseif strcmpi(S.start_mode, 'start')
            S.start_time_sec = S.region_start_sec;
        else
            error('startMode must be "peak" or "start".');
        end

    otherwise
        error('Unsupported regionMode: %s', regionMode);
end

if ~isfinite(S.end_time_sec)
    S.end_time_sec = S.region_start_sec + targetLaps * 0.02;
end
end

function visualize_region_selection_local(RegionSelection, targetLaps, saveFigures, figureDir)
if saveFigures && exist(figureDir, 'dir') ~= 7
    mkdir(figureDir);
end

startSec = RegionSelection.start_time_sec;
regionStart = RegionSelection.region_start_sec;
regionEnd = RegionSelection.end_time_sec;
if ~isfinite(regionEnd)
    regionEnd = regionStart + targetLaps * 0.02;
end

fig = figure('Name', 'Step03 high-speed region selection', 'Color', 'w');
hold on;
plot([regionStart regionEnd], [1 1], 'LineWidth', 8, 'Color', [0.2 0.55 0.85]);
xline(startSec, '--k', 'Start', 'LabelOrientation', 'horizontal');
if isfinite(RegionSelection.peak_time_sec)
    xline(RegionSelection.peak_time_sec, ':r', 'Peak', 'LabelOrientation', 'horizontal');
end
ylim([0.5 1.5]);
yticks(1);
yticklabels({'selected region'});
xlabel('Time (s)');
title(sprintf('High-speed region: %s, sensors %s', ...
    RegionSelection.selection_source, mat2str(RegionSelection.analysis_sensors)), ...
    'Interpreter', 'none');
grid on; box on;

if saveFigures
    exportgraphics(fig, fullfile(figureDir, 'Step03_HighSpeedRegionSelection_20241106.png'), ...
        'Resolution', 300);
end
end

function label = time_label_local(tSec)
if ~isfinite(tSec)
    label = 'TUnknown';
    return;
end
tMs = round(1000 * tSec);
if mod(tMs, 1000) == 0
    label = sprintf('T%03ds', tMs / 1000);
else
    label = sprintf('T%03dp%03ds', floor(tMs / 1000), mod(tMs, 1000));
end
end
