%% Step02A_SelectHighSpeedRegion_20241106
% Build one reusable high-speed region-selection file for the current run.
% Step02 can load this file directly, so switching time regions does not
% require editing the numbering logic itself.

clear; close all; clc;

routeDir = fileparts(mfilename('fullpath'));

%% Parameters
P.mode = 'manual_start';   % 'manual_start' | 'explicit_range' | 'region_plan'
P.analysisSensors = [2 3 5 7];
P.targetBlade = 6;
P.targetLaps = 20;
P.analysisStartTime = 75.0;
P.regionTimeRange = [74.998, 75.369];
P.regionId = 2;
P.regionStartMode = 'peak';  % only used for region_plan: 'peak' or 'start'
P.regionPlanFile = fullfile(routeDir, 'outputs', 'Step05_BTT_WindowPlan_20241106.csv');
P.regionTag = '';
P.saveCsv = true;

if isempty(P.regionTag)
    switch lower(strtrim(P.mode))
        case 'explicit_range'
            P.regionTag = sprintf('%s_to_%s', ...
                make_time_label_local(P.regionTimeRange(1)), ...
                make_time_label_local(P.regionTimeRange(2)));
        otherwise
            P.regionTag = make_time_label_local(P.analysisStartTime);
    end
end

outputDir = fullfile(routeDir, 'output', 'region_selection');
if exist(outputDir, 'dir') ~= 7
    mkdir(outputDir);
end
matFile = fullfile(outputDir, sprintf('HighSpeedRegionSelection_%s_20241106.mat', P.regionTag));
csvFile = fullfile(outputDir, sprintf('HighSpeedRegionSelection_%s_20241106.csv', P.regionTag));

RegionSelection = build_region_selection_local(P);
RegionSelection.region_selection_file = matFile;
save(matFile, 'RegionSelection');

if P.saveCsv
    T = struct2table(RegionSelection);
    writetable(T, csvFile);
end

fprintf('\n=== Step02A: high-speed region selection ===\n');
fprintf('Mode: %s\n', RegionSelection.selection_source);
fprintf('Sensors: %s\n', mat2str(RegionSelection.analysis_sensors));
fprintf('Target blade: %d\n', RegionSelection.target_blade);
fprintf('Target laps: %d\n', RegionSelection.target_laps);
fprintf('Start time: %.6f s\n', RegionSelection.start_time_sec);
fprintf('Region range: [%.6f, %.6f] s\n', RegionSelection.region_start_sec, RegionSelection.end_time_sec);
fprintf('Saved MAT: %s\n', matFile);
if P.saveCsv
    fprintf('Saved CSV: %s\n', csvFile);
end

function RegionSelection = build_region_selection_local(P)
modeText = lower(strtrim(P.mode));
RegionSelection = struct( ...
    'use_saved_selection', true, ...
    'use_region_plan', false, ...
    'selection_source', modeText, ...
    'region_id', NaN, ...
    'start_mode', 'manual', ...
    'target_blade', P.targetBlade, ...
    'region_selected_blade', NaN, ...
    'start_time_sec', P.analysisStartTime, ...
    'region_start_sec', P.analysisStartTime, ...
    'end_time_sec', NaN, ...
    'peak_time_sec', NaN, ...
    'dominant_order', NaN, ...
    'dominant_freq_hz', NaN, ...
    'target_laps', P.targetLaps, ...
    'analysis_sensors', P.analysisSensors, ...
    'region_selection_file', '', ...
    'region_plan_file', P.regionPlanFile);

switch modeText
    case 'manual_start'
        RegionSelection.start_mode = 'manual';
        RegionSelection.start_time_sec = P.analysisStartTime;
        RegionSelection.region_start_sec = P.analysisStartTime;

    case 'explicit_range'
        if numel(P.regionTimeRange) ~= 2 || any(~isfinite(P.regionTimeRange))
            error('P.regionTimeRange must be a finite [t0, t1] pair.');
        end
        RegionSelection.selection_source = 'explicit_range';
        RegionSelection.start_mode = 'start';
        RegionSelection.start_time_sec = P.regionTimeRange(1);
        RegionSelection.region_start_sec = P.regionTimeRange(1);
        RegionSelection.end_time_sec = P.regionTimeRange(2);

    case 'region_plan'
        if exist(P.regionPlanFile, 'file') ~= 2
            error('Region plan file not found: %s', P.regionPlanFile);
        end
        T = readtable(P.regionPlanFile);
        row = T(T.regionId == P.regionId, :);
        if height(row) ~= 1
            error('Region plan must contain exactly one row for regionId=%d.', P.regionId);
        end
        RegionSelection.use_region_plan = true;
        RegionSelection.region_id = P.regionId;
        RegionSelection.start_mode = lower(strtrim(P.regionStartMode));
        RegionSelection.region_selected_blade = row.selectedBladeSlot(1);
        RegionSelection.region_start_sec = row.bttStartSec(1);
        RegionSelection.end_time_sec = row.bttEndSec(1);
        RegionSelection.peak_time_sec = row.peakBttTimeSec(1);
        RegionSelection.dominant_order = row.dominantOrder(1);
        RegionSelection.dominant_freq_hz = row.dominantFreqHz(1);
        if strcmpi(RegionSelection.start_mode, 'peak')
            RegionSelection.start_time_sec = RegionSelection.peak_time_sec;
        elseif strcmpi(RegionSelection.start_mode, 'start')
            RegionSelection.start_time_sec = RegionSelection.region_start_sec;
        else
            error('P.regionStartMode must be "peak" or "start".');
        end

    otherwise
        error('Unsupported P.mode: %s', P.mode);
end
end

function label = make_time_label_local(tSec)
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
