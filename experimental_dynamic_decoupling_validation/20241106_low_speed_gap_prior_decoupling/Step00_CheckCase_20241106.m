%% Step00_CheckCase: print the active 20241106 projection-flow case
% This script does not run identification. Use it before Step00_Run when
% changing blade, sensor, time, or window parameters.

clear; clc;

thisDir = fileparts(mfilename('fullpath'));
oldDir = pwd;
cleanupDir = onCleanup(@() cd(oldDir));
cd(thisDir);

P = ProjectionFlow_Config_20241106();
outDir = fullfile(thisDir, 'outputs');
rootDir = fileparts(thisDir);
rotDir = fullfile(rootDir, '20241106_low_speed_rotating_calibration');

targetBlade = parse_scalar_env_local('BLADE_CASE_BLADE_ID', P.identification.targetBlade);
analysisSensors = parse_int_env_local('BLADE_CASE_SENSOR_IDS', P.identification.analysisSensors);
gapSensors = parse_int_env_local('STEP07J_GAP_SENSORS', ...
    parse_int_env_local('STEP06G_GAP_SENSORS', P.identification.gapSensors));
analysisSensors = reshape(analysisSensors, 1, []);
gapSensors = intersect(gapSensors(:).', analysisSensors, 'stable');
directOnlySensors = setdiff(analysisSensors, gapSensors, 'stable');
sensorTag = ['S', sprintf('%d', analysisSensors)];
gapSensorTag = ['S', sprintf('%d', gapSensors)];
P.identification.analysisStartTimeSec = parse_float_env_local('STEP06G_START_TIME_SEC', P.identification.analysisStartTimeSec);
P.identification.targetBladePasses = parse_scalar_env_local('STEP06G_TARGET_LAPS', P.identification.targetBladePasses);
P.identification.windowBladePasses = parse_scalar_env_local('STEP06G_WINDOW_LAPS', P.identification.windowBladePasses);
P.identification.slidingStepBladePasses = parse_scalar_env_local('STEP06G_SLIDING_STEP_LAPS', P.identification.slidingStepBladePasses);
timeTag = time_label_local(P.identification.analysisStartTimeSec);
plannedWindows = planned_window_count_local(P.identification.targetBladePasses, ...
    P.identification.windowBladePasses, P.identification.slidingStepBladePasses);
runWindows = min(plannedWindows, P.identification.maxRunWindows);

fprintf('\n=== Step00_CheckCase 20241106 ===\n');
fprintf('Run mode      : %s\n', P.run.mode);
fprintf('Target        : B%d, sensors %s, gap sensors %s, direct-only %s\n', ...
    targetBlade, mat2str(analysisSensors), mat2str(gapSensors), mat2str(directOnlySensors));
fprintf('Start/window  : %.6f s (%s), target laps %d, window %d, step %d\n', ...
    P.identification.analysisStartTimeSec, timeTag, ...
    P.identification.targetBladePasses, P.identification.windowBladePasses, ...
    P.identification.slidingStepBladePasses);
fprintf('Windows       : planned %d, Step07J run limit %s, expected run %d\n', ...
    plannedWindows, numeric_limit_text_local(P.identification.maxRunWindows), runWindows);
fprintf('Bundle policy : direct bundle first; compact_case reconstruction fallback; DynamicMap fallback disabled unless explicitly enabled for debugging\n');

fprintf('\nKey files:\n');
responseOk = print_file_status_local('Response surface', fullfile(outDir, ...
    'Step05I_OffsetTilt_Shared_Response_Surface_20241106.mat'));
templateOk = print_file_status_local('Template bank', fullfile(outDir, P.calibration.templateBankFile));
gapBankOk = print_file_status_local('Gap bank', fullfile(outDir, P.calibration.gapBankFile));

dynamicMapFile = fullfile(outDir, sprintf( ...
    'Step06_BuildGapAwareDynamicMap_20241106_B%d_%s.mat', targetBlade, sensorTag));
dynamicMapOk = print_file_status_local('Step06 DynamicMap', dynamicMapFile);
dynamicMapFresh = false;
if isfile(dynamicMapFile)
    S = load(dynamicMapFile, 'DynamicMap');
    if isfield(S, 'DynamicMap') && isfield(S.DynamicMap, 'Window')
        fprintf('    DynamicMap windows: %d\n', numel(S.DynamicMap.Window));
        if isfield(S.DynamicMap, 'SelectionInfo') && isfield(S.DynamicMap.SelectionInfo, 'start_time_sec')
            startDelta = abs(S.DynamicMap.SelectionInfo.start_time_sec - ...
                P.identification.analysisStartTimeSec);
            if startDelta <= 1e-6
                fprintf('    DynamicMap start  : %.6f s\n', S.DynamicMap.SelectionInfo.start_time_sec);
                dynamicMapFresh = true;
            else
                fprintf('    [STALE] DynamicMap start %.6f s does not match config %.6f s; rerun Step06.\n', ...
                    S.DynamicMap.SelectionInfo.start_time_sec, P.identification.analysisStartTimeSec);
            end
        else
            dynamicMapFresh = true;
        end
    end
end

directDir = fullfile(rotDir, 'output', 'new_flow', '06_identification', sensorTag);
directPatterns = {
    sprintf('IdentificationResult_%s_20241106_B%d_%s_L03_eta200.mat', timeTag, targetBlade, sensorTag)
    sprintf('IdentificationResult_%s_20241106_B%d_%s_L03_eta000.mat', timeTag, targetBlade, sensorTag)
    sprintf('IdentificationResult_%s_20241106_B%d_%s_*.mat', timeTag, targetBlade, sensorTag)
    sprintf('IdentificationResult_%s_20241106_B%d_only.mat', timeTag, targetBlade)
    };
[directOk, directFile] = print_patterns_status_local('Direct reference', directDir, directPatterns);

compactDir = fullfile(rotDir, 'output', 'new_flow', '05_waveform_library', sensorTag);
compactPatterns = {
    sprintf('Step06_CompactWaveformCase_B%d_%s_20241106.mat', targetBlade, timeTag)
    };
[compactOk, compactFile] = print_patterns_status_local('Compact waveform', compactDir, compactPatterns);

debugDynamicAllowed = strcmpi(strtrim(getenv('STEP07J_LOCAL_BUNDLE_SOURCE')), 'dynamic_map') && ...
    parse_logical_env_local('STEP07J_ALLOW_DYNAMIC_MAP_FALLBACK', false);
baseFilesOk = responseOk && templateOk && gapBankOk && dynamicMapOk && dynamicMapFresh;
officialBundleOk = directOk || compactOk;

fprintf('\nStep07J readiness:\n');
if baseFilesOk && officialBundleOk
    if directOk
        fprintf('  [READY]   Official route can run from direct bundle: %s\n', directFile);
    else
        fprintf('  [READY]   Official route can run from compact waveform: %s\n', compactFile);
    end
elseif baseFilesOk && debugDynamicAllowed
    fprintf('  [DEBUG]   DynamicMap fallback is explicitly enabled. Use only for diagnosis, not final results.\n');
elseif baseFilesOk
    fprintf('  [BLOCKED] Missing same-time direct reference and compact waveform case; rerun old direct Step06 to generate them.\n');
else
    fprintf('  [BLOCKED] Missing or stale base files above; fix those before Step07J.\n');
end

fprintf('\nEnvironment overrides to watch:\n');
print_env_local('FLOW_RUN_MODE');
print_env_local('BLADE_CASE_BLADE_ID');
print_env_local('BLADE_CASE_SENSOR_IDS');
print_env_local('STEP06G_START_TIME_SEC');
print_env_local('STEP06G_TARGET_LAPS');
print_env_local('STEP06G_WINDOW_LAPS');
print_env_local('STEP06G_SLIDING_STEP_LAPS');
print_env_local('STEP07J_MAX_WINDOWS', 'ignored by current Step07J; use P.identification.maxRunWindows');
print_env_local('STEP07J_LOCAL_BUNDLE_SOURCE');
print_env_local('STEP07J_ALLOW_DYNAMIC_MAP_FALLBACK', 'must be 1 to allow debug-only DynamicMap fallback');

fprintf('\nStep00_CheckCase complete.\n');

function n = planned_window_count_local(targetLaps, windowLaps, stepLaps)
if targetLaps < windowLaps
    error('targetBladePasses must be >= windowBladePasses.');
end
n = floor((targetLaps - windowLaps) / stepLaps) + 1;
end

function label = time_label_local(tSec)
if abs(tSec - round(tSec)) < 1e-9
    label = sprintf('T%03ds', round(tSec));
else
    label = sprintf('T%07.3fs', tSec);
    label = strrep(label, '.', 'p');
end
end

function txt = numeric_limit_text_local(value)
if isinf(value)
    txt = 'inf';
else
    txt = sprintf('%g', value);
end
end

function value = parse_scalar_env_local(name, defaultValue)
raw = strtrim(getenv(name));
if isempty(raw)
    value = defaultValue;
    return;
end
value = str2double(raw);
if ~isfinite(value) || value < 1 || abs(value - round(value)) > eps(value)
    error('%s must be a positive integer.', name);
end
value = round(value);
end

function value = parse_float_env_local(name, defaultValue)
raw = strtrim(getenv(name));
if isempty(raw)
    value = defaultValue;
    return;
end
value = str2double(raw);
if ~isfinite(value) || value <= 0
    error('%s must be a positive number.', name);
end
end

function values = parse_int_env_local(name, defaultValue)
raw = strtrim(getenv(name));
if isempty(raw)
    values = defaultValue;
    return;
end
values = sscanf(raw, '%d').';
if isempty(values) || any(values < 1)
    error('%s must contain positive integer sensor IDs separated by spaces.', name);
end
end

function existsFlag = print_file_status_local(label, fileName)
existsFlag = isfile(fileName);
if existsFlag
    fprintf('  [OK]      %-18s %s\n', label, fileName);
else
    fprintf('  [MISSING] %-18s %s\n', label, fileName);
end
end

function [foundFlag, foundFile] = print_patterns_status_local(label, folderName, patterns)
foundFlag = false;
foundFile = '';
for ip = 1:numel(patterns)
    files = dir(fullfile(folderName, patterns{ip}));
    if ~isempty(files)
        [~, idx] = max([files.datenum]);
        foundFlag = true;
        foundFile = fullfile(files(idx).folder, files(idx).name);
        fprintf('  [OK]      %-18s %s\n', label, foundFile);
        if ip > 1
            fprintf('           %-18s matched fallback pattern: %s\n', '', patterns{ip});
        end
        return;
    end
end
fprintf('  [MISSING] %-18s %s\n', label, fullfile(folderName, patterns{1}));
end

function value = parse_logical_env_local(name, defaultValue)
txt = lower(strtrim(getenv(name)));
if isempty(txt)
    value = defaultValue;
elseif ismember(txt, {'1','true','yes','on'})
    value = true;
elseif ismember(txt, {'0','false','no','off'})
    value = false;
else
    value = defaultValue;
end
end

function print_env_local(name, note)
if nargin < 2
    note = '';
end
value = strtrim(getenv(name));
if isempty(value)
    value = '<empty>';
end
if isempty(note) || strcmp(value, '<empty>')
    fprintf('  %-30s %s\n', name, value);
else
    fprintf('  %-30s %s  (%s)\n', name, value, note);
end
end
