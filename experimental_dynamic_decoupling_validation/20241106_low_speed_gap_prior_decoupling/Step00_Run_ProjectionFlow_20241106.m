%% Step00: canonical 20241106 projection-flow runner
% This is the recommended entry for the reorganized 20241106 route:
%   1) build one all-blade low-speed template bank for CH2/CH3/CH5/CH7,
%   2) build one all-blade gap calibration bank for capacitive CH5/CH7,
%   3) select a target blade/time/sensor set and run identification.
%
% Defaults live in ProjectionFlow_Config_20241106.m. Temporary overrides:
%   FLOW_RUN_MODE              all | calibration | identify | compare
%   BLADE_CASE_BLADE_ID        target blade, for example 4
%   BLADE_CASE_SENSOR_IDS      analysis sensors, for example "2 5 7"
%   STEP06G_GAP_SENSORS        capacitive gap sensors, normally "5 7"
%   STEP07J_GAP_SENSORS        capacitive gap sensors, normally "5 7"
%   STEP07J_RESULT_SUFFIX      output suffix, for example "_hybrid_s2_bank"
%
% Core case defaults should be edited in ProjectionFlow_Config_20241106.m.

clear; clc;

thisDir = fileparts(mfilename('fullpath'));
oldDir = pwd;
cleanupDir = onCleanup(@() cd(oldDir));
cd(thisDir);

P = ProjectionFlow_Config_20241106();
targetBlade = parse_scalar_env_local('BLADE_CASE_BLADE_ID', P.identification.targetBlade);
analysisSensors = parse_int_env_local('BLADE_CASE_SENSOR_IDS', P.identification.analysisSensors);
gapSensors = parse_int_env_local('STEP07J_GAP_SENSORS', ...
    parse_int_env_local('STEP06G_GAP_SENSORS', P.identification.gapSensors));
P.identification.analysisStartTimeSec = parse_float_env_local('STEP06G_START_TIME_SEC', P.identification.analysisStartTimeSec);
P.identification.targetBladePasses = parse_scalar_env_local('STEP06G_TARGET_LAPS', P.identification.targetBladePasses);
P.identification.windowBladePasses = parse_scalar_env_local('STEP06G_WINDOW_LAPS', P.identification.windowBladePasses);
P.identification.slidingStepBladePasses = parse_scalar_env_local('STEP06G_SLIDING_STEP_LAPS', P.identification.slidingStepBladePasses);
gapSensors = intersect(gapSensors(:).', analysisSensors(:).', 'stable');
assert_capacitive_gap_sensors_local(gapSensors);

runMode = lower(strtrim(getenv('FLOW_RUN_MODE')));
if isempty(runMode)
    runMode = P.run.mode;
end
validModes = {'all','calibration','identify','compare'};
if ~ismember(runMode, validModes)
    error('FLOW_RUN_MODE must be one of: %s.', strjoin(validModes, ', '));
end

fprintf('\n=== Step00: 20241106 projection flow ===\n');
fprintf('Mode: %s\n', runMode);
fprintf('Default calibration blades: %s\n', mat2str(P.calibration.buildBladeIds));
fprintf('Low-speed template sensors: %s\n', mat2str(P.calibration.lowSpeedTemplateSensors));
fprintf('Gap-calibration sensors: %s\n', mat2str(P.calibration.gapSensorIds));
fprintf('Selected target: B%d, analysis sensors %s, gap sensors %s\n', ...
    targetBlade, mat2str(analysisSensors), mat2str(gapSensors));
fprintf('Start %.6f s, laps %d, window %d, step %d, max run windows %s\n', ...
    P.identification.analysisStartTimeSec, P.identification.targetBladePasses, ...
    P.identification.windowBladePasses, P.identification.slidingStepBladePasses, ...
    numeric_limit_text_local(P.identification.maxRunWindows));

if isempty(strtrim(getenv('STEP07J_RESULT_SUFFIX'))) && ...
        isfield(P, 'run') && isfield(P.run, 'resultSuffix') && ~isempty(P.run.resultSuffix)
    setenv('STEP07J_RESULT_SUFFIX', P.run.resultSuffix);
end

if any(strcmp(runMode, {'all','calibration'}))
    run_script_local('Step06A_BuildLowSpeedTemplateBank_20241106.m');
    run_script_local('Step06I_Calibrate_AllBladeGapLibrary_20241106.m');
end

if any(strcmp(runMode, {'all','identify','compare'}))
    sync_case_env_local(P, targetBlade, analysisSensors, gapSensors);
    ensure_step07j_direct_inputs_local(thisDir, P, targetBlade, analysisSensors);
    run_script_local('Step06_BuildGapAwareDynamicMap_20241106.m');
end

if any(strcmp(runMode, {'all','identify'}))
    sync_case_env_local(P, targetBlade, analysisSensors, gapSensors);
    ensure_step07j_direct_inputs_local(thisDir, P, targetBlade, analysisSensors);
    run_script_local('Step07J_NestedStaticWarp_VPFullWave_20241106.m');
end

if any(strcmp(runMode, {'all','compare'}))
    oldRunMode = getenv('STEP07J_RUN_MODE');
    oldSuffix = getenv('STEP07J_RESULT_SUFFIX');
    oldCompareSuffix = getenv('STEP07J_COMPARE_SUFFIX');
    cleanupEnv = onCleanup(@() restore_env_local(oldRunMode, oldSuffix, oldCompareSuffix));

    setenv('STEP07J_RUN_MODE', 'comparison');
    resultSuffix = strtrim(getenv('STEP07J_RESULT_SUFFIX'));
    if isempty(resultSuffix)
        resultSuffix = '_bank_compare';
        setenv('STEP07J_RESULT_SUFFIX', resultSuffix);
    end
    sync_case_env_local(P, targetBlade, analysisSensors, gapSensors);
    ensure_step07j_direct_inputs_local(thisDir, P, targetBlade, analysisSensors);
    run_script_local('Step07J_NestedStaticWarp_VPFullWave_20241106.m');

    setenv('STEP07J_COMPARE_SUFFIX', strip_leading_underscore_local(resultSuffix));
    Step07K_Compare_NestedStaticWarp_20241106();
    clear cleanupEnv;
    restore_env_local(oldRunMode, oldSuffix, oldCompareSuffix);
end

fprintf('\nStep00 complete.\n');

function sync_case_env_local(P, targetBlade, analysisSensors, gapSensors)
sensorText = int_vector_to_env_text_local(analysisSensors);
gapText = int_vector_to_env_text_local(gapSensors);
setenv('BLADE_CASE_BLADE_ID', sprintf('%d', targetBlade));
setenv('BLADE_CASE_SENSOR_IDS', sensorText);
setenv('STEP06G_TARGET_BLADE', sprintf('%d', targetBlade));
setenv('STEP06G_ANALYSIS_SENSORS', sensorText);
setenv('STEP06G_GAP_SENSORS', gapText);
setenv('STEP06G_START_TIME_SEC', sprintf('%.15g', P.identification.analysisStartTimeSec));
setenv('STEP06G_TARGET_LAPS', sprintf('%d', P.identification.targetBladePasses));
setenv('STEP06G_WINDOW_LAPS', sprintf('%d', P.identification.windowBladePasses));
setenv('STEP06G_SLIDING_STEP_LAPS', sprintf('%d', P.identification.slidingStepBladePasses));
setenv('STEP07J_ANALYSIS_SENSORS', sensorText);
setenv('STEP07J_GAP_SENSORS', gapText);
end

function run_script_local(scriptName)
STEP06_SKIP_CLEAR_FOR_DRIVER = true; %#ok<NASGU>
run(scriptName);
end

function ensure_step07j_direct_inputs_local(gapRouteDir, P, targetBlade, analysisSensors)
rotDir = fullfile(fileparts(gapRouteDir), '20241106_low_speed_rotating_calibration');
sensorTag = ['S', sprintf('%d', analysisSensors)];
timeTag = step07j_time_label_local(P.identification.analysisStartTimeSec);
compactFile = fullfile(rotDir, 'output', 'new_flow', '05_waveform_library', sensorTag, ...
    sprintf('Step06_CompactWaveformCase_B%d_%s_20241106.mat', targetBlade, timeTag));
directDir = fullfile(rotDir, 'output', 'new_flow', '06_identification', sensorTag);
directPatterns = {
    sprintf('IdentificationResult_%s_20241106_B%d_%s_L03_eta200.mat', timeTag, targetBlade, sensorTag)
    sprintf('IdentificationResult_%s_20241106_B%d_%s_L03_eta000.mat', timeTag, targetBlade, sensorTag)
    sprintf('IdentificationResult_%s_20241106_B%d_%s_*.mat', timeTag, targetBlade, sensorTag)
    sprintf('IdentificationResult_%s_20241106_B%d_only.mat', timeTag, targetBlade)
    };
if isfile(compactFile) && has_matching_file_local(directDir, directPatterns)
    return;
end
setenv('STEP06_TARGET_BLADES', sprintf('%d', targetBlade));
setenv('STEP06_ANALYSIS_SENSORS', strtrim(sprintf('%d ', analysisSensors)));
setenv('STEP06_START_TIME_SEC', sprintf('%.15g', P.identification.analysisStartTimeSec));
setenv('STEP06_TARGET_LAPS', sprintf('%d', P.identification.targetBladePasses));
setenv('STEP06_WINDOW_LAPS', sprintf('%d', P.identification.windowBladePasses));
setenv('STEP06_SLIDING_STEP_LAPS', sprintf('%d', P.identification.slidingStepBladePasses));
setenv('STEP06_OUTPUT_LABEL', sprintf('B%d_only', targetBlade));
run_script_local('..\20241106_low_speed_rotating_calibration\Step06_RunIdentificationByBlade_20241106.m');
end

function tf = has_matching_file_local(folderName, patterns)
tf = false;
for ip = 1:numel(patterns)
    if ~isempty(dir(fullfile(folderName, patterns{ip})))
        tf = true;
        return;
    end
end
end

function label = step07j_time_label_local(tSec)
if abs(tSec - round(tSec)) < 1e-9
    label = sprintf('T%03ds', round(tSec));
else
    label = sprintf('T%07.3fs', tSec);
    label = strrep(label, '.', 'p');
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

function txt = int_vector_to_env_text_local(values)
txt = strtrim(sprintf('%d ', values(:).'));
end

function suffix = strip_leading_underscore_local(suffix)
while startsWith(suffix, '_')
    suffix = extractAfter(suffix, 1);
end
suffix = char(suffix);
end

function restore_env_local(oldRunMode, oldSuffix, oldCompareSuffix)
restore_one_env_local('STEP07J_RUN_MODE', oldRunMode);
restore_one_env_local('STEP07J_RESULT_SUFFIX', oldSuffix);
restore_one_env_local('STEP07J_COMPARE_SUFFIX', oldCompareSuffix);
end

function restore_one_env_local(name, value)
if isempty(value)
    setenv(name, '');
else
    setenv(name, value);
end
end

function txt = numeric_limit_text_local(value)
if isinf(value)
    txt = 'inf';
else
    txt = sprintf('%g', value);
end
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

function assert_capacitive_gap_sensors_local(sensorIds)
invalid = setdiff(sensorIds(:).', [5 7]);
if ~isempty(invalid)
    error(['Only capacitive sensors CH5/CH7 can use the gap library. ' ...
        'Invalid gap sensor(s): %s'], mat2str(invalid));
end
end
