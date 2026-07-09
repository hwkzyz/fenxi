function P = ProjectionFlow_Config_20250527()
%PROJECTIONFLOW_CONFIG_20250527 Central config for the projection-limited route.
%
% This file defines the calibration-first workflow without changing the
% numerical defaults used by the current Step06/Step07J route. The intended
% flow is:
%   1) precompute blade-sensor calibration entries,
%   2) select target high-speed waveforms by blade/time/window parameters,
%   3) identify vibration with the Delta g_proj-limited gap_tilt model.

P = struct();
P.dataset = '20250527';

P.machine.bladeCount = 6;
P.machine.oprPulsesPerRev = 6;
P.machine.oprNote = ['This dataset stores six OPR timing pulses per rotor ' ...
    'revolution. For one fixed target blade, one selected blade pass ' ...
    'corresponds to one rotor revolution.'];

P.calibration.buildBladeIds = 1:P.machine.bladeCount;
P.calibration.candidateSensorIds = [1 3 6 7 8];
P.calibration.mainSensorIds = [1 3 6];
P.calibration.libraryDirName = 'gap_calibration_library';
P.calibration.filePrefix = 'GapCalib';
P.calibration.note = ['Calibration should be prepared per blade-sensor pair. ' ...
    'Identification loads only the entries required by the selected blade ' ...
    'and sensor set.'];

P.identification.targetBlade = 1;
P.identification.analysisSensors = P.calibration.mainSensorIds;
P.identification.dynamicCaseName = read_env_text_or_default_local( ...
    'BLADE_DYNAMIC_CASE_NAME', '20250526_2500-3500_t400');
P.identification.analysisStartMode = 'manual_same_as_low_speed_template_route';
P.identification.analysisStartTimeSec = read_env_double_or_default_local( ...
    'BLADE_ANALYSIS_START_TIME_SEC', 1.5, 0);
P.identification.targetBladePasses = read_env_integer_or_default_local( ...
    'BLADE_TARGET_LAPS', 20, 1);
P.identification.windowBladePasses = read_env_integer_or_default_local( ...
    'BLADE_WINDOW_LAPS', 3, 1);
P.identification.slidingStepBladePasses = read_env_integer_or_default_local( ...
    'BLADE_SLIDING_STEP_LAPS', 1, 1);
P.identification.referenceOrder = 14;
P.identification.pulseWindowSec = read_env_double_or_default_local( ...
    'BLADE_PULSE_WINDOW_SEC', 6e-4, 0);

P.projection.useDeltaGapProjectionLimit = true;
P.projection.deltaGapLimitMm = 0.25;
P.projection.alpha = 2.5;
P.projection.minLimitMm = 0.02;
P.projection.maxLimitMm = P.projection.deltaGapLimitMm;
P.projection.fallbackLimitMm = P.projection.deltaGapLimitMm;
P.projection.derivativeStepMm = 1e-4;
P.projection.sensitivityFloorMvPerMm = 1e-6;
P.projection.deltaMuLimit = 0.080;
P.projection.deltaTauLimitMm = 0.16;
end


function value = read_env_text_or_default_local(name, defaultValue)
raw = strtrim(getenv(name));
if isempty(raw)
    value = defaultValue;
else
    value = raw;
end
end


function value = read_env_double_or_default_local(name, defaultValue, minValue)
raw = strtrim(getenv(name));
if isempty(raw)
    value = defaultValue;
    return;
end
value = str2double(raw);
if ~isfinite(value) || value < minValue
    error('%s must be a finite number >= %.6g.', name, minValue);
end
end


function value = read_env_integer_or_default_local(name, defaultValue, minValue)
raw = strtrim(getenv(name));
if isempty(raw)
    value = defaultValue;
    return;
end
value = str2double(raw);
if ~isfinite(value) || abs(value - round(value)) > eps(value) || value < minValue
    error('%s must be an integer >= %d.', name, minValue);
end
value = round(value);
end
