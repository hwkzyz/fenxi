function P = ProjectionFlow_Config_20251222()
%PROJECTIONFLOW_CONFIG_20251222 Central config for the projection-limited route.
%
% This file defines the calibration-first workflow without changing the
% numerical defaults used by the current Step06/Step07J route. The intended
% flow is:
%   1) precompute blade-sensor calibration entries,
%   2) select target high-speed waveforms by blade/time/window parameters,
%   3) identify vibration with the Delta g_proj-limited gap_tilt model.

P = struct();
P.dataset = '20251222';

P.machine.bladeCount = 6;
P.machine.oprPulsesPerRev = 6;
P.machine.oprNote = ['This dataset stores six OPR timing pulses per rotor ' ...
    'revolution. For one fixed target blade, one selected blade pass ' ...
    'corresponds to one rotor revolution.'];

P.calibration.buildBladeIds = 1:P.machine.bladeCount;
P.calibration.candidateSensorIds = [1 2 3];
P.calibration.mainSensorIds = [1 2 3];
P.calibration.libraryDirName = 'gap_calibration_library';
P.calibration.filePrefix = 'GapCalib';
P.calibration.note = ['Calibration should be prepared per blade-sensor pair. ' ...
    'Identification loads only the entries required by the selected blade ' ...
    'and sensor set.'];

P.identification.analysisSensors = P.calibration.mainSensorIds;
P.identification.dynamicCaseName = '1000_2500_3500';
P.identification.analysisStartMode = 'fixed_resonance_region_catalog';
[resonanceSelection, resonanceCatalog] = ResonanceRegionCatalog_20251222();
P.identification.targetBlade = resonanceSelection.representativeBladeId;
P.identification.resonanceRegionId = resonanceSelection.regionId;
P.identification.resonanceRegionTag = resonanceSelection.tag;
P.identification.resonanceRegionShortTag = resonanceSelection.shortTag;
P.identification.resonanceRegionStartSec = resonanceSelection.regionStartSec;
P.identification.resonanceRegionEndSec = resonanceSelection.regionEndSec;
P.identification.resonanceRegionCatalog = resonanceCatalog;
P.identification.analysisStartTimeSec = resonanceSelection.analysisStartTimeSec;
P.identification.targetBladePasses = resonanceSelection.targetBladePasses;
P.identification.windowBladePasses = resonanceSelection.windowBladePasses;
P.identification.slidingStepBladePasses = resonanceSelection.slidingStepBladePasses;
P.identification.referenceOrder = resonanceSelection.dominantOrder;
P.identification.pulseWindowSec = 6e-4;

[P.identification.analysisStartTimeSec, hasStartOverride] = parse_numeric_override_local( ...
    'BLADE_ANALYSIS_START_TIME_SEC', P.identification.analysisStartTimeSec, true);
P.identification.targetBladePasses = parse_integer_override_local( ...
    'BLADE_TARGET_LAPS', P.identification.targetBladePasses);
P.identification.windowBladePasses = parse_integer_override_local( ...
    'BLADE_WINDOW_LAPS', P.identification.windowBladePasses);
P.identification.slidingStepBladePasses = parse_integer_override_local( ...
    'BLADE_SLIDING_STEP_LAPS', P.identification.slidingStepBladePasses);
P.identification.pulseWindowSec = parse_numeric_override_local( ...
    'BLADE_PULSE_WINDOW_SEC', P.identification.pulseWindowSec, true);
if hasStartOverride
    P.identification.analysisStartMode = 'runtime_start_time_override';
end

P.projection.useDeltaGapProjectionLimit = true;
P.projection.deltaGapLimitMm = 0.25;
P.projection.alpha = 2.5;
P.projection.minLimitMm = 0.02;
P.projection.maxLimitMm = P.projection.deltaGapLimitMm;
P.projection.fallbackLimitMm = P.projection.deltaGapLimitMm;
P.projection.derivativeStepMm = 1e-4;
P.projection.sensitivityFloorMvPerMm = 1e-6;
P.projection.deltaMuLimit = 0.030;
P.projection.deltaTauLimitMm = 0.16;
end

function value = parse_integer_override_local(name, defaultValue)
raw = strtrim(getenv(name));
if isempty(raw)
    value = defaultValue;
    return;
end
value = str2double(raw);
if ~isfinite(value) || value <= 0 || abs(value - round(value)) > eps(value)
    error('%s must be a positive integer.', name);
end
value = round(value);
end

function [value, hasOverride] = parse_numeric_override_local(name, defaultValue, mustBePositive)
raw = strtrim(getenv(name));
hasOverride = ~isempty(raw);
if ~hasOverride
    value = defaultValue;
    return;
end
value = str2double(raw);
if ~isfinite(value) || (mustBePositive && value <= 0)
    error('%s must be a finite positive numeric value.', name);
end
end
