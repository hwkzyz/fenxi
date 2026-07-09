%% Step00_BuildLowSpeedSensorConfig_20251222
% Standardize and audit the 20251222 low-speed Sensor_Config artifact.
%
% This step is the source of the NewFlow low-speed coordinate/numbering
% chain. The mature 20251222 Sensor_Config is produced by the legacy
% preprocessing route, but this script does more than forward that file:
% it validates required fields, extracts the active sensor/blade table, and
% writes a local NewFlow copy plus CSV summaries.

clear; close all; clc;

%% Parameters to tune
P = NewFlow_Config_20251222();
analysisSensors = P.sensors.analysis;
bladeCount = P.machine.bladeCount;

%% Paths
sourceFile = P.data.sensorConfigFile;
outFile = P.files.sensorConfig;
summaryFile = strrep(outFile, '.mat', '_Summary.csv');
angleFile = strrep(outFile, '.mat', '_OPRCenterAngles.csv');

ensure_parent_dir_local(outFile);
require_file_local(sourceFile, 'source Sensor_Config_20251222.mat');

loaded = load(sourceFile, 'Sensor_Config');
if ~isfield(loaded, 'Sensor_Config')
    error('Source file does not contain variable Sensor_Config:\n  %s', sourceFile);
end
Sensor_Config = loaded.Sensor_Config;

validate_sensor_config_local(Sensor_Config, analysisSensors, bladeCount, sourceFile);
[SensorConfigSummary, OPRCenterAngleTable] = build_sensor_config_tables_local( ...
    Sensor_Config, analysisSensors, bladeCount);

Step00Meta = struct();
Step00Meta.dataset = P.dataset;
Step00Meta.source_file = sourceFile;
Step00Meta.analysis_sensors = analysisSensors;
Step00Meta.blade_count = bladeCount;
Step00Meta.opr_pulses_per_rev = P.machine.oprPulsesPerRev;
Step00Meta.output_file = outFile;
Step00Meta.summary_file = summaryFile;
Step00Meta.angle_file = angleFile;
Step00Meta.role = 'NewFlow low-speed sensor configuration and OPR-center angle source';

save(outFile, 'Sensor_Config', 'SensorConfigSummary', 'OPRCenterAngleTable', 'Step00Meta');
writetable(SensorConfigSummary, summaryFile);
writetable(OPRCenterAngleTable, angleFile);

fprintf('\n=== Step00: low-speed Sensor_Config standardization ===\n');
fprintf('Source:  %s\n', sourceFile);
fprintf('Sensors: %s\n', mat2str(analysisSensors));
fprintf('Saved:   %s\n', outFile);
fprintf('Summary: %s\n', summaryFile);
disp(SensorConfigSummary);

function validate_sensor_config_local(Sensor_Config, sensorIds, bladeCount, sourceFile)
requiredFields = {'Sensor_IDs', 'Blades_Num', 'Fingerprints', 'Target_Indices', ...
    'Standard_Relative_Angles_OPRCenter', 'OPRReference'};
for i = 1:numel(requiredFields)
    name = requiredFields{i};
    if ~isfield(Sensor_Config, name)
        error('Sensor_Config is missing field %s:\n  %s', name, sourceFile);
    end
end
if ~isa(Sensor_Config.Fingerprints, 'containers.Map')
    error('Sensor_Config.Fingerprints must be containers.Map.');
end
if ~isa(Sensor_Config.Target_Indices, 'containers.Map')
    error('Sensor_Config.Target_Indices must be containers.Map.');
end
if size(Sensor_Config.Standard_Relative_Angles_OPRCenter, 2) < bladeCount
    error('Standard_Relative_Angles_OPRCenter has fewer than %d blade columns.', bladeCount);
end
for sid = sensorIds(:).'
    if ~isKey(Sensor_Config.Fingerprints, sid)
        error('Sensor_Config.Fingerprints does not contain CH%d.', sid);
    end
    fp = Sensor_Config.Fingerprints(sid);
    if numel(fp) < bladeCount
        error('CH%d fingerprint has %d values; expected at least %d.', ...
            sid, numel(fp), bladeCount);
    end
end
end

function [summaryTable, angleTable] = build_sensor_config_tables_local( ...
        Sensor_Config, sensorIds, bladeCount)
summaryRows = repmat(struct( ...
    'SensorID', NaN, ...
    'HasTargetIndex', false, ...
    'TargetIndexB1', NaN, ...
    'AngleReference', "", ...
    'FingerprintB1', NaN, ...
    'FingerprintB2', NaN, ...
    'FingerprintB3', NaN, ...
    'FingerprintB4', NaN, ...
    'FingerprintB5', NaN, ...
    'FingerprintB6', NaN, ...
    'AngleB1Deg', NaN, ...
    'AngleB2Deg', NaN, ...
    'AngleB3Deg', NaN, ...
    'AngleB4Deg', NaN, ...
    'AngleB5Deg', NaN, ...
    'AngleB6Deg', NaN), numel(sensorIds), 1);
angleRows = repmat(struct( ...
    'SensorID', NaN, ...
    'BladeID', NaN, ...
    'OPRCenterAngleDeg', NaN, ...
    'FingerprintPeak', NaN), numel(sensorIds) * bladeCount, 1);

angleReference = "opr_pulse_center";
if isfield(Sensor_Config, 'Standard_Relative_Angles_Reference')
    angleReference = string(Sensor_Config.Standard_Relative_Angles_Reference);
end

rowIdx = 0;
for i = 1:numel(sensorIds)
    sid = sensorIds(i);
    fp = force_row_vector_local(Sensor_Config.Fingerprints(sid));
    angles = force_row_vector_local(Sensor_Config.Standard_Relative_Angles_OPRCenter(sid, 1:bladeCount));

    summaryRows(i).SensorID = sid;
    summaryRows(i).AngleReference = angleReference;
    if isKey(Sensor_Config.Target_Indices, sid)
        summaryRows(i).HasTargetIndex = true;
        summaryRows(i).TargetIndexB1 = Sensor_Config.Target_Indices(sid);
    end
    for bladeId = 1:bladeCount
        summaryRows(i).(sprintf('FingerprintB%d', bladeId)) = fp(bladeId);
        summaryRows(i).(sprintf('AngleB%dDeg', bladeId)) = angles(bladeId);
        rowIdx = rowIdx + 1;
        angleRows(rowIdx).SensorID = sid;
        angleRows(rowIdx).BladeID = bladeId;
        angleRows(rowIdx).OPRCenterAngleDeg = angles(bladeId);
        angleRows(rowIdx).FingerprintPeak = fp(bladeId);
    end
end
summaryTable = struct2table(summaryRows);
angleTable = struct2table(angleRows);
end

function x = force_row_vector_local(x)
x = x(:).';
end

function require_file_local(pathText, label)
if exist(pathText, 'file') ~= 2
    error('Missing %s:\n  %s', label, pathText);
end
end

function ensure_parent_dir_local(pathText)
folder = fileparts(pathText);
if exist(folder, 'dir') ~= 7
    mkdir(folder);
end
end
