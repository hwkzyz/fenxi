%% Step01_BuildLowSpeedReferenceFingerprint_20251222
% Build the 20251222 low-speed reference fingerprint used by NewFlow.
%
% This step reads the Step00 Sensor_Config, extracts per-sensor six-blade
% fingerprints and OPR-center standard angles, and saves a compact reference
% artifact for low-speed numbering and template construction.

clear; close all; clc;

%% Parameters to tune
P = NewFlow_Config_20251222();
analysisSensors = P.sensors.analysis;
bladeCount = P.machine.bladeCount;
viewEnable = true;
saveFigures = true;

%% Paths
sensorConfigFile = P.files.sensorConfig;
outFile = P.files.lowSpeedReference;
summaryFile = strrep(outFile, '.mat', '_Fingerprint.csv');
figureDir = fullfile(P.outputDir, 'figures', '01_low_speed_reference');

ensure_parent_dir_local(outFile);
require_file_local(sensorConfigFile, 'Step00 Sensor_Config');

loadedCfg = load(sensorConfigFile, 'Sensor_Config', 'SensorConfigSummary', 'OPRCenterAngleTable');
Sensor_Config = loadedCfg.Sensor_Config;

validate_sensor_config_local(Sensor_Config, analysisSensors, bladeCount, sensorConfigFile);
fingerprintTable = build_fingerprint_table_local(Sensor_Config, analysisSensors, bladeCount);

LowSpeedReference = struct();
LowSpeedReference.dataset = P.dataset;
LowSpeedReference.mode = 'sensor_config_fingerprint_and_opr_center_angles';
LowSpeedReference.source_file = sensorConfigFile;
LowSpeedReference.analysis_sensors = analysisSensors;
LowSpeedReference.blade_count = bladeCount;
LowSpeedReference.reference_case = Sensor_Config.ReferenceCase;
LowSpeedReference.opr_reference = Sensor_Config.OPRReference;
LowSpeedReference.fingerprints = Sensor_Config.Fingerprints;
LowSpeedReference.target_indices = Sensor_Config.Target_Indices;
LowSpeedReference.standard_angles_opr_center = Sensor_Config.Standard_Relative_Angles_OPRCenter;
LowSpeedReference.standard_angles_start_edge = get_optional_field_local( ...
    Sensor_Config, 'Standard_Relative_Angles_StartEdge', []);
LowSpeedReference.fingerprint_table = fingerprintTable;
if isfield(loadedCfg, 'SensorConfigSummary')
    LowSpeedReference.sensor_config_summary = loadedCfg.SensorConfigSummary;
end
if isfield(loadedCfg, 'OPRCenterAngleTable')
    LowSpeedReference.opr_center_angle_table = loadedCfg.OPRCenterAngleTable;
end

save(outFile, 'LowSpeedReference');
writetable(fingerprintTable, summaryFile);

fprintf('\n=== Step01: low-speed reference fingerprint ===\n');
fprintf('Source Sensor_Config: %s\n', sensorConfigFile);
fprintf('Sensors: %s\n', mat2str(analysisSensors));
fprintf('Saved:   %s\n', outFile);
fprintf('Summary: %s\n', summaryFile);
disp(fingerprintTable);

if viewEnable
    visualize_reference_fingerprint_local(LowSpeedReference, saveFigures, figureDir);
end

function validate_sensor_config_local(Sensor_Config, sensorIds, bladeCount, sourceFile)
requiredFields = {'ReferenceCase', 'Fingerprints', 'Target_Indices', ...
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

function T = build_fingerprint_table_local(Sensor_Config, sensorIds, bladeCount)
rows = repmat(struct( ...
    'SensorID', NaN, ...
    'TargetIndexB1', NaN, ...
    'FingerprintNormMax', NaN, ...
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
for i = 1:numel(sensorIds)
    sid = sensorIds(i);
    fp = force_row_vector_local(Sensor_Config.Fingerprints(sid));
    angles = force_row_vector_local(Sensor_Config.Standard_Relative_Angles_OPRCenter(sid, 1:bladeCount));
    rows(i).SensorID = sid;
    if isKey(Sensor_Config.Target_Indices, sid)
        rows(i).TargetIndexB1 = Sensor_Config.Target_Indices(sid);
    end
    rows(i).FingerprintNormMax = max(fp(1:bladeCount), [], 'omitnan');
    for bladeId = 1:bladeCount
        rows(i).(sprintf('FingerprintB%d', bladeId)) = fp(bladeId);
        rows(i).(sprintf('AngleB%dDeg', bladeId)) = angles(bladeId);
    end
end
T = struct2table(rows);
end

function visualize_reference_fingerprint_local(LowSpeedReference, saveFigures, figureDir)
if saveFigures && exist(figureDir, 'dir') ~= 7
    mkdir(figureDir);
end
sensorIds = LowSpeedReference.analysis_sensors;
bladeIds = 1:LowSpeedReference.blade_count;
fpMat = nan(numel(sensorIds), numel(bladeIds));
angleMat = nan(numel(sensorIds), numel(bladeIds));
for i = 1:numel(sensorIds)
    sid = sensorIds(i);
    fp = LowSpeedReference.fingerprints(sid);
    fpMat(i, :) = fp(1:numel(bladeIds));
    angleMat(i, :) = LowSpeedReference.standard_angles_opr_center(sid, bladeIds);
end
rowMax = max(fpMat, [], 2, 'omitnan');
rowMax(~isfinite(rowMax) | rowMax == 0) = 1;
fpNorm = fpMat ./ rowMax;

fig = figure('Name', 'Step01 low-speed reference fingerprint 20251222', ...
    'Color', 'w', 'Position', [100 100 1200 520], 'NumberTitle', 'off');
tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
imagesc(bladeIds, sensorIds, fpNorm);
set(gca, 'YDir', 'normal');
xlabel('Physical blade ID');
ylabel('Sensor ID');
title('Normalized low-speed fingerprint');
colorbar; box on;

nexttile;
plot(bladeIds, angleMat.', '-o', 'LineWidth', 1.1);
xlabel('Physical blade ID');
ylabel('OPR-center standard angle (deg)');
title('Sensor/blade standard angles');
legend(compose('CH%d', sensorIds), 'Location', 'best');
grid on; box on;

if saveFigures
    exportgraphics(fig, fullfile(figureDir, 'Step01_LowSpeedReferenceFingerprint_20251222.png'), ...
        'Resolution', 300);
end
end

function value = get_optional_field_local(S, name, defaultValue)
if isstruct(S) && isfield(S, name)
    value = S.(name);
else
    value = defaultValue;
end
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
