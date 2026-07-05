%% Step01_BuildLowSpeedReferenceFingerprint_20241106
% Build the low-speed reference fingerprint artifact for the new numbering
% workflow. This step reads the local Step00 Sensor_Config file and writes a
% compact reference artifact for later numbering/matching steps.

clear; close all; clc;

%% Parameters to tune
analysisSensors = [2 3 5 7];
bladeCount = 6;
lowSpeedCase = '900';

viewEnable = true;
saveFigures = true;

%% Paths
routeDir = fileparts(mfilename('fullpath'));
lowSpeedReferenceDir = fullfile(routeDir, 'output', 'reference', [lowSpeedCase, '_reference']);

sensorConfigFile = fullfile(lowSpeedReferenceDir, 'Sensor_Config_20241106.mat');
outDir = fullfile(routeDir, 'output', 'new_flow', '01_low_speed_reference');
outFile = fullfile(outDir, 'LowSpeedReferenceFingerprint_20241106.mat');
figureDir = fullfile(routeDir, 'output', 'new_flow', 'figures', '01_low_speed_reference');

if exist(outDir, 'dir') ~= 7
    mkdir(outDir);
end

if exist(sensorConfigFile, 'file') ~= 2
    error(['Missing low-speed Sensor_Config input:\n  %s\n\n' ...
        'Run Step00_BuildLowSpeedSensorConfig_20241106 first, then rerun Step01.'], ...
        sensorConfigFile);
end

loaded = load(sensorConfigFile, 'Sensor_Config');
if ~isfield(loaded, 'Sensor_Config')
    error('File does not contain variable Sensor_Config:\n  %s', sensorConfigFile);
end
Sensor_Config = loaded.Sensor_Config;

validate_sensor_config_local(Sensor_Config, analysisSensors, bladeCount, sensorConfigFile);
fingerprintRows = build_fingerprint_table_local(Sensor_Config, analysisSensors, bladeCount);

LowSpeedReference = struct();
LowSpeedReference.source_file = sensorConfigFile;
LowSpeedReference.analysis_sensors = analysisSensors;
LowSpeedReference.blade_count = bladeCount;
LowSpeedReference.reference_sensor_id = Sensor_Config.ReferenceSensorID;
LowSpeedReference.reference_revolution_id = Sensor_Config.ReferenceRevolutionID;
LowSpeedReference.reference_fingerprint = Sensor_Config.ReferenceFingerprint;
LowSpeedReference.fingerprint_table = fingerprintRows;
LowSpeedReference.fingerprints = Sensor_Config.Fingerprints;
LowSpeedReference.target_indices = Sensor_Config.Target_Indices;
LowSpeedReference.revolution_ids = Sensor_Config.Revolution_IDs;
LowSpeedReference.revolution_local_peaks = Sensor_Config.Revolution_Local_Peaks;
LowSpeedReference.revolution_pulse_times = Sensor_Config.Revolution_Pulse_Times;
LowSpeedReference.revolution_physical_peaks = Sensor_Config.Revolution_Physical_Peaks;
LowSpeedReference.revolution_best_shifts = Sensor_Config.Revolution_Best_Shifts;
LowSpeedReference.revolution_best_scores = Sensor_Config.Revolution_Best_Scores;
LowSpeedReference.revolution_score_margins = Sensor_Config.Revolution_Score_Margins;
LowSpeedReference.sensor_dominant_shifts = Sensor_Config.Sensor_Dominant_Shifts;
LowSpeedReference.sensor_local_to_physical = Sensor_Config.Sensor_Local_To_Physical;
LowSpeedReference.sensor_physical_to_local = Sensor_Config.Sensor_Physical_To_Local;
if isfield(Sensor_Config, 'Standard_Relative_Angles_OPRCenter')
    LowSpeedReference.standard_angles_opr_center = Sensor_Config.Standard_Relative_Angles_OPRCenter;
else
    LowSpeedReference.standard_angles_opr_center = Sensor_Config.Standard_Relative_Angles;
end

save(outFile, 'LowSpeedReference');

fprintf('\n=== Step01: low-speed reference fingerprint ===\n');
fprintf('Sensors: %s\n', mat2str(analysisSensors));
fprintf('Input Sensor_Config: %s\n', sensorConfigFile);
fprintf('Saved: %s\n', outFile);
disp(fingerprintRows);

if viewEnable
    visualize_fingerprint_local(LowSpeedReference, saveFigures, figureDir);
end

function validate_sensor_config_local(Sensor_Config, sensorIds, bladeCount, sensorConfigFile)
requiredFields = {'Fingerprints', 'Target_Indices', 'ReferenceSensorID', 'ReferenceRevolutionID', ...
    'Revolution_IDs', 'Revolution_Pulse_Times', 'Revolution_Physical_Peaks', ...
    'Sensor_Dominant_Shifts', 'Sensor_Physical_To_Local'};
for i = 1:numel(requiredFields)
    name = requiredFields{i};
    if ~isfield(Sensor_Config, name)
        error('Sensor_Config is missing field %s:\n  %s', name, sensorConfigFile);
    end
end

if ~isa(Sensor_Config.Fingerprints, 'containers.Map')
    error('Sensor_Config.Fingerprints must be containers.Map.');
end
if ~isa(Sensor_Config.Target_Indices, 'containers.Map')
    error('Sensor_Config.Target_Indices must be containers.Map.');
end

hasAngles = isfield(Sensor_Config, 'Standard_Relative_Angles_OPRCenter') || ...
    isfield(Sensor_Config, 'Standard_Relative_Angles');
if ~hasAngles
    error(['Sensor_Config must contain Standard_Relative_Angles_OPRCenter ' ...
        'or Standard_Relative_Angles.']);
end

for i = 1:numel(sensorIds)
    sid = sensorIds(i);
    if ~isKey(Sensor_Config.Fingerprints, sid)
        error('Sensor_Config.Fingerprints does not contain CH%d.', sid);
    end
    if ~isKey(Sensor_Config.Target_Indices, sid)
        error('Sensor_Config.Target_Indices does not contain CH%d.', sid);
    end
    fp = Sensor_Config.Fingerprints(sid);
    if numel(fp) < bladeCount
        error('Fingerprint for CH%d has only %d values, but bladeCount=%d.', ...
            sid, numel(fp), bladeCount);
    end
end
end

function T = build_fingerprint_table_local(Sensor_Config, sensorIds, bladeCount)
rows = repmat(struct('SensorID', NaN, 'ReferenceRevolutionID', NaN, 'DominantShift', NaN, ...
    'ValidRevolutionCount', NaN, 'Fingerprint', ''), numel(sensorIds), 1);
for i = 1:numel(sensorIds)
    sid = sensorIds(i);
    fp = Sensor_Config.Fingerprints(sid);
    rows(i).SensorID = sid;
    rows(i).ReferenceRevolutionID = Sensor_Config.ReferenceRevolutionID;
    rows(i).DominantShift = Sensor_Config.Sensor_Dominant_Shifts(sid);
    rows(i).ValidRevolutionCount = numel(Sensor_Config.Revolution_IDs(sid));
    rows(i).Fingerprint = mat2str(fp(1:bladeCount), 6);
end
T = struct2table(rows);
end

function visualize_fingerprint_local(LowSpeedReference, saveFigures, figureDir)
if saveFigures && exist(figureDir, 'dir') ~= 7
    mkdir(figureDir);
end

sensorIds = LowSpeedReference.analysis_sensors;
bladeIds = 1:LowSpeedReference.blade_count;
fpMat = nan(numel(sensorIds), numel(bladeIds));
for i = 1:numel(sensorIds)
    sid = sensorIds(i);
    fp = LowSpeedReference.fingerprints(sid);
    fpMat(i, :) = fp(1:numel(bladeIds));
end

rowMax = max(fpMat, [], 2, 'omitnan');
rowMax(~isfinite(rowMax) | rowMax == 0) = 1;
fpNorm = fpMat ./ rowMax;

fig = figure('Name', 'Step01 low-speed reference fingerprint', 'Color', 'w');
tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot(bladeIds, fpMat.', '-o', 'LineWidth', 1.2);
xlabel('Physical blade ID');
ylabel('Low-speed fingerprint peak');
title('Six-peak reference fingerprint');
legend(compose('CH%d', sensorIds), 'Location', 'best');
grid on; box on;

nexttile;
imagesc(bladeIds, sensorIds, fpNorm);
set(gca, 'YDir', 'normal');
xlabel('Physical blade ID');
ylabel('Sensor ID');
title('Normalized fingerprint');
colorbar;
box on;

if saveFigures
    exportgraphics(fig, fullfile(figureDir, 'Step01_LowSpeedReferenceFingerprint_20241106.png'), ...
        'Resolution', 300);
end
end
