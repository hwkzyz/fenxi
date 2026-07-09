%% Step06A: build all-blade low-speed template bank for 20241106
% This step prepares the direct/template baseline for every blade and every
% low-speed sensor used by the route. Eddy-current probes CH2/CH3 are kept
% here as direct-only template channels; only CH5/CH7 are used later by the
% gap-correction bank.

clear; clc;

thisDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(thisDir);
outDir = fullfile(thisDir, 'outputs');
P = ProjectionFlow_Config_20241106();

bankDir = fullfile(outDir, P.calibration.templateBankDirName);
if exist(bankDir, 'dir') ~= 7
    mkdir(bankDir);
end

rotTemplateLibraryDir = fullfile(rootDir, '20241106_low_speed_rotating_calibration', ...
    'output', 'new_flow', 'template_library');
rotLegacySingleDir = fullfile(rootDir, '20241106_low_speed_rotating_calibration', ...
    'output', 'new_flow', 'templates', P.calibration.templateBundleSensorTag);

bladeIds = P.calibration.buildBladeIds(:).';
sensorIds = P.calibration.lowSpeedTemplateSensors(:).';
entry = repmat(struct('bladeId', NaN, 'sensorId', NaN, 'templateFile', "", ...
    'sensorTemplate', struct(), 'sourceTemplateFile', "", 'sourceMode', ""), ...
    numel(bladeIds) * numel(sensorIds), 1);
rows = cell(numel(entry), 1);

fprintf('\n=== Step06A: low-speed template bank 20241106 ===\n');
fprintf('Single-template library dir: %s\n', rotTemplateLibraryDir);
fprintf('Legacy single-template fallback dir: %s\n', rotLegacySingleDir);
fprintf('Blades: %s | sensors: %s\n', mat2str(bladeIds), mat2str(sensorIds));

k = 0;
for ib = 1:numel(bladeIds)
    bladeId = bladeIds(ib);
    for is = 1:numel(sensorIds)
        sid = sensorIds(is);
        [sourceTemplateFile, sourceMode] = find_single_template_file_local( ...
            rotTemplateLibraryDir, rotLegacySingleDir, bladeId, sid, P.calibration.templateSuffixPreference);
        sensor = load_single_template_sensor_local(sourceTemplateFile, sid);
        sensor = normalize_sensor_voltage_units_local(sensor);
        k = k + 1;
        sensorFile = fullfile(bankDir, sprintf( ...
            'LowSpeedTemplate_20241106_B%d_S%d.mat', bladeId, sid));

        LowSpeedTemplateEntry = struct(); %#ok<NASGU>
        LowSpeedTemplateEntry.dataset = P.dataset;
        LowSpeedTemplateEntry.bladeId = bladeId;
        LowSpeedTemplateEntry.sensorId = sid;
        LowSpeedTemplateEntry.sensorTemplate = sensor;
        LowSpeedTemplateEntry.sourceTemplateFile = sourceTemplateFile;
        LowSpeedTemplateEntry.sourceMode = sourceMode;
        LowSpeedTemplateEntry.role = sensor_role_local(sid, P);
        save(sensorFile, 'LowSpeedTemplateEntry', '-v7.3');

        entry(k).bladeId = bladeId;
        entry(k).sensorId = sid;
        entry(k).templateFile = string(sensorFile);
        entry(k).sensorTemplate = sensor;
        entry(k).sourceTemplateFile = string(sourceTemplateFile);
        entry(k).sourceMode = string(sourceMode);
        rows{k} = table(string(P.dataset), bladeId, sid, string(sensorFile), ...
            string(sourceTemplateFile), string(sourceMode), string(LowSpeedTemplateEntry.role), ...
            'VariableNames', {'dataset','bladeId','sensorId','templateFile', ...
            'sourceTemplateFile','sourceMode','sensorRole'});
        fprintf('  B%d CH%d -> %s\n', bladeId, sid, sensorFile);
    end
end

LowSpeedTemplateBank = struct(); %#ok<NASGU>
LowSpeedTemplateBank.dataset = P.dataset;
LowSpeedTemplateBank.method = 'all_blade_single_sensor_low_speed_template_bank';
LowSpeedTemplateBank.description = ['Low-speed direct/template baseline for all blades and sensors. ' ...
    'Entries are assembled from single blade-sensor templates, not TemplateBundle files. ' ...
    'CH2/CH3 are direct-only sensors; CH5/CH7 are also eligible for gap correction.'];
LowSpeedTemplateBank.bladeIds = bladeIds;
LowSpeedTemplateBank.sensorIds = sensorIds;
LowSpeedTemplateBank.gapSensors = P.calibration.gapSensorIds(:).';
LowSpeedTemplateBank.directOnlySensors = P.calibration.directOnlySensorIds(:).';
LowSpeedTemplateBank.entry = entry;
LowSpeedTemplateBank.cfg = P;

indexTable = vertcat(rows{:}); %#ok<NASGU>
bankFile = fullfile(outDir, P.calibration.templateBankFile);
indexFile = fullfile(bankDir, 'LowSpeedTemplateBank_Index_20241106.csv');
save(bankFile, 'LowSpeedTemplateBank', 'indexTable', '-v7.3');
writetable(indexTable, indexFile);

fprintf('\nStep06A complete.\n');
fprintf('Bank : %s\n', bankFile);
fprintf('Index: %s\n', indexFile);
disp(indexTable);

function [templateFile, sourceMode] = find_single_template_file_local( ...
        templateLibraryDir, legacySingleDir, bladeId, sensorId, suffixPreference)
candidate = fullfile(templateLibraryDir, sprintf('Template_20241106_B%d_S%d.mat', bladeId, sensorId));
if isfile(candidate)
    templateFile = candidate;
    sourceMode = 'rotating_calibration_single_template_library';
    return;
end
patterns = {
    sprintf('Template_LowSpeedRotating_B%d_CH%d_%s_20241106.mat', bladeId, sensorId, suffixPreference)
    sprintf('Template_LowSpeedRotating_B%d_CH%d_OPRCenterStdRef_20241106.mat', bladeId, sensorId)
    sprintf('Template_LowSpeedRotating_B%d_CH%d_20241106.mat', bladeId, sensorId)
    };
for i = 1:numel(patterns)
    candidate = fullfile(legacySingleDir, patterns{i});
    if isfile(candidate)
        templateFile = candidate;
        sourceMode = 'rotating_calibration_legacy_single_template';
        return;
    end
end
error('No single template found for B%d CH%d under %s or %s.', ...
    bladeId, sensorId, templateLibraryDir, legacySingleDir);
end

function sensor = load_single_template_sensor_local(templateFile, sid)
S = load(templateFile);
if isfield(S, 'Template')
    Template = S.Template;
elseif isfield(S, 'TemplateSingle')
    Template = S.TemplateSingle;
elseif isfield(S, 'sensor_template')
    Template = S.sensor_template;
else
    error('Unsupported single template file: %s', templateFile);
end
if ~isfield(Template, 'Sensor') || isempty(Template.Sensor)
    error('Single template file has no Sensor field: %s', templateFile);
end
sensor = Template.Sensor;
if numel(sensor) > 1
    sensor = sensor([sensor.sensor_id] == sid);
end
if isempty(sensor)
    error('Single template file %s does not contain CH%d.', templateFile, sid);
end
end

function sensor = normalize_sensor_voltage_units_local(sensor)
scale = infer_voltage_scale_local(sensor);
sensor.v_grid = sensor.v_grid(:) * scale;
if isfield(sensor, 'baseline') && isfinite(sensor.baseline)
    sensor.baseline = sensor.baseline * scale;
else
    sensor.baseline = 0;
end
if isfield(sensor, 'dv_dx') && ~isempty(sensor.dv_dx)
    sensor.dv_dx = sensor.dv_dx(:) * scale;
end
if isfield(sensor, 'threshold') && isfinite(sensor.threshold)
    sensor.threshold = sensor.threshold * scale;
end
sensor.voltageUnit = 'mV';
end

function scale = infer_voltage_scale_local(sensor)
scale = 1;
if isfield(sensor, 'voltageUnit') && ~isempty(sensor.voltageUnit)
    unit = strtrim(char(string(sensor.voltageUnit)));
    if strcmpi(unit, 'v')
        scale = 1000;
    end
    return;
end
v = abs(sensor.v_grid(:));
v = v(isfinite(v));
if ~isempty(v) && median(v, 'omitnan') < 20
    scale = 1000;
end
end

function role = sensor_role_local(sid, P)
if ismember(sid, P.calibration.gapSensorIds)
    role = 'gap_capacitive';
elseif ismember(sid, P.calibration.directOnlySensorIds)
    role = 'direct_only';
else
    role = 'low_speed_template_only';
end
end
