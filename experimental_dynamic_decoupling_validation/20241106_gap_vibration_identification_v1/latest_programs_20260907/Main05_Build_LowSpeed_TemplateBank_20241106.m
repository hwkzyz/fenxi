%% 05: build the all-blade low-speed template bank for 20241106
% This step prepares the direct/template baseline for every blade and every
% low-speed sensor used by the route. Eddy-current probes CH2/CH3 are kept
% here as direct-only template channels; only CH5/CH7 are used later by the
% gap-correction bank.

clear; clc;

thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir, fullfile(thisDir, 'functions', 'gap_aware'), ...
    fullfile(thisDir, 'functions', 'utilities'));
packageCfg = Config_20241106();
rootDir = fileparts(thisDir);
outDir = packageCfg.paths.calibrationRuntime;
P = ProjectionFlow_Config_20241106();

bankDir = fullfile(outDir, P.calibration.templateBankDirName);
if exist(bankDir, 'dir') ~= 7
    mkdir(bankDir);
end

bundledTemplateFile = fullfile(packageCfg.paths.calibrationFoundation, ...
    'step04_low_speed_template', ...
    'Template_OPRCenterStd_LowSpeed_AllBlades_S2357_20241106.mat');
assert(isfile(bundledTemplateFile), ...
    'The formal template source is missing: %s', bundledTemplateFile);
loadedBundled = load(bundledTemplateFile, 'Template');
assert(isfield(loadedBundled, 'Template') && isfield(loadedBundled.Template, 'SensorBlade'), ...
    'Bundled Step04 template has no SensorBlade records: %s', bundledTemplateFile);
bundledTemplate = loadedBundled.Template;

% Do not reach into the historical rotating-calibration directory during a
% formal rebuild.  That directory remains archival; the V1 package is
% reproducible from its bundled Step04 artifact.
rotTemplateLibraryDir = '';
rotLegacySingleDir = '';

availableBladeIds = unique([bundledTemplate.SensorBlade.blade_id]);
assert(ismember(P.identification.targetBlade, availableBladeIds), ...
    'Bundled template does not contain the configured target blade B%d.', P.identification.targetBlade);
% The formal 20241106 package currently contains the validated B4 case.
% Never fabricate entries for blades absent from the bundled source.
bladeIds = P.identification.targetBlade;
sensorIds = P.calibration.lowSpeedTemplateSensors(:).';
entry = repmat(struct('bladeId', NaN, 'sensorId', NaN, 'templateFile', "", ...
    'sensorTemplate', struct(), 'sourceTemplateFile', "", 'sourceMode', ""), ...
    numel(bladeIds) * numel(sensorIds), 1);
rows = cell(numel(entry), 1);

fprintf('\n=== Step06A: low-speed template bank 20241106 ===\n');
fprintf('Bundled common-window template: %s\n', bundledTemplateFile);
fprintf('Blades: %s | sensors: %s\n', mat2str(bladeIds), mat2str(sensorIds));

k = 0;
for ib = 1:numel(bladeIds)
    bladeId = bladeIds(ib);
    for is = 1:numel(sensorIds)
        sid = sensorIds(is);
        idxTemplate = find(arrayfun(@(q) q.blade_id == bladeId && q.sensor_id == sid, ...
            bundledTemplate.SensorBlade), 1, 'first');
        assert(~isempty(idxTemplate), 'Bundled template lacks B%d CH%d.', bladeId, sid);
        sensor = bundledTemplate.SensorBlade(idxTemplate);
        sourceTemplateFile = bundledTemplateFile;
        sourceMode = 'bundled_step04_common_low_speed_window';
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
    'All sensors are extracted from one bundled Step04 template built from a common low-speed lap window. ' ...
    'CH2/CH3 are direct-only sensors; CH5/CH7 are also eligible for gap correction.'];
LowSpeedTemplateBank.commonLowSpeedLapRange = packageCfg.calibration.commonLowSpeedLapRange;
LowSpeedTemplateBank.commonLowSpeedWindowRequired = packageCfg.calibration.requireCommonLowSpeedWindow;
LowSpeedTemplateBank.sourceTemplateFile = bundledTemplateFile;
LowSpeedTemplateBank.sourceTemplateMetadata = get_template_metadata_local(bundledTemplate);
LowSpeedTemplateBank.formalScope = 'configured_target_blade_only';
LowSpeedTemplateBank.requestedBuildBladeIds = P.calibration.buildBladeIds(:).';
LowSpeedTemplateBank.availableBundledBladeIds = availableBladeIds(:).';
LowSpeedTemplateBank.bladeIds = bladeIds;
LowSpeedTemplateBank.sensorIds = sensorIds;
LowSpeedTemplateBank.gapSensors = P.calibration.gapSensorIds(:).';
LowSpeedTemplateBank.directOnlySensors = P.calibration.directOnlySensorIds(:).';
LowSpeedTemplateBank.entry = entry;
LowSpeedTemplateBank.cfg = P;

indexTable = vertcat(rows{:}); %#ok<NASGU>
bankFile = fullfile(outDir, P.calibration.runtimeTemplateBankFile);
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

function metadata = get_template_metadata_local(Template)
metadata = struct();
names = {'Dataset','Case_Name','CreatedBy','CreatedOn', ...
    'LowSpeed_Template_Source_Mode','OPR_Timing_Method', ...
    'OPR_Events_Per_Revolution','Standard_Angle_Source'};
for i = 1:numel(names)
    if isfield(Template, names{i})
        metadata.(names{i}) = Template.(names{i});
    end
end
if isfield(Template, 'Summary_Table')
    metadata.summaryTable = Template.Summary_Table;
end
end
