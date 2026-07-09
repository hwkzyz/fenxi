%% Step06I: calibrate 20241106 sensor gap library from unified low-speed template
% Main mode now reuses the 20241106 rotating-calibration low-speed template.
% The old "rebuild from 900 rpm raw pulses" path is retained as a fallback
% for debugging and comparison.
%
% Run after:
%   Step05_Build_Response_Surface_20241106
%   Step05I_Learn_OffsetTilt_Shared_Response_Surface_20241106
% Optional:
%   Step06_Build_HighMap_20241106

clear; clc; close all;

%% 1. Paths and simple settings
thisDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(thisDir);
outDir = fullfile(thisDir, 'outputs');
figDir = fullfile(outDir, 'figures_step06i_20241106_gap_library');
C0 = CaseConfig();
if ~exist(figDir, 'dir')
    mkdir(figDir);
end

rotDir = fullfile(rootDir, '20241106_low_speed_rotating_calibration');

step02File = fullfile(outDir, 'Step02_BTT_Displacement_20241106.mat');
responseFile = fullfile(outDir, 'Step05I_OffsetTilt_Shared_Response_Surface_20241106.mat');
highMapFile = fullfile(outDir, 'Step06_HighMap_20241106.mat');

if ~isfile(responseFile)
    error('Run Step05I first. Missing file: %s', responseFile);
end

if ~isfile(step02File)
    warning('Legacy Step02 file is missing. Legacy template fallback mode will be unavailable: %s', step02File);
end

Srf = load(responseFile, 'OffsetTiltResponseSurface');
responseSurface = Srf.OffsetTiltResponseSurface;

if isfile(highMapFile)
    H = load(highMapFile, 'highMap');
    highMap = H.highMap;
else
    highMap = struct();
end
defaultSensors = C0.sensorIds;
defaultBlade = C0.bladeId;

cfg = struct();
cfg.templateSourceMode = 'rotating_step05_template'; % rotating_step05_template | legacy_generate_from_900rpm
cfg.targetBlade = defaultBlade;
cfg.analysisSensors = defaultSensors;
cfg.rotTemplateFile = '';
cfg.rotTemplateSuffixPreference = 'GradientXRange030_OPRCenterStd';
cfg.lowPulseCountPerSensor = 40;
cfg.lowPulseHalfWindowSec = 1.0e-3;
cfg.pointsPerPulse = 151;
cfg.templateXCount = 201;
cfg.edgeFractionForBaseline = 0.18;
cfg.minValidPulsePoints = 45;
cfg.blockStepSamples = 1e7;
cfg.blockLabelStep = 1000;
cfg.rawCacheMaxEntries = 8;
cfg.gBoundsMm = [min(responseSurface.gTrainMm(:)), max(responseSurface.gTrainMm(:))];
cfg.tauBoundsMm = [-0.80, 0.80];
cfg.xScaleBounds = [0.65, 1.45];
cfg.muBounds = [-0.35, 0.35];
cfg.overshootPenaltyMv = 800;
cfg.fminOptions = optimset('Display', 'off', 'MaxIter', 450, ...
    'MaxFunEvals', 1600, 'TolX', 1e-6, 'TolFun', 1e-6);

cfg.analysisSensors = parse_int_env_local('STEP06I_ANALYSIS_SENSORS', cfg.analysisSensors);
assert_capacitive_gap_sensors_local(cfg.analysisSensors, 'STEP06I_ALLOW_NONCAP_GAP_SENSORS');
tmpBlade = parse_int_env_local('STEP06I_TARGET_BLADE', cfg.targetBlade);
cfg.targetBlade = tmpBlade(1);
cfg.lowPulseCountPerSensor = round(parse_scalar_env_local( ...
    'STEP06I_LOW_PULSE_COUNT', cfg.lowPulseCountPerSensor));
cfg.templateSourceMode = parse_mode_env_local('STEP06I_TEMPLATE_SOURCE_MODE', cfg.templateSourceMode, ...
    {'rotating_step05_template', 'legacy_generate_from_900rpm'});
cfg.rotTemplateFile = strtrim(getenv('STEP06I_ROTATING_TEMPLATE_FILE'));

sensorTag = ['S', sprintf('%d', cfg.analysisSensors)];

[Template, templateMeta] = load_or_build_low_speed_template_local( ...
    cfg, rotDir, step02File, responseSurface, highMap);

fprintf('\n=== Step06I 20241106 offset-tilt shared gap library calibration ===\n');
fprintf('Response surface: %s\n', responseFile);
fprintf('Template source mode: %s\n', cfg.templateSourceMode);
fprintf('Template source file: %s\n', templateMeta.templateFile);
fprintf('Sensors: %s, blade slot: %d\n', mat2str(cfg.analysisSensors), cfg.targetBlade);

%% 2. Fit per-sensor corrected library parameters
sensorCorr = repmat(struct('sensorId', NaN, 'g0Mm', NaN, 'tauMm', NaN, ...
    'xScale', NaN, 'muGapPerXMm', NaN, 'tiltAngleDeg', NaN, ...
    'voltageGain', NaN, 'voltageOffsetMv', NaN, 'lowFitRmseMv', NaN, ...
    'overshootRmseMm', NaN, 'pointCount', NaN, 'x', [], 'vLowMv', [], ...
    'vFitMv', [], 'xLib', [], 'gEff', [], 'residualMv', [], ...
    'etaMedianMm', NaN, 'etaIqrMm', NaN, 'etaLimitMm', NaN), ...
    numel(cfg.analysisSensors), 1);
rows = cell(numel(cfg.analysisSensors), 1);

for is = 1:numel(cfg.analysisSensors)
    sid = cfg.analysisSensors(is);
    Tsen = get_template_sensor_local(Template, sid);
    validate_template_sensor_fields_local(Tsen, sid);
    fit = fit_low_template_local(Tsen, responseSurface, cfg);
    sensorCorr(is) = fit;
    rows{is} = table(sid, fit.g0Mm, fit.tauMm, fit.xScale, fit.muGapPerXMm, ...
        fit.tiltAngleDeg, fit.voltageGain, fit.voltageOffsetMv, fit.lowFitRmseMv, ...
        fit.overshootRmseMm, fit.pointCount, fit.etaMedianMm, fit.etaIqrMm, fit.etaLimitMm, ...
        'VariableNames', {'sensorId','g0Mm','tauMm','xScale','muGapPerXMm', ...
        'tiltAngleDeg','voltageGain','voltageOffsetMv','lowFitRmseMv', ...
        'overshootRmseMm','pointCount','etaMedianMm','etaIqrMm','etaLimitMm'});

    pulseCount = get_optional_field_local(Tsen, 'pulseCount', NaN);
    fprintf('  P%d: g0 %.4f mm, tau %.4f mm, k %.4f, mu %.4f (%.2f deg), RMSE %.2f mV, pulses %.0f\n', ...
        sid, fit.g0Mm, fit.tauMm, fit.xScale, fit.muGapPerXMm, ...
        fit.tiltAngleDeg, fit.lowFitRmseMv, pulseCount);
    fprintf('       eta median %.4f mm, IQR %.4f mm, limit %.4f mm\n', ...
        fit.etaMedianMm, fit.etaIqrMm, fit.etaLimitMm);
end

calibrationTable = vertcat(rows{:});
highCheckTable = build_highmap_check_local(highMap, sensorCorr, responseSurface, cfg);

CorrectedGapLibrary = struct();
CorrectedGapLibrary.dataset = '20241106';
CorrectedGapLibrary.method = 'unified_low_speed_template_offset_tilt_shared_gap_library';
CorrectedGapLibrary.description = ['The 20241106 rotating_calibration low-speed template is the default baseline. ' ...
    'The gap library supplies only the clearance-change increment needed to reach the current high-speed static condition. ' ...
    'The old 900 rpm template rebuild path is retained only as a fallback/debug mode.'];
CorrectedGapLibrary.responseFile = responseFile;
CorrectedGapLibrary.templateFile = '';
CorrectedGapLibrary.templateSourceMode = cfg.templateSourceMode;
CorrectedGapLibrary.templateSourceFile = string(templateMeta.templateFile);
CorrectedGapLibrary.targetBlade = cfg.targetBlade;
CorrectedGapLibrary.analysisSensors = cfg.analysisSensors(:).';
CorrectedGapLibrary.responseSurface = responseSurface;
CorrectedGapLibrary.lowSpeedTemplate = Template;
CorrectedGapLibrary.sensor = sensorCorr;
CorrectedGapLibrary.formula = 'T_low(x)+a*(F_raw(g+mu*(x-tau),k*(x-tau))-F_raw(g0+mu*(x-tau),k*(x-tau)))';
CorrectedGapLibrary.cfg = cfg;

matFile = fullfile(outDir, sprintf('Step06I_OffsetTiltShared_GapLibrary_20241106_B%d_%s.mat', ...
    cfg.targetBlade, sensorTag));
legacyTemplateFile = fullfile(outDir, sprintf('Step06I_Generated_LowSpeed_Template_20241106_B%d_%s.mat', ...
    cfg.targetBlade, sensorTag));
csvFile = fullfile(outDir, sprintf('Step06I_OffsetTiltShared_GapLibrary_Calibration_20241106_B%d_%s.csv', ...
    cfg.targetBlade, sensorTag));
highCsvFile = fullfile(outDir, sprintf('Step06I_OffsetTiltShared_GapLibrary_HighMap_Check_20241106_B%d_%s.csv', ...
    cfg.targetBlade, sensorTag));
figLow = fullfile(figDir, sprintf('Step06I_LowTemplate_Fit_20241106_B%d_%s.png', ...
    cfg.targetBlade, sensorTag));
figPath = fullfile(figDir, sprintf('Step06I_TiltPath_20241106_B%d_%s.png', ...
    cfg.targetBlade, sensorTag));

if strcmpi(cfg.templateSourceMode, 'legacy_generate_from_900rpm')
    CorrectedGapLibrary.templateFile = legacyTemplateFile;
    save(legacyTemplateFile, 'Template', '-v7.3');
else
    CorrectedGapLibrary.templateFile = string(templateMeta.templateFile);
end
save(matFile, 'CorrectedGapLibrary', 'Template', 'calibrationTable', ...
    'highCheckTable', 'responseFile', 'highMapFile', 'templateMeta', '-v7.3');
writetable(calibrationTable, csvFile);
if ~isempty(highCheckTable)
    writetable(highCheckTable, highCsvFile);
end

plot_low_fit_local(sensorCorr, Template, figLow);
plot_tilt_path_local(sensorCorr, responseSurface, figPath);

fprintf('\nStep06I complete.\n');
disp(calibrationTable);
if ~isempty(highCheckTable)
    disp(highCheckTable);
end
fprintf('Saved:\n  %s\n  %s\n', matFile, csvFile);
if strcmpi(cfg.templateSourceMode, 'legacy_generate_from_900rpm')
    fprintf('Legacy template snapshot:\n  %s\n', legacyTemplateFile);
end

%% Local functions
function [Template, meta] = load_or_build_low_speed_template_local(cfg, rotDir, step02File, responseSurface, highMap)
meta = struct('templateFile', '', 'mode', cfg.templateSourceMode);
switch lower(strtrim(cfg.templateSourceMode))
    case 'rotating_step05_template'
        [Template, meta.templateFile] = load_rotating_step05_template_local(cfg, rotDir);
        meta.mode = 'rotating_step05_template';
    case 'legacy_generate_from_900rpm'
        if ~isfile(step02File)
            error('Legacy template fallback requires Step02 file: %s', step02File);
        end
        Template = build_legacy_low_speed_template_local(cfg, step02File, responseSurface, highMap);
        meta.templateFile = '[legacy_generated_from_900rpm_raw_probe_pulses]';
        meta.mode = 'legacy_generate_from_900rpm';
    otherwise
        error('Unsupported templateSourceMode: %s', cfg.templateSourceMode);
end
Template.analysisSensors = cfg.analysisSensors(:).';
Template.targetBlade = cfg.targetBlade;
end

function [Template, templateFile] = load_rotating_step05_template_local(cfg, rotDir)
if ~isempty(cfg.rotTemplateFile)
    templateFile = cfg.rotTemplateFile;
    if ~isfile(templateFile)
        error('STEP06I_ROTATING_TEMPLATE_FILE does not exist: %s', templateFile);
    end
else
    templateFile = find_rotating_template_file_local(rotDir, cfg.targetBlade, cfg.analysisSensors, cfg.rotTemplateSuffixPreference);
end
Template = load_rotating_template_from_source_local(templateFile, cfg);
Template = enrich_rotating_template_eta_local(Template, cfg, rotDir);
end

function templateFile = find_rotating_template_file_local(rotDir, targetBlade, sensorIds, suffixPreference)
templateDirLegacy = fullfile(rotDir, 'output', 'templates');
sensorSubdirLegacy = fullfile(templateDirLegacy, ['S', sprintf('%d', sensorIds)]);
templateDirNewFlow = fullfile(rotDir, 'output', 'new_flow', 'templates');
sensorSubdirNewFlow = fullfile(templateDirNewFlow, 'S2357');
sensorTag = ['S', sprintf('%d', sensorIds)];
bundlePatterns = {
    sprintf('TemplateBundle_LowSpeedRotating_B%d_S*_%s_20241106.mat', targetBlade, suffixPreference)
    sprintf('TemplateBundle_LowSpeedRotating_B%d_S*_OPRCenterStdRef_20241106.mat', targetBlade)
    sprintf('TemplateBundle_LowSpeedRotating_B%d_S*_20241106.mat', targetBlade)
    };
searchRoots = {sensorSubdirNewFlow, templateDirNewFlow};
for ir = 1:numel(searchRoots)
    if exist(searchRoots{ir}, 'dir') ~= 7
        continue;
    end
    for ip = 1:numel(bundlePatterns)
        files = dir(fullfile(searchRoots{ir}, bundlePatterns{ip}));
        for jf = 1:numel(files)
            candidate = fullfile(files(jf).folder, files(jf).name);
            if bundle_index_contains_sensors_local(candidate, sensorIds)
                templateFile = candidate;
                return;
            end
        end
    end
end
patterns = {
    sprintf('Template_LowSpeedRotating_B%d_%s_%s_20241106.mat', targetBlade, sensorTag, suffixPreference)
    sprintf('Template_LowSpeedRotating_B%d_%s_OPRCenterStdRef_20241106.mat', targetBlade, sensorTag)
    sprintf('Template_LowSpeedRotating_B%d_%s_GradientXRange030_OPRAnchored_20241106.mat', targetBlade, sensorTag)
    sprintf('Template_LowSpeedRotating_B%d_%s_GradientXRange030_20241106.mat', targetBlade, sensorTag)
    sprintf('Template_LowSpeedRotating_B%d_%s_20241106.mat', targetBlade, sensorTag)
    };
for i = 1:numel(patterns)
    searchRoots = {templateDirNewFlow, sensorSubdirNewFlow, templateDirLegacy, sensorSubdirLegacy};
    for ir = 1:numel(searchRoots)
        if exist(searchRoots{ir}, 'dir') ~= 7
            continue;
        end
        files = dir(fullfile(searchRoots{ir}, patterns{i}));
        if ~isempty(files)
            templateFile = fullfile(files(1).folder, files(1).name);
            return;
        end
    end
end
error('No rotating Step05 template file found in %s, %s, %s, or %s for blade %d sensors %s.', ...
    templateDirNewFlow, sensorSubdirNewFlow, templateDirLegacy, sensorSubdirLegacy, targetBlade, mat2str(sensorIds));
end

function tf = bundle_index_contains_sensors_local(bundleFile, sensorIds)
tf = false;
S = load(bundleFile, 'TemplateBundle');
if ~isfield(S, 'TemplateBundle') || ~isfield(S.TemplateBundle, 'SensorFiles') || isempty(S.TemplateBundle.SensorFiles)
    return;
end
bundleSensorIds = S.TemplateBundle.SensorFiles.sensor_id(:);
tf = all(ismember(sensorIds(:), bundleSensorIds));
end

function Template = load_rotating_template_from_source_local(templateFile, cfg)
S = load(templateFile);
if isfield(S, 'Template')
    Template = normalize_rotating_template_local(S.Template, cfg);
    return;
end
if isfield(S, 'sensor_template')
    Template = normalize_rotating_template_local(S.sensor_template, cfg);
    return;
end
if isfield(S, 'TemplateBundle')
    Template = load_rotating_template_from_bundle_index_local(S.TemplateBundle, cfg);
    return;
end
error('Unsupported rotating template source file: %s', templateFile);
end

function Template = load_rotating_template_from_bundle_index_local(bundle, cfg)
if ~isfield(bundle, 'SensorFiles') || isempty(bundle.SensorFiles)
    error('TemplateBundle is missing SensorFiles index.');
end
sensorArray = repmat(struct(), 1, 0);
templateMeta = [];
for sid = cfg.analysisSensors(:).'
    idx = find(bundle.SensorFiles.sensor_id == sid, 1, 'first');
    if isempty(idx)
        error('TemplateBundle does not contain requested sensor %d.', sid);
    end
    sourceFile = char(bundle.SensorFiles.template_file(idx));
    S = load(sourceFile);
    if isfield(S, 'sensor_template')
        sensorTemplate = S.sensor_template;
    elseif isfield(S, 'Template')
        sensorTemplate = S.Template;
    else
        error('Sensor template file is missing supported root variable: %s', sourceFile);
    end
    sensorTemplate = normalize_rotating_template_local(sensorTemplate, struct('analysisSensors', sid));
    if isempty(templateMeta)
        templateMeta = sensorTemplate;
    end
    sensorArray = append_sensor_struct_local(sensorArray, sensorTemplate.Sensor(1));
end
Template = templateMeta;
Template.Sensor = sensorArray;
Template.SensorIDs = cfg.analysisSensors(:).';
Template.SensorTag = ['S', sprintf('%d', cfg.analysisSensors)];
if ~isfield(Template, 'xGrid') || isempty(Template.xGrid)
    Template.xGrid = Template.Sensor(1).x_grid(:);
end
end

function sensorArray = append_sensor_struct_local(sensorArray, sensorNow)
if isempty(sensorArray)
    sensorArray = sensorNow;
    return;
end
allFields = unique([fieldnames(sensorArray); fieldnames(sensorNow)], 'stable');
for i = 1:numel(sensorArray)
    sensorArray(i) = ensure_struct_fields_local(sensorArray(i), allFields);
end
sensorNow = ensure_struct_fields_local(sensorNow, allFields);
sensorArray(end + 1) = sensorNow;
end

function s = ensure_struct_fields_local(s, fields)
for i = 1:numel(fields)
    name = fields{i};
    if ~isfield(s, name)
        s.(name) = [];
    end
end
s = orderfields(s, fields);
end

function Template = enrich_rotating_template_eta_local(Template, cfg, rotDir)
[summaryFile, auditFile] = find_rotating_eta_artifacts_local(rotDir);
summaryTable = [];
if ~isempty(summaryFile) && isfile(summaryFile)
    summaryTable = readtable(summaryFile);
end
auditStruct = [];
if ~isempty(auditFile) && isfile(auditFile)
    Sa = load(auditFile, 'LowSpeedTemplateAudit');
    if isfield(Sa, 'LowSpeedTemplateAudit')
        auditStruct = Sa.LowSpeedTemplateAudit;
    end
end
for is = 1:numel(Template.Sensor)
    sid = Template.Sensor(is).sensor_id;
    Template.Sensor(is) = backfill_sensor_eta_local(Template.Sensor(is), cfg.targetBlade, sid, summaryTable, auditStruct);
end
end

function [summaryFile, auditFile] = find_rotating_eta_artifacts_local(rotDir)
summaryFile = '';
auditFile = '';
summaryHits = dir(fullfile(rotDir, 'output', 'new_flow', '05A_low_speed_template_library', '**', 'LowSpeedTemplateLibrary_20241106.csv'));
if isempty(summaryHits)
    summaryHits = dir(fullfile(rotDir, 'output', '05A_low_speed_template_library', '**', 'LowSpeedTemplateLibrary_20241106.csv'));
end
if ~isempty(summaryHits)
    [~, idx] = max([summaryHits.datenum]);
    summaryFile = fullfile(summaryHits(idx).folder, summaryHits(idx).name);
end
auditHits = dir(fullfile(rotDir, 'output', 'new_flow', '05A_low_speed_template_library', '**', 'LowSpeedTemplateAudit_20241106.mat'));
if isempty(auditHits)
    auditHits = dir(fullfile(rotDir, 'output', '05A_low_speed_template_library', '**', 'LowSpeedTemplateAudit_20241106.mat'));
end
if ~isempty(auditHits)
    [~, idx] = max([auditHits.datenum]);
    auditFile = fullfile(auditHits(idx).folder, auditHits(idx).name);
end
end

function sensorStruct = backfill_sensor_eta_local(sensorStruct, targetBlade, sid, summaryTable, auditStruct)
if ~isfield(sensorStruct, 'eta_app_mm') || isempty(sensorStruct.eta_app_mm)
    sensorStruct.eta_app_mm = [];
end
if ~isfield(sensorStruct, 'eta_median_mm') || ~isfinite(sensorStruct.eta_median_mm)
    sensorStruct.eta_median_mm = NaN;
end
if ~isfield(sensorStruct, 'eta_iqr_mm') || ~isfinite(sensorStruct.eta_iqr_mm)
    sensorStruct.eta_iqr_mm = NaN;
end
if ~isfield(sensorStruct, 'eta_limit_mm') || ~isfinite(sensorStruct.eta_limit_mm)
    sensorStruct.eta_limit_mm = NaN;
end
if ~isempty(summaryTable) && all(ismember({'BladeID','SensorID'}, summaryTable.Properties.VariableNames))
    mask = summaryTable.BladeID == targetBlade & summaryTable.SensorID == sid;
    if any(mask)
        row = summaryTable(find(mask, 1, 'first'), :);
        if ismember('EtaMedianMM', summaryTable.Properties.VariableNames) && ~isfinite(sensorStruct.eta_median_mm)
            sensorStruct.eta_median_mm = row.EtaMedianMM;
        end
        if ismember('EtaIQRMM', summaryTable.Properties.VariableNames) && ~isfinite(sensorStruct.eta_iqr_mm)
            sensorStruct.eta_iqr_mm = row.EtaIQRMM;
        end
        if ismember('EtaLimitMM', summaryTable.Properties.VariableNames) && ~isfinite(sensorStruct.eta_limit_mm)
            sensorStruct.eta_limit_mm = row.EtaLimitMM;
        end
    end
end
if isempty(sensorStruct.eta_app_mm) && ~isempty(auditStruct) && isfield(auditStruct, 'blade')
    for ib = 1:numel(auditStruct.blade)
        bladeNow = auditStruct.blade(ib);
        if ~isfield(bladeNow, 'blade_id') || bladeNow.blade_id ~= targetBlade || ~isfield(bladeNow, 'sensor')
            continue;
        end
        for js = 1:numel(bladeNow.sensor)
            sensorNow = bladeNow.sensor(js);
            if isfield(sensorNow, 'sensor_id') && sensorNow.sensor_id == sid
                if isfield(sensorNow, 'eta_app_mm') && ~isempty(sensorNow.eta_app_mm)
                    sensorStruct.eta_app_mm = sensorNow.eta_app_mm(:);
                end
                if ~isfinite(sensorStruct.eta_median_mm) && isfield(sensorNow, 'eta_median_mm')
                    sensorStruct.eta_median_mm = sensorNow.eta_median_mm;
                end
                if ~isfinite(sensorStruct.eta_iqr_mm) && isfield(sensorNow, 'eta_iqr_mm')
                    sensorStruct.eta_iqr_mm = sensorNow.eta_iqr_mm;
                end
                if ~isfinite(sensorStruct.eta_limit_mm) && isfield(sensorNow, 'eta_limit_mm')
                    sensorStruct.eta_limit_mm = sensorNow.eta_limit_mm;
                end
                return;
            end
        end
    end
end
end

function Template = normalize_rotating_template_local(rotTemplate, cfg)
if isfield(rotTemplate, 'Sensor')
    Template = rotTemplate;
else
    error('Unsupported rotating template structure: missing Sensor field.');
end
Template.dataset = '20241106';
Template.method = 'adapted_from_rotating_step05_template';
if ~isfield(Template, 'xGrid') || isempty(Template.xGrid)
    idx = find([Template.Sensor.sensor_id] == cfg.analysisSensors(1), 1, 'first');
    if isempty(idx)
        error('Rotating template does not contain requested sensor %d.', cfg.analysisSensors(1));
    end
    Template.xGrid = Template.Sensor(idx).x_grid(:);
end
Template.Sensor = filter_template_sensors_local(Template.Sensor, cfg.analysisSensors);
for is = 1:numel(Template.Sensor)
    if ~isfield(Template.Sensor(is), 'baseline') || ~isfinite(Template.Sensor(is).baseline)
        Template.Sensor(is).baseline = 0;
    end
    if ~isfield(Template.Sensor(is), 'xc') || ~isfinite(Template.Sensor(is).xc)
        Template.Sensor(is).xc = 0;
    end
    if ~isfield(Template.Sensor(is), 'pulseCount')
        Template.Sensor(is).pulseCount = get_optional_field_local(Template.Sensor(is), 'lap_count', NaN);
    end
    if ~isfield(Template.Sensor(is), 'pulsePeakMv')
        Template.Sensor(is).pulsePeakMv = [];
    end
    if ~isfield(Template.Sensor(is), 'eventTime')
        Template.Sensor(is).eventTime = [];
    end
end
Template = normalize_template_voltage_units_local(Template);
end

function sensors = filter_template_sensors_local(sensorArray, analysisSensors)
sensors = repmat(sensorArray(1), 1, 0);
for sid = analysisSensors(:).'
    idx = find([sensorArray.sensor_id] == sid, 1, 'first');
    if isempty(idx)
        error('Rotating template is missing requested sensor %d.', sid);
    end
    sensors(end + 1) = sensorArray(idx); %#ok<AGROW>
end
end

function Template = build_legacy_low_speed_template_local(cfg, step02File, responseSurface, highMap)
S02 = load(step02File);
if isstruct(highMap) && isfield(highMap, 'analysisSensors') && ~isempty(highMap.analysisSensors)
    defaultSensors = highMap.analysisSensors;
else
    defaultSensors = intersect([2 3 5 7], S02.sensorIds, 'stable');
end
if isempty(cfg.analysisSensors)
    cfg.analysisSensors = defaultSensors;
end

xGrid = responseSurface.xGrid(:);
if isfield(responseSurface, 'effectiveWindow') && any(responseSurface.effectiveWindow)
    xMin = min(xGrid(responseSurface.effectiveWindow));
    xMax = max(xGrid(responseSurface.effectiveWindow));
else
    xMin = min(xGrid);
    xMax = max(xGrid);
end
templateXGrid = linspace(xMin, xMax, cfg.templateXCount).';

lowOpr = load_opr_local(S02.lowDir, S02.sampleRateHz);
rawCache = containers.Map('KeyType', 'char', 'ValueType', 'any');

Template = struct();
Template.dataset = '20241106';
Template.method = 'generated_from_900rpm_raw_probe_pulses';
Template.lowDir = S02.lowDir;
Template.targetBlade = cfg.targetBlade;
Template.analysisSensors = cfg.analysisSensors;
Template.xGrid = templateXGrid;
Template.Sensor = repmat(struct('sensor_id', NaN, 'x_grid', [], 'v_grid', [], ...
    'baseline', 0, 'xc', 0, 'pulseCount', 0, 'pulsePeakMv', [], 'eventTime', [], ...
    'eta_app_mm', [], 'eta_median_mm', NaN, 'eta_iqr_mm', NaN, 'eta_limit_mm', NaN), ...
    1, numel(cfg.analysisSensors));

for is = 1:numel(cfg.analysisSensors)
    sid = cfg.analysisSensors(is);
    sensorIdx = find([S02.lowCheck.sensorId] == sid, 1, 'first');
    if isempty(sensorIdx)
        error('No low-speed Step02 data for sensor %d.', sid);
    end
    slot = S02.lowCheck(sensorIdx).slot(cfg.targetBlade);
    eventTimesAll = slot.time(:);
    eventTimesAll = eventTimesAll(isfinite(eventTimesAll));
    pick = pick_evenly_local(numel(eventTimesAll), cfg.lowPulseCountPerSensor);
    eventTimes = eventTimesAll(pick);

    pulseMatrix = NaN(numel(templateXGrid), numel(eventTimes));
    peakMv = NaN(numel(eventTimes), 1);
    usedEventTimes = NaN(numel(eventTimes), 1);

    for ie = 1:numel(eventTimes)
        tEvent = eventTimes(ie);
        rotFreqHz = local_rot_freq_local(lowOpr.tCenter, tEvent);
        tipSpeedMmS = 2 * pi * S02.rTipMm * rotFreqHz;
        [tq, vq, ok] = read_raw_pulse_local(S02.lowDir, sid, tEvent, ...
            cfg.lowPulseHalfWindowSec, cfg.pointsPerPulse, S02.sampleRateHz, ...
            cfg.blockStepSamples, cfg.blockLabelStep, rawCache, cfg.rawCacheMaxEntries);
        if ~ok
            continue;
        end

        xPulse = (tq - tEvent) * tipSpeedMmS;
        edgeCount = max(5, round(cfg.edgeFractionForBaseline * numel(vq)));
        baselineV = median([vq(1:edgeCount); vq(end-edgeCount+1:end)], 'omitnan');
        vCorrMv = 1000 * (vq - baselineV);
        keep = isfinite(xPulse) & isfinite(vCorrMv) & ...
            xPulse >= min(templateXGrid) & xPulse <= max(templateXGrid);
        if nnz(keep) < cfg.minValidPulsePoints
            continue;
        end

        pulseMatrix(:, ie) = interp1(xPulse(keep), vCorrMv(keep), templateXGrid, 'linear', NaN);
        peakMv(ie) = max(vCorrMv(keep), [], 'omitnan');
        usedEventTimes(ie) = tEvent;
    end

    good = sum(isfinite(pulseMatrix), 1) >= 0.6 * numel(templateXGrid);
    if nnz(good) < 5
        error('Too few usable low-speed pulses for sensor %d.', sid);
    end

    vTemplate = median(pulseMatrix(:, good), 2, 'omitnan');
    Template.Sensor(is).sensor_id = sid;
    Template.Sensor(is).x_grid = templateXGrid;
    Template.Sensor(is).v_grid = vTemplate;
    Template.Sensor(is).baseline = 0;
    Template.Sensor(is).xc = 0;
    Template.Sensor(is).pulseCount = nnz(good);
    Template.Sensor(is).pulsePeakMv = peakMv(good);
    Template.Sensor(is).eventTime = usedEventTimes(good);
end
Template = normalize_template_voltage_units_local(Template);
end

function validate_template_sensor_fields_local(Tsen, sid)
requiredFields = {'x_grid', 'v_grid', 'baseline', 'xc'};
for i = 1:numel(requiredFields)
    name = requiredFields{i};
    if ~isfield(Tsen, name)
        error('Template sensor CH%d is missing required field %s.', sid, name);
    end
end
if isempty(Tsen.x_grid) || isempty(Tsen.v_grid)
    error('Template sensor CH%d has empty x_grid/v_grid.', sid);
end
end

function values = parse_int_env_local(name, defaultValues)
raw = strtrim(getenv(name));
if isempty(raw)
    values = defaultValues;
    return;
end
values = sscanf(raw, '%d').';
if isempty(values)
    values = defaultValues;
end
end

function value = parse_scalar_env_local(name, defaultValue)
raw = strtrim(getenv(name));
if isempty(raw)
    value = defaultValue;
    return;
end
tmp = str2double(raw);
if isfinite(tmp)
    value = tmp;
else
    value = defaultValue;
end
end

function mode = parse_mode_env_local(name, defaultMode, allowedModes)
raw = lower(strtrim(getenv(name)));
if isempty(raw)
    mode = defaultMode;
    return;
end
if ~ismember(raw, allowedModes)
    error('%s must be one of: %s', name, strjoin(allowedModes, ', '));
end
mode = raw;
end

function value = get_optional_field_local(S, fieldName, defaultValue)
if isstruct(S) && isfield(S, fieldName)
    value = S.(fieldName);
else
    value = defaultValue;
end
end

function opr = load_opr_local(caseDir, sampleRateHz)
f = fullfile(caseDir, 'jiluOPR.mat');
if ~isfile(f)
    error('Missing OPR file: %s', f);
end
S = load(f, 'jiluOPR');
raw = S.jiluOPR;
opr.raw = raw;
opr.tCenter = mean(raw(:, 1:2), 2) / sampleRateHz;
end

function idx = pick_evenly_local(n, maxCount)
maxCount = max(1, round(maxCount));
if n <= maxCount
    idx = (1:n).';
else
    idx = unique(round(linspace(1, n, maxCount))).';
end
end

function f = local_rot_freq_local(oprTime, t)
idx = find(oprTime <= t, 1, 'last');
if isempty(idx)
    idx = 1;
end
idx = min(idx, numel(oprTime) - 1);
dt = oprTime(idx + 1) - oprTime(idx);
if ~(isfinite(dt) && dt > 0)
    dt = median(diff(oprTime), 'omitnan');
end
f = 1 / dt;
end

function [tq, vq, ok] = read_raw_pulse_local(caseDir, sensorId, tCenter, halfWindowSec, ...
    pointsPerPulse, sampleRateHz, blockStepSamples, blockLabelStep, rawCache, rawCacheMaxEntries)
tq = linspace(tCenter - halfWindowSec, tCenter + halfWindowSec, pointsPerPulse).';
sampleStart = floor((tCenter - halfWindowSec) * sampleRateHz);
sampleEnd = ceil((tCenter + halfWindowSec) * sampleRateHz);

blockStart = sample_to_block_label_local(sampleStart, blockStepSamples, blockLabelStep);
blockEnd = sample_to_block_label_local(sampleEnd, blockStepSamples, blockLabelStep);
blockLabels = blockStart:blockLabelStep:blockEnd;

rawAll = zeros(0, 2);
for blockLabel = blockLabels
    raw = load_raw_block_cached_local(caseDir, sensorId, blockLabel, rawCache, rawCacheMaxEntries);
    if ~isempty(raw)
        rawAll = [rawAll; raw]; %#ok<AGROW>
    end
end

if isempty(rawAll)
    vq = NaN(size(tq));
    ok = false;
    return;
end

rawAll(rawAll(:, 1) == 0, :) = [];
rawAll = sortrows(rawAll, 1);
rawT = rawAll(:, 1) / sampleRateHz;
rawV = rawAll(:, 2);
keep = rawT >= tq(1) & rawT <= tq(end) & isfinite(rawV);
if nnz(keep) < 4
    vq = NaN(size(tq));
    ok = false;
    return;
end

vq = interp1(rawT(keep), rawV(keep), tq, 'linear', NaN);
ok = nnz(isfinite(vq)) >= 0.8 * numel(vq);
end

function blockLabel = sample_to_block_label_local(sample, blockStepSamples, blockLabelStep)
blockIndex = ceil((sample - 10000) / blockStepSamples);
blockIndex = max(1, blockIndex);
blockLabel = blockIndex * blockLabelStep;
end

function raw = load_raw_block_cached_local(caseDir, sensorId, blockLabel, rawCache, rawCacheMaxEntries)
key = sprintf('S%d_B%d', sensorId, blockLabel);
if isKey(rawCache, key)
    raw = rawCache(key);
    return;
end

filePath = fullfile(caseDir, sprintf('4-%d-%d.mat', sensorId, blockLabel));
if ~isfile(filePath)
    raw = zeros(0, 2);
    rawCache(key) = raw;
    return;
end

varName = sprintf('jilu0%d', sensorId);
S = load(filePath, varName);
if isfield(S, varName)
    raw = S.(varName);
else
    raw = zeros(0, 2);
end
rawCache(key) = raw;
if rawCache.Count > rawCacheMaxEntries
    cacheKeys = keys(rawCache);
    remove(rawCache, cacheKeys(1:max(1, numel(cacheKeys) - rawCacheMaxEntries)));
end
end

function fit = fit_low_template_local(Tsen, responseSurface, cfg)
x = Tsen.x_grid(:);
vLow = Tsen.v_grid(:) - Tsen.baseline;
mask = isfinite(x) & isfinite(vLow);
x = x(mask);
v = vLow(mask);
if numel(x) < 50
    error('Low-speed template for sensor %d has too few points.', Tsen.sensor_id);
end

gInit = median(responseSurface.gTrainMm(:), 'omitnan');
p0 = [gInit, 0, 1, 0];
obj = @(p) low_template_objective_local(p, x, v, responseSurface, cfg);
p = fminsearch(obj, p0, cfg.fminOptions);
p = clip_params_local(p, cfg);

[model, over] = eval_raw_tilt_path_local(responseSurface, p(1), p(2), p(3), p(4), x);
[gain, offset, vFit] = solve_gain_offset_local(model, v);
residual = v - vFit;

fit = struct();
fit.sensorId = Tsen.sensor_id;
fit.g0Mm = p(1);
fit.tauMm = p(2);
fit.xScale = p(3);
fit.muGapPerXMm = p(4);
fit.tiltAngleDeg = atan(p(4)) * 180 / pi;
fit.voltageGain = gain;
fit.voltageOffsetMv = offset;
fit.lowFitRmseMv = sqrt(mean(residual.^2, 'omitnan'));
fit.overshootRmseMm = sqrt(mean(over.^2, 'omitnan'));
fit.pointCount = numel(x);
fit.x = x;
fit.vLowMv = v;
fit.vFitMv = vFit;
fit.xLib = p(3) .* (x - p(2));
fit.gEff = p(1) + p(4) .* (x - p(2));
fit.residualMv = residual;
fit.etaMedianMm = get_eta_scalar_local(Tsen, 'eta_median_mm');
fit.etaIqrMm = get_eta_scalar_local(Tsen, 'eta_iqr_mm');
fit.etaLimitMm = get_eta_scalar_local(Tsen, 'eta_limit_mm');
end

function value = low_template_objective_local(p, x, v, responseSurface, cfg)
if any(~isfinite(p))
    value = inf;
    return;
end
pClip = clip_params_local(p, cfg);
penalty = sum((p - pClip).^2) * cfg.overshootPenaltyMv;
[model, over] = eval_raw_tilt_path_local(responseSurface, pClip(1), pClip(2), pClip(3), pClip(4), x);
[~, ~, vFit] = solve_gain_offset_local(model, v);
residual = v - vFit;
value = sqrt(mean(residual.^2, 'omitnan')) + penalty + cfg.overshootPenaltyMv * sqrt(mean(over.^2, 'omitnan'));
end

function p = clip_params_local(p, cfg)
p(1) = min(max(p(1), cfg.gBoundsMm(1)), cfg.gBoundsMm(2));
p(2) = min(max(p(2), cfg.tauBoundsMm(1)), cfg.tauBoundsMm(2));
p(3) = min(max(p(3), cfg.xScaleBounds(1)), cfg.xScaleBounds(2));
p(4) = min(max(p(4), cfg.muBounds(1)), cfg.muBounds(2));
end

function [model, overshoot] = eval_raw_tilt_path_local(responseSurface, g0, tau, k, mu, xOpr)
xLib = k .* (xOpr - tau);
gEff = g0 + mu .* (xOpr - tau);
xMin = min(responseSurface.xGrid(:));
xMax = max(responseSurface.xGrid(:));
gMin = min(responseSurface.gTrainMm(:));
gMax = max(responseSurface.gTrainMm(:));
overX = max(0, xMin - xLib) + max(0, xLib - xMax);
overG = max(0, gMin - gEff) + max(0, gEff - gMax);
overshoot = hypot(overX, overG);
xEval = min(max(xLib, xMin), xMax);
gEval = min(max(gEff, gMin), gMax);
model = eval_response_surface_local(responseSurface, gEval, xEval);
end

function F = eval_response_surface_local(responseSurface, g, xq)
B0 = interp1(responseSurface.xGrid, responseSurface.coeff(:, 1), xq, 'linear', NaN);
B1 = interp1(responseSurface.xGrid, responseSurface.coeff(:, 2), xq, 'linear', NaN);
B2 = interp1(responseSurface.xGrid, responseSurface.coeff(:, 3), xq, 'linear', NaN);
F = B0 + B1 ./ g + B2 .* log(g ./ responseSurface.g0Mm);
end

function [gain, offset, yFit] = solve_gain_offset_local(model, y)
mask = isfinite(model) & isfinite(y);
if nnz(mask) < 4
    gain = 1;
    offset = 0;
    yFit = NaN(size(y));
    return;
end
A = [model(mask), ones(nnz(mask), 1)];
ab = A \ y(mask);
gain = ab(1);
offset = ab(2);
yFit = gain .* model + offset;
end

function highCheckTable = build_highmap_check_local(highMap, sensorCorr, responseSurface, cfg)
if ~isstruct(highMap) || ~isfield(highMap, 'S_v') || isempty(highMap.S_v)
    highCheckTable = table();
    return;
end

rows = {};
for is = 1:numel(sensorCorr)
    sid = sensorCorr(is).sensorId;
    mask = highMap.S_v == sid & isfinite(highMap.x_v) & isfinite(highMap.V_a);
    if ~any(mask)
        continue;
    end
    x = highMap.x_v(mask);
    v = highMap.V_a(mask);
    keep = isfinite(x) & isfinite(v);
    x = x(keep);
    v = v(keep);
    if numel(x) < 20
        continue;
    end

    [model, over] = eval_raw_tilt_path_local( ...
        responseSurface, sensorCorr(is).g0Mm, sensorCorr(is).tauMm, ...
        sensorCorr(is).xScale, sensorCorr(is).muGapPerXMm, x);
    [gain, offset, vFit] = solve_gain_offset_local(model, v);
    residual = v - vFit;
    rows{end + 1, 1} = table(sid, gain, offset, ...
        sqrt(mean(residual.^2, 'omitnan')), sqrt(mean(over.^2, 'omitnan')), numel(x), ...
        'VariableNames', {'sensorId','highMapGain','highMapOffsetMv', ...
        'weightedRmseMv','overshootRmseMm','pointCount'}); %#ok<SAGROW>
end

if isempty(rows)
    highCheckTable = table();
else
    highCheckTable = vertcat(rows{:});
end
end

function plot_low_fit_local(sensorCorr, Template, figFile)
fig = figure('Name', '20241106 Step06I low-template fit', 'Color', 'w', ...
    'Position', [100, 90, 1360, 760], 'NumberTitle', 'off');
tiledlayout(fig, 2, max(2, numel(sensorCorr)), 'TileSpacing', 'compact', 'Padding', 'compact');

for is = 1:numel(sensorCorr)
    nexttile;
    hold on; grid on; box on;
    plot(sensorCorr(is).x, sensorCorr(is).vLowMv, 'k-', 'LineWidth', 1.0, 'DisplayName', 'Template');
    plot(sensorCorr(is).x, sensorCorr(is).vFitMv, 'r--', 'LineWidth', 1.2, 'DisplayName', 'Fit');
    xlabel('x (mm)');
    ylabel('Voltage (mV)');
    title(sprintf('CH%d low fit', sensorCorr(is).sensorId));
    legend('Location', 'best');

    nexttile;
    hold on; grid on; box on;
    plot(sensorCorr(is).x, sensorCorr(is).residualMv, 'b-', 'LineWidth', 1.0);
    yline(0, ':k');
    xlabel('x (mm)');
    ylabel('Residual (mV)');
    title(sprintf('CH%d residual', sensorCorr(is).sensorId));
end

exportgraphics(fig, figFile, 'Resolution', 300);
end

function plot_tilt_path_local(sensorCorr, responseSurface, figFile)
fig = figure('Name', '20241106 Step06I tilt path', 'Color', 'w', ...
    'Position', [120, 100, 1280, 720], 'NumberTitle', 'off');
tiledlayout(fig, 1, max(1, numel(sensorCorr)), 'TileSpacing', 'compact', 'Padding', 'compact');

for is = 1:numel(sensorCorr)
    nexttile;
    hold on; grid on; box on;
    x = sensorCorr(is).x;
    plot(sensorCorr(is).xLib, sensorCorr(is).gEff, 'r-', 'LineWidth', 1.2);
    xline(min(responseSurface.xGrid(:)), ':k');
    xline(max(responseSurface.xGrid(:)), ':k');
    yline(min(responseSurface.gTrainMm(:)), ':k');
    yline(max(responseSurface.gTrainMm(:)), ':k');
    xlabel('x_{lib} (mm)');
    ylabel('g_{eff} (mm)');
    title(sprintf('CH%d path, mu=%.4f', sensorCorr(is).sensorId, sensorCorr(is).muGapPerXMm));
    if ~isempty(x)
        text(sensorCorr(is).xLib(1), sensorCorr(is).gEff(1), 'start', 'Color', [0.1 0.4 0.1]);
        text(sensorCorr(is).xLib(end), sensorCorr(is).gEff(end), 'end', 'Color', [0.7 0.1 0.1]);
    end
end

exportgraphics(fig, figFile, 'Resolution', 300);
end

function value = get_eta_scalar_local(Tsen, fieldName)
value = NaN;
if ~isfield(Tsen, fieldName)
    return;
end
raw = Tsen.(fieldName);
if isempty(raw)
    return;
end
raw = raw(isfinite(raw));
if isempty(raw)
    return;
end
value = raw(1);
end

function Tsen = get_template_sensor_local(Template, sid)
idx = find([Template.Sensor.sensor_id] == sid, 1, 'first');
if isempty(idx)
    error('Template does not contain CH%d.', sid);
end
Tsen = Template.Sensor(idx);
end

function Template = normalize_template_voltage_units_local(Template)
for is = 1:numel(Template.Sensor)
    scale = infer_template_voltage_scale_local(Template.Sensor(is));
    Template.Sensor(is).v_grid = Template.Sensor(is).v_grid(:) * scale;
    if isfield(Template.Sensor(is), 'baseline') && isfinite(Template.Sensor(is).baseline)
        Template.Sensor(is).baseline = Template.Sensor(is).baseline * scale;
    else
        Template.Sensor(is).baseline = 0;
    end
    if isfield(Template.Sensor(is), 'dv_dx') && ~isempty(Template.Sensor(is).dv_dx)
        Template.Sensor(is).dv_dx = Template.Sensor(is).dv_dx(:) * scale;
    end
    if isfield(Template.Sensor(is), 'threshold') && isfinite(Template.Sensor(is).threshold)
        Template.Sensor(is).threshold = Template.Sensor(is).threshold * scale;
    end
    if ~isfield(Template.Sensor(is), 'eta_median_mm')
        Template.Sensor(is).eta_median_mm = NaN;
    end
    if ~isfield(Template.Sensor(is), 'eta_iqr_mm')
        Template.Sensor(is).eta_iqr_mm = NaN;
    end
    if ~isfield(Template.Sensor(is), 'eta_limit_mm')
        Template.Sensor(is).eta_limit_mm = NaN;
    end
    Template.Sensor(is).voltageUnit = 'mV';
end
Template.voltageUnit = 'mV';
end

function scale = infer_template_voltage_scale_local(Tsen)
scale = 1;
if isfield(Tsen, 'voltageUnit') && ~isempty(Tsen.voltageUnit)
    unit = lower(strtrim(char(string(Tsen.voltageUnit))));
    if strcmp(unit, 'v')
        scale = 1000;
        return;
    elseif strcmp(unit, 'mv')
        scale = 1;
        return;
    end
end
v = abs(Tsen.v_grid(:));
v = v(isfinite(v));
if isempty(v)
    return;
end
if median(v, 'omitnan') < 20
    scale = 1000;
end
end

function assert_capacitive_gap_sensors_local(sensorIds, overrideEnv)
allowed = [5 7];
if all(ismember(sensorIds(:).', allowed))
    return;
end
if parse_logical_env_local(overrideEnv, false)
    warning('Using non-capacitive sensors %s in the gap library because %s is enabled.', ...
        mat2str(sensorIds), overrideEnv);
    return;
end
error(['The 20241106 gap library is calibrated for capacitive sensors [5 7]. ' ...
    'Requested sensors %s include non-capacitive probes. Set %s=1 only for diagnostics.'], ...
    mat2str(sensorIds), overrideEnv);
end

function value = parse_logical_env_local(name, defaultValue)
raw = lower(strtrim(getenv(name)));
if isempty(raw)
    value = logical(defaultValue);
    return;
end
if ismember(raw, {'1','true','yes','on'})
    value = true;
elseif ismember(raw, {'0','false','no','off'})
    value = false;
else
    error('%s must be one of on/off, true/false, yes/no, or 1/0.', name);
end
end
