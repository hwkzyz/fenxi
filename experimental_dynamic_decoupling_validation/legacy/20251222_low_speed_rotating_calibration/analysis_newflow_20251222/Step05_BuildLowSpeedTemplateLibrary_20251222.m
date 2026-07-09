%% Step05_BuildLowSpeedTemplateLibrary_20251222
% Build the NewFlow low-speed template library for 20251222.
%
% 20251222 already has a mature OPRCenterStd low-speed template builder.
% This step resolves the accepted template, validates the required sensors,
% exports a NewFlow library artifact, and writes a per-sensor audit table so
% Step06 can identify directly from an explicit low-speed template library.

clear; close all; clc;

%% Parameters to tune
P = NewFlow_Config_20251222();
analysisSensors = P.sensors.analysis;
targetBlade = P.case.bladeId;
templateSuffix = P.template.suffix;

%% Paths
lowSpeedNumberingFile = P.files.lowSpeedNumbering;
outFile = P.files.lowSpeedTemplateLibrary;
summaryCsv = strrep(outFile, '.mat', '_Summary.csv');

ensure_parent_dir_local(outFile);
require_file_local(lowSpeedNumberingFile, 'Step02 low-speed numbering');

[templateFile, templateSourceMode] = resolveTemplateFile_OPRCenterStd_20251222( ...
    P.routeDir, P.case, templateSuffix);
require_file_local(templateFile, 'resolved OPRCenterStd template');

loaded = load(templateFile, 'Template');
if ~isfield(loaded, 'Template')
    error('Resolved template file does not contain variable Template:\n  %s', templateFile);
end
Template = loaded.Template;
validate_template_local(Template, analysisSensors, targetBlade, templateFile);

TemplateSummary = build_template_summary_local(Template, analysisSensors, targetBlade, templateFile);

LowSpeedTemplateLibrary = struct();
LowSpeedTemplateLibrary.dataset = P.dataset;
LowSpeedTemplateLibrary.mode = 'newflow_template_library_from_oprcenterstd_template';
LowSpeedTemplateLibrary.source_template_file = templateFile;
LowSpeedTemplateLibrary.source_template_mode = templateSourceMode;
LowSpeedTemplateLibrary.source_low_speed_numbering_file = lowSpeedNumberingFile;
LowSpeedTemplateLibrary.analysis_sensors = analysisSensors;
LowSpeedTemplateLibrary.sensor_tag = ['S', sprintf('%d', analysisSensors)];
LowSpeedTemplateLibrary.target_blade = targetBlade;
LowSpeedTemplateLibrary.template_suffix = templateSuffix;
LowSpeedTemplateLibrary.Template = Template;
LowSpeedTemplateLibrary.TemplateSummary = TemplateSummary;
LowSpeedTemplateLibrary.note = ['Step05 is the explicit NewFlow low-speed template library. ' ...
    'The template waveform is the accepted 20251222 OPRCenterStd template; ' ...
    'this artifact validates and indexes it for Step06 inline identification.'];

save(outFile, 'LowSpeedTemplateLibrary', 'Template', 'TemplateSummary', '-v7.3');
writetable(TemplateSummary, summaryCsv);

fprintf('\n=== Step05: low-speed template library ===\n');
fprintf('Template source: %s [%s]\n', templateFile, templateSourceMode);
fprintf('Target blade: B%d, sensors: %s\n', targetBlade, mat2str(analysisSensors));
fprintf('Saved:   %s\n', outFile);
fprintf('Summary: %s\n', summaryCsv);
disp(TemplateSummary);

function validate_template_local(Template, sensorIds, targetBlade, templateFile)
requiredFields = {'Sensor', 'SensorIDs'};
for i = 1:numel(requiredFields)
    name = requiredFields{i};
    if ~isfield(Template, name)
        error('Template is missing field %s:\n  %s', name, templateFile);
    end
end
available = [Template.Sensor.sensor_id];
for sid = sensorIds(:).'
    idx = find(available == sid, 1, 'first');
    if isempty(idx)
        error('Template for B%d does not contain CH%d:\n  %s', targetBlade, sid, templateFile);
    end
    S = Template.Sensor(idx);
    sensorRequired = {'x_grid', 'v_grid', 'dv_dx', 'x_domain', 'xc'};
    for k = 1:numel(sensorRequired)
        name = sensorRequired{k};
        if ~isfield(S, name) || isempty(S.(name))
            error('Template CH%d is missing nonempty field %s:\n  %s', sid, name, templateFile);
        end
    end
end
end

function T = build_template_summary_local(Template, sensorIds, targetBlade, templateFile)
available = [Template.Sensor.sensor_id];
rows = repmat(struct( ...
    'BladeID', NaN, ...
    'SensorID', NaN, ...
    'PointCount', NaN, ...
    'XDomainLeftMM', NaN, ...
    'XDomainRightMM', NaN, ...
    'XGridMinMM', NaN, ...
    'XGridMaxMM', NaN, ...
    'XcMM', NaN, ...
    'ThresholdV', NaN, ...
    'EtaLimitMM', NaN, ...
    'TemplateFile', ""), numel(sensorIds), 1);
for i = 1:numel(sensorIds)
    sid = sensorIds(i);
    S = Template.Sensor(available == sid);
    rows(i).BladeID = targetBlade;
    rows(i).SensorID = sid;
    rows(i).PointCount = numel(S.x_grid);
    rows(i).XDomainLeftMM = S.x_domain(1);
    rows(i).XDomainRightMM = S.x_domain(2);
    rows(i).XGridMinMM = min(S.x_grid(:), [], 'omitnan');
    rows(i).XGridMaxMM = max(S.x_grid(:), [], 'omitnan');
    rows(i).XcMM = S.xc;
    rows(i).ThresholdV = get_optional_field_local(S, 'threshold', NaN);
    rows(i).EtaLimitMM = get_optional_field_local(S, 'eta_limit_mm', NaN);
    rows(i).TemplateFile = string(templateFile);
end
T = struct2table(rows);
end

function value = get_optional_field_local(S, fieldName, defaultValue)
if isstruct(S) && isfield(S, fieldName)
    value = S.(fieldName);
else
    value = defaultValue;
end
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
