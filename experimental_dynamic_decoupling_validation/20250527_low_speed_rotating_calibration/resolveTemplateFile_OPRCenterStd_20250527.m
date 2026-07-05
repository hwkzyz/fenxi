function [templateFile, sourceMode] = resolveTemplateFile_OPRCenterStd_20250527(routeDir, C, templateSuffix)
% Resolve a low-speed template file, preferring the single blade-sensor library.

if nargin < 3
    templateSuffix = 'GradientXRange030_OPRCenterStd';
end

if ~strcmpi(strtrim(templateSuffix), 'OPRCenterStdRef')
    [templateFile, ok] = assemble_from_single_library_local(routeDir, C);
    if ok
        sourceMode = 'single_blade_sensor_library';
        return;
    end
end

templateDir = fullfile(routeDir, 'output', 'templates');
templateCandidates = {
    fullfile(templateDir, sprintf('Template_LowSpeedRotating_%s_%s_%s.mat', C.caseTag, templateSuffix, C.dataset))
    fullfile(templateDir, sprintf('Template_LowSpeedRotating_B%d_%s_%s_%s.mat', C.bladeId, C.sensorTag, templateSuffix, C.dataset))
    fullfile(templateDir, sprintf('Template_LowSpeedRotating_B%d_%s_%s.mat', C.bladeId, C.sensorTag, C.dataset))
    };

for i = 1:numel(templateCandidates)
    if exist(templateCandidates{i}, 'file') == 2
        templateFile = templateCandidates{i};
        sourceMode = 'legacy_combined_template';
        return;
    end
end

templateFile = find_legacy_template_for_sensors_local(templateDir, C.bladeId, C.sensorIds);
sourceMode = 'legacy_combined_template_fallback';
end

function [assembledFile, ok] = assemble_from_single_library_local(routeDir, C)
ok = false;
assembledFile = '';
libraryDir = fullfile(routeDir, 'output', 'template_library');
sensorIds = C.sensorIds(:).';
singleFiles = cell(numel(sensorIds), 1);
for i = 1:numel(sensorIds)
    sid = sensorIds(i);
    singleFiles{i} = fullfile(libraryDir, sprintf('Template_%s_B%d_S%d.mat', ...
        C.dataset, C.bladeId, sid));
    if exist(singleFiles{i}, 'file') ~= 2
        return;
    end
end

loaded = load(singleFiles{1}, 'Template');
Template = loaded.Template;
sensorList = repmat(Template.Sensor(1), 1, numel(sensorIds));
for i = 1:numel(sensorIds)
    loaded = load(singleFiles{i}, 'Template');
    singleTemplate = loaded.Template;
    if ~isfield(singleTemplate, 'Sensor') || isempty(singleTemplate.Sensor)
        error('Single template has no Sensor entry: %s', singleFiles{i});
    end
    sensorList(i) = singleTemplate.Sensor(1);
end

Template.Sensor = sensorList;
Template.SensorIDs = sensorIds;
Template.SensorTag = ['S', sprintf('%d', sensorIds)];
Template.TargetBlade = C.bladeId;
Template.SourceTemplateMode = 'single_blade_sensor_library';
Template.SourceTemplateLibraryFiles = singleFiles;

assembledDir = fullfile(libraryDir, '_assembled');
if exist(assembledDir, 'dir') ~= 7
    mkdir(assembledDir);
end
assembledFile = fullfile(assembledDir, sprintf('Template_%s_B%d_%s_assembled.mat', ...
    C.dataset, C.bladeId, Template.SensorTag));
save(assembledFile, 'Template', '-v7.3');
ok = true;
end

function templateFile = find_legacy_template_for_sensors_local(templateDir, targetBlade, sensorIds)
sensorTag = ['S', sprintf('%d', sensorIds)];
caseTag = sprintf('B%d_%s', targetBlade, sensorTag);
candidates = {
    sprintf('Template_LowSpeedRotating_%s_GradientXRange030_OPRCenterStd_20250527.mat', caseTag)
    sprintf('Template_LowSpeedRotating_%s_GradientXRange030_OPRAnchored_20250527.mat', caseTag)
    sprintf('Template_LowSpeedRotating_%s_GradientXRange030_20250527.mat', caseTag)
    sprintf('Template_LowSpeedRotating_%s_20250527.mat', caseTag)
    };
for i = 1:numel(candidates)
    candidate = fullfile(templateDir, candidates{i});
    if exist(candidate, 'file') == 2
        templateFile = candidate;
        return;
    end
end

files = dir(fullfile(templateDir, sprintf('Template_LowSpeedRotating_%s*.mat', caseTag)));
if isempty(files)
    files = dir(fullfile(templateDir, sprintf('Template_LowSpeedRotating_B%d_*.mat', targetBlade)));
end
if isempty(files)
    error('Template file not found for blade %d sensors %s in %s.', ...
        targetBlade, mat2str(sensorIds), templateDir);
end
[~, idx] = max([files.datenum]);
templateFile = fullfile(files(idx).folder, files(idx).name);
end
