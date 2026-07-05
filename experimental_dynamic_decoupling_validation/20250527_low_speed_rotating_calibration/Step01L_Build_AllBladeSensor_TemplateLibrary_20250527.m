%% Step01L_Build_AllBladeSensor_TemplateLibrary_20250527
% Build/publish a single blade-sensor low-speed template library.
% The underlying template construction still uses Step01_Main unchanged.

clear; clc;

routeDir = fileparts(mfilename('fullpath'));
dataset = '20250527';
templateSuffix = 'GradientXRange030_OPRCenterStd';
bladeIds = 1:6;
sensorIds = [1 3 6];
buildMissingCombined = parse_logical_env_local('STEP01L_BUILD_MISSING', true);
strictMissing = parse_logical_env_local('STEP01L_STRICT_MISSING', false);

bladeOverride = strtrim(getenv('STEP01L_BLADE_IDS'));
if ~isempty(bladeOverride)
    bladeIds = sscanf(bladeOverride, '%d').';
end
sensorOverride = strtrim(getenv('STEP01L_SENSOR_IDS'));
if ~isempty(sensorOverride)
    sensorIds = sscanf(sensorOverride, '%d').';
end

libraryDir = fullfile(routeDir, 'output', 'template_library');
if exist(libraryDir, 'dir') ~= 7
    mkdir(libraryDir);
end

oldBladeEnv = getenv('BLADE_CASE_BLADE_ID');
oldSensorEnv = getenv('BLADE_CASE_SENSOR_IDS');
cleanupObj = onCleanup(@() restore_env_local(oldBladeEnv, oldSensorEnv));

maxRows = numel(bladeIds) * numel(sensorIds);
indexRows = cell(maxRows, 10);
indexCount = 0;
missingRows = cell(maxRows, 7);
missingCount = 0;
fprintf('\n=== Step01L: single blade-sensor template library (%s) ===\n', dataset);
fprintf('Blades: %s, sensors: %s, build missing combined templates: %d, strict missing: %d\n', ...
    mat2str(bladeIds), mat2str(sensorIds), buildMissingCombined, strictMissing);

for ib = 1:numel(bladeIds)
    bladeId = bladeIds(ib);
    setenv('BLADE_CASE_BLADE_ID', num2str(bladeId));
    setenv('BLADE_CASE_SENSOR_IDS', sprintf('%d ', sensorIds));
    C = CaseConfig();
    stablePlanFile = fullfile(routeDir, 'output', 'stable_window_plan', sprintf( ...
        'Step01_CoverageFirstStableWindowPlan_%s_%s.csv', C.caseTag, dataset));
    referenceTemplateFile = fullfile(routeDir, 'output', 'templates', sprintf( ...
        'Template_LowSpeedRotating_%s_OPRCenterStdRef_%s.mat', C.caseTag, dataset));
    combinedFile = fullfile(routeDir, 'output', 'templates', sprintf( ...
        'Template_LowSpeedRotating_%s_%s_%s.mat', C.caseTag, templateSuffix, dataset));

    if exist(combinedFile, 'file') ~= 2 && buildMissingCombined
        try
            if exist(stablePlanFile, 'file') ~= 2
                if exist(referenceTemplateFile, 'file') ~= 2
                    fprintf('\nBuild missing reference template: B%d %s\n', bladeId, C.sensorTag);
                    run_script_in_child_matlab_local(routeDir, ...
                        'Step01_Main_Build_OPRCenterStd_Template_20250527.m', bladeId, sensorIds, ...
                        {'STEP01_ALLOW_REFERENCE_ONLY', '1'});
                end
                if ~has_dynamic_map_for_case_local(routeDir, bladeId, sensorIds)
                    fprintf('\nBuild missing bootstrap dynamic map: B%d %s\n', bladeId, C.sensorTag);
                    run_script_in_child_matlab_local(routeDir, ...
                        'Step02_Main_Build_OPRCenterStd_DynamicMap_20250527.m', bladeId, sensorIds, ...
                        {'STEP02_TEMPLATE_SUFFIX', 'OPRCenterStdRef'; ...
                         'STEP02_DYNAMIC_SUFFIX', 'BootstrapFromRef'});
                end
                fprintf('\nBuild missing stable-window plan: B%d %s\n', bladeId, C.sensorTag);
                run_script_in_child_matlab_local(routeDir, ...
                    'Step00_Build_CoverageFirstStableWindowPlan_20250527.m', bladeId, sensorIds);
            end
            fprintf('\nBuild missing combined template: B%d %s\n', bladeId, C.sensorTag);
            run_script_in_child_matlab_local(routeDir, ...
                'Step01_Main_Build_OPRCenterStd_Template_20250527.m', bladeId, sensorIds);
        catch ME
            if strictMissing
                rethrow(ME);
            end
            warning('Skip B%d after failed upstream build: %s', bladeId, ME.message);
            missingCount = missingCount + 1;
            missingRows(missingCount, :) = {dataset, bladeId, ['S', sprintf('%d', sensorIds)], ...
                combinedFile, stablePlanFile, 'upstream_build_failed', ME.message}; %#ok<AGROW>
            continue;
        end
    end

    if exist(combinedFile, 'file') ~= 2
        warning('Skip B%d: combined template is missing: %s', bladeId, combinedFile);
        if exist(stablePlanFile, 'file') ~= 2
            reason = 'stable_plan_and_combined_template_missing';
        else
            reason = 'combined_template_missing';
        end
        missingCount = missingCount + 1;
        missingRows(missingCount, :) = {dataset, bladeId, ['S', sprintf('%d', sensorIds)], ...
            combinedFile, stablePlanFile, reason, ''}; %#ok<AGROW>
        continue;
    end

    fprintf('\nPublish from combined template: %s\n', combinedFile);
    rows = publish_single_templates_local(combinedFile, libraryDir, dataset, bladeId, sensorIds);
    for ir = 1:size(rows, 1)
        indexCount = indexCount + 1;
        indexRows(indexCount, :) = rows(ir, :);
    end
end

indexFile = fullfile(libraryDir, sprintf('Template_Index_%s_All.csv', dataset));
write_index_csv_local(indexFile, indexRows(1:indexCount, :));
fprintf('\nTemplate index: %s\n', indexFile);
missingFile = fullfile(libraryDir, sprintf('Template_Missing_%s_All.csv', dataset));
write_missing_csv_local(missingFile, missingRows(1:missingCount, :));
fprintf('Template missing report: %s\n', missingFile);

function rows = publish_single_templates_local(combinedFile, libraryDir, dataset, bladeId, sensorIds)
loaded = load(combinedFile, 'Template');
CombinedTemplate = loaded.Template;
rows = {};
for i = 1:numel(sensorIds)
    sid = sensorIds(i);
    sensorIndex = find([CombinedTemplate.Sensor.sensor_id] == sid, 1, 'first');
    if isempty(sensorIndex)
        warning('Skip B%d S%d: sensor not present in %s', bladeId, sid, combinedFile);
        continue;
    end

    Template = CombinedTemplate;
    Template.Sensor = CombinedTemplate.Sensor(sensorIndex);
    Template.SensorIDs = sid;
    Template.SensorTag = sprintf('S%d', sid);
    Template.TargetBlade = bladeId;
    Template.SourceTemplateMode = 'single_blade_sensor_library';
    Template.SourceCombinedTemplateFile = combinedFile;
    Template.LibraryCreatedBy = mfilename;

    templateFile = fullfile(libraryDir, sprintf('Template_%s_B%d_S%d.mat', dataset, bladeId, sid));
    save(templateFile, 'Template', '-v7.3');

    S = Template.Sensor;
    rows(end+1, :) = {dataset, bladeId, sid, templateFile, combinedFile, ...
        get_numeric_field_local(S, 'point_count'), get_numeric_field_local(S, 'lap_count'), ...
        get_numeric_field_local(S, 'eta_limit_mm'), min(S.x_grid(:)), max(S.x_grid(:))}; %#ok<AGROW>
    fprintf('  B%d S%d -> %s\n', bladeId, sid, templateFile);
end
end

function write_index_csv_local(indexFile, rows)
fid = fopen(indexFile, 'w', 'n', 'UTF-8');
if fid < 0
    error('Cannot write template index: %s', indexFile);
end
cleanupObj = onCleanup(@() fclose(fid));
fprintf(fid, 'dataset,bladeId,sensorId,templateFile,sourceCombinedTemplateFile,pointCount,lapCount,etaLimitMm,xMinMm,xMaxMm\n');
for i = 1:size(rows, 1)
    fprintf(fid, '%s,%d,%d,"%s","%s",%.15g,%.15g,%.15g,%.15g,%.15g\n', rows{i, :});
end
end

function write_missing_csv_local(missingFile, rows)
fid = fopen(missingFile, 'w', 'n', 'UTF-8');
if fid < 0
    error('Cannot write missing-template report: %s', missingFile);
end
cleanupObj = onCleanup(@() fclose(fid));
fprintf(fid, 'dataset,bladeId,sensorTag,combinedTemplateFile,stablePlanFile,reason,message\n');
for i = 1:size(rows, 1)
    fprintf(fid, '%s,%d,%s,"%s","%s",%s,"%s"\n', rows{i, 1}, rows{i, 2}, ...
        rows{i, 3}, rows{i, 4}, rows{i, 5}, rows{i, 6}, sanitize_csv_text_local(rows{i, 7}));
end
end

function txt = sanitize_csv_text_local(txt)
txt = strrep(char(txt), '"', '""');
txt = strrep(txt, newline, ' ');
end

function value = get_numeric_field_local(S, name)
if isfield(S, name) && ~isempty(S.(name))
    value = double(S.(name));
else
    value = NaN;
end
end

function value = parse_logical_env_local(name, defaultValue)
txt = lower(strtrim(getenv(name)));
if isempty(txt)
    value = defaultValue;
elseif ismember(txt, {'1','true','yes','on'})
    value = true;
elseif ismember(txt, {'0','false','no','off'})
    value = false;
else
    error('%s must be logical.', name);
end
end

function tf = has_dynamic_map_for_case_local(routeDir, bladeId, sensorIds)
dynamicDir = fullfile(routeDir, 'output', 'dynamic_maps');
sensorTag = ['S', sprintf('%d', sensorIds)];
files = dir(fullfile(dynamicDir, sprintf('DynamicMap_B%d_%s_*BootstrapFromRef_%s.mat', ...
    bladeId, sensorTag, '20250527')));
tf = ~isempty(files);
end

function run_script_in_child_matlab_local(routeDir, scriptName, bladeId, sensorIds, extraEnvPairs)
if nargin < 5
    extraEnvPairs = {};
end
matlabExe = fullfile(matlabroot, 'bin', 'matlab');
if ispc
    matlabExe = [matlabExe, '.exe'];
end
sensorText = strtrim(sprintf('%d ', sensorIds));
batchCmd = sprintf(['cd(''%s''); setenv(''BLADE_CASE_BLADE_ID'',''%d''); ' ...
    'setenv(''BLADE_CASE_SENSOR_IDS'',''%s''); run(''%s'');'], ...
    escape_single_quote_local(routeDir), bladeId, sensorText, ...
    escape_single_quote_local(fullfile(routeDir, scriptName)));
for i = 1:size(extraEnvPairs, 1)
    batchCmd = sprintf('setenv(''%s'',''%s''); %s', ...
        extraEnvPairs{i, 1}, extraEnvPairs{i, 2}, batchCmd);
end
cmd = sprintf('"%s" -batch "%s"', matlabExe, batchCmd);
[status, outputText] = system(cmd);
fprintf('%s\n', outputText);
if status ~= 0
    error('Step01 child MATLAB failed for B%d. Command: %s', bladeId, cmd);
end
end

function txt = escape_single_quote_local(txt)
txt = strrep(txt, '''', '''''');
end

function restore_env_local(oldBladeEnv, oldSensorEnv)
setenv('BLADE_CASE_BLADE_ID', oldBladeEnv);
setenv('BLADE_CASE_SENSOR_IDS', oldSensorEnv);
end
