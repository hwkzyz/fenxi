function Summary = Run_ParameterStudy_20241106(selection, executeRuns)
%RUN_PARAMETERSTUDY_20241106 Run isolated one-factor formal identifications.
%
% List available cases without running:
%   Run_ParameterStudy_20241106
%
% Preview one group:
%   Run_ParameterStudy_20241106('vp_top_k', false)
%
% Execute one group or one named case:
%   Run_ParameterStudy_20241106('vp_top_k', true)
%   Run_ParameterStudy_20241106('window_5lap', true)

if nargin < 1 || isempty(selection)
    selection = 'list';
end
if nargin < 2
    executeRuns = false;
end

studyDir = fileparts(mfilename('fullpath'));
packageRoot = fileparts(fileparts(studyDir));
addpath(packageRoot, studyDir, ...
    fullfile(packageRoot, 'functions', 'foundation'), ...
    fullfile(packageRoot, 'functions', 'preparation'), ...
    fullfile(packageRoot, 'functions', 'gap_aware'), ...
    fullfile(packageRoot, 'functions', 'utilities'));
Study = Config_ParameterStudy_20241106();

if strcmpi(string(selection), "list")
    Summary = case_catalog_local(Study.Cases);
    disp(Summary);
    fprintf('\nRun one group with, for example:\n');
    fprintf('  Run_ParameterStudy_20241106(''vp_top_k'', true)\n');
    return;
end

selectedCases = select_cases_local(Study.Cases, selection);
Plan = case_catalog_local(selectedCases);
disp(Plan);
if ~executeRuns
    fprintf('\nPreview only. Pass true as the second argument to execute.\n');
    Summary = Plan;
    return;
end

if ~exist(Study.outputRoot, 'dir')
    mkdir(Study.outputRoot);
end
rows = cell(numel(selectedCases), 1);
batchTag = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
for iCase = 1:numel(selectedCases)
    fprintf('\n[%d/%d] Parameter case: %s\n', ...
        iCase, numel(selectedCases), selectedCases(iCase).name);
    rows{iCase} = run_case_local(packageRoot, Study, selectedCases(iCase));
    CurrentSummary = vertcat(rows{1:iCase});
    writetable(CurrentSummary, fullfile(Study.outputRoot, ...
        ['ParameterStudy_Summary_', batchTag, '.csv']));
    writetable(CurrentSummary, fullfile(Study.outputRoot, ...
        'ParameterStudy_Summary_Latest.csv'));
end
Summary = vertcat(rows{:});
fprintf('\nParameter study complete: %s\n', Study.outputRoot);
end


function row = run_case_local(packageRoot, Study, C)
caseIndex = find(strcmp({Study.Cases.name}, C.name), 1, 'first');
assert(~isempty(caseIndex), 'Parameter case is not registered: %s', C.name);
caseDir = fullfile(Study.outputRoot, sprintf('C%02d', caseIndex));
if ~exist(caseDir, 'dir')
    mkdir(caseDir);
end

StudyOverride = C.override;
StudyOverride.paths.foundationResults = fullfile(caseDir, 'f');
StudyOverride.paths.gapResults = fullfile(caseDir, 'g');
StudyOverride.paths.comparison = fullfile(caseDir, 'c');
StudyOverride.paths.figures = fullfile(caseDir, 'fig');
StudyMetadata = struct('name', C.name, 'group', C.group, ...
    'description', C.description, 'environment', C.environment, ...
    'createdAt', char(datetime('now')));
overrideFile = fullfile(caseDir, 'StudyOverride_20241106.mat');
save(overrideFile, 'StudyOverride', 'StudyMetadata');

envNames = [{'BTT_PARAMETER_STUDY_OVERRIDE_FILE'}, fieldnames(C.environment).'];
envValues = cellfun(@getenv, envNames, 'UniformOutput', false);
cleanupEnv = onCleanup(@() restore_environment_local(envNames, envValues));
setenv('BTT_PARAMETER_STUDY_OVERRIDE_FILE', overrideFile);
environmentNames = fieldnames(C.environment);
for i = 1:numel(environmentNames)
    setenv(environmentNames{i}, C.environment.(environmentNames{i}));
end

started = tic;
status = "pass";
message = "";
trendFile = "";
try
    for programNumber = Study.formalPrograms
        programFile = find_program_local(packageRoot, programNumber);
        evalin('base', sprintf('run(''%s'')', ...
            strrep(programFile, '''', '''''')));
    end
    trendFile = find_gap_trend_local(StudyOverride.paths.gapResults);
    metrics = summarize_trend_local(trendFile);
catch ME
    status = "failed";
    message = string(getReport(ME, 'basic', 'hyperlinks', 'off'));
    metrics = empty_metrics_local();
end
elapsedSec = toc(started);

EffectiveConfig = Config_20241106();
save(fullfile(caseDir, 'StudyManifest_20241106.mat'), ...
    'C', 'StudyMetadata', 'StudyOverride', 'EffectiveConfig', ...
    'status', 'message', 'trendFile', 'elapsedSec');
row = table(string(C.name), string(C.group), string(C.description), ...
    status, message, elapsedSec, string(trendFile), ...
    EffectiveConfig.case.analysisStartTimeSec, ...
    EffectiveConfig.window.windowBladePasses, ...
    EffectiveConfig.window.slidingStepBladePasses, ...
    EffectiveConfig.frequency.searchHz(1), ...
    EffectiveConfig.frequency.searchHz(2), ...
    EffectiveConfig.frequency.vpTopK, ...
    EffectiveConfig.frequency.refineHalfWidthHz, ...
    EffectiveConfig.model.deltaGapLimitMm, ...
    metrics.windowCount, metrics.EO12Count, metrics.dominantEO, ...
    metrics.meanFrequencyHz, metrics.meanAmplitudeMm, ...
    metrics.meanPlainRmseMv, metrics.meanAbsDeltaGapMm, ...
    metrics.meanEoMarginPercent, metrics.frequencyBoundFraction, ...
    'VariableNames', {'case_name','group','description','status', ...
    'message','elapsed_seconds','trend_file','analysis_time_s', ...
    'window_laps','step_laps','frequency_min_hz','frequency_max_hz', ...
    'vp_top_k','refine_half_width_hz','delta_gap_limit_mm', ...
    'window_count','EO12_count','dominant_EO','mean_frequency_hz', ...
    'mean_amplitude_mm','mean_plain_rmse_mV','mean_abs_delta_gap_mm', ...
    'mean_EO_margin_percent','frequency_bound_fraction'});
clear cleanupEnv;
restore_environment_local(envNames, envValues);
end


function metrics = summarize_trend_local(trendFile)
T = readtable(trendFile, 'TextType', 'string');
required = {'gap_EO','gap_frequency_hz','gap_amplitude_mm', ...
    'gap_plain_rmse_mV','gap_mean_delta_gap_mm', ...
    'gap_EO_margin_percent','gap_frequency_offset_at_bound'};
assert(all(ismember(required, T.Properties.VariableNames)), ...
    'Gap trend is missing required parameter-study columns: %s', trendFile);
valid = isfinite(T.gap_EO) & isfinite(T.gap_amplitude_mm);
assert(any(valid), 'Gap trend contains no finite identified windows: %s', trendFile);
metrics = empty_metrics_local();
metrics.windowCount = height(T);
metrics.EO12Count = nnz(T.gap_EO(valid) == 12);
metrics.dominantEO = mode(T.gap_EO(valid));
metrics.meanFrequencyHz = mean(T.gap_frequency_hz(valid), 'omitnan');
metrics.meanAmplitudeMm = mean(T.gap_amplitude_mm(valid), 'omitnan');
metrics.meanPlainRmseMv = mean(T.gap_plain_rmse_mV(valid), 'omitnan');
metrics.meanAbsDeltaGapMm = mean(abs(T.gap_mean_delta_gap_mm(valid)), 'omitnan');
metrics.meanEoMarginPercent = mean(T.gap_EO_margin_percent(valid), 'omitnan');
metrics.frequencyBoundFraction = mean(T.gap_frequency_offset_at_bound(valid) ~= 0, 'omitnan');
end


function metrics = empty_metrics_local()
metrics = struct('windowCount', NaN, 'EO12Count', NaN, 'dominantEO', NaN, ...
    'meanFrequencyHz', NaN, 'meanAmplitudeMm', NaN, ...
    'meanPlainRmseMv', NaN, 'meanAbsDeltaGapMm', NaN, ...
    'meanEoMarginPercent', NaN, 'frequencyBoundFraction', NaN);
end


function file = find_program_local(packageRoot, number)
matches = dir(fullfile(packageRoot, sprintf('Main%02d_*.m', number)));
assert(isscalar(matches), ...
    'Expected exactly one Main%02d program in %s.', number, packageRoot);
file = fullfile(matches(1).folder, matches(1).name);
end


function file = find_gap_trend_local(resultDir)
matches = dir(fullfile(resultDir, ...
    'Step07J_NestedStaticWarp_VPFullWave_Trend_*.csv'));
assert(~isempty(matches), 'No GapAware trend was generated in %s.', resultDir);
[~, order] = sort([matches.datenum], 'descend');
file = fullfile(matches(order(1)).folder, matches(order(1)).name);
end


function selected = select_cases_local(Cases, selection)
requested = string(selection);
if isscalar(requested) && strcmpi(requested, "all")
    selected = Cases;
    return;
end
mask = false(size(Cases));
for i = 1:numel(Cases)
    mask(i) = any(strcmpi(requested, string(Cases(i).name))) || ...
        any(strcmpi(requested, string(Cases(i).group)));
end
assert(any(mask), 'Unknown parameter-study case or group: %s', ...
    strjoin(requested, ', '));
selected = Cases(mask);
end


function T = case_catalog_local(Cases)
T = table(string({Cases.name}).', string({Cases.group}).', ...
    string({Cases.description}).', ...
    'VariableNames', {'case_name','group','description'});
end


function restore_environment_local(names, values)
for i = 1:numel(names)
    setenv(names{i}, values{i});
end
end
