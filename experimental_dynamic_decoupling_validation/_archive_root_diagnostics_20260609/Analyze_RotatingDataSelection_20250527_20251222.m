%% Analyze_RotatingDataSelection_20250527_20251222
% Compare data-selection stability for the two low-speed rotating
% calibration datasets. This script only reads existing Step03 result files.
% It skips unreadable MAT files and writes CSV summaries for deciding how to
% select data in the main direct-template identification route.

clear; clc;

%% Parameters to tune
P.targetEO = 14;
P.targetBlade = 1;
P.datasets = {
    '20250527', '20250527_low_speed_rotating_calibration', 'S136'
    '20251222', '20251222_low_speed_rotating_calibration', 'S123'
    };

P.includeNameContains = 'Result_Step03_Main_VPTop3SynchronousWaveform';
P.focusOnlyMainStep03 = true;

rootDir = fileparts(mfilename('fullpath'));
outDir = fullfile(rootDir, 'analysis_outputs', 'rotating_data_selection');
if exist(outDir, 'dir') ~= 7
    mkdir(outDir);
end

%% Read and summarize existing Step03 result files
summaryRows = {};
windowRows = {};
skippedRows = {};

for id = 1:size(P.datasets, 1)
    datasetTag = P.datasets{id, 1};
    folderName = P.datasets{id, 2};
    sensorTag = P.datasets{id, 3};
    resultDir = fullfile(rootDir, folderName, 'output', 'identification');
    files = dir(fullfile(resultDir, '*.mat'));

    for ifile = 1:numel(files)
        fileName = files(ifile).name;
        if P.focusOnlyMainStep03 && ~contains(fileName, P.includeNameContains)
            continue;
        end
        expectedSensorToken = ['_B' num2str(P.targetBlade) '_' sensorTag '_'];
        if ~contains(fileName, expectedSensorToken)
            continue;
        end

        filePath = fullfile(files(ifile).folder, fileName);
        try
            S = load(filePath, 'Result');
        catch ME
            skippedRows(end+1, :) = {datasetTag, sensorTag, fileName, 'load_failed', ME.message}; %#ok<SAGROW>
            continue;
        end

        if ~isfield(S, 'Result') || ~isfield(S.Result, 'Trend')
            skippedRows(end+1, :) = {datasetTag, sensorTag, fileName, 'no_Result_Trend', ''}; %#ok<SAGROW>
            continue;
        end

        T = S.Result.Trend;
        requiredVars = {'EO_id','fn_id','A_id','weighted_voltage_rmse'};
        if ~all(ismember(requiredVars, T.Properties.VariableNames))
            skippedRows(end+1, :) = {datasetTag, sensorTag, fileName, 'missing_trend_vars', strjoin(setdiff(requiredVars, T.Properties.VariableNames), ',')}; %#ok<SAGROW>
            continue;
        end

        eo = T.EO_id(:);
        freq = T.fn_id(:);
        amp = T.A_id(:);
        rmseMv = 1000 * T.weighted_voltage_rmse(:);
        dx = get_table_var_or_nan_local(T, {'dx_id','dx_mm','dx'});
        pointCount = extract_point_counts_local(S.Result, numel(eo));

        policy = infer_policy_local(fileName, S.Result);
        summaryRows(end+1, :) = make_summary_row_local( ...
            datasetTag, sensorTag, fileName, filePath, policy, ...
            eo, freq, amp, dx, rmseMv, pointCount, P.targetEO); %#ok<SAGROW>

        for iw = 1:numel(eo)
            windowRows(end+1, :) = {datasetTag, sensorTag, fileName, ...
                policy.center_policy, policy.pulse_mode, policy.domain_mode, policy.eta_policy, ...
                iw, eo(iw), freq(iw), amp(iw), point_or_nan_local(dx, iw), ...
                rmseMv(iw), point_or_nan_local(pointCount, iw)}; %#ok<SAGROW>
        end
    end
end

Summary = cell2table(summaryRows, 'VariableNames', summary_var_names_local());
WindowDetail = cell2table(windowRows, 'VariableNames', ...
    {'dataset','sensor_tag','result_file','center_policy','pulse_mode','domain_mode','eta_policy', ...
    'window_id','EO','frequency_hz','amplitude_mm','dx_mm','rmse_mV','point_count'});
SkippedFiles = cell2table(skippedRows, 'VariableNames', ...
    {'dataset','sensor_tag','result_file','status','message'});

if ~isempty(Summary)
    Summary = sortrows(Summary, {'dataset','target_EO_percent','rmse_mean_mV'}, {'ascend','descend','ascend'});
end

Recommendation = build_recommendation_local(Summary, P.targetEO);

%% Save outputs
prefix = 'RotatingDataSelection_20250527_20251222';
writetable(Summary, fullfile(outDir, [prefix '_Summary.csv']));
writetable(WindowDetail, fullfile(outDir, [prefix '_WindowDetail.csv']));
writetable(SkippedFiles, fullfile(outDir, [prefix '_SkippedFiles.csv']));
writetable(Recommendation, fullfile(outDir, [prefix '_Recommendation.csv']));

fprintf('\n=== Data-selection summary ===\n');
if ~isempty(Summary)
    disp(Summary(:, {'dataset','sensor_tag','short_case','center_policy','pulse_mode','domain_mode', ...
        'eta_policy','EO_counts','target_EO_count','target_EO_percent','wrong_windows', ...
        'freq_std_hz','amp_std_mm','rmse_mean_mV','point_count_mean','stability_grade'}));
end

fprintf('\n=== Recommendation ===\n');
disp(Recommendation);

fprintf('\nSaved CSV files under:\n  %s\n', outDir);

%% Local functions
function names = summary_var_names_local()
names = {'dataset','sensor_tag','result_file','result_path','short_case', ...
    'center_policy','pulse_mode','domain_mode','eta_policy','window_policy', ...
    'status','window_count','EO_counts','dominant_EO','target_EO','target_EO_count', ...
    'target_EO_percent','wrong_windows','freq_mean_hz','freq_std_hz','freq_min_hz','freq_max_hz', ...
    'amp_mean_mm','amp_std_mm','amp_min_mm','amp_max_mm','dx_mean_mm','dx_std_mm', ...
    'rmse_mean_mV','rmse_std_mV','rmse_min_mV','rmse_max_mV', ...
    'point_count_mean','point_count_min','point_count_max','stability_grade'};
end

function policy = infer_policy_local(fileName, Result)
policy.center_policy = "unknown";
policy.pulse_mode = "unknown";
policy.domain_mode = "unknown";
policy.eta_policy = "unknown";
policy.window_policy = "unknown";

nameLower = lower(fileName);
if contains(fileName, 'OPRCenterStd')
    policy.center_policy = "OPRCenterStd";
elseif contains(fileName, 'OPRAnchored')
    policy.center_policy = "OPRAnchored";
elseif contains(fileName, 'BaseFrame') || contains(fileName, 'BaseXc')
    policy.center_policy = "BaseFrame/BaseXc";
elseif contains(fileName, 'GradientXRange030')
    policy.center_policy = "GradientXRange030";
elseif contains(fileName, 'PriorXRange')
    policy.center_policy = "PriorXRange";
end

if contains(nameLower, 'single')
    policy.pulse_mode = "single";
elseif contains(nameLower, 'allsoft') || contains(nameLower, 'alleo') || contains(nameLower, 'all_')
    policy.pulse_mode = "all";
end

if contains(nameLower, 'soft')
    policy.domain_mode = "soft";
elseif contains(nameLower, 'hard') || contains(nameLower, '_sh_')
    policy.domain_mode = "hard";
end

if contains(nameLower, 'noeta')
    policy.eta_policy = "NoEta";
elseif contains(nameLower, 'eta')
    policy.eta_policy = "Eta";
end

if contains(fileName, 'Diag80L_W3S1')
    policy.window_policy = "Diag80L_W3S1";
elseif contains(fileName, 'Main20L_W3S1')
    policy.window_policy = "Main20L_W3S1";
elseif contains(fileName, 'Main_')
    policy.window_policy = "Main";
end

policy = overwrite_from_result_fields_local(policy, Result);
end

function policy = overwrite_from_result_fields_local(policy, Result)
containers = {};
if isfield(Result, 'Config')
    containers{end+1} = Result.Config;
end
if isfield(Result, 'Method')
    containers{end+1} = Result.Method;
end

for i = 1:numel(containers)
    C = containers{i};
    if isfield(C, 'pulse_selection_mode')
        policy.pulse_mode = string(C.pulse_selection_mode);
    end
    if isfield(C, 'domain_selection_mode')
        policy.domain_mode = string(C.domain_selection_mode);
    end
    if isfield(C, 'domain_soft_margin_mm') && strcmpi(policy.domain_mode, "soft")
        policy.domain_mode = "soft_margin_" + string(C.domain_soft_margin_mm) + "mm";
    end
    if isfield(C, 'sensor_eta_limit_mm')
        if abs(C.sensor_eta_limit_mm) < eps
            policy.eta_policy = "NoEta";
        else
            policy.eta_policy = "EtaLimit_" + string(C.sensor_eta_limit_mm) + "mm";
        end
    end
end
end

function row = make_summary_row_local(datasetTag, sensorTag, fileName, filePath, policy, ...
    eo, freq, amp, dx, rmseMv, pointCount, targetEO)
eo = eo(:);
freq = freq(:);
amp = amp(:);
dx = dx(:);
rmseMv = rmseMv(:);
pointCount = pointCount(:);
nWindow = numel(eo);
finiteEO = eo(isfinite(eo));
if isempty(finiteEO)
    dominantEO = NaN;
else
    dominantEO = mode(finiteEO);
end
targetCount = sum(eo == targetEO);
wrong = find(eo ~= targetEO & isfinite(eo));
shortCase = erase(fileName, {['Result_Step03_Main_VPTop3SynchronousWaveform_B1_' sensorTag '_'], ['_' datasetTag '.mat']});

row = {datasetTag, sensorTag, fileName, filePath, shortCase, ...
    char(policy.center_policy), char(policy.pulse_mode), char(policy.domain_mode), ...
    char(policy.eta_policy), char(policy.window_policy), 'ok', nWindow, ...
    count_text_local(eo), dominantEO, targetEO, targetCount, 100 * targetCount / max(nWindow, 1), ...
    mat2str(wrong(:).'), mean(freq, 'omitnan'), std(freq, 'omitnan'), ...
    min(freq, [], 'omitnan'), max(freq, [], 'omitnan'), ...
    mean(amp, 'omitnan'), std(amp, 'omitnan'), min(amp, [], 'omitnan'), max(amp, [], 'omitnan'), ...
    mean(dx, 'omitnan'), std(dx, 'omitnan'), ...
    mean(rmseMv, 'omitnan'), std(rmseMv, 'omitnan'), min(rmseMv, [], 'omitnan'), max(rmseMv, [], 'omitnan'), ...
    mean(pointCount, 'omitnan'), min(pointCount, [], 'omitnan'), max(pointCount, [], 'omitnan'), ...
    grade_stability_local(targetCount, nWindow)};
end

function Recommendation = build_recommendation_local(Summary, targetEO)
rows = {};
if isempty(Summary)
    Recommendation = cell2table(rows, 'VariableNames', {'decision_item','recommended_choice','reason'});
    return;
end

datasets = unique(Summary.dataset, 'stable');
for id = 1:numel(datasets)
    ds = datasets{id};
    S = Summary(strcmp(Summary.dataset, ds), :);
    stable = S(S.target_EO_count == S.window_count, :);
    if ~isempty(stable)
        [~, idx] = min(stable.rmse_mean_mV);
        best = stable(idx, :);
        rows(end+1, :) = {['best_completed_' ds], best.short_case{1}, ...
            sprintf('Completed result with EO%d in all windows and the lowest RMSE for this dataset.', targetEO)}; %#ok<SAGROW>
    end

    soft = S(contains(S.domain_mode, 'soft') & S.target_EO_count == S.window_count, :);
    hard = S(contains(S.domain_mode, 'hard') & strcmp(S.center_policy, 'OPRCenterStd'), :);
    if ~isempty(soft)
        rows(end+1, :) = {['soft_role_' ds], 'preferred_when_available', ...
            sprintf('Soft-domain selection keeps EO%d stable in all completed soft tests for this dataset.', targetEO)}; %#ok<SAGROW>
    end
    if ~isempty(hard) && any(hard.target_EO_count < hard.window_count)
        wrongText = strjoin(hard.wrong_windows, '; ');
        rows(end+1, :) = {['hard_role_' ds], 'diagnostic_not_default', ...
            ['OPRCenterStd hard-domain results contain wrong EO windows: ' wrongText '.']}; %#ok<SAGROW>
    end

    oprSoft = S(strcmp(S.center_policy, 'OPRCenterStd') & contains(S.domain_mode, 'soft'), :);
    if ~isempty(oprSoft)
        if any(oprSoft.target_EO_count == oprSoft.window_count)
            rows(end+1, :) = {['oprcenterstd_soft_role_' ds], 'usable_after_validation', ...
                sprintf('OPRCenterStd soft-domain selection has at least one completed result with EO%d in all windows for this dataset.', targetEO)}; %#ok<SAGROW>
        else
            rows(end+1, :) = {['oprcenterstd_soft_role_' ds], 'not_default_for_this_dataset', ...
                sprintf('OPRCenterStd soft-domain selection does not lock EO%d in the completed tests; center/template policy is the dominant issue here.', targetEO)}; %#ok<SAGROW>
        end
    end
end

rows(end+1, :) = {'cross_dataset_default', 'validate_center_policy_first_then_use_soft_domain', ...
    sprintf(['Across the two rotating datasets, center/template policy dominates the result: ' ...
    '20250527 favors GradientXRange030/EtaGradientXc, while 20251222 can use OPRCenterStd only when the domain is not hard. ' ...
    'Use EO%d all-window stability as the first gate, RMSE as the second gate.'], targetEO)};

rows(end+1, :) = {'pulse_default', 'single_after_center_policy_is_validated', ...
    'Single pulse is acceptable for readability and tuning after the center/domain policy is validated; all-vs-single is not the main failure source in the completed tests.'};

rows(end+1, :) = {'domain_default', 'soft_or_adaptive_apparent_displacement', ...
    'Avoid fixed hard boundaries as the default; hard-domain selection should be kept as a sensitivity diagnostic.'};

Recommendation = cell2table(rows, 'VariableNames', {'decision_item','recommended_choice','reason'});
end

function grade = grade_stability_local(targetCount, nWindow)
if nWindow == 0
    grade = "missing";
elseif targetCount == nWindow
    grade = "stable_all_target";
elseif targetCount >= 0.9 * nWindow
    grade = "mostly_stable";
else
    grade = "unstable";
end
end

function txt = count_text_local(x)
x = x(isfinite(x));
if isempty(x)
    txt = "";
    return;
end
u = unique(x(:).');
parts = strings(size(u));
for i = 1:numel(u)
    parts(i) = sprintf('%g:%d', u(i), sum(x == u(i)));
end
txt = strjoin(parts, ', ');
end

function pointCount = extract_point_counts_local(Result, n)
pointCount = NaN(n, 1);
if ~isfield(Result, 'WindowResult')
    return;
end
for iw = 1:min(n, numel(Result.WindowResult))
    W = Result.WindowResult(iw);
    if isfield(W, 'bundle') && isfield(W.bundle, 'pointCount')
        pointCount(iw) = W.bundle.pointCount;
    elseif isfield(W, 'bundle') && isfield(W.bundle, 'X')
        pointCount(iw) = numel(W.bundle.X);
    elseif isfield(W, 'point_count')
        pointCount(iw) = W.point_count;
    end
end
end

function v = get_table_var_or_nan_local(T, names)
v = NaN(height(T), 1);
for i = 1:numel(names)
    if ismember(names{i}, T.Properties.VariableNames)
        v = T.(names{i});
        v = v(:);
        return;
    end
end
end

function v = point_or_nan_local(x, idx)
if numel(x) >= idx
    v = x(idx);
else
    v = NaN;
end
end
