%% Analyze_DataSelection_Stability_20251222
% Read the 20251222 rotating-calibration and gap-prior result folders, then
% compare which data-selection/model route gives stable vibration
% identification. This script only reads existing result files and writes
% diagnostic CSV tables. It does not rerun identification.

clear; clc;

%% Parameters to tune
C0 = CaseConfig();
P.targetEO = 14;
P.targetBlade = C0.bladeId;
P.sensorTag = C0.sensorTag;
P.caseTag = C0.caseTag;
P.dataset = C0.dataset;

thisDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(thisDir);
rotDir = fullfile(rootDir, '20251222_low_speed_rotating_calibration');
gapDir = thisDir;
rotResultDir = fullfile(rotDir, 'output', 'identification');
gapOutDir = fullfile(gapDir, 'outputs');
analysisDir = fullfile(gapOutDir, 'data_selection_stability_analysis');
if exist(analysisDir, 'dir') ~= 7
    mkdir(analysisDir);
end

%% Existing result cases
directCases = {
    'old_gradient_xrange030',       'GradientXRange030',  'mixed/old', 'mixed/old', 'legacy', ...
        sprintf('Result_Step03_Main_VPTop3SynchronousWaveform_B%d_%s_Main_GradientXRange030_%s.mat', P.targetBlade, P.sensorTag, P.dataset)
    'oprcenterstd_all_soft_noeta',  'OPRCenterStd',       'all',       'soft',      'NoEta', ...
        sprintf('Result_Step03_Main_VPTop3SynchronousWaveform_B%d_%s_Diag_AllSoft_OPRCenterStd_NoEta_%s.mat', P.targetBlade, P.sensorTag, P.dataset)
    'oprcenterstd_single_soft_noeta','OPRCenterStd',      'single',    'soft',      'NoEta', ...
        sprintf('Result_Step03_Main_VPTop3SynchronousWaveform_B%d_%s_Diag_SingleSoft_OPRCenterStd_NoEta_%s.mat', P.targetBlade, P.sensorTag, P.dataset)
    'oprcenterstd_all_hard_noeta',  'OPRCenterStd',       'all',       'hard',      'NoEta', ...
        sprintf('Result_Step03_Main_VPTop3SynchronousWaveform_B%d_%s_Diag_AllHard_OPRCenterStd_NoEta_%s.mat', P.targetBlade, P.sensorTag, P.dataset)
    'oprcenterstd_single_hard_noeta','OPRCenterStd',      'single',    'hard',      'NoEta', ...
        sprintf('Result_Step03_Main_VPTop3SynchronousWaveform_B%d_%s_Main_DirectTemplate_OPRCenterStd_Dx035_NoEta_%s.mat', P.targetBlade, P.sensorTag, P.dataset)
    };

gapCases = {
    'current_main_gaptilt',    'current main route: gap_tilt only', ...
        sprintf('Step07J_NestedStaticWarp_VPFullWave_%s_%s_main_gaptilt.mat', P.dataset, P.caseTag)
    'current_method_compare',  'current comparison route: fixed/gap_only/gap_tilt', ...
        sprintf('Step07J_NestedStaticWarp_VPFullWave_%s_%s_method_compare.mat', P.dataset, P.caseTag)
    };

%% 1. Direct-template data-selection stability
directRows = {};
directWindowRows = {};
for ic = 1:size(directCases, 1)
    caseName = directCases{ic, 1};
    templatePolicy = directCases{ic, 2};
    pulseMode = directCases{ic, 3};
    domainMode = directCases{ic, 4};
    etaPolicy = directCases{ic, 5};
    fileName = directCases{ic, 6};
    filePath = fullfile(rotResultDir, fileName);

    if exist(filePath, 'file') ~= 2
        directRows(end+1, :) = missing_direct_row_local( ...
            caseName, templatePolicy, pulseMode, domainMode, etaPolicy, filePath); %#ok<SAGROW>
        continue;
    end

    S = load(filePath, 'Result');
    T = S.Result.Trend;
    eo = T.EO_id(:);
    freq = T.fn_id(:);
    amp = T.A_id(:);
    rmseMv = 1000 * T.weighted_voltage_rmse(:);
    pointCount = extract_step03_point_counts_local(S.Result);
    directRows(end+1, :) = make_stability_row_local( ...
        caseName, 'direct_template', templatePolicy, pulseMode, domainMode, etaPolicy, ...
        filePath, eo, freq, amp, NaN(size(eo)), rmseMv, pointCount, P.targetEO); %#ok<SAGROW>

    for iw = 1:numel(eo)
        directWindowRows(end+1, :) = {caseName, iw, eo(iw), freq(iw), amp(iw), NaN, rmseMv(iw), point_or_nan_local(pointCount, iw)}; %#ok<SAGROW>
    end
end
DirectStability = cell2table(directRows, 'VariableNames', stability_var_names_local());
DirectWindow = cell2table(directWindowRows, 'VariableNames', ...
    {'case_name','window_id','EO','frequency_hz','amplitude_mm','dx_mm','rmse_mV','point_count'});

%% 2. Gap-prior model stability
modelNames = {'direct','fixed','gap_only','gap_tilt'};
gapRows = {};
gapWindowRows = {};
gapRowsForRecommendation = table();
gapPointRows = {};
for ic = 1:size(gapCases, 1)
    caseName = gapCases{ic, 1};
    routePolicy = gapCases{ic, 2};
    filePath = fullfile(gapOutDir, gapCases{ic, 3});
    if exist(filePath, 'file') ~= 2
        filePath = fullfile(gapOutDir, strrep(gapCases{ic, 3}, ...
            ['_' P.caseTag '_'], sprintf('_B%d_%s_', P.targetBlade, P.sensorTag)));
    end
    if exist(filePath, 'file') ~= 2
        for im = 1:numel(modelNames)
            gapRows(end+1, :) = missing_stability_row_local(caseName, modelNames{im}, routePolicy, filePath); %#ok<SAGROW>
        end
        continue;
    end

    S = load(filePath, 'Result');
    R = S.Result;
    nWindow = numel(R.WindowResult);
    pointCount = extract_step07j_point_counts_local(R);
    for im = 1:numel(modelNames)
        model = modelNames{im};
        [eo, freq, amp, dx, rmseMv] = extract_step07j_model_series_local(R, model);
        gapRows(end+1, :) = make_stability_row_local( ...
            caseName, model, routePolicy, model, '', '', filePath, ...
            eo, freq, amp, dx, rmseMv, pointCount, P.targetEO); %#ok<SAGROW>
        for iw = 1:nWindow
            gapWindowRows(end+1, :) = {caseName, model, iw, eo(iw), freq(iw), amp(iw), dx(iw), rmseMv(iw), point_or_nan_local(pointCount, iw)}; %#ok<SAGROW>
        end
    end

    if isfield(R.WindowResult(1).modelFits, 'gap_only')
        for iw = 1:nWindow
            fit = R.WindowResult(iw).modelFits.gap_only;
            sensorIds = R.WindowResult(iw).bundle.sensorIds;
            for is = 1:numel(sensorIds)
                gapPointRows(end+1, :) = {caseName, iw, sensorIds(is), fit.deltaGapMm(is), fit.EO, fit.freqHz, fit.amplitudeMm, fit.dxMm, fit.weightedRmseMv}; %#ok<SAGROW>
            end
        end
    end
end
GapModelStability = cell2table(gapRows, 'VariableNames', stability_var_names_local());
GapWindow = cell2table(gapWindowRows, 'VariableNames', ...
    {'case_name','model','window_id','EO','frequency_hz','amplitude_mm','dx_mm','rmse_mV','point_count'});
GapOnlyDeltaGapWindow = cell2table(gapPointRows, 'VariableNames', ...
    {'case_name','window_id','sensor_id','delta_gap_mm','EO','frequency_hz','amplitude_mm','dx_mm','rmse_mV'});
if ~isempty(GapOnlyDeltaGapWindow)
    GapOnlyDeltaGapSensorSummary = groupsummary(GapOnlyDeltaGapWindow, {'case_name','sensor_id'}, ...
        {'mean','std'}, {'delta_gap_mm'});
else
    GapOnlyDeltaGapSensorSummary = table();
end

%% 3. Recommendation table
Recommendation = build_recommendation_table_local(DirectStability, GapModelStability, P.targetEO);

%% 4. Save outputs
prefix = sprintf('StabilityAnalysis_%s_B%d_%s', P.dataset, P.targetBlade, P.sensorTag);
writetable(DirectStability, fullfile(analysisDir, [prefix '_DirectTemplate_Summary.csv']));
writetable(DirectWindow, fullfile(analysisDir, [prefix '_DirectTemplate_Window.csv']));
writetable(GapModelStability, fullfile(analysisDir, [prefix '_GapPrior_ModelSummary.csv']));
writetable(GapWindow, fullfile(analysisDir, [prefix '_GapPrior_Window.csv']));
writetable(GapOnlyDeltaGapWindow, fullfile(analysisDir, [prefix '_GapOnly_DeltaGap_WindowSensor.csv']));
writetable(GapOnlyDeltaGapSensorSummary, fullfile(analysisDir, [prefix '_GapOnly_DeltaGap_SensorSummary.csv']));
writetable(Recommendation, fullfile(analysisDir, [prefix '_Recommendation.csv']));

fprintf('\n=== Direct-template stability ===\n');
disp(DirectStability(:, {'case_name','EO_counts','target_EO_count','target_EO_percent','wrong_windows','freq_std_hz','amp_std_mm','rmse_mean_mV'}));
fprintf('\n=== Gap-prior model stability ===\n');
disp(GapModelStability(:, {'case_name','model_name','EO_counts','target_EO_count','target_EO_percent','wrong_windows','freq_std_hz','amp_std_mm','rmse_mean_mV'}));
fprintf('\n=== Recommendation ===\n');
disp(Recommendation);
fprintf('\nSaved analysis tables under:\n  %s\n', analysisDir);

%% Local functions
function names = stability_var_names_local()
names = {'case_name','model_name','template_or_route','pulse_or_model','domain_mode', ...
    'eta_policy','result_file','status','window_count','EO_counts','dominant_EO', ...
    'target_EO','target_EO_count','target_EO_percent','wrong_windows', ...
    'freq_mean_hz','freq_std_hz','freq_min_hz','freq_max_hz', ...
    'amp_mean_mm','amp_std_mm','amp_min_mm','amp_max_mm', ...
    'dx_mean_mm','dx_std_mm','dx_min_mm','dx_max_mm', ...
    'rmse_mean_mV','rmse_std_mV','rmse_min_mV','rmse_max_mV', ...
    'point_count_mean','point_count_min','point_count_max','stability_grade'};
end

function row = make_stability_row_local(caseName, modelName, routeText, pulseText, domainText, etaText, ...
    resultFile, eo, freq, amp, dx, rmseMv, pointCount, targetEO)
eo = eo(:); freq = freq(:); amp = amp(:); dx = dx(:); rmseMv = rmseMv(:);
finiteEO = eo(isfinite(eo));
if isempty(finiteEO)
    dominantEO = NaN;
else
    dominantEO = mode(finiteEO);
end
targetCount = sum(eo == targetEO);
nWindow = numel(eo);
wrong = find(eo ~= targetEO & isfinite(eo));
grade = grade_stability_local(targetCount, nWindow);
row = {caseName, modelName, routeText, pulseText, domainText, etaText, resultFile, 'ok', ...
    nWindow, count_text_local(eo), dominantEO, targetEO, targetCount, ...
    100 * targetCount / max(nWindow, 1), mat2str(wrong(:).'), ...
    mean(freq, 'omitnan'), std(freq, 'omitnan'), min(freq, [], 'omitnan'), max(freq, [], 'omitnan'), ...
    mean(amp, 'omitnan'), std(amp, 'omitnan'), min(amp, [], 'omitnan'), max(amp, [], 'omitnan'), ...
    mean(dx, 'omitnan'), std(dx, 'omitnan'), min(dx, [], 'omitnan'), max(dx, [], 'omitnan'), ...
    mean(rmseMv, 'omitnan'), std(rmseMv, 'omitnan'), min(rmseMv, [], 'omitnan'), max(rmseMv, [], 'omitnan'), ...
    mean(pointCount, 'omitnan'), min(pointCount, [], 'omitnan'), max(pointCount, [], 'omitnan'), grade};
end

function row = missing_direct_row_local(caseName, templatePolicy, pulseMode, domainMode, etaPolicy, filePath)
row = {caseName, 'direct_template', templatePolicy, pulseMode, domainMode, etaPolicy, filePath, ...
    'missing', 0, '', NaN, NaN, 0, NaN, '', NaN, NaN, NaN, NaN, ...
    NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, ...
    NaN, NaN, NaN, 'missing'};
end

function row = missing_stability_row_local(caseName, modelName, routePolicy, filePath)
row = {caseName, modelName, routePolicy, modelName, '', '', filePath, ...
    'missing', 0, '', NaN, NaN, 0, NaN, '', NaN, NaN, NaN, NaN, ...
    NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, ...
    NaN, NaN, NaN, 'missing'};
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

function pointCount = extract_step03_point_counts_local(Result)
n = height(Result.Trend);
pointCount = NaN(n, 1);
if ~isfield(Result, 'WindowResult')
    return;
end
for iw = 1:min(n, numel(Result.WindowResult))
    W = Result.WindowResult(iw);
    if isfield(W, 'bundle') && isfield(W.bundle, 'X')
        pointCount(iw) = numel(W.bundle.X);
    elseif isfield(W, 'point_count')
        pointCount(iw) = W.point_count;
    end
end
end

function pointCount = extract_step07j_point_counts_local(Result)
n = numel(Result.WindowResult);
pointCount = NaN(n, 1);
for iw = 1:n
    W = Result.WindowResult(iw);
    if isfield(W, 'bundle') && isfield(W.bundle, 'pointCount')
        pointCount(iw) = W.bundle.pointCount;
    elseif isfield(W, 'bundle') && isfield(W.bundle, 'X')
        pointCount(iw) = numel(W.bundle.X);
    end
end
end

function [eo, freq, amp, dx, rmseMv] = extract_step07j_model_series_local(Result, model)
n = numel(Result.WindowResult);
eo = NaN(n, 1); freq = NaN(n, 1); amp = NaN(n, 1); dx = NaN(n, 1); rmseMv = NaN(n, 1);
for iw = 1:n
    if strcmpi(model, 'direct')
        T = Result.Trend;
        eo(iw) = T.direct_EO(iw);
        freq(iw) = T.direct_frequency_hz(iw);
        amp(iw) = T.direct_amplitude_mm(iw);
        rmseMv(iw) = T.direct_rmse_mV(iw);
        continue;
    end
    fit = Result.WindowResult(iw).modelFits.(model);
    eo(iw) = fit.EO;
    freq(iw) = fit.freqHz;
    amp(iw) = fit.amplitudeMm;
    dx(iw) = fit.dxMm;
    rmseMv(iw) = fit.weightedRmseMv;
end
end

function v = point_or_nan_local(pointCount, iw)
if numel(pointCount) >= iw
    v = pointCount(iw);
else
    v = NaN;
end
end

function Recommendation = build_recommendation_table_local(DirectStability, GapModelStability, targetEO)
rows = {};
okDirect = DirectStability(strcmp(DirectStability.status, 'ok'), :);
if ~isempty(okDirect)
    stableDirect = okDirect(okDirect.target_EO_count == okDirect.window_count, :);
    if ~isempty(stableDirect)
        [~, idx] = min(stableDirect.rmse_mean_mV);
        bestDirect = stableDirect(idx, :);
        rows(end+1, :) = {'lowest_rmse_direct_reference', bestDirect.case_name{1}, ...
            sprintf('Reference only: stable EO%d in all windows with the lowest RMSE among completed direct-template routes.', targetEO)};
    end

    cleanSoft = okDirect(strcmp(okDirect.case_name, 'oprcenterstd_single_soft_noeta'), :);
    if ~isempty(cleanSoft) && cleanSoft.target_EO_count(1) == cleanSoft.window_count(1)
        rows(end+1, :) = {'clean_direct_template_default', cleanSoft.case_name{1}, ...
            sprintf('Recommended for the new OPRCenterStd route: single pulse keeps the data simple, soft domain keeps EO%d stable in all windows.', targetEO)};
    elseif ~isempty(stableDirect)
        stableClean = stableDirect(strcmp(stableDirect.template_or_route, 'OPRCenterStd'), :);
        if ~isempty(stableClean)
            [~, idx] = min(stableClean.rmse_mean_mV);
            bestClean = stableClean(idx, :);
            rows(end+1, :) = {'clean_direct_template_default', bestClean.case_name{1}, ...
                sprintf('Best stable OPRCenterStd direct-template route found in the completed cases for EO%d.', targetEO)};
        end
    end
end

gap = GapModelStability(strcmp(GapModelStability.status, 'ok') & strcmp(GapModelStability.model_name, 'gap_only'), :);
if ~isempty(gap)
    stableGap = gap(gap.target_EO_count == gap.window_count, :);
    if ~isempty(stableGap)
        [~, idx] = min(stableGap.rmse_mean_mV);
        bestGap = stableGap(idx, :);
        rows(end+1, :) = {'lowest_rmse_gap_prior_completed', bestGap.case_name{1}, ...
            sprintf('Reference only: gap_only keeps EO%d stable in all windows and has the lowest RMSE among completed Step07J outputs.', targetEO)};
    end

    cleanHardGap = gap(strcmp(gap.case_name, 'new_oprcenterstd_hard'), :);
    if ~isempty(cleanHardGap)
        rows(end+1, :) = {'completed_clean_hard_gap_prior', cleanHardGap.case_name{1}, ...
            sprintf('The new OPRCenterStd hard-domain Step07J gap_only result is stable for EO%d, but its direct/fixed reference still has wrong EO windows.', targetEO)};
    end
end

hardDirect = okDirect(contains(okDirect.domain_mode, 'hard'), :);
if ~isempty(hardDirect) && any(hardDirect.target_EO_count < hardDirect.window_count)
    rows(end+1, :) = {'hard_domain_role', 'diagnostic_only', ...
        'Hard domain exposes direct/fixed EO instability and should be used as a sensitivity test, not as the default direct-template route.'};
end

rows(end+1, :) = {'next_gap_prior_test', 'rerun_step07j_with_oprcenterstd_single_soft_noeta', ...
    'This is the missing clean comparison: keep the new OPRCenterStd center-time/template idea, but use single+soft direct data selection before judging the final gap-prior route.'};

Recommendation = cell2table(rows, 'VariableNames', {'decision_item','recommended_choice','reason'});
end
