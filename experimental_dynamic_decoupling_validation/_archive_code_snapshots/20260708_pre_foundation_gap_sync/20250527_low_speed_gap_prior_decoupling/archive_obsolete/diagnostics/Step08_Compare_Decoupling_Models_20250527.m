%% Step08: All-window comparison of identification routes
% This script does not re-run nonlinear optimization. It compares all saved
% sliding windows from:
%   1) the current hybrid low-speed-template/gap-library Step07 route,
%   2) the copied 20250527 gap-library-only Step07 route when available,
%   3) the gap-library static-template full-wave Step07A route when available,
%   4) the low-speed rotating-template identification route when available.

clear; clc;

thisDir = fileparts(mfilename('fullpath'));
outDir = fullfile(thisDir, 'outputs');
figDir = fullfile(outDir, 'figures_step08');
if ~exist(figDir, 'dir')
    mkdir(figDir);
end

cfg = struct();
cfg.targetBlade = 1;
cfg.analysisSensors = [1, 3, 6];
sensorOverride = strtrim(getenv('STEP07_ANALYSIS_SENSORS'));
if ~isempty(sensorOverride)
    parsedSensors = sscanf(sensorOverride, '%d').';
    if isempty(parsedSensors)
        error('STEP07_ANALYSIS_SENSORS must contain integer sensor IDs, for example "1 3 6".');
    end
    cfg.analysisSensors = parsedSensors;
end
sensorTag = ['S', sprintf('%d', cfg.analysisSensors)];

currentFile = fullfile(outDir, sprintf( ...
    'Step07_Main_Decoupled_Identification_20250527_B%d_%s.mat', ...
    cfg.targetBlade, sensorTag));
gapFullWaveFile = fullfile(outDir, sprintf( ...
    'Step07A_GapLibrary_StaticTemplate_FullWave_20250527_B%d_%s.mat', ...
    cfg.targetBlade, sensorTag));
oldGapFile = fullfile(thisDir, '..', '20250527', 'outputs', sprintf( ...
    'Step07_Main_Decoupled_Identification_20250527_B%d_%s.mat', ...
    cfg.targetBlade, sensorTag));
lowSpeedDir = fullfile(thisDir, '..', '20250527_low_speed_rotating_calibration', ...
    'output', 'identification');
lowSpeedMainFile = first_existing_file_local({
    fullfile(lowSpeedDir, sprintf( ...
        'Result_Step03_Main_VPTop3SynchronousWaveform_B%d_%s_Main_GradientXRange030_20250527.mat', ...
        cfg.targetBlade, sensorTag))
    fullfile(lowSpeedDir, sprintf( ...
        'Result_Step03_04_FirstOrderVPAdaptive_B%d_%s_Main_GradientXRange030_20250527.mat', ...
        cfg.targetBlade, sensorTag))
    fullfile(lowSpeedDir, sprintf( ...
        'Result_Step03_01_DirectLowSpeedWaveform_B%d_%s_Main_GradientXRange030_20250527.mat', ...
        cfg.targetBlade, sensorTag))
    fullfile(lowSpeedDir, sprintf( ...
        'Result_ProposedFixedTemplateSliding_B%d_%s_20250527.mat', ...
        cfg.targetBlade, sensorTag))
    });
lowSpeedBaselineFile = first_existing_file_local({
    fullfile(lowSpeedDir, sprintf( ...
        'Result_Step03_Baseline_AllEOFullWaveform_B%d_%s_Main_GradientXRange030_20250527.mat', ...
        cfg.targetBlade, sensorTag))
    fullfile(lowSpeedDir, sprintf( ...
        'Result_Step03_01_DirectLowSpeedWaveform_B%d_%s_Main_GradientXRange030_20250527.mat', ...
        cfg.targetBlade, sensorTag))
    });

allRows = {};
summaryRows = {};
diagnosticRows = {};
if isfile(currentFile)
    S = load(currentFile, 'result');
    T = normalize_decoupled_trend_local(S.result.Trend, ...
        'hybrid_low_speed_template_plus_gap_library', currentFile, true);
    allRows{end+1} = T;
    summaryRows(end+1, :) = summarize_windows_local(T);
    Tprior = normalize_low_prior_trend_local(S.result.Trend, ...
        'low_speed_prior_gap_library_corrected_template', currentFile);
    if ~isempty(Tprior)
        allRows{end+1} = Tprior;
        summaryRows(end+1, :) = summarize_windows_local(Tprior);
    end
    D = normalize_step07_reference_diagnostics_local(S.result.Trend, currentFile);
    if ~isempty(D)
        diagnosticRows{end+1} = D;
    end
else
    warning('Missing current Step07 result: %s', currentFile);
end
if isfile(oldGapFile)
    S = load(oldGapFile, 'result');
    T = normalize_decoupled_trend_local(S.result.Trend, ...
        'gap_library_only_original_20250527', oldGapFile, true);
    allRows{end+1} = T;
    summaryRows(end+1, :) = summarize_windows_local(T);
else
    warning('Missing original gap-only result: %s', oldGapFile);
end
if isfile(gapFullWaveFile)
    S = load(gapFullWaveFile, 'result');
    T = normalize_decoupled_trend_local(S.result.Trend, ...
        'gap_library_static_template_full_wave', gapFullWaveFile, true);
    allRows{end+1} = T;
    summaryRows(end+1, :) = summarize_windows_local(T);
else
    warning('Missing Step07A gap-library full-wave result: %s', gapFullWaveFile);
end
if isfile(lowSpeedMainFile)
    S = load(lowSpeedMainFile, 'Result');
    T = normalize_low_speed_trend_local(S.Result.Trend, ...
        'low_speed_template_main_vp_top3_synchronous', lowSpeedMainFile);
    allRows{end+1} = T;
    summaryRows(end+1, :) = summarize_windows_local(T);
else
    warning('Missing low-speed rotating-template main result.');
end
if isfile(lowSpeedBaselineFile) && ~strcmpi(lowSpeedBaselineFile, lowSpeedMainFile)
    S = load(lowSpeedBaselineFile, 'Result');
    T = normalize_low_speed_trend_local(S.Result.Trend, ...
        'low_speed_template_baseline_all_eo_full_waveform', lowSpeedBaselineFile);
    allRows{end+1} = T;
    summaryRows(end+1, :) = summarize_windows_local(T);
end

if isempty(allRows)
    error('No saved result files were found for comparison.');
end

allWindowTable = vertcat(allRows{:});
summaryTable = cell2table(summaryRows, 'VariableNames', { ...
    'method', 'windowCount', 'frequencyMeanHz', 'frequencyStdHz', ...
    'frequencyRangeHz', 'eoMean', 'eoStd', 'amplitudeMeanMm', ...
    'amplitudeStdMm', 'amplitudeRangeMm', 'meanGapMeanMm', ...
    'meanGapStdMm', 'weightedRmseMeanMv', 'weightedRmseStdMv', ...
    'weightedRmseMinMv', 'weightedRmseMaxMv'});

csvAllFile = fullfile(outDir, sprintf('Step08_AllWindow_Method_Comparison_20250527_B%d_%s.csv', ...
    cfg.targetBlade, sensorTag));
csvSummaryFile = fullfile(outDir, sprintf('Step08_AllWindow_Method_Summary_20250527_B%d_%s.csv', ...
    cfg.targetBlade, sensorTag));
csvDiagnosticFile = fullfile(outDir, sprintf('Step08_Step07_Internal_Reference_Diagnostics_20250527_B%d_%s.csv', ...
    cfg.targetBlade, sensorTag));
matFile = fullfile(outDir, sprintf('Step08_AllWindow_Method_Comparison_20250527_B%d_%s.mat', ...
    cfg.targetBlade, sensorTag));
figFile = fullfile(figDir, sprintf('Step08_AllWindow_Method_Comparison_20250527_B%d_%s.png', ...
    cfg.targetBlade, sensorTag));

writetable(allWindowTable, csvAllFile);
writetable(summaryTable, csvSummaryFile);
if isempty(diagnosticRows)
    referenceDiagnosticTable = table();
else
    referenceDiagnosticTable = vertcat(diagnosticRows{:});
    writetable(referenceDiagnosticTable, csvDiagnosticFile);
end
save(matFile, 'allWindowTable', 'summaryTable', 'cfg', 'currentFile', ...
    'oldGapFile', 'lowSpeedMainFile', 'lowSpeedBaselineFile', ...
    'referenceDiagnosticTable', '-v7.3');
plot_all_window_comparison_local(allWindowTable, figFile);

fprintf('\nStep08 all-window comparison complete.\n');
disp(summaryTable);
fprintf('Saved:\n  %s\n  %s\n  %s\n  %s\n', csvAllFile, csvSummaryFile, matFile, figFile);
if exist('referenceDiagnosticTable', 'var') && ~isempty(referenceDiagnosticTable)
    fprintf('Step07 internal reference diagnostics:\n  %s\n', csvDiagnosticFile);
end

function T = normalize_decoupled_trend_local(trendTable, methodName, sourceFile, hasGap)
n = height(trendTable);
if hasGap
    meanGap = trendTable.mean_gap_mm;
else
    meanGap = nan(n, 1);
end
T = table(repmat(string(methodName), n, 1), repmat(string(sourceFile), n, 1), ...
    trendTable.window_id, trendTable.lap_start, trendTable.lap_end, ...
    trendTable.window_center_time, trendTable.frequency_hz, trendTable.EO, ...
    trendTable.amplitude_mm, meanGap, trendTable.weighted_rmse_mV, ...
    trendTable.point_count, ...
    'VariableNames', {'method', 'sourceFile', 'window_id', 'lap_start', ...
    'lap_end', 'window_center_time', 'frequencyHz', 'EO', 'amplitudeMm', ...
    'meanGapMm', 'weightedRmseMv', 'pointCount'});
end

function filePath = first_existing_file_local(candidates)
filePath = '';
for i = 1:numel(candidates)
    if isfile(candidates{i})
        filePath = candidates{i};
        return;
    end
end
if ~isempty(candidates)
    filePath = candidates{1};
end
end

function T = normalize_low_speed_trend_local(trendTable, methodName, sourceFile)
n = height(trendTable);
if ismember('frequency_hz', trendTable.Properties.VariableNames)
    freqHz = trendTable.frequency_hz;
elseif ismember('fn_id', trendTable.Properties.VariableNames)
    freqHz = trendTable.fn_id;
else
    error('Low-speed trend table does not contain frequency_hz or fn_id.');
end
if ismember('EO', trendTable.Properties.VariableNames)
    eo = trendTable.EO;
elseif ismember('EO_id', trendTable.Properties.VariableNames)
    eo = trendTable.EO_id;
else
    error('Low-speed trend table does not contain EO or EO_id.');
end
if ismember('amplitude_mm', trendTable.Properties.VariableNames)
    ampMm = trendTable.amplitude_mm;
elseif ismember('A_id', trendTable.Properties.VariableNames)
    ampMm = trendTable.A_id;
else
    error('Low-speed trend table does not contain amplitude_mm or A_id.');
end
if ismember('weighted_rmse_V', trendTable.Properties.VariableNames)
    rmseMv = 1000 * trendTable.weighted_rmse_V;
elseif ismember('weighted_voltage_rmse', trendTable.Properties.VariableNames)
    rmseMv = 1000 * trendTable.weighted_voltage_rmse;
else
    error('Low-speed trend table does not contain weighted_rmse_V or weighted_voltage_rmse.');
end
T = table(repmat(string(methodName), n, 1), repmat(string(sourceFile), n, 1), ...
    trendTable.window_id, trendTable.lap_start, trendTable.lap_end, ...
    trendTable.window_center_time, freqHz, eo, ...
    ampMm, nan(n, 1), rmseMv, ...
    trendTable.point_count, ...
    'VariableNames', {'method', 'sourceFile', 'window_id', 'lap_start', ...
    'lap_end', 'window_center_time', 'frequencyHz', 'EO', 'amplitudeMm', ...
    'meanGapMm', 'weightedRmseMv', 'pointCount'});
end

function T = normalize_low_prior_trend_local(trendTable, methodName, sourceFile)
required = {'low_prior_frequency_hz', 'low_prior_EO', 'low_prior_amplitude_mm', ...
    'low_prior_weighted_rmse_mV', 'low_prior_delta_gap_mm'};
if ~all(ismember(required, trendTable.Properties.VariableNames))
    T = table();
    return;
end
n = height(trendTable);
T = table(repmat(string(methodName), n, 1), repmat(string(sourceFile), n, 1), ...
    trendTable.window_id, trendTable.lap_start, trendTable.lap_end, ...
    trendTable.window_center_time, trendTable.low_prior_frequency_hz, ...
    trendTable.low_prior_EO, trendTable.low_prior_amplitude_mm, ...
    low_prior_mean_gap_local(trendTable), trendTable.low_prior_weighted_rmse_mV, ...
    trendTable.point_count, ...
    'VariableNames', {'method', 'sourceFile', 'window_id', 'lap_start', ...
    'lap_end', 'window_center_time', 'frequencyHz', 'EO', 'amplitudeMm', ...
    'meanGapMm', 'weightedRmseMv', 'pointCount'});
end

function meanGap = low_prior_mean_gap_local(trendTable)
if ismember('low_prior_mean_gap_mm', trendTable.Properties.VariableNames)
    meanGap = trendTable.low_prior_mean_gap_mm;
elseif ismember('mean_gap_mm', trendTable.Properties.VariableNames)
    meanGap = trendTable.mean_gap_mm;
else
    meanGap = nan(height(trendTable), 1);
end
end

function T = normalize_step07_reference_diagnostics_local(trendTable, sourceFile)
required = {'template_only_rmse_mV', 'fixed_gap_hybrid_rmse_mV', ...
    'template_to_free_rmse_drop_mV', 'fixed_gap_to_free_rmse_drop_mV'};
if ~all(ismember(required, trendTable.Properties.VariableNames))
    T = table();
    return;
end
n = height(trendTable);
T = table(repmat(string(sourceFile), n, 1), trendTable.window_id, ...
    trendTable.frequency_hz, trendTable.EO, trendTable.weighted_rmse_mV, ...
    trendTable.template_only_rmse_mV, trendTable.fixed_gap_hybrid_rmse_mV, ...
    trendTable.template_to_free_rmse_drop_mV, trendTable.fixed_gap_to_free_rmse_drop_mV, ...
    get_optional_column_local(trendTable, 'low_prior_delta_gap_mm'), ...
    get_optional_column_local(trendTable, 'low_prior_mean_gap_mm'), ...
    get_optional_column_local(trendTable, 'low_prior_frequency_hz'), ...
    get_optional_column_local(trendTable, 'low_prior_EO'), ...
    get_optional_column_local(trendTable, 'low_prior_amplitude_mm'), ...
    get_optional_column_local(trendTable, 'low_prior_weighted_rmse_mV'), ...
    'VariableNames', {'sourceFile', 'window_id', 'frequencyHz', 'EO', ...
    'freeGapHybridRmseMv', 'templateOnlyRmseMv', 'fixedGapHybridRmseMv', ...
    'templateToFreeRmseDropMv', 'fixedGapToFreeRmseDropMv', ...
    'lowPriorDeltaGapMm', 'lowPriorMeanGapMm', 'lowPriorFrequencyHz', 'lowPriorEO', ...
    'lowPriorAmplitudeMm', 'lowPriorWeightedRmseMv'});
end

function x = get_optional_column_local(T, name)
if ismember(name, T.Properties.VariableNames)
    x = T.(name);
else
    x = nan(height(T), 1);
end
end

function row = summarize_windows_local(T)
freqRange = max(T.frequencyHz) - min(T.frequencyHz);
ampRange = max(T.amplitudeMm) - min(T.amplitudeMm);
row = {T.method(1), height(T), mean(T.frequencyHz, 'omitnan'), ...
    std(T.frequencyHz, 'omitnan'), freqRange, mean(T.EO, 'omitnan'), ...
    std(T.EO, 'omitnan'), mean(T.amplitudeMm, 'omitnan'), ...
    std(T.amplitudeMm, 'omitnan'), ampRange, mean(T.meanGapMm, 'omitnan'), ...
    std(T.meanGapMm, 'omitnan'), mean(T.weightedRmseMv, 'omitnan'), ...
    std(T.weightedRmseMv, 'omitnan'), min(T.weightedRmseMv), ...
    max(T.weightedRmseMv)};
end

function plot_all_window_comparison_local(T, figFile)
methods = unique(T.method, 'stable');
colors = lines(numel(methods));

fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 17, 15]);
tiledlayout(fig, 4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

legendHandles = plot_metric_local(T, methods, colors, 'frequencyHz', 'f (Hz)', 'Frequency');
plot_metric_local(T, methods, colors, 'amplitudeMm', 'A (mm)', 'Amplitude');
plot_metric_local(T, methods, colors, 'meanGapMm', 'g (mm)', 'Mean clearance');
plot_metric_local(T, methods, colors, 'weightedRmseMv', 'RMSE (mV)', 'Voltage residual');

legend(legendHandles, cellstr(methods), 'Location', 'southoutside', 'Orientation', 'horizontal', ...
    'Interpreter', 'none', 'FontName', 'Times New Roman', 'FontSize', 8);

exportgraphics(fig, figFile, 'Resolution', 300);
end

function handles = plot_metric_local(T, methods, colors, fieldName, yLabelText, titleText)
nexttile;
hold on;
handles = gobjects(numel(methods), 1);
for i = 1:numel(methods)
    idx = T.method == methods(i);
    y = T.(fieldName)(idx);
    h = plot(T.window_id(idx), y, '-o', 'LineWidth', 1.0, ...
        'MarkerSize', 3.5, 'Color', colors(i, :));
    handles(i) = h;
end
hold off;
box on;
set(gca, 'TickDir', 'in', 'FontName', 'Times New Roman', 'FontSize', 8);
xlabel('Window index', 'FontName', 'Times New Roman', 'FontSize', 9);
ylabel(yLabelText, 'FontName', 'Times New Roman', 'FontSize', 9);
title(titleText, 'FontName', 'Times New Roman', 'FontSize', 9, 'FontWeight', 'normal');
end
