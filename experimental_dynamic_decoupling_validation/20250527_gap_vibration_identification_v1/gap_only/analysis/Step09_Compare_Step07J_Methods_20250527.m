function Step09_Compare_Step07J_Methods_20250527
%% Compare Step07J methods, 20250527
% This is the separate comparison entry. The main Step07J program computes
% only gap_only by default; this script uses STEP07J_RUN_MODE=comparison to
% generate/read a separate result containing fixed, gap_only and gap_tilt.

clc; close all;

thisDir = fileparts(mfilename('fullpath'));
C0 = CaseConfig();
oldDir = pwd;
cleanupDir = onCleanup(@() cd(oldDir));
cd(thisDir);
outDir = fullfile(thisDir, 'outputs');
figDir = fullfile(outDir, 'figures_step07j_method_comparison_current');
if ~exist(figDir, 'dir')
    mkdir(figDir);
end

dataset = C0.dataset;
bladeTag = C0.caseTag;
defaultSuffix = 'method_compare';
resultSuffix = strtrim(getenv('STEP07J_COMPARE_SUFFIX'));
if isempty(resultSuffix)
    resultSuffix = defaultSuffix;
end

resultFile = strtrim(getenv('STEP07J_COMPARE_RESULT_FILE'));
if isempty(resultFile)
    resultFile = fullfile(outDir, sprintf( ...
        'Step07J_NestedStaticWarp_VPFullWave_%s_%s_%s.mat', ...
        dataset, bladeTag, resultSuffix));
end
if ~isfile(resultFile)
    fprintf('Comparison result not found. Running Step07J in comparison mode...\n');
    previousRunMode = getenv('STEP07J_RUN_MODE');
    previousSuffix = getenv('STEP07J_RESULT_SUFFIX');
    cleanupEnv = onCleanup(@() restore_step07j_env_local(previousRunMode, previousSuffix));
    setenv('STEP07J_RUN_MODE', 'comparison');
    setenv('STEP07J_RESULT_SUFFIX', resultSuffix);
    evalin('base', 'Step07J_NestedStaticWarp_VPFullWave_20250527');
    delete(cleanupEnv);
    restore_step07j_env_local(previousRunMode, previousSuffix);
end
if ~isfile(resultFile)
    error('Missing comparison Step07J result file after rerun: %s', resultFile);
end

S = load(resultFile, 'Result');
Result = S.Result;
Trend = Result.Trend;

fprintf('\n=== Step07J current method comparison: %s %s ===\n', dataset, bladeTag);
fprintf('Result: %s\n', resultFile);
fprintf('Main model recorded in Result: %s\n', string(Result.cfg.mainModel));
if ~isfield(Result.cfg, 'runMode') || ~strcmpi(Result.cfg.runMode, 'comparison')
    error('This is not a comparison result. Use STEP07J_RUN_MODE=comparison or delete/regenerate: %s', resultFile);
end

Summary = build_method_summary_local(Trend);
summaryFile = fullfile(outDir, sprintf( ...
    'Compare_Step07J_CurrentMethods_Summary_%s_%s_%s.csv', ...
    dataset, bladeTag, resultSuffix));
writetable(Summary, summaryFile);

windowFile = fullfile(outDir, sprintf( ...
    'Compare_Step07J_CurrentMethods_Window_%s_%s_%s.csv', ...
    dataset, bladeTag, resultSuffix));
windowVars = {'window_id','direct_EO','direct_frequency_hz','direct_amplitude_mm','direct_rmse_mV', ...
    'fixed_EO','fixed_frequency_hz','fixed_amplitude_mm','fixed_rmse_mV', ...
    'gap_EO','gap_frequency_hz','gap_amplitude_mm','gap_rmse_mV','gap_mean_delta_gap_mm', ...
    'tilt_EO','tilt_frequency_hz','tilt_amplitude_mm','tilt_rmse_mV','tilt_mean_delta_gap_mm','tilt_mean_delta_mu'};
windowVars = windowVars(ismember(windowVars, Trend.Properties.VariableNames));
writetable(Trend(:, windowVars), windowFile);

figFile = fullfile(figDir, sprintf( ...
    'Compare_Step07J_CurrentMethods_Trend_%s_%s_%s.png', ...
    dataset, bladeTag, resultSuffix));
plot_method_trend_local(Trend, figFile, dataset, bladeTag);

disp(Summary);
fprintf('Saved method summary:\n  %s\n', summaryFile);
fprintf('Saved window comparison:\n  %s\n', windowFile);
fprintf('Saved method trend figure:\n  %s\n', figFile);

makeWaveformFigures = parse_logical_env_local('STEP07J_COMPARE_WAVEFORM_FIGURES', true);
if makeWaveformFigures
    fprintf('\nGenerating separate waveform figures for fixed, gap_only and gap_tilt...\n');
    previousMethodResultFile = getenv('STEP07J_METHOD_RESULT_FILE');
    setenv('STEP07J_METHOD_RESULT_FILE', resultFile);
    Step10_Visualize_Step07J_MethodSeparatePipeline_20250527;
    setenv('STEP07J_METHOD_RESULT_FILE', previousMethodResultFile);
end

delete(cleanupDir);
cd(oldDir);
end

function Summary = build_method_summary_local(T)
spec = struct( ...
    'method', {'direct_low_template_main','fixed','gap_only','gap_tilt'}, ...
    'eo', {'direct_EO','fixed_EO','gap_EO','tilt_EO'}, ...
    'freq', {'direct_frequency_hz','fixed_frequency_hz','gap_frequency_hz','tilt_frequency_hz'}, ...
    'amp', {'direct_amplitude_mm','fixed_amplitude_mm','gap_amplitude_mm','tilt_amplitude_mm'}, ...
    'rmse', {'direct_rmse_mV','fixed_rmse_mV','gap_rmse_mV','tilt_rmse_mV'});

n = numel(spec);
method = strings(n, 1);
windowCount = zeros(n, 1);
eoDistribution = strings(n, 1);
dominantEO = NaN(n, 1);
meanFrequencyHz = NaN(n, 1);
stdFrequencyHz = NaN(n, 1);
meanAmplitudeMm = NaN(n, 1);
meanRmseMv = NaN(n, 1);
medianRmseMv = NaN(n, 1);

for i = 1:n
    method(i) = string(spec(i).method);
    eo = T.(spec(i).eo);
    freq = T.(spec(i).freq);
    amp = T.(spec(i).amp);
    rmse = T.(spec(i).rmse);
    windowCount(i) = sum(isfinite(eo));
    eoDistribution(i) = eo_distribution_local(eo);
    dominantEO(i) = mode_finite_local(eo);
    meanFrequencyHz(i) = mean(freq, 'omitnan');
    stdFrequencyHz(i) = std(freq, 'omitnan');
    meanAmplitudeMm(i) = mean(amp, 'omitnan');
    meanRmseMv(i) = mean(rmse, 'omitnan');
    medianRmseMv(i) = median(rmse, 'omitnan');
end

Summary = table(method, windowCount, eoDistribution, dominantEO, ...
    meanFrequencyHz, stdFrequencyHz, meanAmplitudeMm, meanRmseMv, medianRmseMv, ...
    'VariableNames', {'method','window_count','EO_distribution','dominant_EO', ...
    'mean_frequency_hz','std_frequency_hz','mean_amplitude_mm','mean_rmse_mV','median_rmse_mV'});
end

function plot_method_trend_local(T, figFile, dataset, bladeTag)
fig = figure('Name', sprintf('%s %s Step07J method comparison', dataset, bladeTag), ...
    'Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 18, 14]);
tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
colors = struct('direct', [0.10, 0.10, 0.10], 'fixed', [0.10, 0.32, 0.72], ...
    'gap', [0.82, 0.22, 0.18], 'tilt', [0.10, 0.55, 0.25]);

nexttile; hold on; grid on; box on;
plot(T.window_id, T.direct_EO, '-.', 'Color', colors.direct, 'LineWidth', 1.0, 'DisplayName', 'direct');
plot(T.window_id, T.fixed_EO, '-o', 'Color', colors.fixed, 'LineWidth', 1.0, 'MarkerSize', 3.8, 'DisplayName', 'fixed');
plot(T.window_id, T.gap_EO, '-s', 'Color', colors.gap, 'LineWidth', 1.0, 'MarkerSize', 3.8, 'DisplayName', 'gap only');
plot(T.window_id, T.tilt_EO, '-^', 'Color', colors.tilt, 'LineWidth', 1.0, 'MarkerSize', 3.8, 'DisplayName', 'gap + tilt');
ylabel('EO'); title(sprintf('%s %s order comparison', dataset, bladeTag)); legend('Location', 'best');

nexttile; hold on; grid on; box on;
plot(T.window_id, T.direct_amplitude_mm, '-.', 'Color', colors.direct, 'LineWidth', 1.0, 'DisplayName', 'direct');
plot(T.window_id, T.fixed_amplitude_mm, '-o', 'Color', colors.fixed, 'LineWidth', 1.0, 'MarkerSize', 3.8, 'DisplayName', 'fixed');
plot(T.window_id, T.gap_amplitude_mm, '-s', 'Color', colors.gap, 'LineWidth', 1.0, 'MarkerSize', 3.8, 'DisplayName', 'gap only');
plot(T.window_id, T.tilt_amplitude_mm, '-^', 'Color', colors.tilt, 'LineWidth', 1.0, 'MarkerSize', 3.8, 'DisplayName', 'gap + tilt');
ylabel('Amplitude (mm)'); title('Identified vibration amplitude');

nexttile; hold on; grid on; box on;
plot(T.window_id, T.direct_rmse_mV, '-.', 'Color', colors.direct, 'LineWidth', 1.0, 'DisplayName', 'direct');
plot(T.window_id, T.fixed_rmse_mV, '-o', 'Color', colors.fixed, 'LineWidth', 1.0, 'MarkerSize', 3.8, 'DisplayName', 'fixed');
plot(T.window_id, T.gap_rmse_mV, '-s', 'Color', colors.gap, 'LineWidth', 1.0, 'MarkerSize', 3.8, 'DisplayName', 'gap only');
plot(T.window_id, T.tilt_rmse_mV, '-^', 'Color', colors.tilt, 'LineWidth', 1.0, 'MarkerSize', 3.8, 'DisplayName', 'gap + tilt');
xlabel('Step07J window'); ylabel('RMSE (mV)'); title('Final waveform objective');

exportgraphics(fig, figFile, 'Resolution', 300);
end

function text = eo_distribution_local(v)
v = v(isfinite(v));
if isempty(v)
    text = "";
    return;
end
vals = unique(v(:)).';
parts = strings(1, numel(vals));
for i = 1:numel(vals)
    parts(i) = sprintf('EO%d=%d', vals(i), sum(v == vals(i)));
end
text = strjoin(parts, ', ');
end

function y = mode_finite_local(v)
v = v(isfinite(v));
if isempty(v)
    y = NaN;
    return;
end
vals = unique(v(:));
counts = zeros(size(vals));
for i = 1:numel(vals)
    counts(i) = sum(v == vals(i));
end
[~, idx] = max(counts);
y = vals(idx);
end

function tf = parse_logical_env_local(name, defaultValue)
txt = lower(strtrim(getenv(name)));
if isempty(txt)
    tf = defaultValue;
elseif any(strcmp(txt, {'1','true','yes','on'}))
    tf = true;
elseif any(strcmp(txt, {'0','false','no','off'}))
    tf = false;
else
    error('%s must be true/false.', name);
end
end

function restore_step07j_env_local(previousRunMode, previousSuffix)
setenv('STEP07J_RUN_MODE', previousRunMode);
setenv('STEP07J_RESULT_SUFFIX', previousSuffix);
end
