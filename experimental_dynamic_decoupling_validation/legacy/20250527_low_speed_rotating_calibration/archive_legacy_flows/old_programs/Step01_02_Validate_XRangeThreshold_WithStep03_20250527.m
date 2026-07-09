function Step01_02_Validate_XRangeThreshold_WithStep03_20250527()
%% Step01_02_Validate_XRangeThreshold_WithStep03_20250527
% Validate x-range gradient thresholds by rebuilding templates and running
% Step03_04 identification on each candidate threshold.

clc; close all;

route_dir = fileparts(mfilename('fullpath'));
plan_file = fullfile(route_dir, 'output', 'stable_window_plan', ...
    'Step01_CoverageFirstStableWindowPlan_B1_S136_20250527.csv');
if exist(plan_file, 'file') ~= 2
    error('Required 30-lap plan not found: %s', plan_file);
end

gradient_ratios = parse_numeric_list_env_local('STEP01_02_GRADIENT_RATIOS', [0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50]);
amplitude_ratio = parse_nonnegative_numeric_env_local('STEP01_02_AMPLITUDE_RATIO', 0);
max_windows = parse_positive_integer_env_local('STEP01_02_MAX_WINDOWS', inf);
if isfinite(max_windows)
    max_windows_text = sprintf('%d', max_windows);
else
    max_windows_text = '';
end

output_dir = fullfile(route_dir, 'output');
summary_dir = fullfile(output_dir, 'stable_window_plan');
figure_dir = fullfile(output_dir, 'figures');
result_summary_csv = fullfile(summary_dir, 'Step01_02_XRangeThresholdStep03Summary_B1_S136_20250527.csv');

fprintf('\n=== Step01_02: x-range threshold validation with Step03 ===\n');
fprintf('Gradient ratios: %s\n', mat2str(gradient_ratios));
fprintf('Amplitude ratio: %.3f\n', amplitude_ratio);

rows = table();
cleanup_obj = onCleanup(@clear_candidate_env_local);
for gr = gradient_ratios(:).'
    tag = sprintf('XThr_g%03d_a%03d', round(1000 * gr), round(1000 * amplitude_ratio));
    fprintf('\n--- Candidate %s ---\n', tag);

    setenv('STEP01_STABLE_WINDOW_PLAN', plan_file);
    setenv('STEP01_TEMPLATE_SUFFIX', tag);
    setenv('STEP01_CENTER_MODE', 'sgfit');
    setenv('STEP01_XRANGE_MODE', 'threshold');
    setenv('STEP01_XRANGE_GRADIENT_MIN_RATIO', sprintf('%.8g', gr));
    setenv('STEP01_XRANGE_AMPLITUDE_MIN_RATIO', sprintf('%.8g', amplitude_ratio));
    setenv('STEP01_XRANGE_MIN_HALF_WIDTH_MM', '2.5');
    setenv('STEP01_XRANGE_MAX_HALF_WIDTH_MM', '4.2');
    run_step01_script_local();

    template_file = fullfile(route_dir, 'output', 'templates', ...
        sprintf('Template_LowSpeedRotating_B1_S136_%s_20250527.mat', tag));

    setenv('STEP03_TEMPLATE_FILE', template_file);
    setenv('STEP03_RESULT_SUFFIX', tag);
    setenv('STEP03D_DYNAMIC_EFFECTIVE_MODE', 'gradient');
    setenv('STEP03D_DOMAIN_SELECTION_MODE', 'hard');
    setenv('STEP03_04_MAX_WINDOWS', max_windows_text);
    run_step03_script_local();

    result_file = fullfile(route_dir, 'output', 'identification', ...
        sprintf('Result_Step03_04_FirstOrderVPAdaptive_B1_S13_%s_20250527.mat', tag));
    loaded_template = load(template_file, 'Template');
    loaded_result = load(result_file, 'Result');
    row = summarize_candidate_local(loaded_template.Template, loaded_result.Result, gr, amplitude_ratio, tag, template_file, result_file);
    rows = [rows; row]; %#ok<AGROW>
end

clear_candidate_env_local();
rows = score_summary_local(rows);
writetable(rows, result_summary_csv);
fprintf('\nSaved summary: %s\n', result_summary_csv);
disp(rows(:, {'Tag','GradientRatio','MeanHalfWidthMM','DominantEO','EOConsistency', ...
    'MeanFreqHz','MeanRMSE','MeanClampFraction','MeanPointCount','Score'}));

plot_threshold_summary_local(rows, figure_dir);
end

%% Local functions
function run_step01_script_local()
evalin('base', 'Step01_Build_Rotating_Template_20250527;');
end

function run_step03_script_local()
evalin('base', 'Step03_04_Run_FirstOrderVPAdaptive_20250527;');
end

function vals = parse_numeric_list_env_local(name, default_value)
raw = strtrim(getenv(name));
if isempty(raw)
    vals = default_value;
    return;
end
vals = sscanf(raw, '%f').';
if isempty(vals)
    error('%s must contain numeric values separated by spaces.', name);
end
end

function v = parse_nonnegative_numeric_env_local(name, default_value)
raw = strtrim(getenv(name));
if isempty(raw)
    v = default_value;
    return;
end
v = str2double(raw);
if ~isfinite(v) || v < 0
    error('%s must be a nonnegative numeric value.', name);
end
end

function v = parse_positive_integer_env_local(name, default_value)
raw = strtrim(getenv(name));
if isempty(raw)
    v = default_value;
    return;
end
tmp = str2double(raw);
if ~isfinite(tmp) || tmp <= 0
    error('%s must be a positive integer.', name);
end
v = max(1, floor(tmp));
end

function row = summarize_candidate_local(Template, Result, gr, ar, tag, template_file, result_file)
domains = vertcat(Template.Sensor.x_domain);
half_widths = max(abs(domains), [], 2);
trend = Result.Trend;
eo = trend.EO_id;
dominant_eo = mode(eo(isfinite(eo)));
eo_consistency = mean(eo == dominant_eo, 'omitnan');
clamp = nan(height(trend), 1);
point_count = nan(height(trend), 1);
effective_count = nan(height(trend), 1);
for iw = 1:numel(Result.WindowResult)
    if ~isfield(Result.WindowResult(iw), 'bundle')
        continue;
    end
    meta = Result.WindowResult(iw).bundle.window_meta;
    point_count(iw) = sum([meta.valid_points]);
    if isfield(meta, 'dynamic_effective_points')
        effective_count(iw) = sum([meta.dynamic_effective_points]);
    end
    if isfield(Result.WindowResult(iw).Result, 'Coverage')
        clamp(iw) = Result.WindowResult(iw).Result.Coverage.clamp_fraction;
    end
end
row = table();
row.Tag = string(tag);
row.GradientRatio = gr;
row.AmplitudeRatio = ar;
row.MeanHalfWidthMM = mean(half_widths, 'omitnan');
row.MinHalfWidthMM = min(half_widths, [], 'omitnan');
row.MaxHalfWidthMM = max(half_widths, [], 'omitnan');
row.DominantEO = dominant_eo;
row.EOConsistency = eo_consistency;
row.MeanFreqHz = mean(trend.fn_id, 'omitnan');
row.StdFreqHz = std(trend.fn_id, 'omitnan');
row.MeanRMSE = mean(trend.weighted_voltage_rmse, 'omitnan');
row.MedianRMSE = median(trend.weighted_voltage_rmse, 'omitnan');
row.MeanClampFraction = mean(clamp, 'omitnan');
row.MeanPointCount = mean(point_count, 'omitnan');
row.MeanDynamicEffectiveCount = mean(effective_count, 'omitnan');
row.TemplateFile = string(template_file);
row.ResultFile = string(result_file);
end

function rows = score_summary_local(rows)
rmse_ref = max(median(rows.MeanRMSE, 'omitnan'), 0.02);
freq_ref = max(median(rows.StdFreqHz, 'omitnan'), 1e-6);
width_target = 3.7;
rows.Score = 3.0 * (1 - rows.EOConsistency) + ...
    1.5 * rows.MeanRMSE ./ rmse_ref + ...
    0.8 * rows.StdFreqHz ./ freq_ref + ...
    1.0 * rows.MeanClampFraction + ...
    0.4 * abs(rows.MeanHalfWidthMM - width_target) ./ width_target;
end

function clear_candidate_env_local()
vars = {'STEP01_STABLE_WINDOW_PLAN','STEP01_TEMPLATE_SUFFIX','STEP01_CENTER_MODE', ...
    'STEP01_XRANGE_MODE','STEP01_XRANGE_GRADIENT_MIN_RATIO','STEP01_XRANGE_AMPLITUDE_MIN_RATIO', ...
    'STEP01_XRANGE_MIN_HALF_WIDTH_MM','STEP01_XRANGE_MAX_HALF_WIDTH_MM', ...
    'STEP03_TEMPLATE_FILE','STEP03_RESULT_SUFFIX','STEP03D_DYNAMIC_EFFECTIVE_MODE', ...
    'STEP03D_DOMAIN_SELECTION_MODE','STEP03_04_MAX_WINDOWS'};
for i = 1:numel(vars)
    setenv(vars{i}, '');
end
end

function plot_threshold_summary_local(rows, figure_dir)
fig = figure('Name', 'X-range threshold Step03 validation', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 18, 11]);
tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot(rows.GradientRatio, rows.MeanHalfWidthMM, 'o-', 'LineWidth', 1.2);
xlabel('Gradient threshold');
ylabel('Mean half width (mm)');
title('Template width');
style_axes_local();

nexttile;
plot(rows.GradientRatio, rows.EOConsistency, 'o-', 'LineWidth', 1.2);
xlabel('Gradient threshold');
ylabel('EO consistency');
title('EO stability');
style_axes_local();

nexttile;
yyaxis left;
plot(rows.GradientRatio, rows.MeanRMSE, 'o-', 'LineWidth', 1.2);
ylabel('Mean RMSE (V)');
yyaxis right;
plot(rows.GradientRatio, rows.MeanPointCount, 's-', 'LineWidth', 1.2);
ylabel('Mean valid point count');
xlabel('Gradient threshold');
title('Fit and information');
style_axes_local();

nexttile;
plot(rows.GradientRatio, rows.Score, 'o-', 'LineWidth', 1.2); hold on;
[~, best_idx] = min(rows.Score);
plot(rows.GradientRatio(best_idx), rows.Score(best_idx), 'kp', 'MarkerFaceColor', 'y', 'MarkerSize', 12);
xlabel('Gradient threshold');
ylabel('Score');
title('Combined score');
style_axes_local();

exportgraphics(fig, fullfile(figure_dir, 'Step01_02_XRangeThresholdStep03Summary_B1_S136_20250527.png'), 'Resolution', 300);
exportgraphics(fig, fullfile(figure_dir, 'Step01_02_XRangeThresholdStep03Summary_B1_S136_20250527.pdf'), 'ContentType', 'vector');
end

function style_axes_local()
set(gca, 'FontName', 'Times New Roman', 'FontSize', 8.5, 'TickDir', 'in', 'Box', 'on');
end
