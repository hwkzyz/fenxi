function Analyze_PriorStaticXRangeSelection_20250527_20251222
%% Analyze_PriorStaticXRangeSelection_20250527_20251222
% Prior static x-domain selection without using final identified vibration
% parameters. The pre-identification displacement envelope is estimated as
% u_app = x_dynamic - T^{-1}(V_dynamic), then used to score gradient-ratio
% candidate static domains before Step03 waveform identification.

clc;

root_dir = fileparts(mfilename('fullpath'));

cfg = struct();
cfg.gradient_ratios = parse_numeric_list_env_local('STEP01_PRIOR_XRANGE_RATIOS', ...
    [0.02, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50]);
cfg.u_app_quantile = parse_nonnegative_numeric_env_local('STEP01_PRIOR_U_APP_QUANTILE', 98);
cfg.u_app_margin_mm = parse_nonnegative_numeric_env_local('STEP01_PRIOR_U_APP_MARGIN_MM', 0.05);
cfg.u_app_max_abs_mm = parse_nonnegative_numeric_env_local('STEP01_PRIOR_U_APP_MAX_ABS_MM', 1.50);
cfg.u_app_template_gradient_ratio = parse_nonnegative_numeric_env_local('STEP01_PRIOR_U_APP_TEMPLATE_GRADIENT_RATIO', 0.08);
cfg.edge_margin_mm = parse_nonnegative_numeric_env_local('STEP01_PRIOR_EDGE_MARGIN_MM', 0.20);
cfg.coverage_accept_mm = parse_nonnegative_numeric_env_local('STEP01_PRIOR_COVERAGE_ACCEPT_MM', 0.12);
cfg.edge_accept_fraction = parse_nonnegative_numeric_env_local('STEP01_PRIOR_EDGE_ACCEPT_FRACTION', 0.05);
cfg.run_step03 = parse_logical_env_local('STEP01_PRIOR_RUN_STEP03', true);

case_list = build_case_list_local(root_dir);
all_rows = table();
selected_rows = table();
validation_rows = table();

fprintf('\n=== Prior static x-domain selection ===\n');
fprintf('Candidate gradient ratios: %s\n', mat2str(cfg.gradient_ratios));
fprintf('u_app envelope: P%.1f(abs(u_app - median)) + %.3f mm\n', ...
    cfg.u_app_quantile, cfg.u_app_margin_mm);
fprintf('u_app prior mask: |u_app| <= %.3f mm, template gradient >= %.3f of max\n', ...
    cfg.u_app_max_abs_mm, cfg.u_app_template_gradient_ratio);

for ic = 1:numel(case_list)
    case_cfg = case_list(ic);
    fprintf('\n--- %s ---\n', case_cfg.tag);
    ensure_candidate_templates_local(case_cfg, cfg.gradient_ratios);

    envelope_template_file = candidate_template_file_local(case_cfg, min(cfg.gradient_ratios));
    envelope = estimate_prior_envelope_local(case_cfg, envelope_template_file, cfg);
    rows = score_case_candidates_local(case_cfg, cfg, envelope);
    rows = choose_candidate_local(rows, cfg);
    all_rows = [all_rows; rows]; %#ok<AGROW>

    chosen = rows(rows.IsSelected, :);
    selected_rows = [selected_rows; chosen]; %#ok<AGROW>
    fprintf('Selected ratio %.3f (%s): score %.4g, coverage %.4f mm, edge %.3f%%, mean U_pre %.4f mm\n', ...
        chosen.GradientRatio, chosen.SelectionReason, chosen.Score, ...
        chosen.CoverageRMSMM, 100 * chosen.EdgeRiskFraction, chosen.MeanUPreMM);

    if cfg.run_step03
        validation_rows = [validation_rows; run_selected_step03_local(case_cfg, chosen)]; %#ok<AGROW>
    end
end

output_csv = fullfile(root_dir, 'PriorStaticXRangeSelection_Candidates.csv');
selected_csv = fullfile(root_dir, 'PriorStaticXRangeSelection_Selected.csv');
writetable(all_rows, output_csv);
writetable(selected_rows, selected_csv);

if ~isempty(validation_rows)
    validation_csv = fullfile(root_dir, 'PriorStaticXRangeSelection_Step03Validation.csv');
    writetable(validation_rows, validation_csv);
else
    validation_csv = "";
end

figure_dir = fullfile(root_dir, 'figures');
if exist(figure_dir, 'dir') ~= 7; mkdir(figure_dir); end
figure_file = fullfile(figure_dir, 'PriorStaticXRangeSelection_Candidates.png');
plot_prior_selection_local(all_rows, figure_file, cfg);

fprintf('\nSaved candidate table: %s\n', output_csv);
fprintf('Saved selected table:  %s\n', selected_csv);
if strlength(validation_csv) > 0
    fprintf('Saved Step03 validation table: %s\n', validation_csv);
end
fprintf('Saved figure: %s\n', figure_file);
if ~isempty(validation_rows)
    disp(validation_rows);
end

end

function case_list = build_case_list_local(root_dir)
case_list = struct([]);

case_list(1).tag = "20250527";
case_list(1).folder = fullfile(root_dir, '20250527_low_speed_rotating_calibration');
case_list(1).target_blade = 1;
case_list(1).sensor_tag = 'S136';
case_list(1).analysis_sensors = '1 3 6';
case_list(1).stable_plan_file = fullfile(case_list(1).folder, 'output', 'stable_window_plan', ...
    'Step01_CoverageFirstStableWindowPlan_B1_S136_20250527.csv');
case_list(1).step01_script = fullfile(case_list(1).folder, 'Step01_Build_Rotating_Template_20250527.m');
case_list(1).step03_wrapper = fullfile(case_list(1).folder, 'Step03_Main_Run_GradientXRange030_Identification_20250527.m');
case_list(1).dynamic_map_file = fullfile(case_list(1).folder, 'output', 'dynamic_maps', ...
    'DynamicMap_B1_S136_SlidingWindows_20250527.mat');
case_list(1).date_tag = '20250527';

case_list(2).tag = "20251222";
case_list(2).folder = fullfile(root_dir, '20251222_low_speed_rotating_calibration');
case_list(2).target_blade = 1;
case_list(2).sensor_tag = 'S123';
case_list(2).analysis_sensors = '1 2 3';
case_list(2).stable_plan_file = fullfile(case_list(2).folder, 'output', 'stable_window_plan', ...
    'Step01_CoverageFirstStableWindowPlan_B1_S123_20251222.csv');
case_list(2).step01_script = fullfile(case_list(2).folder, 'Step01_Build_Rotating_Template_20251222.m');
case_list(2).step03_wrapper = fullfile(case_list(2).folder, 'Step03_Main_Run_GradientXRange030_Identification_20251222.m');
case_list(2).dynamic_map_file = fullfile(case_list(2).folder, 'output', 'dynamic_maps', ...
    'DynamicMap_B1_S123_SlidingWindows_20251222.mat');
case_list(2).date_tag = '20251222';
end

function ensure_candidate_templates_local(case_cfg, ratios)
for ratio = ratios(:).'
    template_file = candidate_template_file_local(case_cfg, ratio);
    if exist(template_file, 'file') == 2
        fprintf('Reuse template ratio %.3f: %s\n', ratio, template_file);
        continue;
    end

    suffix = candidate_suffix_local(ratio);
    fprintf('Build template ratio %.3f: %s\n', ratio, suffix);
    old_env = capture_env_local({'STEP01_STABLE_WINDOW_PLAN','STEP01_TEMPLATE_SUFFIX', ...
        'STEP01_CENTER_MODE','STEP01_XRANGE_MODE','STEP01_XRANGE_GRADIENT_MIN_RATIO', ...
        'STEP01_XRANGE_AMPLITUDE_MIN_RATIO','STEP01_XRANGE_MIN_HALF_WIDTH_MM', ...
        'STEP01_XRANGE_MAX_HALF_WIDTH_MM'});
    cleaner = onCleanup(@() restore_env_local(old_env));

    setenv('STEP01_STABLE_WINDOW_PLAN', case_cfg.stable_plan_file);
    setenv('STEP01_TEMPLATE_SUFFIX', suffix);
    setenv('STEP01_CENTER_MODE', 'sgfit');
    setenv('STEP01_XRANGE_MODE', 'threshold');
    setenv('STEP01_XRANGE_GRADIENT_MIN_RATIO', sprintf('%.8g', ratio));
    setenv('STEP01_XRANGE_AMPLITUDE_MIN_RATIO', '0');
    setenv('STEP01_XRANGE_MIN_HALF_WIDTH_MM', '2.5');
    setenv('STEP01_XRANGE_MAX_HALF_WIDTH_MM', '4.2');
    run_script_in_base_local(case_cfg.step01_script);
    clear cleaner;
end
end

function envelope = estimate_prior_envelope_local(case_cfg, template_file, cfg)
loaded_template = load(template_file, 'Template');
loaded_dynamic = load(case_cfg.dynamic_map_file, 'DynamicMap');
Template = loaded_template.Template;
DynamicMap = loaded_dynamic.DynamicMap;

sensor_ids = [Template.Sensor.sensor_id];
stats = repmat(struct('sensor_id', NaN, 'u_pre_mm', NaN, 'u_app_median_mm', NaN, ...
    'u_app_q_mm', NaN, 'point_count', NaN, 'x_p01', NaN, 'x_p99', NaN, ...
    'x_eff', [], 'w_eff', []), numel(sensor_ids), 1);

for is = 1:numel(sensor_ids)
    sid = sensor_ids(is);
    Tpl = Template.Sensor(is);
    all_u = [];
    all_w = [];
    all_x = [];
    for iw = 1:numel(DynamicMap.Window)
        S = DynamicMap.Window(iw).Sensor([DynamicMap.Window(iw).Sensor.sensor_id] == sid);
        if isempty(S) || isempty(S.x_rel)
            continue;
        end
        x_dyn = S.x_rel(:);
        v_dyn = S.V(:);
        x_inv = inverse_template_voltage_local(Tpl, v_dyn, x_dyn);
        g_inv = interp1(Tpl.x_grid(:), abs(Tpl.dv_dx(:)), x_inv, 'linear', 0);
        v_domain = Tpl.x_grid(:) >= Tpl.x_domain(1) & Tpl.x_grid(:) <= Tpl.x_domain(2);
        v_min = min(Tpl.v_grid(v_domain), [], 'omitnan');
        v_max = max(Tpl.v_grid(v_domain), [], 'omitnan');
        valid = isfinite(x_dyn) & isfinite(x_inv) & isfinite(v_dyn);
        if isfield(S, 'W') && ~isempty(S.W)
            w = S.W(:);
        else
            w = ones(size(x_dyn));
        end
        u_app = x_dyn - x_inv;
        valid = valid & isfinite(w) & w > 0 & ...
            v_dyn >= v_min & v_dyn <= v_max & ...
            g_inv >= cfg.u_app_template_gradient_ratio * max(abs(Tpl.dv_dx(:)), [], 'omitnan') & ...
            abs(u_app) <= cfg.u_app_max_abs_mm;
        all_u = [all_u; u_app(valid)]; %#ok<AGROW>
        all_w = [all_w; w(valid)]; %#ok<AGROW>
        all_x = [all_x; x_dyn(valid)]; %#ok<AGROW>
    end
    u_med = weighted_quantile_local(all_u, all_w, 0.50);
    u_abs = abs(all_u - u_med);
    u_q = weighted_quantile_local(u_abs, all_w, cfg.u_app_quantile / 100);
    stats(is).sensor_id = sid;
    stats(is).u_pre_mm = u_q + cfg.u_app_margin_mm;
    stats(is).u_app_median_mm = u_med;
    stats(is).u_app_q_mm = u_q;
    stats(is).point_count = numel(all_u);
    stats(is).x_p01 = prctile(all_x, 1);
    stats(is).x_p99 = prctile(all_x, 99);
    stats(is).x_eff = all_x;
    stats(is).w_eff = all_w;
end

envelope = struct();
envelope.Sensor = stats;
envelope.TemplateFile = template_file;
envelope.DynamicMapFile = case_cfg.dynamic_map_file;
end

function rows = score_case_candidates_local(case_cfg, cfg, envelope)
rows = table();

for ratio = cfg.gradient_ratios(:).'
    template_file = candidate_template_file_local(case_cfg, ratio);
    loaded_template = load(template_file, 'Template');
    Template = loaded_template.Template;

    coverage_num = 0;
    coverage_den = 0;
    edge_num = 0;
    edge_den = 0;
    widths = nan(numel(Template.Sensor), 1);
    mean_grad = nan(numel(Template.Sensor), 1);
    grad_balance = nan(numel(Template.Sensor), 1);
    u_pre = nan(numel(Template.Sensor), 1);
    point_count = 0;

    for is = 1:numel(Template.Sensor)
        Tpl = Template.Sensor(is);
        sid = Tpl.sensor_id;
        env_idx = find([envelope.Sensor.sensor_id] == sid, 1);
        if isempty(env_idx)
            continue;
        end
        U_pre = envelope.Sensor(env_idx).u_pre_mm;
        u_pre(is) = U_pre;
        x_domain = Tpl.x_domain(:).';
        widths(is) = diff(x_domain);

        in_domain = Tpl.x_grid >= x_domain(1) & Tpl.x_grid <= x_domain(2);
        g = abs(Tpl.dv_dx(:));
        mean_grad(is) = mean(g(in_domain), 'omitnan') / max(mean(g, 'omitnan'), eps);
        left_g = sum(g(in_domain & Tpl.x_grid(:) <= mean(x_domain)), 'omitnan');
        right_g = sum(g(in_domain & Tpl.x_grid(:) > mean(x_domain)), 'omitnan');
        grad_balance(is) = abs(left_g - right_g) / max(left_g + right_g, eps);

        x_dyn = envelope.Sensor(env_idx).x_eff(:);
        w = envelope.Sensor(env_idx).w_eff(:);
        left_over = max(0, x_domain(1) - (x_dyn - U_pre));
        right_over = max(0, (x_dyn + U_pre) - x_domain(2));
        coverage_num = coverage_num + sum(w .* (left_over.^2 + right_over.^2), 'omitnan');
        coverage_den = coverage_den + sum(w, 'omitnan');

        dist_to_edge = min(x_dyn - x_domain(1), x_domain(2) - x_dyn);
        edge_num = edge_num + sum(w .* (dist_to_edge < cfg.edge_margin_mm), 'omitnan');
        edge_den = edge_den + sum(w, 'omitnan');
        point_count = point_count + numel(x_dyn);
    end

    coverage_rms = sqrt(coverage_num / max(coverage_den, eps));
    edge_fraction = edge_num / max(edge_den, eps);
    row = table(string(case_cfg.tag), ratio, string(candidate_suffix_local(ratio)), ...
        mean(widths, 'omitnan'), min(widths), max(widths), ...
        mean(u_pre, 'omitnan'), coverage_rms, edge_fraction, ...
        mean(mean_grad, 'omitnan'), mean(grad_balance, 'omitnan'), point_count, ...
        string(template_file), ...
        'VariableNames', {'CaseTag','GradientRatio','CandidateSuffix', ...
        'MeanWidthMM','MinWidthMM','MaxWidthMM','MeanUPreMM','CoverageRMSMM', ...
        'EdgeRiskFraction','MeanGradientGain','GradientImbalance','DynamicPointCount', ...
        'TemplateFile'});
    rows = [rows; row]; %#ok<AGROW>
end

rows.Score = score_candidates_local(rows);
end

function rows = choose_candidate_local(rows, cfg)
rows.IsSelected = false(height(rows), 1);
rows.SelectionReason = strings(height(rows), 1);
acceptable = rows.CoverageRMSMM <= cfg.coverage_accept_mm & ...
    rows.EdgeRiskFraction <= cfg.edge_accept_fraction;

if any(acceptable)
    sub = rows(acceptable, :);
    max_ratio = max(sub.GradientRatio, [], 'omitnan');
    sub2 = sub(sub.GradientRatio == max_ratio, :);
    [~, idx_sub2] = min(sub2.Score);
    idx = find(acceptable);
    idx_max = idx(rows.GradientRatio(acceptable) == max_ratio);
    idx = idx_max(idx_sub2);
    rows.SelectionReason(idx) = "accepted_highest_gradient_ratio";
else
    [~, idx] = min(rows.Score);
    rows.SelectionReason(idx) = "fallback_min_score";
end
rows.IsSelected(idx) = true;
end

function validation = run_selected_step03_local(case_cfg, chosen)
ratio = chosen.GradientRatio;
suffix = char(chosen.CandidateSuffix);
result_suffix = ['Main_', suffix];

fprintf('Run Step03 validation for %s ratio %.3f...\n', case_cfg.tag, ratio);
old_env = capture_env_local({'STEP01_CURRENT_TEMPLATE_SUFFIX','STEP03_CURRENT_RESULT_SUFFIX', ...
    'STEP03_TEMPLATE_FILE','STEP03_RESULT_SUFFIX','STEP03_ANALYSIS_SENSORS', ...
    'STEP03_MAIN_TOP_K_EO','STEP03D_DYNAMIC_EFFECTIVE_MODE','STEP03D_PULSE_MODE', ...
    'STEP03D_DOMAIN_SELECTION_MODE','STEP03D_DOMAIN_SOFT_MARGIN_MM', ...
    'STEP03D_OVERSHOOT_PENALTY_WEIGHT'});
cleaner = onCleanup(@() restore_env_local(old_env));

setenv('STEP01_CURRENT_TEMPLATE_SUFFIX', suffix);
setenv('STEP03_CURRENT_RESULT_SUFFIX', result_suffix);
run_script_in_base_local(case_cfg.step03_wrapper);

clear cleaner;

result_file = fullfile(case_cfg.folder, 'output', 'identification', sprintf( ...
    'Result_Step03_Main_VPTop3SynchronousWaveform_B%d_%s_%s_%s.mat', ...
    case_cfg.target_blade, case_cfg.sensor_tag, result_suffix, case_cfg.date_tag));
loaded = load(result_file, 'Result');
T = loaded.Result.Trend;
dominant_eo = mode(T.EO_id);
validation = table(string(case_cfg.tag), ratio, string(suffix), string(result_file), ...
    dominant_eo, mean(T.EO_id == dominant_eo, 'omitnan'), ...
    mean(T.fn_id, 'omitnan'), median(T.fn_id, 'omitnan'), ...
    mean(T.weighted_voltage_rmse, 'omitnan'), median(T.weighted_voltage_rmse, 'omitnan'), ...
    mean(T.A_id, 'omitnan'), ...
    'VariableNames', {'CaseTag','GradientRatio','CandidateSuffix','ResultFile', ...
    'DominantEO','EOConsistency','MeanFrequencyHz','MedianFrequencyHz', ...
    'MeanRMSE','MedianRMSE','MeanAmplitudeMM'});
end

function score = score_candidates_local(rows)
coverage = normalize_robust_local(rows.CoverageRMSMM);
edge = normalize_robust_local(rows.EdgeRiskFraction);
width = normalize_robust_local(rows.MeanWidthMM);
grad_loss = normalize_robust_local(1 ./ max(rows.MeanGradientGain, eps));
imbalance = normalize_robust_local(rows.GradientImbalance);
score = 1.40 * coverage + 0.80 * edge + 0.25 * width + 0.35 * grad_loss + 0.20 * imbalance;
end

function y = normalize_robust_local(x)
x = x(:);
lo = min(x, [], 'omitnan');
hi = max(x, [], 'omitnan');
if ~isfinite(lo) || ~isfinite(hi) || hi <= lo
    y = zeros(size(x));
else
    y = (x - lo) ./ (hi - lo);
end
end

function x_inv = inverse_template_voltage_local(Tpl, v_query, x_hint)
x = Tpl.x_grid(:);
v = Tpl.v_grid(:);
valid = isfinite(x) & isfinite(v);
x = x(valid);
v = v(valid);
[~, peak_idx] = max(abs(v - Tpl.baseline));
x_peak = x(peak_idx);

left = x <= x_peak;
right = x >= x_peak;
x_inv_left = inverse_branch_local(x(left), v(left), v_query);
x_inv_right = inverse_branch_local(x(right), v(right), v_query);
x_inv = x_inv_right;
x_inv(x_hint <= x_peak) = x_inv_left(x_hint <= x_peak);

bad = ~isfinite(x_inv);
if any(bad)
    x_inv_all = nearest_voltage_x_local(x, v, v_query(bad));
    x_inv(bad) = x_inv_all;
end
end

function xq = inverse_branch_local(x, v, vq)
[v_sort, order] = sort(v(:));
x_sort = x(order);
[v_unique, ia] = unique(v_sort, 'stable');
x_unique = x_sort(ia);
if numel(v_unique) < 2 || range(v_unique) <= eps
    xq = nan(size(vq));
    return;
end
xq = interp1(v_unique, x_unique, vq(:), 'linear', NaN);
end

function x_near = nearest_voltage_x_local(x, v, vq)
x_near = nan(size(vq));
for i = 1:numel(vq)
    [~, idx] = min(abs(v - vq(i)));
    if ~isempty(idx)
        x_near(i) = x(idx);
    end
end
end

function q = weighted_quantile_local(x, w, p)
valid = isfinite(x) & isfinite(w) & w > 0;
x = x(valid);
w = w(valid);
if isempty(x)
    q = NaN;
    return;
end
[x, order] = sort(x(:));
w = w(order);
cw = cumsum(w) ./ sum(w);
q = interp1(cw, x, min(max(p, 0), 1), 'linear', 'extrap');
end

function template_file = candidate_template_file_local(case_cfg, ratio)
template_file = fullfile(case_cfg.folder, 'output', 'templates', sprintf( ...
    'Template_LowSpeedRotating_B%d_%s_%s_%s.mat', ...
    case_cfg.target_blade, case_cfg.sensor_tag, candidate_suffix_local(ratio), case_cfg.date_tag));
end

function suffix = candidate_suffix_local(ratio)
suffix = sprintf('PriorXRange_g%03d_a000', round(1000 * ratio));
end

function plot_prior_selection_local(rows, output_png, cfg)
fig = figure('Name', 'Prior static x-domain candidate selection', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 17, 13]);
tiledlayout(fig, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
cases = unique(rows.CaseTag, 'stable');
colors = lines(numel(cases));

nexttile; hold on;
for ic = 1:numel(cases)
    M = rows(rows.CaseTag == cases(ic), :);
    plot(M.GradientRatio, M.CoverageRMSMM, 'o-', 'Color', colors(ic,:), ...
        'LineWidth', 1.2, 'MarkerSize', 4.5, 'DisplayName', cases(ic));
    scatter(M.GradientRatio(M.IsSelected), M.CoverageRMSMM(M.IsSelected), ...
        55, colors(ic,:), 'filled', 'HandleVisibility', 'off');
end
yline(cfg.coverage_accept_mm, 'k:', 'LineWidth', 1.0, 'DisplayName', 'accept limit');
ylabel('Coverage RMS (mm)');
title('Prior envelope coverage risk');
style_axes_local();
legend('Location', 'best');

nexttile; hold on;
for ic = 1:numel(cases)
    M = rows(rows.CaseTag == cases(ic), :);
    plot(M.GradientRatio, 100 * M.EdgeRiskFraction, 's-', 'Color', colors(ic,:), ...
        'LineWidth', 1.2, 'MarkerSize', 4.5, 'DisplayName', cases(ic));
    scatter(M.GradientRatio(M.IsSelected), 100 * M.EdgeRiskFraction(M.IsSelected), ...
        55, colors(ic,:), 'filled', 'HandleVisibility', 'off');
end
yline(100 * cfg.edge_accept_fraction, 'k:', 'LineWidth', 1.0, 'DisplayName', 'accept limit');
ylabel('Edge risk (%)');
title('Dynamic points near static-domain edge');
style_axes_local();

nexttile; hold on;
for ic = 1:numel(cases)
    M = rows(rows.CaseTag == cases(ic), :);
    plot(M.GradientRatio, M.Score, '^-', 'Color', colors(ic,:), ...
        'LineWidth', 1.2, 'MarkerSize', 4.5, 'DisplayName', cases(ic));
    scatter(M.GradientRatio(M.IsSelected), M.Score(M.IsSelected), ...
        55, colors(ic,:), 'filled', 'HandleVisibility', 'off');
end
xlabel('Gradient threshold ratio');
ylabel('Score');
title('Fallback score');
style_axes_local();

exportgraphics(fig, output_png, 'Resolution', 300);
end

function style_axes_local()
box on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 8.5, 'TickDir', 'in');
end

function values = parse_numeric_list_env_local(name, default_values)
raw = strtrim(getenv(name));
if isempty(raw)
    values = default_values;
    return;
end
values = sscanf(regexprep(raw, '[,;]', ' '), '%f').';
if isempty(values) || any(~isfinite(values))
    error('%s must be a numeric list.', name);
end
end

function value = parse_nonnegative_numeric_env_local(name, default_value)
raw = strtrim(getenv(name));
if isempty(raw)
    value = default_value;
    return;
end
value = str2double(raw);
if ~isfinite(value) || value < 0
    error('%s must be nonnegative numeric.', name);
end
end

function value = parse_logical_env_local(name, default_value)
raw = lower(strtrim(getenv(name)));
if isempty(raw)
    value = default_value;
elseif ismember(raw, {'1','true','yes','on'})
    value = true;
elseif ismember(raw, {'0','false','no','off'})
    value = false;
else
    error('%s must be logical.', name);
end
end

function old_env = capture_env_local(names)
old_env = struct();
old_env.names = names;
old_env.values = cell(size(names));
for i = 1:numel(names)
    old_env.values{i} = getenv(names{i});
end
end

function restore_env_local(old_env)
for i = 1:numel(old_env.names)
    setenv(old_env.names{i}, old_env.values{i});
end
end

function run_script_in_base_local(script_file)
evalin('base', sprintf('run(''%s'');', strrep(script_file, '''', '''''')));
end
