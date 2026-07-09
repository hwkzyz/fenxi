function Analyze_Step04_Calibration_Methods_20250527(action)
%ANALYZE_STEP04_CALIBRATION_METHODS_20250527 Compare Step04 calibration policies.
%
% This utility keeps the 20250527 data, OPR semantics, Step05 windowing, and
% identification settings fixed.  Only the low-speed calibration template is
% changed, so the Step05 trend comparison isolates the effect of Step04.
%
% Usage:
%   Analyze_Step04_Calibration_Methods_20250527('build')
%   % run Step05 once per generated template with STEP05_TEMPLATE_OVERRIDE_FILE
%   Analyze_Step04_Calibration_Methods_20250527('summarize')

if nargin < 1 || isempty(action)
    action = 'build';
end
action = lower(strtrim(string(action)));

cfg = BTTProjectConfig_20250527();
cfg = apply_compare_defaults_local(cfg);
Analysis = analysis_settings_local(cfg);
cmp_dir = fullfile(cfg.output_root, 'step04_calibration_method_compare');
fig_dir = fullfile(cfg.figure_root, 'step04_calibration_method_compare');
ensure_dir_local(cmp_dir);
ensure_dir_local(fig_dir);

switch action
    case "build"
        build_template_variants_local(cfg, cmp_dir, fig_dir, Analysis);
    case "summarize"
        summarize_step05_comparison_local(cfg, cmp_dir, fig_dir, Analysis);
    case "all"
        build_template_variants_local(cfg, cmp_dir, fig_dir, Analysis);
        summarize_step05_comparison_local(cfg, cmp_dir, fig_dir, Analysis);
    otherwise
        error('Unknown action: %s', action);
end
end


function build_template_variants_local(cfg, cmp_dir, fig_dir, Analysis)
base_template_file = fullfile(cfg.step04_output_dir, ...
    sprintf('Template_OPRCenterStd_LowSpeed_AllBlades_S%s_%s.mat', ...
    sensor_tag_local(cfg.sensor_ids), cfg.dataset));
source_file = fullfile(cfg.step04_output_dir, ...
    sprintf('LowSpeed_Template_SourceData_%s.mat', cfg.dataset));

if ~isfile(base_template_file)
    error('Missing base Step04 template. Run Step04 first:\n  %s', base_template_file);
end
if ~isfile(source_file)
    error('Missing Step04 source point cloud. Run Step04 first:\n  %s', source_file);
end

loaded = load(base_template_file, 'Template');
src = load(source_file, 'point_cloud');
TemplateBase = loaded.Template;
point_cloud = src.point_cloud;

methods = method_defs_local(cfg, cmp_dir);
summary_rows = [];

for im = 1:numel(methods)
    M = methods(im);
    T = TemplateBase;
    T.Calibration_Method_Tag = M.tag;
    T.Calibration_Method_Label = M.label;
    T.LowSpeed_Template_Source_Mode = M.source_mode;
    T.Method_Comparison_Source_Template = base_template_file;
    T.Method_Comparison_Source_PointCloud = source_file;
    T.Method_Comparison_CreatedOn = datestr(now, 31);
    if isfield(T, 'Metadata')
        T.Metadata.calibration_method_tag = M.tag;
        T.Metadata.calibration_method_label = M.label;
    end

    for it = 1:numel(T.SensorBlade)
        sid = T.SensorBlade(it).sensor_id;
        bid = T.SensorBlade(it).blade_id;
        pc = find_point_cloud_local(point_cloud, sid, bid);
        base_entry = T.SensorBlade(it);

        switch M.policy
            case "wide_absolute"
                entry = build_wide_absolute_entry_local(cfg, base_entry, pc);
            case "wide_relative"
                entry = build_wide_relative_entry_local(cfg, base_entry);
            case "trusted_buffer"
                entry = build_trusted_buffer_entry_local(cfg, base_entry, pc, M.buffer_mm);
            case "trusted_strict"
                entry = base_entry;
                entry.comparison_policy = 'trusted_strict_current_20250527';
            otherwise
                error('Unsupported policy: %s', M.policy);
        end

        entry.calibration_method_tag = M.tag;
        entry.calibration_method_label = M.label;
        entry = harmonize_entry_fields_local(T.SensorBlade(it), entry);
        T.SensorBlade = ensure_struct_fields_local(T.SensorBlade, fieldnames(entry));
        entry = ensure_entry_fields_local(entry, fieldnames(T.SensorBlade));
        T.SensorBlade(it) = entry;
    end

    Template = T; %#ok<NASGU>
    save(M.template_file, 'Template', '-v7.3');
    rows = template_summary_rows_local(Template, M);
    summary_rows = [summary_rows; rows]; %#ok<AGROW>
    fprintf('Saved %s template:\n  %s\n', M.tag, M.template_file);
end

Summary = struct2table(summary_rows, 'AsArray', true);
writetable(Summary, fullfile(cmp_dir, ...
    sprintf('Step04_CalibrationMethod_TemplateSummary_%s.csv', cfg.dataset)));
plot_template_method_overview_local(Summary, fig_dir, cfg);
    write_step05_run_commands_local(methods, cfg, cmp_dir, Analysis);
end


function Analysis = analysis_settings_local(cfg)
Analysis = struct();
Analysis.case_name = cfg.dynamic_cases{1};
Analysis.target_blade = 1;
Analysis.sensors = cfg.sensor_ids;
Analysis.sensor_tag = sensor_tag_local(Analysis.sensors);
Analysis.extra_step05_env = {};
end


function cfg = apply_compare_defaults_local(cfg)
if ~isfield(cfg, 'step04_x_collect_abs_limit_mm') || isempty(cfg.step04_x_collect_abs_limit_mm)
    cfg.step04_x_collect_abs_limit_mm = 12.0;
end
if ~isfield(cfg, 'step04_grid_dx_mm') || isempty(cfg.step04_grid_dx_mm)
    cfg.step04_grid_dx_mm = 0.02;
end
if ~isfield(cfg, 'step04_min_bin_count') || isempty(cfg.step04_min_bin_count)
    cfg.step04_min_bin_count = 5;
end
if ~isfield(cfg, 'step04_smooth_span_bins') || isempty(cfg.step04_smooth_span_bins)
    cfg.step04_smooth_span_bins = 9;
end
if mod(cfg.step04_smooth_span_bins, 2) == 0
    cfg.step04_smooth_span_bins = cfg.step04_smooth_span_bins + 1;
end
if ~isfield(cfg, 'step04_gradient_min_ratio') || isempty(cfg.step04_gradient_min_ratio)
    cfg.step04_gradient_min_ratio = 0.30;
end
if ~isfield(cfg, 'step04_amplitude_min_ratio') || isempty(cfg.step04_amplitude_min_ratio)
    cfg.step04_amplitude_min_ratio = 0.02;
end
if ~isfield(cfg, 'step04_min_half_width_mm') || isempty(cfg.step04_min_half_width_mm)
    cfg.step04_min_half_width_mm = 2.5;
end
if ~isfield(cfg, 'step04_max_half_width_mm') || isempty(cfg.step04_max_half_width_mm)
    cfg.step04_max_half_width_mm = 4.2;
end
if ~isfield(cfg, 'step04_query_safe_margin_mm') || isempty(cfg.step04_query_safe_margin_mm)
    cfg.step04_query_safe_margin_mm = 0.20;
end
if ~isfield(cfg, 'step04_min_template_points') || isempty(cfg.step04_min_template_points)
    cfg.step04_min_template_points = 1500;
end
if ~isfield(cfg, 'step04_min_domain_width_mm') || isempty(cfg.step04_min_domain_width_mm)
    cfg.step04_min_domain_width_mm = 2.0 * cfg.step04_min_half_width_mm;
end
end


function methods = method_defs_local(cfg, cmp_dir)
base = fullfile(cmp_dir, 'templates');
ensure_dir_local(base);

methods = struct([]);
methods(1).tag = 'm20250527_wide_abs';
methods(1).label = '20250527 style: wide cloud direct template, absolute x';
methods(1).policy = "wide_absolute";
methods(1).source_mode = 'method_compare_20250527_wide_absolute';
methods(1).buffer_mm = 0;

methods(2).tag = 'm20251222_trusted_buffer';
methods(2).label = '20251222 style: xc + trusted-domain recalibration with small buffer';
methods(2).policy = "trusted_buffer";
methods(2).source_mode = 'method_compare_20251222_trusted_buffer';
methods(2).buffer_mm = 0.20;

methods(3).tag = 'm20250527_trusted_strict';
methods(3).label = '20250527 current foundation template';
methods(3).policy = "trusted_strict";
methods(3).source_mode = 'method_compare_20250527_current_foundation';
methods(3).buffer_mm = 0;

for i = 1:numel(methods)
    methods(i).template_file = fullfile(base, sprintf( ...
        'Template_MethodCompare_%s_S%s_%s.mat', ...
        methods(i).tag, sensor_tag_local(cfg.sensor_ids), cfg.dataset));
end
methods(1).step05_output_tag = 'c25wide';
methods(2).step05_output_tag = 'c22buf';
methods(3).step05_output_tag = 'c25curr';
end


function entry = build_wide_absolute_entry_local(cfg, base_entry, pc)
entry = base_entry;
if isempty(pc)
    entry.quality_status = 'missing_point_cloud';
    return;
end

x = pc.x_mm(:);
v = pc.v(:);
entry = aggregate_points_for_entry_local(cfg, entry, x, v, [], ...
    'wide_absolute_all_points', false);
entry.xc = NaN;
entry.x_coordinate_mode = 'x_abs in OPRCenterStd frame';
entry.source_point_policy = '20250527-style wide point cloud directly builds final template';
entry.final_template_point_policy = 'wide_point_cloud_direct';
end


function entry = build_wide_relative_entry_local(cfg, base_entry)
entry = base_entry;
if ~isfield(base_entry, 'wide_x_grid') || isempty(base_entry.wide_x_grid)
    entry.comparison_policy = 'wide_relative_unavailable';
    return;
end
entry.x_grid = base_entry.wide_x_grid(:);
entry.v_grid = base_entry.wide_v_grid(:);
entry.v_grid_raw = base_entry.wide_v_grid_raw(:);
entry.dv_dx = base_entry.wide_dv_dx(:);
entry.weight_grid = base_entry.wide_weight_grid(:);
entry.count_grid = base_entry.wide_count_grid(:);
entry.valid_grid_mask = base_entry.wide_valid_grid_mask(:);
entry.domain_mask = base_entry.wide_domain_mask(:);
entry.domain_effective_mask = base_entry.wide_domain_effective_mask(:);
entry.x_query_safe_domain = shrink_domain_local(entry.x_domain, cfg.step04_query_safe_margin_mm);
entry.comparison_policy = 'wide_relative_all_points';
entry.source_point_policy = 'wide relative point cloud directly builds final template';
entry.final_template_point_policy = 'wide_point_cloud_direct';
end


function entry = build_trusted_buffer_entry_local(cfg, base_entry, pc, buffer_mm)
entry = base_entry;
if isempty(pc)
    entry.quality_status = 'missing_point_cloud';
    return;
end
if ~isfield(base_entry, 'xc') || ~isfinite(base_entry.xc)
    entry.quality_status = 'missing_xc';
    return;
end
if ~isfield(base_entry, 'x_domain') || any(~isfinite(base_entry.x_domain))
    entry.quality_status = 'missing_domain';
    return;
end

x = pc.x_mm(:) - base_entry.xc;
v = pc.v(:);
fit_domain = [base_entry.x_domain(1) - buffer_mm, base_entry.x_domain(2) + buffer_mm];
fit_mask = isfinite(x) & isfinite(v) & x >= fit_domain(1) & x <= fit_domain(2);
entry = aggregate_points_for_entry_local(cfg, entry, x(fit_mask), v(fit_mask), ...
    base_entry.x_domain, 'trusted_domain_buffer_recalibration', true, fit_domain);
entry.xc = base_entry.xc;
entry.x_coordinate_mode = 'x_rel = x_abs - xc';
entry.x_fit_domain = fit_domain;
entry.x_query_safe_domain = shrink_domain_local(base_entry.x_domain, cfg.step04_query_safe_margin_mm);
entry.source_point_policy = sprintf('trusted-domain points plus %.3g mm buffer', buffer_mm);
entry.final_template_point_policy = 'trusted_domain_buffer_recalibration';
end


function entry = aggregate_points_for_entry_local(cfg, entry, x, v, forced_domain, policy, use_forced_domain, grid_range)
x = x(:);
v = v(:);
valid = isfinite(x) & isfinite(v);
x = x(valid);
v = v(valid);

entry.point_count = numel(x);
entry.comparison_policy = policy;
entry.used_for_final_template_count = numel(x);

if numel(x) < cfg.step04_min_template_points
    entry.quality_status = 'insufficient_points';
    return;
end

if nargin >= 8 && ~isempty(grid_range) && all(isfinite(grid_range))
    x_lo = grid_range(1);
    x_hi = grid_range(2);
elseif nargin < 6 || ~use_forced_domain || isempty(forced_domain) || any(~isfinite(forced_domain))
    x_lo = max(prctile(x, 0.5), -cfg.step04_x_collect_abs_limit_mm);
    x_hi = min(prctile(x, 99.5), cfg.step04_x_collect_abs_limit_mm);
else
    x_lo = forced_domain(1);
    x_hi = forced_domain(2);
end
if x_hi <= x_lo
    entry.quality_status = 'degenerate_x_range';
    return;
end

[x_grid, v_raw, v_smooth, count_grid, valid_mask] = ...
    grid_median_smooth_local(cfg, x, v, [x_lo, x_hi]);
if nnz(valid_mask) < 5
    entry.quality_status = 'insufficient_coverage';
    return;
end

baseline = estimate_baseline_local(x, v, v_smooth);
amplitude = max(v_smooth, [], 'omitnan') - baseline;
if ~isfinite(amplitude) || amplitude <= 0
    amplitude = max(v_smooth, [], 'omitnan') - min(v_smooth, [], 'omitnan');
end
if ~isfinite(amplitude)
    amplitude = 0;
end

dv_dx = gradient(v_smooth, x_grid);
max_abs_grad = max(abs(dv_dx(valid_mask)), [], 'omitnan');
if ~isfinite(max_abs_grad)
    max_abs_grad = 0;
end

if nargin >= 5 && ~isempty(forced_domain) && all(isfinite(forced_domain))
    x_domain = forced_domain;
    domain_mask = x_grid >= x_domain(1) & x_grid <= x_domain(2) & valid_mask;
else
    [x_domain, domain_mask] = determine_domain_simple_local(cfg, x_grid, ...
        v_smooth, dv_dx, valid_mask, baseline, amplitude, max_abs_grad);
end

weight_grid = build_weight_grid_local(count_grid, valid_mask, domain_mask, dv_dx, max_abs_grad);

entry.x_grid = x_grid(:);
entry.v_grid = v_smooth(:);
entry.v_grid_raw = v_raw(:);
entry.v_grid_baseline_removed = v_smooth(:) - baseline;
entry.dv_dx = dv_dx(:);
entry.weight_grid = weight_grid(:);
entry.count_grid = count_grid(:);
entry.valid_grid_mask = valid_mask(:);
entry.domain_mask = domain_mask(:);
entry.domain_effective_mask = domain_mask(:);
entry.x_domain = x_domain;
entry.x_fit_domain = [x_lo, x_hi];
entry.x_query_safe_domain = shrink_domain_local(x_domain, cfg.step04_query_safe_margin_mm);
entry.baseline = baseline;
entry.amplitude = amplitude;
entry.max_abs_gradient_v_per_mm = max_abs_grad;
entry.quality_status = quality_status_local(cfg, numel(x), x_domain, valid_mask);
end


function [x_grid, v_raw, v_smooth, count_grid, valid_mask] = ...
    grid_median_smooth_local(cfg, x, v, x_range)
edges = x_range(1):cfg.step04_grid_dx_mm:x_range(2);
if numel(edges) < 5
    edges = linspace(x_range(1), x_range(2), 20);
end
if edges(end) < x_range(2)
    edges(end + 1) = x_range(2); %#ok<AGROW>
end
x_grid = 0.5 * (edges(1:end-1) + edges(2:end));
x_grid = x_grid(:);
bin = discretize(x, edges);
n = numel(x_grid);
v_raw = nan(n, 1);
count_grid = zeros(n, 1);
for i = 1:n
    mask = bin == i;
    count_grid(i) = nnz(mask);
    if count_grid(i) > 0
        v_raw(i) = median(v(mask), 'omitnan');
    end
end
valid_mask = count_grid >= cfg.step04_min_bin_count & isfinite(v_raw);
v_interp = v_raw;
missing = ~isfinite(v_interp);
if any(missing) && nnz(valid_mask) >= 2
    v_interp(missing) = interp1(x_grid(valid_mask), v_raw(valid_mask), ...
        x_grid(missing), 'linear', 'extrap');
end
v_smooth = smoothdata(v_interp, 'movmedian', cfg.step04_smooth_span_bins, 'omitnan');
v_smooth = smoothdata(v_smooth, 'movmean', cfg.step04_smooth_span_bins, 'omitnan');
end


function [x_domain, domain_mask] = determine_domain_simple_local(cfg, x_grid, v_grid, dv_dx, valid_mask, baseline, amplitude, max_abs_grad)
if ~isfinite(max_abs_grad) || max_abs_grad <= 0
    x_domain = [NaN NaN];
    domain_mask = false(size(x_grid));
    return;
end
amp_norm = abs(v_grid - baseline) ./ max(abs(amplitude), eps);
mask = valid_mask & abs(dv_dx) >= cfg.step04_gradient_min_ratio * max_abs_grad;
mask = mask & amp_norm >= cfg.step04_amplitude_min_ratio;
center_x = 0;
if nnz(mask) < 2
    mask = valid_mask;
end
left_candidates = find(x_grid < center_x & mask);
right_candidates = find(x_grid > center_x & mask);
if isempty(left_candidates) || isempty(right_candidates)
    valid_x = x_grid(valid_mask);
    if isempty(valid_x)
        x_domain = [NaN NaN];
    else
        x_domain = [min(valid_x), max(valid_x)];
    end
else
    left_L = max(center_x - x_grid(left_candidates(1)), cfg.step04_min_half_width_mm);
    right_L = max(x_grid(right_candidates(end)) - center_x, cfg.step04_min_half_width_mm);
    left_L = min(left_L, cfg.step04_max_half_width_mm);
    right_L = min(right_L, cfg.step04_max_half_width_mm);
    x_domain = [center_x - left_L, center_x + right_L];
end
domain_mask = x_grid >= x_domain(1) & x_grid <= x_domain(2) & valid_mask;
end


function weight = build_weight_grid_local(count_grid, valid_mask, domain_mask, dv_dx, max_abs_grad)
count_grid = count_grid(:);
valid_mask = valid_mask(:);
domain_mask = domain_mask(:);
weight = zeros(size(count_grid));
if any(valid_mask)
    c = sqrt(max(count_grid, 0));
    c = c ./ max(c(valid_mask), [], 'omitnan');
    weight(valid_mask) = c(valid_mask);
end
if isfinite(max_abs_grad) && max_abs_grad > 0
    g = abs(dv_dx(:)) ./ max_abs_grad;
    weight = weight .* max(0.25, min(1, g));
end
weight(~domain_mask) = 0;
weight(~isfinite(weight)) = 0;
end


function baseline = estimate_baseline_local(x, v, v_grid)
edge = abs(x) >= prctile(abs(x), 70);
if nnz(edge) >= 20
    baseline = median(v(edge), 'omitnan');
else
    baseline = prctile(v, 10);
end
if ~isfinite(baseline)
    baseline = min(v_grid, [], 'omitnan');
end
end


function status = quality_status_local(cfg, point_count, x_domain, valid_mask)
status = 'good';
if point_count < cfg.step04_min_template_points
    status = 'insufficient_points';
elseif any(~isfinite(x_domain)) || x_domain(2) <= x_domain(1)
    status = 'no_trusted_domain';
elseif diff(x_domain) < cfg.step04_min_domain_width_mm
    status = 'narrow_domain';
elseif nnz(valid_mask) < 5
    status = 'insufficient_coverage';
end
end


function pc = find_point_cloud_local(point_cloud, sid, bid)
pc = [];
for i = 1:numel(point_cloud)
    if point_cloud(i).sensor_id == sid && point_cloud(i).blade_id == bid
        pc = point_cloud(i);
        return;
    end
end
end


function out = harmonize_entry_fields_local(template_entry, new_entry)
out = template_entry;
names = fieldnames(new_entry);
for i = 1:numel(names)
    out.(names{i}) = new_entry.(names{i});
end
end


function S = ensure_struct_fields_local(S, names)
current = fieldnames(S);
for i = 1:numel(names)
    if ~ismember(names{i}, current)
        [S.(names{i})] = deal([]);
    end
end
end


function entry = ensure_entry_fields_local(entry, names)
current = fieldnames(entry);
for i = 1:numel(names)
    if ~ismember(names{i}, current)
        entry.(names{i}) = [];
    end
end
end


function rows = template_summary_rows_local(Template, M)
rows = struct([]);
for i = 1:numel(Template.SensorBlade)
    e = Template.SensorBlade(i);
    row = struct();
    row.method_tag = M.tag;
    row.method_label = M.label;
    row.sensor_id = e.sensor_id;
    row.blade_id = e.blade_id;
    row.quality_status = char(string(get_field_default_local(e, 'quality_status', 'unknown')));
    row.xc_mm = get_field_default_local(e, 'xc', NaN);
    row.domain_left_mm = e.x_domain(1);
    row.domain_right_mm = e.x_domain(2);
    row.domain_width_mm = diff(e.x_domain);
    row.grid_count = numel(e.x_grid);
    row.valid_bin_count = nnz(e.valid_grid_mask);
    row.final_point_count = get_field_default_local(e, 'used_for_final_template_count', ...
        get_field_default_local(e, 'trusted_point_count', NaN));
    row.amplitude_v = get_field_default_local(e, 'amplitude', NaN);
    row.max_abs_gradient_v_per_mm = max(abs(e.dv_dx), [], 'omitnan');
    rows = [rows; row]; %#ok<AGROW>
end
end


function summarize_step05_comparison_local(cfg, cmp_dir, fig_dir, Analysis)
methods = method_defs_local(cfg, cmp_dir);
rows = struct([]);
for i = 1:numel(methods)
    M = methods(i);
    trend_file = fullfile(cfg.output_root, ...
        ['step05_single_sync_direct_template_', M.step05_output_tag], ...
        Analysis.case_name, ...
        sprintf('Trend_Step05_SingleSyncDirectTemplate_B%d_S%s_%s.csv', ...
        Analysis.target_blade, Analysis.sensor_tag, cfg.dataset));
    row = summarize_one_trend_local(M, trend_file);
    rows = [rows; row]; %#ok<AGROW>
end

T = struct2table(rows, 'AsArray', true);
out_csv = fullfile(cmp_dir, sprintf('Step05_CalibrationMethod_Comparison_%s.csv', cfg.dataset));
writetable(T, out_csv);
plot_step05_method_comparison_local(T, fig_dir, cfg);
fprintf('Saved Step05 comparison summary:\n  %s\n', out_csv);
disp(T);
end


function row = summarize_one_trend_local(M, trend_file)
row = struct();
row.method_tag = M.tag;
row.method_label = M.label;
row.trend_file = trend_file;
row.file_exists = isfile(trend_file);
row.window_count = 0;
row.ok_count = 0;
row.eo12_count = 0;
row.eo_mode = NaN;
row.eo_unique = "";
row.mean_freq_hz = NaN;
row.std_freq_hz = NaN;
row.mean_rmse_v = NaN;
row.median_rmse_v = NaN;
row.mean_amplitude_mm = NaN;
row.std_amplitude_mm = NaN;
row.mean_dx_mm = NaN;
row.std_dx_mm = NaN;
row.mean_point_count = NaN;
row.min_point_count = NaN;
row.mean_core_point_count = NaN;
row.max_eta_mm = NaN;
row.correctness_score = NaN;
if ~isfile(trend_file)
    return;
end
T = readtable(trend_file);
if isempty(T)
    return;
end
ok = strcmpi(string(T.status), "ok") & isfinite(T.EO_id);
row.window_count = height(T);
row.ok_count = nnz(ok);
row.eo12_count = nnz(ok & T.EO_id == 12);
if any(ok)
    eos = T.EO_id(ok);
    row.eo_mode = mode(eos);
    row.eo_unique = strjoin(string(unique(eos(:).')), ',');
    row.mean_freq_hz = mean(T.fn_id(ok), 'omitnan');
    row.std_freq_hz = std(T.fn_id(ok), 'omitnan');
    row.mean_rmse_v = mean(T.weighted_voltage_rmse(ok), 'omitnan');
    row.median_rmse_v = median(T.weighted_voltage_rmse(ok), 'omitnan');
    row.mean_amplitude_mm = mean(T.A_id(ok), 'omitnan');
    row.std_amplitude_mm = std(T.A_id(ok), 'omitnan');
    row.mean_dx_mm = mean(T.dx_c_id(ok), 'omitnan');
    row.std_dx_mm = std(T.dx_c_id(ok), 'omitnan');
    row.mean_point_count = mean(T.point_count(ok), 'omitnan');
    row.min_point_count = min(T.point_count(ok), [], 'omitnan');
    row.mean_core_point_count = mean(T.core_point_count(ok), 'omitnan');
    row.max_eta_mm = max(T.sensor_eta_max_abs_mm(ok), [], 'omitnan');
    eo_fraction = row.eo12_count / max(row.window_count, 1);
    stability = 1 / (1 + row.std_freq_hz);
    residual = 1 / (1 + 100 * row.mean_rmse_v);
    coverage = min(row.ok_count / max(row.window_count, 1), 1);
    row.correctness_score = 0.45 * eo_fraction + 0.25 * coverage + ...
        0.15 * stability + 0.15 * residual;
end
end


function plot_template_method_overview_local(Summary, fig_dir, cfg)
fig = figure('Color', 'w', 'Position', [80, 80, 1200, 760], 'Visible', 'off');
tiledlayout(fig, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

methods = unique(string(Summary.method_tag), 'stable');
colors = lines(numel(methods));

nexttile;
hold on; grid on;
for i = 1:numel(methods)
    mask = string(Summary.method_tag) == methods(i);
    scatter(Summary.domain_width_mm(mask), Summary.valid_bin_count(mask), ...
        42, colors(i, :), 'filled', 'DisplayName', methods(i));
end
xlabel('domain width (mm)');
ylabel('valid bin count');
title('Template coverage');
legend('Location', 'best', 'Interpreter', 'none');

nexttile;
boxchart(categorical(Summary.method_tag), Summary.max_abs_gradient_v_per_mm);
ylabel('max |dV/dx| (V/mm)');
title('Template slope');
grid on;

nexttile;
boxchart(categorical(Summary.method_tag), Summary.amplitude_v);
ylabel('amplitude (V)');
title('Template amplitude');
grid on;

nexttile;
status = categorical(Summary.quality_status);
histogram(status);
title('Quality status count');
grid on;

save_figure_local(fig, fig_dir, sprintf('Step04_CalibrationMethod_TemplateOverview_%s', cfg.dataset));
close(fig);
end


function plot_step05_method_comparison_local(T, fig_dir, cfg)
fig = figure('Color', 'w', 'Position', [80, 80, 1200, 760], 'Visible', 'off');
tiledlayout(fig, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
labels = method_short_labels_local(string(T.method_tag));
x = categorical(labels, labels, 'Ordinal', true);

nexttile;
bar(x, T.eo12_count ./ max(T.window_count, 1));
ylabel('EO12 fraction');
ylim([0, 1.05]);
title('EO consistency');
grid on;

nexttile;
bar(x, T.mean_rmse_v * 1000);
ylabel('mean weighted RMSE (mV)');
title('Residual');
grid on;

nexttile;
bar(x, T.std_freq_hz * 1000);
ylabel('frequency std (mHz)');
title('Frequency stability');
grid on;

nexttile;
bar(x, T.correctness_score);
ylabel('score');
ylim([0, 1.05]);
title('Overall diagnostic score');
grid on;

save_figure_local(fig, fig_dir, sprintf('Step05_CalibrationMethod_Comparison_%s', cfg.dataset));
close(fig);

fig2 = figure('Color', 'w', 'Position', [90, 90, 1200, 420], 'Visible', 'off');
tiledlayout(fig2, 1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
nexttile;
bar(x, T.mean_amplitude_mm);
ylabel('mean A (mm)');
title('Identified amplitude');
grid on;
nexttile;
bar(x, T.mean_dx_mm);
ylabel('mean dx_c (mm)');
title('Static offset absorbed by fit');
grid on;
nexttile;
bar(x, T.max_eta_mm);
ylabel('max |eta| (mm)');
title('Sensor offset usage');
grid on;
save_figure_local(fig2, fig_dir, sprintf('Step05_CalibrationMethod_AmplitudeDxEta_%s', cfg.dataset));
close(fig2);
end


function write_step05_run_commands_local(methods, cfg, cmp_dir, Analysis)
cmd_file = fullfile(cmp_dir, sprintf('Run_Step05_CalibrationMethod_Commands_%s.txt', cfg.dataset));
fid = fopen(cmd_file, 'w');
cleanup = onCleanup(@() fclose(fid));
project_dir = cfg.route_dir;
fprintf(fid, 'Run these commands from PowerShell, one at a time:\n\n');
for i = 1:numel(methods)
    M = methods(i);
    fprintf(fid, ['matlab -batch "cd(''%s''); setenv(''STEP05_OUTPUT_TAG'',''%s''); ', ...
        'setenv(''STEP05_TEMPLATE_OVERRIDE_FILE'',''%s''); ', ...
        'setenv(''STEP05_FORCE_REBUILD'',''true''); setenv(''STEP05_SAVE_FIGURES'',''true''); ', ...
        'setenv(''STEP05_SHOW_PLOTS'',''false''); Step05_SingleSync_DirectTemplate_Identification_20250527"\n'], ...
        project_dir, M.step05_output_tag, M.template_file);
end
fprintf(fid, '\nThen run:\n');
fprintf(fid, 'matlab -batch "cd(''%s''); Analyze_Step04_Calibration_Methods_20250527(''summarize'')"\n', project_dir);
fprintf('Saved Step05 run commands:\n  %s\n', cmd_file);
end


function safe_domain = shrink_domain_local(x_domain, margin_mm)
safe_domain = [NaN NaN];
if numel(x_domain) ~= 2 || any(~isfinite(x_domain)) || x_domain(2) <= x_domain(1)
    return;
end
safe_domain = [x_domain(1) + margin_mm, x_domain(2) - margin_mm];
if safe_domain(2) <= safe_domain(1)
    c = mean(x_domain);
    safe_domain = [c, c];
end
end


function val = get_field_default_local(s, name, default_val)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    val = s.(name);
else
    val = default_val;
end
end


function labels = method_short_labels_local(tags)
labels = strings(size(tags));
for i = 1:numel(tags)
    switch char(tags(i))
        case 'm20250527_wide_abs'
            labels(i) = "20250527-wide";
        case 'm20251222_trusted_buffer'
            labels(i) = "20251222-buffer";
        case 'm20250527_trusted_strict'
            labels(i) = "20250527-current";
        otherwise
            labels(i) = tags(i);
    end
end
end


function tag = sensor_tag_local(sensor_ids)
tag = sprintf('%d', sensor_ids(:).');
end


function ensure_dir_local(p)
if ~exist(p, 'dir')
    mkdir(p);
end
end


function save_figure_local(fig, fig_dir, name)
ensure_dir_local(fig_dir);
exportgraphics(fig, fullfile(fig_dir, [name, '.png']), 'Resolution', 300);
exportgraphics(fig, fullfile(fig_dir, [name, '.pdf']), 'ContentType', 'vector');
end
