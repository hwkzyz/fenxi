%% Step01_01_Validate_StaticXRange_20250527
% Validate the static template x-range with fixed 30-lap calibration windows.
%
% This script does not scan LapCount. It uses the current 30-lap window plan
% and compares candidate x-range definitions by static fit quality, gradient
% sensitivity, edge/background contamination, and dynamic effective coverage.

clear; clc; close all;

%% Settings
route_dir = fileparts(mfilename('fullpath'));
validation_root = fileparts(route_dir);
legacy_dir = fullfile(validation_root, '20250527', 'legacy');
if exist(legacy_dir, 'dir') ~= 7
    error('Legacy 20250527 helper folder not found: %s', legacy_dir);
end
addpath(legacy_dir);

legacy_cfg = Get_20250527_BTT_Config();
target_blade = 1;
analysis_sensors = [1, 3, 6];
sensor_tag = ['S', sprintf('%d', analysis_sensors)];
low_speed_data_dir = fullfile(legacy_cfg.dataset_root, legacy_cfg.low_speed_case);
sensor_config_file = fullfile(legacy_cfg.reference_output_dir, 'Sensor_Config_20250527.mat');

output_dir = fullfile(route_dir, 'output');
plan_dir = fullfile(output_dir, 'stable_window_plan');
dynamic_dir = fullfile(output_dir, 'dynamic_maps');
figure_dir = fullfile(output_dir, 'figures');
if exist(plan_dir, 'dir') ~= 7; mkdir(plan_dir); end
if exist(figure_dir, 'dir') ~= 7; mkdir(figure_dir); end

plan_file = strtrim(getenv('STEP01_XRANGE_PLAN_FILE'));
if isempty(plan_file)
    plan_file = fullfile(plan_dir, sprintf('Step01_CoverageFirstStableWindowPlan_B%d_%s_20250527.csv', ...
        target_blade, sensor_tag));
end
if exist(plan_file, 'file') ~= 2
    error('Window plan not found: %s', plan_file);
end

gradient_ratios = parse_numeric_list_env_local('STEP01_XRANGE_GRADIENT_RATIOS', [0.04, 0.06, 0.08, 0.10, 0.12, 0.15, 0.20]);
amplitude_ratios = parse_numeric_list_env_local('STEP01_XRANGE_AMPLITUDE_RATIOS', [0.06, 0.08, 0.10, 0.12, 0.15]);
fixed_half_widths = parse_numeric_list_env_local('STEP01_XRANGE_FIXED_HALF_WIDTHS_MM', [2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0]);
dynamic_query_guard_mm = parse_nonnegative_numeric_env_local('STEP01_XRANGE_QUERY_GUARD_MM', 0.75);
dynamic_quantile = [1, 99];
weight_floor = 0.05;
segment_expand_factor = 0.30;

result_csv = fullfile(plan_dir, sprintf('Step01_StaticXRangeValidation_B%d_%s_20250527.csv', ...
    target_blade, sensor_tag));
summary_csv = fullfile(plan_dir, sprintf('Step01_StaticXRangeRecommendation_B%d_%s_20250527.csv', ...
    target_blade, sensor_tag));

fprintf('\n=== Step01_01: static x-range validation ===\n');
fprintf('Plan source: %s\n', plan_file);
fprintf('Sensors: %s\n', mat2str(analysis_sensors));

loaded_cfg = load(sensor_config_file, 'Sensor_Config');
Sensor_Config = loaded_cfg.Sensor_Config;
window_plan = readtable(plan_file);
[raw_stream, opr_times] = load_low_speed_streams_local(low_speed_data_dir, legacy_cfg, analysis_sensors);
F_omega_deg = build_phase_speed_local(opr_times, legacy_cfg.blades_num);
dynamic_file = find_dynamic_file_for_sensors_local(dynamic_dir, target_blade, analysis_sensors);
loaded_dynamic = load(dynamic_file, 'DynamicMap');
DynamicMap = loaded_dynamic.DynamicMap;

all_results = table();
plot_cache = struct('sensor_id', {}, 'x_rel', {}, 'v', {}, 'g_norm', {}, ...
    'amp_norm', {}, 'dynamic_core', {}, 'dynamic_query', {}, 'chosen', {});

for is = 1:numel(analysis_sensors)
    sid = analysis_sensors(is);
    plan_row = window_plan(window_plan.SensorID == sid, :);
    if height(plan_row) ~= 1
        error('Plan must contain exactly one row for CH%d.', sid);
    end

    threshold = legacy_cfg.sensor_thresholds(sid);
    theta_std = Sensor_Config.Standard_Relative_Angles(sid, target_blade);
    R = raw_stream(sid);
    segments = extract_sensor_segments_local(R.T, R.V, threshold, legacy_cfg.gap_points);
    baseline = estimate_background_baseline_local(R.V, segments);
    if isKey(Sensor_Config.Target_Indices, sid)
        start_idx = Sensor_Config.Target_Indices(sid);
    else
        error('Sensor_Config.Target_Indices does not contain CH%d.', sid);
    end

    lap_count = round(plan_row.LapCount(1));
    start_offset_lap = round(plan_row.StartOffsetLap(1));
    selected_idx = start_idx + start_offset_lap * legacy_cfg.blades_num + ...
        (0:(lap_count - 1)) * legacy_cfg.blades_num + (target_blade - 1);
    [x_wide, v_wide, pulse_id] = restore_low_speed_point_cloud_with_pulse_id_local( ...
        R, segments, selected_idx, opr_times, F_omega_deg, theta_std, ...
        legacy_cfg.r_tip_mm, segment_expand_factor);

    P = build_static_cloud_metrics_local(x_wide, v_wide, pulse_id, baseline, threshold, weight_floor);
    [dyn_x_abs, dyn_v, dyn_t] = collect_dynamic_points_local(DynamicMap, sid);
    dyn_mask = build_dynamic_effective_mask_local(dyn_x_abs - P.xc, dyn_v, dyn_t, P, threshold);
    dyn_rel = dyn_x_abs(dyn_mask) - P.xc;
    dynamic_core = prctile(dyn_rel, dynamic_quantile);
    dynamic_query = [dynamic_core(1) - dynamic_query_guard_mm, dynamic_core(2) + dynamic_query_guard_mm];

    sensor_results = table();
    for gr = gradient_ratios(:).'
        for ar = amplitude_ratios(:).'
            mask = P.valid & P.g_norm >= gr & P.amp_norm >= ar;
            row = evaluate_xrange_candidate_local(P, dyn_rel, dynamic_query, mask, ...
                sprintf('thr_g%.3f_a%.3f', gr, ar), 'threshold', gr, ar, NaN, dynamic_query_guard_mm);
            sensor_results = [sensor_results; row]; %#ok<AGROW>
        end
    end
    for hw = fixed_half_widths(:).'
        mask = P.valid & abs(P.x_rel) <= hw & P.amp_norm >= 0.06;
        row = evaluate_xrange_candidate_local(P, dyn_rel, dynamic_query, mask, ...
            sprintf('fixed_hw%.2f', hw), 'fixed_half_width', NaN, 0.06, hw, dynamic_query_guard_mm);
        sensor_results = [sensor_results; row]; %#ok<AGROW>
    end
    sensor_results.BladeID(:) = target_blade;
    sensor_results.SensorID(:) = sid;
    sensor_results.LapCount(:) = lap_count;
    sensor_results.StartOffsetLap(:) = start_offset_lap;
    sensor_results = score_xrange_candidates_local(sensor_results);
    chosen = choose_xrange_candidate_local(sensor_results);
    all_results = [all_results; sensor_results]; %#ok<AGROW>

    pc = struct();
    pc.sensor_id = sid;
    pc.x_rel = P.x_rel;
    pc.v = P.v;
    pc.g_norm = P.g_norm;
    pc.amp_norm = P.amp_norm;
    pc.dynamic_core = dynamic_core;
    pc.dynamic_query = dynamic_query;
    pc.chosen = chosen;
    plot_cache(end + 1) = pc; %#ok<SAGROW>

    fprintf('CH%d recommended %s: x=[%.3f %.3f] mm, safe=%.3f, rmse=%.4f, edge=%.3f, score=%.3f\n', ...
        sid, char(chosen.CandidateName), chosen.x_left, chosen.x_right, ...
        chosen.safe_frac, chosen.static_rmse, chosen.edge_contamination, chosen.Score);
end

recommendation = table();
for sid = analysis_sensors
    rows = all_results(all_results.SensorID == sid, :);
    recommendation = [recommendation; choose_xrange_candidate_local(rows)]; %#ok<AGROW>
end
writetable(all_results, result_csv);
writetable(recommendation, summary_csv);
fprintf('Saved validation table: %s\n', result_csv);
fprintf('Saved recommendation: %s\n', summary_csv);

plot_xrange_tradeoff_local(all_results, recommendation, figure_dir, target_blade, sensor_tag);
plot_xrange_selected_clouds_local(plot_cache, figure_dir, target_blade, sensor_tag);

%% Local functions
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

function dynamic_file = find_dynamic_file_for_sensors_local(dynamic_dir, target_blade, analysis_sensors)
sensor_tag = ['S', sprintf('%d', analysis_sensors)];
candidate = fullfile(dynamic_dir, sprintf('DynamicMap_B%d_%s_SlidingWindows_20250527.mat', target_blade, sensor_tag));
if isfile(candidate)
    dynamic_file = candidate;
    return;
end
files = dir(fullfile(dynamic_dir, sprintf('DynamicMap_B%d_*.mat', target_blade)));
for i = 1:numel(files)
    candidate = fullfile(files(i).folder, files(i).name);
    loaded = load(candidate, 'DynamicMap');
    if isfield(loaded.DynamicMap, 'SensorIDs') && all(ismember(analysis_sensors, loaded.DynamicMap.SensorIDs))
        dynamic_file = candidate;
        return;
    end
end
error('DynamicMap for sensors %s was not found. Run Step02 first.', sensor_tag);
end

function [raw_stream, opr_times] = load_low_speed_streams_local(case_dir, cfg, sensor_ids)
file_ids = list_case_file_ids_local(case_dir, cfg.opr_id);
max_sid = max([sensor_ids(:); cfg.opr_id]);
raw_stream(max_sid) = struct('T', [], 'V', []);
opr_T = [];
opr_V = [];
last_end = [];
for i = 1:numel(file_ids)
    file_id = file_ids(i);
    [t_opr, v_opr] = load_raw_channel_local(case_dir, cfg.opr_id, file_id, cfg.pinlv);
    if isempty(t_opr)
        continue;
    end
    if isempty(last_end)
        offset = 0;
    else
        offset = last_end + 1 / cfg.pinlv - t_opr(1);
    end
    t_opr = t_opr(:) + offset;
    last_end = t_opr(end);
    if i == 1 && cfg.initial_trim_points > 0 && numel(t_opr) > cfg.initial_trim_points
        keep = (cfg.initial_trim_points + 1):numel(t_opr);
        t_opr = t_opr(keep);
        v_opr = v_opr(keep);
    end
    opr_T = [opr_T; t_opr(:)]; %#ok<AGROW>
    opr_V = [opr_V; v_opr(:)]; %#ok<AGROW>
    for sid = sensor_ids
        if ~isfile(fullfile(case_dir, sprintf('4-%d-%d.mat', sid, file_id)))
            continue;
        end
        [t_local, v_local] = load_raw_channel_local(case_dir, sid, file_id, cfg.pinlv);
        t_local = t_local(:) + offset;
        if i == 1 && cfg.initial_trim_points > 0 && numel(t_local) > cfg.initial_trim_points
            keep = (cfg.initial_trim_points + 1):numel(t_local);
            t_local = t_local(keep);
            v_local = v_local(keep);
        end
        raw_stream(sid).T = [raw_stream(sid).T; t_local(:)]; %#ok<AGROW>
        raw_stream(sid).V = [raw_stream(sid).V; v_local(:)]; %#ok<AGROW>
    end
end
opr_segments = extract_opr_segments_local(opr_T, opr_V, cfg.opr_threshold, cfg.gap_points);
opr_times = opr_segments.arrival_time(:);
end

function file_ids = list_case_file_ids_local(case_dir, channel_id)
d = dir(fullfile(case_dir, sprintf('4-%d-*.mat', channel_id)));
file_ids = nan(numel(d), 1);
for i = 1:numel(d)
    tok = regexp(d(i).name, ['4-' num2str(channel_id) '-(\d+)\.mat'], 'tokens', 'once');
    if ~isempty(tok)
        file_ids(i) = str2double(tok{1});
    end
end
file_ids = sort(unique(file_ids(~isnan(file_ids))));
end

function [t_sec, v] = load_raw_channel_local(case_dir, channel_id, file_id, fs)
filepath = fullfile(case_dir, sprintf('4-%d-%d.mat', channel_id, file_id));
if ~isfile(filepath)
    t_sec = [];
    v = [];
    return;
end
loaded = load(filepath);
fn = fieldnames(loaded);
raw = loaded.(fn{1});
raw(raw(:, 1) == 0, :) = [];
t_sec = raw(:, 1) / fs;
v = raw(:, 2);
end

function segments = extract_opr_segments_local(t, v, threshold, gap_points)
idx = find(v(:) > threshold);
segments = split_segments_to_struct_local(t, v, idx, gap_points, threshold, false);
end

function segments = extract_sensor_segments_local(t, v, threshold, gap_points)
if numel(v) >= 21
    v_smooth = sgolayfilt(v(:), 3, 21);
elseif numel(v) >= 5
    v_smooth = smoothdata(v(:), 'movmean', max(3, 2 * floor(numel(v) / 4) + 1));
else
    v_smooth = v(:);
end
idx = find(v_smooth > threshold);
segments = split_segments_to_struct_local(t, v, idx, gap_points, threshold, true);
end

function segments = split_segments_to_struct_local(t, v, idx, gap_points, threshold, use_centroid)
segments = struct('start_idx', [], 'end_idx', [], 'start_time', [], ...
    'end_time', [], 'arrival_time', [], 'peak_value', []);
if isempty(idx)
    return;
end
jumps = find(diff(idx) > gap_points);
seg_start_pos = [1; jumps(:) + 1];
seg_end_pos = [jumps(:); numel(idx)];
n = numel(seg_start_pos);
segments.start_idx = zeros(n, 1);
segments.end_idx = zeros(n, 1);
segments.start_time = zeros(n, 1);
segments.end_time = zeros(n, 1);
segments.arrival_time = zeros(n, 1);
segments.peak_value = zeros(n, 1);
for k = 1:n
    a = idx(seg_start_pos(k));
    b = idx(seg_end_pos(k));
    segments.start_idx(k) = a;
    segments.end_idx(k) = b;
    segments.start_time(k) = t(a);
    segments.end_time(k) = t(b);
    [pk, ipk] = max(v(a:b));
    segments.peak_value(k) = pk;
    if use_centroid
        tt = t(a:b);
        vv = max(v(a:b) - threshold, 0);
        if sum(vv) > eps
            segments.arrival_time(k) = sum(tt(:) .* vv(:)) / sum(vv);
        else
            segments.arrival_time(k) = t(a + ipk - 1);
        end
    else
        segments.arrival_time(k) = t(a);
    end
end
end

function baseline = estimate_background_baseline_local(v, segments)
mask = true(size(v));
for k = 1:numel(segments.start_idx)
    a = segments.start_idx(k);
    b = segments.end_idx(k);
    pad = max(10, round(1.5 * (b - a + 1)));
    aa = max(1, a - pad);
    bb = min(numel(v), b + pad);
    mask(aa:bb) = false;
end
if nnz(mask) > 100
    baseline = median(v(mask), 'omitnan');
else
    baseline = median(v(v <= prctile(v, 30)), 'omitnan');
end
end

function F_omega_deg = build_phase_speed_local(opr_times, blades_num)
spd_t = opr_times(1:(end - blades_num));
spd_v = 360 ./ max(opr_times((blades_num + 1):end) - opr_times(1:(end - blades_num)), eps);
F_omega_deg = griddedInterpolant(spd_t, spd_v, 'linear', 'nearest');
end

function theta_points = map_segment_to_relative_angle_local(t_ref, t_seg, F_omega_deg)
theta_points = zeros(size(t_seg));
for i = 1:numel(t_seg)
    t_grid = linspace(t_ref, t_seg(i), 10);
    theta_points(i) = trapz(t_grid, F_omega_deg(t_grid));
end
end

function [x_wide, v_wide, pulse_id] = restore_low_speed_point_cloud_with_pulse_id_local( ...
    R, segments, selected_idx, opr_times, F_omega_deg, theta_std, r_tip_mm, expand_factor)
x_wide = [];
v_wide = [];
pulse_id = [];
for ii = 1:numel(selected_idx)
    pidx = selected_idx(ii);
    a = segments.start_idx(pidx);
    b = segments.end_idx(pidx);
    pad = max(8, round(expand_factor * (b - a + 1)));
    aa = max(1, a - pad);
    bb = min(numel(R.T), b + pad);
    t_seg = R.T(aa:bb);
    v_seg = R.V(aa:bb);
    t_arrival = segments.arrival_time(pidx);
    idx_prev = find(opr_times < t_arrival, 1, 'last');
    if isempty(idx_prev)
        continue;
    end
    theta_points_deg = map_segment_to_relative_angle_local(opr_times(idx_prev), t_seg, F_omega_deg);
    theta_diff_deg = mod(theta_points_deg - theta_std + 180, 360) - 180;
    x_points_mm = theta_diff_deg * (pi / 180) * r_tip_mm;
    x_wide = [x_wide; x_points_mm(:)]; %#ok<AGROW>
    v_wide = [v_wide; v_seg(:)]; %#ok<AGROW>
    pulse_id = [pulse_id; repmat(ii, numel(v_seg), 1)]; %#ok<AGROW>
end
end

function P = build_static_cloud_metrics_local(x_wide, v_wide, pulse_id, baseline, threshold, weight_floor)
finite = isfinite(x_wide) & isfinite(v_wide) & isfinite(pulse_id);
x_wide = x_wide(finite);
v_wide = v_wide(finite);
pulse_id = pulse_id(finite);
v_zero = max(v_wide - baseline, 0);
high = v_zero >= prctile(v_zero, 85) & v_zero > 0;
if nnz(high) >= 5 && sum(v_zero(high)) > eps
    xc = sum(x_wide(high) .* v_zero(high)) / sum(v_zero(high));
else
    [~, imax] = max(v_zero);
    xc = x_wide(imax);
end
xc = refine_center_by_sg_fit_local(x_wide(v_wide >= threshold), v_wide(v_wide >= threshold), baseline, xc);
x_rel = x_wide - xc;
[g_norm, amp_norm] = build_profile_metrics_local(x_rel, v_wide, baseline);
P = struct();
P.xc = xc;
P.x_rel = x_rel(:);
P.v = v_wide(:);
P.pulse_id = pulse_id(:);
P.g_norm = g_norm(:);
P.amp_norm = amp_norm(:);
P.valid = isfinite(x_rel(:)) & isfinite(v_wide(:)) & v_wide(:) >= 0.5 * threshold;
P.weight_floor = weight_floor;
P.local_xc_mad = estimate_local_center_mad_local(x_wide, v_wide, pulse_id, baseline, threshold);
end

function [g_norm, amp_norm] = build_profile_metrics_local(x_rel, v, baseline)
[x_sort, idx] = sort(x_rel(:));
v_zero = max(v(idx) - baseline, 0);
grid_n = min(700, max(120, round(numel(x_sort) / 150)));
edges = linspace(min(x_sort), max(x_sort), grid_n + 1).';
x_grid = 0.5 * (edges(1:end-1) + edges(2:end));
bin_id = discretize(x_sort, edges);
valid = ~isnan(bin_id);
v_med = accumarray(bin_id(valid), v_zero(valid), [grid_n, 1], @median, NaN);
v_fill = fillmissing(v_med, 'linear', 'EndValues', 'nearest');
span = min(41, 2 * floor((numel(v_fill) - 1) / 2) + 1);
if span >= 5
    v_smooth = smoothdata(v_fill, 'sgolay', span);
else
    v_smooth = v_fill;
end
g = abs(gradient(v_smooth, x_grid));
g_norm_grid = normalize01_local(g);
amp_norm_grid = normalize01_local(v_smooth);
g_norm = interp1(x_grid, g_norm_grid, x_rel(:), 'linear', 0);
amp_norm = interp1(x_grid, amp_norm_grid, x_rel(:), 'linear', 0);
end

function xc = refine_center_by_sg_fit_local(x, v, baseline, xc_initial)
finite = isfinite(x) & isfinite(v);
x = x(finite);
v = v(finite);
if numel(x) < 50
    xc = xc_initial;
    return;
end
[peak_val, idx_peak] = max(v);
B0 = max(peak_val - baseline, 0.1);
w0 = max(std(x), 0.2);
n0 = 3.0;
xc0 = xc_initial;
if ~isfinite(xc0)
    xc0 = x(idx_peak);
end
sg_model = @(p, xx) p(1) .* exp(-abs((xx - p(4)) ./ max(p(2), 1e-6)).^max(p(3), 1e-6)) + baseline;
obj_fun = @(p) sum((v - sg_model(p, x)).^2);
p0 = [B0, w0, n0, xc0];
lb = [0.05, 0.05, 1.2, min(x) - 0.5];
ub = [max(10, 2 * B0 + 0.5), 6.0, 10.0, max(x) + 0.5];
try
    if exist('fmincon', 'file') == 2
        opts = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp');
        p_opt = fmincon(obj_fun, p0, [], [], [], [], lb, ub, [], opts);
    else
        penalty_obj = @(p) obj_fun(min(max(p, lb), ub));
        opts = optimset('Display', 'off', 'MaxIter', 300, 'MaxFunEvals', 900);
        p_opt = fminsearch(penalty_obj, p0, opts);
        p_opt = min(max(p_opt, lb), ub);
    end
    xc = p_opt(4);
    if ~isfinite(xc)
        xc = xc_initial;
    end
catch
    xc = xc_initial;
end
end

function y = normalize01_local(x)
x = x(:);
xmin = min(x, [], 'omitnan');
xmax = max(x, [], 'omitnan');
if ~isfinite(xmin) || ~isfinite(xmax) || xmax <= xmin
    y = zeros(size(x));
else
    y = (x - xmin) ./ (xmax - xmin);
end
end

function local_xc_mad = estimate_local_center_mad_local(x, v, pulse_id, baseline, threshold)
ids = unique(pulse_id(:));
xc = nan(numel(ids), 1);
for i = 1:numel(ids)
    mask = pulse_id == ids(i);
    xv = x(mask);
    vv = v(mask);
    vz = max(vv - baseline, 0);
    high = vv >= threshold;
    if nnz(high) < 5
        high = vz >= prctile(vz, 75);
    end
    if nnz(high) >= 3 && sum(vz(high)) > eps
        xc(i) = sum(xv(high) .* vz(high)) / sum(vz(high));
    end
end
local_xc_mad = mad(xc, 1);
end

function [x_abs, v, t] = collect_dynamic_points_local(DynamicMap, sid)
x_abs = [];
v = [];
t = [];
for iw = 1:numel(DynamicMap.Window)
    D = DynamicMap.Window(iw).Sensor([DynamicMap.Window(iw).Sensor.sensor_id] == sid);
    if isempty(D)
        continue;
    end
    x_abs = [x_abs; D.x_abs(:)]; %#ok<AGROW>
    v = [v; D.V(:)]; %#ok<AGROW>
    t = [t; D.t(:)]; %#ok<AGROW>
end
finite = isfinite(x_abs) & isfinite(v) & isfinite(t);
x_abs = x_abs(finite);
v = v(finite);
t = t(finite);
end

function mask = build_dynamic_effective_mask_local(x_rel, v, t, P, threshold)
finite = isfinite(x_rel) & isfinite(v) & isfinite(t);
g_static = interp1(P.x_rel, P.g_norm, x_rel(:), 'linear', 0);
mask_static = g_static >= 0.08;
g_time = zeros(size(v));
if nnz(finite) >= 5 && range(t(finite)) > 0
    g_time(finite) = abs(gradient(v(finite), t(finite)));
end
g_time = normalize01_local(g_time);
mask_time = g_time >= 0.15;
mask_peak = v >= max(threshold, prctile(v(finite), 85));
mask = finite & (mask_static | mask_time | mask_peak) & v >= 0.5 * threshold;
if nnz(mask) < 30
    mask = finite & v >= threshold;
end
end

function row = evaluate_xrange_candidate_local(P, dyn_rel, dynamic_query, mask, name, mode, gr, ar, hw, query_guard_mm)
row = table();
row.CandidateName = string(name);
row.Mode = string(mode);
row.GradientRatio = gr;
row.AmplitudeRatio = ar;
row.HalfWidthMM = hw;
if nnz(mask) < 50
    row.PointCount = nnz(mask);
    row.x_left = NaN;
    row.x_right = NaN;
    row.WidthMM = NaN;
    row.static_rmse = NaN;
    row.edge_contamination = NaN;
    row.gradient_balance = NaN;
    row.safe_frac = NaN;
    row.coverage_deficit = inf;
    row.mean_gradient = NaN;
    row.Score = inf;
    return;
end
x = P.x_rel(mask);
v = P.v(mask);
w = max(P.weight_floor, P.g_norm(mask));
row.PointCount = nnz(mask);
row.x_left = min(x);
row.x_right = max(x);
row.WidthMM = row.x_right - row.x_left;
row.static_rmse = estimate_static_rmse_local(x, v, w);
edge_band = max(0.04, 0.03 * row.WidthMM);
row.edge_contamination = mean(x <= row.x_left + edge_band | x >= row.x_right - edge_band);
left_g = mean(P.g_norm(mask & P.x_rel < 0), 'omitnan');
right_g = mean(P.g_norm(mask & P.x_rel > 0), 'omitnan');
row.gradient_balance = abs(left_g - right_g);
row.safe_frac = mean(dyn_rel >= row.x_left + query_guard_mm & dyn_rel <= row.x_right - query_guard_mm);
row.coverage_deficit = max(row.x_left - dynamic_query(1), 0) + max(dynamic_query(2) - row.x_right, 0);
row.mean_gradient = mean(P.g_norm(mask), 'omitnan');
end

function rmse = estimate_static_rmse_local(x, v, w)
if numel(x) < 50 || range(x) <= 0
    rmse = NaN;
    return;
end
grid_n = min(401, max(81, round(numel(x) / 80)));
edges = linspace(min(x), max(x), grid_n + 1).';
x_grid = 0.5 * (edges(1:end-1) + edges(2:end));
bin_id = discretize(x, edges);
valid = ~isnan(bin_id);
v_med = accumarray(bin_id(valid), v(valid), [grid_n, 1], @median, NaN);
valid_grid = isfinite(v_med);
if nnz(valid_grid) < 20
    rmse = NaN;
    return;
end
v_fit = interp1(x_grid(valid_grid), v_med(valid_grid), x, 'linear', 'extrap');
rmse = sqrt(sum(w(:) .* (v(:) - v_fit(:)).^2) / max(sum(w), eps));
end

function T = score_xrange_candidates_local(T)
rmse_ref = max(median(T.static_rmse, 'omitnan'), 0.02);
edge_ref = max(median(T.edge_contamination, 'omitnan'), 0.02);
width_target = 6.0;
width_penalty = abs(T.WidthMM - width_target) / width_target;
T.Score = 4.0 * T.coverage_deficit + ...
    3.0 * (1 - T.safe_frac) + ...
    1.0 * T.static_rmse ./ rmse_ref + ...
    0.8 * T.edge_contamination ./ edge_ref + ...
    0.7 * T.gradient_balance + ...
    0.6 * width_penalty - ...
    0.5 * T.mean_gradient;
T.Score(~isfinite(T.Score)) = inf;
end

function chosen = choose_xrange_candidate_local(T)
T = T(isfinite(T.Score), :);
if isempty(T)
    error('No finite x-range candidates.');
end
T = sortrows(T, {'Score', 'coverage_deficit', 'static_rmse'}, {'ascend', 'ascend', 'ascend'});
chosen = T(1, :);
end

function plot_xrange_tradeoff_local(T, Rec, figure_dir, target_blade, sensor_tag)
fig = figure('Name', 'Static x-range tradeoff', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 18, 11]);
tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
scatter(T.WidthMM, T.static_rmse, 22, T.Score, 'filled'); hold on;
scatter(Rec.WidthMM, Rec.static_rmse, 80, 'kp', 'filled');
xlabel('Width (mm)');
ylabel('Static RMSE (V)');
title('Fit vs. width');
colorbar;
style_axes_local();

nexttile;
scatter(T.WidthMM, T.safe_frac, 22, T.Score, 'filled'); hold on;
scatter(Rec.WidthMM, Rec.safe_frac, 80, 'kp', 'filled');
xlabel('Width (mm)');
ylabel('Safe fraction');
title('Dynamic coverage');
colorbar;
style_axes_local();

nexttile;
scatter(T.edge_contamination, T.mean_gradient, 22, T.Score, 'filled'); hold on;
scatter(Rec.edge_contamination, Rec.mean_gradient, 80, 'kp', 'filled');
xlabel('Edge contamination');
ylabel('Mean gradient');
title('Sensitivity vs. edge');
colorbar;
style_axes_local();

nexttile;
scatter(T.GradientRatio, T.AmplitudeRatio, 22, T.Score, 'filled'); hold on;
mask_thr = Rec.Mode == "threshold";
scatter(Rec.GradientRatio(mask_thr), Rec.AmplitudeRatio(mask_thr), 80, 'kp', 'filled');
xlabel('Gradient threshold');
ylabel('Amplitude threshold');
title('Threshold candidates');
colorbar;
style_axes_local();

exportgraphics(fig, fullfile(figure_dir, sprintf('Step01_01_StaticXRangeTradeoff_B%d_%s_20250527.png', ...
    target_blade, sensor_tag)), 'Resolution', 300);
exportgraphics(fig, fullfile(figure_dir, sprintf('Step01_01_StaticXRangeTradeoff_B%d_%s_20250527.pdf', ...
    target_blade, sensor_tag)), 'ContentType', 'vector');
end

function plot_xrange_selected_clouds_local(plot_cache, figure_dir, target_blade, sensor_tag)
fig = figure('Name', 'Recommended static x-range', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 18, 4 + 4 * numel(plot_cache)]);
tiledlayout(fig, numel(plot_cache), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:numel(plot_cache)
    P = plot_cache(i);
    nexttile;
    [x_ds, v_ds, g_ds] = downsample_three_local(P.x_rel, P.v, P.g_norm, 30000);
    scatter(x_ds, v_ds, 4, g_ds, 'filled'); hold on;
    xline(P.chosen.x_left, 'k-.', 'LineWidth', 1.1);
    xline(P.chosen.x_right, 'k-.', 'LineWidth', 1.1);
    xline(P.dynamic_query(1), '-', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.0);
    xline(P.dynamic_query(2), '-', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.0);
    xlabel('x relative to center (mm)');
    ylabel('Voltage (V)');
    title(sprintf('CH%d: %s', P.sensor_id, char(P.chosen.CandidateName)));
    cb = colorbar;
    cb.Label.String = 'Normalized |dV/dx|';
    style_axes_local();
end
exportgraphics(fig, fullfile(figure_dir, sprintf('Step01_01_StaticXRangeSelected_B%d_%s_20250527.png', ...
    target_blade, sensor_tag)), 'Resolution', 300);
exportgraphics(fig, fullfile(figure_dir, sprintf('Step01_01_StaticXRangeSelected_B%d_%s_20250527.pdf', ...
    target_blade, sensor_tag)), 'ContentType', 'vector');
end

function [x_ds, v_ds, g_ds] = downsample_three_local(x, v, g, max_points)
keep = isfinite(x) & isfinite(v) & isfinite(g);
x = x(keep);
v = v(keep);
g = g(keep);
if numel(x) > max_points
    idx = round(linspace(1, numel(x), max_points));
else
    idx = 1:numel(x);
end
x_ds = x(idx);
v_ds = v(idx);
g_ds = g(idx);
end

function style_axes_local()
set(gca, 'FontName', 'Times New Roman', 'FontSize', 8.5, 'TickDir', 'in', 'Box', 'on');
end
