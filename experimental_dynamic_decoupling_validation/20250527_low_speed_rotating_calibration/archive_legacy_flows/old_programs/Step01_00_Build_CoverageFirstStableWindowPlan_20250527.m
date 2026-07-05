%% Step01_00_Build_CoverageFirstStableWindowPlan_20250527
% Build a coverage-first stable-window plan for Step01.
%
% The scan treats static-window stability as a constraint and prioritizes
% long windows whose gradient-informed template domain covers the dynamic
% query region used later by Step03.

clear; clc; close all;

%% Settings
script_dir = fileparts(mfilename('fullpath'));
route_dir = fileparts(fileparts(script_dir));
validation_root = fileparts(route_dir);
legacy_dir = fullfile(validation_root, '20250527_low_speed_gap_prior_decoupling', 'legacy');
if exist(legacy_dir, 'dir') ~= 7
    error('Legacy 20250527 helper folder not found: %s', legacy_dir);
end
addpath(legacy_dir);

legacy_cfg = Get_20250527_BTT_Config();
low_speed_case = legacy_cfg.low_speed_case;
low_speed_data_dir = fullfile(legacy_cfg.dataset_root, low_speed_case);
sensor_config_file = fullfile(legacy_cfg.reference_output_dir, 'Sensor_Config_20250527.mat');

target_blade = parse_positive_integer_env_local('STEP01_TARGET_BLADE', 1);
analysis_sensors = [1, 3, 6];
sensor_override = strtrim(getenv('STEP01_ANALYSIS_SENSORS'));
if ~isempty(sensor_override)
    parsed_sensors = sscanf(sensor_override, '%d').';
    if isempty(parsed_sensors)
        error('STEP01_ANALYSIS_SENSORS must contain integer sensor IDs, for example "1 3 6".');
    end
    analysis_sensors = parsed_sensors;
end
sensor_tag = ['S', sprintf('%d', analysis_sensors)];

candidate_lap_counts = parse_integer_list_env_local('STEP01_PLAN_LAP_COUNTS', [30, 40, 50]);
candidate_lap_counts = unique(candidate_lap_counts(candidate_lap_counts >= 10), 'stable');
scan_step_laps = parse_positive_integer_env_local('STEP01_PLAN_SCAN_STEP_LAPS', 2);
max_offset_lap = parse_nonnegative_numeric_env_local('STEP01_PLAN_MAX_OFFSET_LAP', inf);
scan_suffix = sanitize_scan_suffix_local(strtrim(getenv('STEP01_PLAN_SCAN_SUFFIX')));
query_guard_mm = parse_nonnegative_numeric_env_local('STEP01_PLAN_QUERY_GUARD_MM', 0.75);
dynamic_quantile = [1, 99];
gradient_energy_quantile = parse_fraction_env_local('STEP01_PLAN_GRADIENT_ENERGY_QUANTILE', 0.995);
gradient_edge_margin_mm = parse_nonnegative_numeric_env_local('STEP01_PLAN_GRADIENT_EDGE_MARGIN_MM', 0.02);
weight_floor = 0.05;
min_dynamic_effective_points = 30;

output_dir = fullfile(route_dir, 'output');
plan_dir = fullfile(output_dir, 'stable_window_plan');
figure_dir = fullfile(output_dir, 'figures');
dynamic_dir = fullfile(output_dir, 'dynamic_maps');
if exist(plan_dir, 'dir') ~= 7; mkdir(plan_dir); end
if exist(figure_dir, 'dir') ~= 7; mkdir(figure_dir); end

plan_csv = fullfile(plan_dir, sprintf('Step01_CoverageFirstStableWindowPlan_B%d_%s%s_20250527.csv', ...
    target_blade, sensor_tag, scan_suffix));
plan_mat = fullfile(plan_dir, sprintf('Step01_CoverageFirstStableWindowPlan_B%d_%s%s_20250527.mat', ...
    target_blade, sensor_tag, scan_suffix));
candidate_csv = fullfile(plan_dir, sprintf('Step01_CoverageFirstStableWindowCandidates_B%d_%s%s_20250527.csv', ...
    target_blade, sensor_tag, scan_suffix));

fprintf('\n=== Step01_00: coverage-first stable-window scan ===\n');
fprintf('Blade: %d\n', target_blade);
fprintf('Sensors: %s\n', mat2str(analysis_sensors));
fprintf('Candidate lap counts: %s, scan step: %d lap(s)\n', mat2str(candidate_lap_counts), scan_step_laps);
if isfinite(max_offset_lap)
    fprintf('Max start offset: %d lap(s)\n', floor(max_offset_lap));
end
fprintf('Dynamic query guard: %.3f mm; gradient energy quantile: %.4f\n', ...
    query_guard_mm, gradient_energy_quantile);

if ~isfolder(low_speed_data_dir)
    error('Low-speed data folder not found: %s', low_speed_data_dir);
end
if ~isfile(sensor_config_file)
    error('Sensor_Config not found. Run legacy Step1 once: %s', sensor_config_file);
end

loaded_cfg = load(sensor_config_file, 'Sensor_Config');
Sensor_Config = loaded_cfg.Sensor_Config;
[raw_stream, opr_times] = load_low_speed_streams_local(low_speed_data_dir, legacy_cfg, analysis_sensors);
F_omega_deg = build_phase_speed_local(opr_times, legacy_cfg.blades_num);

dynamic_file = find_dynamic_file_for_sensors_local(dynamic_dir, target_blade, analysis_sensors);
loaded_dynamic = load(dynamic_file, 'DynamicMap');
DynamicMap = loaded_dynamic.DynamicMap;
fprintf('Dynamic map source: %s\n', dynamic_file);

all_candidates = table();
selected_rows = table();
plot_cache = struct('sensor_id', {}, 'x_selected', {}, 'v_selected', {}, ...
    'x_wide', {}, 'v_wide', {}, 'weight_selected', {}, 'x_domain', {}, ...
    'dynamic_core', {}, 'dynamic_query', {}, 'threshold', {}, 'xc', {});

for is = 1:numel(analysis_sensors)
    sid = analysis_sensors(is);
    threshold = legacy_cfg.sensor_thresholds(sid);
    theta_std = Sensor_Config.Standard_Relative_Angles(sid, target_blade);
    if ~isfinite(theta_std)
        error('Missing standard relative angle for CH%d blade %d.', sid, target_blade);
    end

    R = raw_stream(sid);
    segments = extract_sensor_segments_local(R.T, R.V, threshold, legacy_cfg.gap_points);
    if isempty(segments.start_idx)
        error('No low-speed pulses detected for CH%d.', sid);
    end
    baseline = estimate_background_baseline_local(R.V, segments);
    if isKey(Sensor_Config.Target_Indices, sid)
        start_idx = Sensor_Config.Target_Indices(sid);
    else
        error('Sensor_Config.Target_Indices does not contain CH%d.', sid);
    end

    [dyn_x_abs, dyn_v, dyn_t] = collect_dynamic_points_local(DynamicMap, sid);
    dyn_effective = build_dynamic_effective_mask_local(dyn_x_abs, dyn_v, dyn_t, threshold, min_dynamic_effective_points);
    dyn_x_abs_eff = dyn_x_abs(dyn_effective);
    if numel(dyn_x_abs_eff) < min_dynamic_effective_points
        error('Too few effective dynamic points for CH%d.', sid);
    end

    max_target_pulses = floor((numel(segments.start_idx) - start_idx - target_blade + 1) / legacy_cfg.blades_num) + 1;
    pulse_cloud = build_target_pulse_cloud_cache_local( ...
        R, segments, start_idx, max_target_pulses, target_blade, ...
        legacy_cfg.blades_num, opr_times, F_omega_deg, theta_std, ...
        legacy_cfg.r_tip_mm, 0.30);
    fprintf('CH%d: cached %d target-blade low-speed pulses.\n', sid, numel(pulse_cloud));

    sensor_candidates = table();
    for lap_count = candidate_lap_counts(:).'
        max_offset = max_target_pulses - lap_count;
        if max_offset < 0
            continue;
        end
        max_offset_scan = min(max_offset, floor(max_offset_lap));
        for start_offset_lap = 0:scan_step_laps:max_offset_scan
            selected_lap_ids = start_offset_lap + (1:lap_count);
            if selected_lap_ids(end) > numel(pulse_cloud)
                continue;
            end

            [x_wide, v_wide, pulse_id] = concat_pulse_clouds_local(pulse_cloud, selected_lap_ids);

            select_cfg = struct('threshold', threshold, 'baseline', baseline, ...
                'weight_floor', weight_floor, ...
                'gradient_energy_quantile', gradient_energy_quantile, ...
                'gradient_edge_margin_mm', gradient_edge_margin_mm);
            P = select_gradient_template_points_local(x_wide, v_wide, pulse_id, select_cfg);
            if numel(P.x_selected) < 50
                continue;
            end

            dyn_rel = dyn_x_abs_eff - P.xc;
            dynamic_core = prctile(dyn_rel, dynamic_quantile);
            dynamic_query = [dynamic_core(1) - query_guard_mm, dynamic_core(2) + query_guard_mm];
            static_domain = [min(P.x_selected - P.xc), max(P.x_selected - P.xc)];
            safe_frac = mean(dyn_rel >= static_domain(1) + query_guard_mm & ...
                             dyn_rel <= static_domain(2) - query_guard_mm);
            coverage_deficit = max(static_domain(1) - dynamic_query(1), 0) + ...
                               max(dynamic_query(2) - static_domain(2), 0);
            left_margin = dynamic_query(1) - static_domain(1);
            right_margin = static_domain(2) - dynamic_query(2);

            candidate = table();
            candidate.BladeID = target_blade;
            candidate.SensorID = sid;
            candidate.LapCount = numel(selected_lap_ids);
            candidate.StartOffsetLap = start_offset_lap;
            candidate.StartTime = pulse_cloud(selected_lap_ids(1)).arrival_time;
            candidate.EndTime = pulse_cloud(selected_lap_ids(end)).arrival_time;
            candidate.PointCount = numel(P.x_selected);
            candidate.WidePointCount = numel(P.x_wide);
            candidate.xc = P.xc;
            candidate.local_xc_mad = P.local_xc_mad;
            candidate.rmse = P.rmse;
            candidate.bound_ratio = P.bound_ratio;
            candidate.left_grad_ratio = P.left_grad_ratio;
            candidate.right_grad_ratio = P.right_grad_ratio;
            candidate.static_x_left = static_domain(1);
            candidate.static_x_right = static_domain(2);
            candidate.dynamic_core_left = dynamic_core(1);
            candidate.dynamic_core_right = dynamic_core(2);
            candidate.dynamic_query_left = dynamic_query(1);
            candidate.dynamic_query_right = dynamic_query(2);
            candidate.query_guard_mm = query_guard_mm;
            candidate.safe_frac = safe_frac;
            candidate.coverage_deficit = coverage_deficit;
            candidate.left_margin = left_margin;
            candidate.right_margin = right_margin;
            sensor_candidates = [sensor_candidates; candidate]; %#ok<AGROW>
        end
    end

    if isempty(sensor_candidates)
        error('No valid stable-window candidates were built for CH%d.', sid);
    end

    sensor_candidates = score_candidates_local(sensor_candidates);
    chosen = choose_sensor_candidate_local(sensor_candidates);
    selected_rows = [selected_rows; chosen]; %#ok<AGROW>
    all_candidates = [all_candidates; sensor_candidates]; %#ok<AGROW>

    selected_lap_ids = chosen.StartOffsetLap + (1:chosen.LapCount);
    [x_wide, v_wide, pulse_id] = concat_pulse_clouds_local(pulse_cloud, selected_lap_ids);
    select_cfg = struct('threshold', threshold, 'baseline', baseline, ...
        'weight_floor', weight_floor, ...
        'gradient_energy_quantile', gradient_energy_quantile, ...
        'gradient_edge_margin_mm', gradient_edge_margin_mm);
    P = select_gradient_template_points_local(x_wide, v_wide, pulse_id, select_cfg);

    pc = struct();
    pc.sensor_id = sid;
    pc.x_selected = P.x_selected(:) - P.xc;
    pc.v_selected = P.v_selected(:);
    pc.x_wide = P.x_wide(:) - P.xc;
    pc.v_wide = P.v_wide(:);
    pc.weight_selected = P.weight_selected(:);
    pc.x_domain = [min(pc.x_selected), max(pc.x_selected)];
    pc.dynamic_core = [chosen.dynamic_core_left, chosen.dynamic_core_right];
    pc.dynamic_query = [chosen.dynamic_query_left, chosen.dynamic_query_right];
    pc.threshold = threshold;
    pc.xc = P.xc;
    plot_cache(end + 1) = pc; %#ok<SAGROW>

    fprintf(['CH%d selected: LapCount=%d, StartOffset=%d, score=%.4f, ' ...
        'safe=%.3f, deficit=%.3f mm, domain=[%.3f %.3f] mm, query=[%.3f %.3f] mm\n'], ...
        sid, chosen.LapCount, chosen.StartOffsetLap, chosen.ChosenScore, ...
        chosen.safe_frac, chosen.coverage_deficit, chosen.static_x_left, ...
        chosen.static_x_right, chosen.dynamic_query_left, chosen.dynamic_query_right);
end

window_plan = selected_rows(:, {'BladeID', 'SensorID', 'LapCount', 'StartOffsetLap'});
window_plan.StartTime = selected_rows.StartTime;
window_plan.EndTime = selected_rows.EndTime;
window_plan.ChosenScore = selected_rows.ChosenScore;
window_plan.xc = selected_rows.xc;
window_plan.local_xc_mad = selected_rows.local_xc_mad;
window_plan.rmse = selected_rows.rmse;
window_plan.bound_ratio = selected_rows.bound_ratio;
window_plan.safe_frac = selected_rows.safe_frac;
window_plan.coverage_deficit = selected_rows.coverage_deficit;
window_plan.static_x_left = selected_rows.static_x_left;
window_plan.static_x_right = selected_rows.static_x_right;
window_plan.dynamic_core_left = selected_rows.dynamic_core_left;
window_plan.dynamic_core_right = selected_rows.dynamic_core_right;
window_plan.dynamic_query_left = selected_rows.dynamic_query_left;
window_plan.dynamic_query_right = selected_rows.dynamic_query_right;
window_plan.SelectionRule = repmat("coverage_first_long_gradient_window", height(window_plan), 1);

writetable(all_candidates, candidate_csv);
writetable(window_plan, plan_csv);
save(plan_mat, 'window_plan', 'all_candidates', 'plot_cache', 'dynamic_file', '-v7.3');
fprintf('Saved candidates: %s\n', candidate_csv);
fprintf('Saved plan: %s\n', plan_csv);

plot_candidate_score_map_local(all_candidates, window_plan, figure_dir, target_blade, sensor_tag);
plot_domain_coverage_local(window_plan, figure_dir, target_blade, sensor_tag);
plot_selected_clouds_local(plot_cache, figure_dir, target_blade, sensor_tag);

%% Local functions
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

function vals = parse_integer_list_env_local(name, default_value)
raw = strtrim(getenv(name));
if isempty(raw)
    vals = default_value;
    return;
end
vals = sscanf(raw, '%d').';
if isempty(vals)
    error('%s must contain integer values separated by spaces.', name);
end
end

function suffix = sanitize_scan_suffix_local(raw)
if isempty(raw)
    suffix = '';
    return;
end
suffix = regexprep(raw, '[^A-Za-z0-9_\\-]', '_');
if suffix(1) ~= '_'
    suffix = ['_', suffix];
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

function v = parse_fraction_env_local(name, default_value)
v = parse_nonnegative_numeric_env_local(name, default_value);
if v <= 0 || v >= 1
    error('%s must be between 0 and 1.', name);
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

function mask = build_dynamic_effective_mask_local(x_abs, v, t, threshold, min_points)
mask = false(size(v));
if isempty(v)
    return;
end
if numel(v) >= 5 && range(t) > 0
    g = abs(gradient(v(:), t(:)));
    g_thr = prctile(g(isfinite(g)), 65);
    mask = (v(:) >= threshold) | (g(:) >= g_thr & v(:) >= prctile(v, 45));
else
    mask = v(:) >= threshold;
end
mask = mask & isfinite(x_abs(:)) & isfinite(v(:));
if nnz(mask) < min_points
    mask = v(:) >= prctile(v(:), 65);
end
end

function [raw_stream, opr_times] = load_low_speed_streams_local(case_dir, cfg, sensor_ids)
file_ids = list_case_file_ids_local(case_dir, cfg.opr_id);
if isempty(file_ids)
    error('No OPR files found in %s.', case_dir);
end
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

function pulse_cloud = build_target_pulse_cloud_cache_local( ...
    R, segments, start_idx, max_target_pulses, target_blade, blades_num, ...
    opr_times, F_omega_deg, theta_std, r_tip_mm, expand_factor)
pulse_cloud = repmat(struct('x', [], 'v', [], 'arrival_time', NaN, 'segment_idx', NaN), max_target_pulses, 1);
for lap_id = 1:max_target_pulses
    pidx = start_idx + (lap_id - 1) * blades_num + (target_blade - 1);
    pulse_cloud(lap_id).segment_idx = pidx;
    if pidx < 1 || pidx > numel(segments.start_idx)
        continue;
    end
    a = segments.start_idx(pidx);
    b = segments.end_idx(pidx);
    pad = max(8, round(expand_factor * (b - a + 1)));
    aa = max(1, a - pad);
    bb = min(numel(R.T), b + pad);
    t_seg = R.T(aa:bb);
    v_seg = R.V(aa:bb);
    t_arrival = segments.arrival_time(pidx);
    idx_prev = find(opr_times < t_arrival, 1, 'last');
    pulse_cloud(lap_id).arrival_time = t_arrival;
    if isempty(idx_prev)
        continue;
    end
    theta_points_deg = map_segment_to_relative_angle_local(opr_times(idx_prev), t_seg, F_omega_deg);
    theta_diff_deg = mod(theta_points_deg - theta_std + 180, 360) - 180;
    pulse_cloud(lap_id).x = theta_diff_deg(:) * (pi / 180) * r_tip_mm;
    pulse_cloud(lap_id).v = v_seg(:);
end
end

function [x_wide, v_wide, pulse_id] = concat_pulse_clouds_local(pulse_cloud, selected_lap_ids)
x_wide = [];
v_wide = [];
pulse_id = [];
for ii = 1:numel(selected_lap_ids)
    lap_id = selected_lap_ids(ii);
    if lap_id < 1 || lap_id > numel(pulse_cloud) || isempty(pulse_cloud(lap_id).x)
        continue;
    end
    x = pulse_cloud(lap_id).x(:);
    v = pulse_cloud(lap_id).v(:);
    x_wide = [x_wide; x]; %#ok<AGROW>
    v_wide = [v_wide; v]; %#ok<AGROW>
    pulse_id = [pulse_id; repmat(ii, numel(x), 1)]; %#ok<AGROW>
end
end

function P = select_gradient_template_points_local(x_wide, v_wide, pulse_id, cfg)
finite = isfinite(x_wide) & isfinite(v_wide) & isfinite(pulse_id);
x_wide = x_wide(finite);
v_wide = v_wide(finite);
pulse_id = pulse_id(finite);
v_zero = max(v_wide - cfg.baseline, 0);
if numel(x_wide) < 50
    P = empty_template_points_local();
    return;
end
high_level = prctile(v_zero, 85);
high_mask = v_zero >= high_level & v_zero > 0;
if nnz(high_mask) >= 5 && sum(v_zero(high_mask)) > eps
    xc = sum(x_wide(high_mask) .* v_zero(high_mask)) / sum(v_zero(high_mask));
else
    [~, imax] = max(v_zero);
    xc = x_wide(imax);
end
xc = refine_center_by_sg_fit_local(x_wide(v_wide >= cfg.threshold), ...
    v_wide(v_wide >= cfg.threshold), cfg.baseline, xc);

[left_L, right_L, gradient_ok] = estimate_gradient_energy_span_local( ...
    x_wide, v_wide, cfg.baseline, xc, cfg.gradient_energy_quantile, cfg.gradient_edge_margin_mm);
strict_mask = v_wide >= cfg.threshold;
if gradient_ok
    trust_mask = strict_mask & x_wide >= xc - left_L & x_wide <= xc + right_L;
else
    dx = x_wide(strict_mask) - xc;
    if isempty(dx)
        trust_mask = v_zero >= prctile(v_zero, 60);
    else
        L = max(prctile(abs(dx), 99), 0.1);
        trust_mask = strict_mask & abs(x_wide - xc) <= L;
    end
end
if nnz(trust_mask) < 30
    trust_mask = strict_mask;
end

x_selected = x_wide(trust_mask);
v_selected = v_wide(trust_mask);
pid_selected = pulse_id(trust_mask);
weight_selected = build_gradient_weight_from_wide_cloud_local( ...
    x_selected, x_wide, v_wide, cfg.baseline, cfg.weight_floor);
[rmse, bound_ratio] = estimate_static_fit_quality_local(x_selected - xc, v_selected, weight_selected);
local_xc_mad = estimate_local_center_mad_local(x_wide, v_wide, pulse_id, cfg.baseline, cfg.threshold);
[left_grad_ratio, right_grad_ratio] = estimate_gradient_side_balance_local( ...
    x_selected, weight_selected, xc);

P = struct();
P.xc = xc;
P.x_wide = x_wide;
P.v_wide = v_wide;
P.x_selected = x_selected;
P.v_selected = v_selected;
P.weight_selected = weight_selected;
P.pid_selected = pid_selected;
P.rmse = rmse;
P.bound_ratio = bound_ratio;
P.local_xc_mad = local_xc_mad;
P.left_grad_ratio = left_grad_ratio;
P.right_grad_ratio = right_grad_ratio;
end

function P = empty_template_points_local()
P = struct('xc', NaN, 'x_wide', [], 'v_wide', [], 'x_selected', [], ...
    'v_selected', [], 'weight_selected', [], 'pid_selected', [], 'rmse', NaN, ...
    'bound_ratio', NaN, 'local_xc_mad', NaN, 'left_grad_ratio', NaN, ...
    'right_grad_ratio', NaN);
end

function [left_L, right_L, ok] = estimate_gradient_energy_span_local( ...
    x, v, baseline, xc, energy_quantile, edge_margin_mm)
left_L = NaN;
right_L = NaN;
ok = false;
finite = isfinite(x) & isfinite(v);
x = x(finite);
v_zero = max(v(finite) - baseline, 0);
if numel(x) < 50 || max(v_zero) <= 0
    return;
end
[x_sort, idx] = sort(x);
v_sort = v_zero(idx);
grid_n = min(700, max(120, round(numel(x_sort) / 120)));
edges = linspace(min(x_sort), max(x_sort), grid_n + 1).';
x_grid = 0.5 * (edges(1:end-1) + edges(2:end));
bin_id = discretize(x_sort, edges);
valid = ~isnan(bin_id);
v_med = accumarray(bin_id(valid), v_sort(valid), [grid_n, 1], @median, NaN);
if nnz(isfinite(v_med)) < 15
    return;
end
v_fill = fillmissing(v_med, 'linear', 'EndValues', 'nearest');
span = min(41, 2 * floor((numel(v_fill) - 1) / 2) + 1);
if span >= 5
    v_smooth = smoothdata(v_fill, 'sgolay', span);
else
    v_smooth = v_fill;
end
g = abs(gradient(v_smooth, x_grid));
signal_gate = v_smooth >= prctile(v_smooth, 40);
energy = g(:) .* signal_gate(:) + 0.05 * max(g(:)) * normalize01_local(v_smooth(:));
if sum(energy) <= eps
    return;
end
tail = (1 - energy_quantile) / 2;
cum = cumsum(energy) ./ sum(energy);
[cum_unique, ia] = unique(cum, 'stable');
x_unique = x_grid(ia);
if numel(cum_unique) < 2
    return;
end
x_lo = interp1(cum_unique, x_unique, tail, 'linear', 'extrap');
x_hi = interp1(cum_unique, x_unique, 1 - tail, 'linear', 'extrap');
if ~isfinite(x_lo) || ~isfinite(x_hi) || x_hi <= x_lo || xc <= x_lo || xc >= x_hi
    return;
end
left_L = max(xc - x_lo - edge_margin_mm, 0.1);
right_L = max(x_hi - xc - edge_margin_mm, 0.1);
ok = true;
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

function W = build_gradient_weight_from_wide_cloud_local(x_query, x_wide, v_wide, baseline, weight_floor)
finite = isfinite(x_wide) & isfinite(v_wide);
x_wide = x_wide(finite);
v_zero = max(v_wide(finite) - baseline, 0);
if numel(x_wide) < 20 || all(v_zero == 0)
    W = ones(size(x_query));
    return;
end
[x_sort, idx] = sort(x_wide);
v_sort = v_zero(idx);
grid_n = min(500, max(80, round(numel(x_sort) / 200)));
edges = linspace(min(x_sort), max(x_sort), grid_n + 1).';
x_grid = 0.5 * (edges(1:end-1) + edges(2:end));
bin_id = discretize(x_sort, edges);
valid = ~isnan(bin_id);
v_med = accumarray(bin_id(valid), v_sort(valid), [grid_n, 1], @median, NaN);
valid_grid = isfinite(v_med);
if nnz(valid_grid) < 10
    W = ones(size(x_query));
    return;
end
v_fill = fillmissing(v_med, 'linear', 'EndValues', 'nearest');
span = min(41, 2 * floor((numel(v_fill) - 1) / 2) + 1);
if span >= 5
    v_smooth = smoothdata(v_fill, 'sgolay', span);
else
    v_smooth = v_fill;
end
dVdx = abs(gradient(v_smooth, x_grid));
w_grid = weight_floor + (1 - weight_floor) * dVdx ./ max(dVdx + eps);
W = interp1(x_grid, w_grid, x_query(:), 'linear', weight_floor);
W(~isfinite(W)) = weight_floor;
W = max(W, weight_floor);
end

function [rmse, bound_ratio] = estimate_static_fit_quality_local(x_rel, v, w)
finite = isfinite(x_rel) & isfinite(v) & isfinite(w);
x_rel = x_rel(finite);
v = v(finite);
w = w(finite);
if numel(x_rel) < 30 || range(x_rel) <= 0
    rmse = NaN;
    bound_ratio = NaN;
    return;
end
grid_n = min(401, max(81, round(numel(x_rel) / 60)));
edges = linspace(min(x_rel), max(x_rel), grid_n + 1).';
x_grid = 0.5 * (edges(1:end-1) + edges(2:end));
bin_id = discretize(x_rel, edges);
valid = ~isnan(bin_id);
v_med = accumarray(bin_id(valid), v(valid), [grid_n, 1], @median, NaN);
valid_grid = isfinite(v_med);
if nnz(valid_grid) < 20
    rmse = NaN;
else
    v_fit = interp1(x_grid(valid_grid), v_med(valid_grid), x_rel, 'linear', 'extrap');
    rmse = sqrt(sum(w(:) .* (v(:) - v_fit(:)).^2) / max(sum(w), eps));
end
edge_band = max(0.04, 0.03 * range(x_rel));
bound_ratio = mean(x_rel <= min(x_rel) + edge_band | x_rel >= max(x_rel) - edge_band);
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

function [left_ratio, right_ratio] = estimate_gradient_side_balance_local(x_selected, weight_selected, xc)
dx = x_selected(:) - xc;
w = max(weight_selected(:), 0);
left_sum = sum(w(dx < 0), 'omitnan');
right_sum = sum(w(dx > 0), 'omitnan');
total = left_sum + right_sum;
if total <= eps
    left_ratio = NaN;
    right_ratio = NaN;
else
    left_ratio = left_sum / total;
    right_ratio = right_sum / total;
end
end

function T = score_candidates_local(T)
rmse_ref = max(median(T.rmse, 'omitnan'), 0.02);
xc_ref = max(median(T.local_xc_mad, 'omitnan'), 0.002);
bound_ref = max(median(T.bound_ratio, 'omitnan'), 0.02);
lap_span = max(T.LapCount) - min(T.LapCount);
if lap_span <= 0
    lap_reward = zeros(height(T), 1);
else
    lap_reward = (T.LapCount - min(T.LapCount)) ./ lap_span;
end
balance_penalty = abs(T.left_grad_ratio - T.right_grad_ratio);
T.IsStableCandidate = T.local_xc_mad <= max(3 * xc_ref, 0.012) & ...
    T.rmse <= max(2.5 * rmse_ref, 0.10) & ...
    T.bound_ratio <= max(2.5 * bound_ref, 0.20);
score = 80 * T.coverage_deficit + ...
    6 * (1 - T.safe_frac) + ...
    1.5 * T.rmse ./ rmse_ref + ...
    1.2 * T.local_xc_mad ./ xc_ref + ...
    0.8 * T.bound_ratio ./ bound_ref + ...
    0.8 * balance_penalty - ...
    1.0 * lap_reward;
score(~T.IsStableCandidate) = score(~T.IsStableCandidate) + 25;
T.ChosenScore = score;
end

function chosen = choose_sensor_candidate_local(T)
stable = T(T.IsStableCandidate, :);
if isempty(stable)
    stable = T;
end
stable = sortrows(stable, {'coverage_deficit', 'ChosenScore', 'LapCount'}, {'ascend', 'ascend', 'descend'});
chosen = stable(1, :);
end

function plot_candidate_score_map_local(C, Plan, figure_dir, target_blade, sensor_tag)
fig = figure('Name', 'Coverage-first stable-window candidates', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 18, 5.8 + 3.2 * numel(unique(C.SensorID))]);
tiledlayout(fig, numel(unique(C.SensorID)), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
sids = unique(C.SensorID).';
for sid = sids
    nexttile;
    rows = C(C.SensorID == sid, :);
    scatter(rows.StartOffsetLap, rows.LapCount, 28, rows.ChosenScore, 'filled'); hold on;
    chosen = Plan(Plan.SensorID == sid, :);
    plot(chosen.StartOffsetLap, chosen.LapCount, 'kp', 'MarkerSize', 13, ...
        'MarkerFaceColor', 'y', 'DisplayName', 'selected');
    xlabel('Start offset (lap)');
    ylabel('Lap count');
    title(sprintf('CH%d candidate score', sid));
    cb = colorbar;
    cb.Label.String = 'Score';
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 8.5, 'TickDir', 'in', 'Box', 'on');
end
exportgraphics(fig, fullfile(figure_dir, sprintf('Step01_00_CandidateScore_B%d_%s_20250527.png', ...
    target_blade, sensor_tag)), 'Resolution', 300);
exportgraphics(fig, fullfile(figure_dir, sprintf('Step01_00_CandidateScore_B%d_%s_20250527.pdf', ...
    target_blade, sensor_tag)), 'ContentType', 'vector');
end

function plot_domain_coverage_local(Plan, figure_dir, target_blade, sensor_tag)
fig = figure('Name', 'Static-domain and dynamic-query coverage', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 18, 7.5]);
hold on;
for i = 1:height(Plan)
    y = i;
    plot([Plan.static_x_left(i), Plan.static_x_right(i)], [y, y], 'k-', 'LineWidth', 4);
    plot([Plan.dynamic_query_left(i), Plan.dynamic_query_right(i)], [y, y], ...
        '-', 'Color', [0.85 0.33 0.10], 'LineWidth', 2);
    plot([Plan.dynamic_core_left(i), Plan.dynamic_core_right(i)], [y, y], ...
        '-', 'Color', [0 0.45 0.74], 'LineWidth', 2);
    text(Plan.static_x_right(i) + 0.08, y, sprintf('safe=%.2f', Plan.safe_frac(i)), ...
        'FontName', 'Times New Roman', 'FontSize', 8, 'VerticalAlignment', 'middle');
end
ylim([0.5, height(Plan) + 0.5]);
yticks(1:height(Plan));
yticklabels(compose('CH%d', Plan.SensorID));
xlabel('x relative to candidate center (mm)');
title('Static template domain vs. dynamic query interval');
legend({'static domain', 'dynamic query', 'dynamic core'}, 'Location', 'best');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 8.5, 'TickDir', 'in', 'Box', 'on');
exportgraphics(fig, fullfile(figure_dir, sprintf('Step01_00_DomainCoverage_B%d_%s_20250527.png', ...
    target_blade, sensor_tag)), 'Resolution', 300);
exportgraphics(fig, fullfile(figure_dir, sprintf('Step01_00_DomainCoverage_B%d_%s_20250527.pdf', ...
    target_blade, sensor_tag)), 'ContentType', 'vector');
end

function plot_selected_clouds_local(plot_cache, figure_dir, target_blade, sensor_tag)
fig = figure('Name', 'Selected gradient-gated static clouds', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 18, 5 + 4 * numel(plot_cache)]);
tiledlayout(fig, numel(plot_cache), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:numel(plot_cache)
    P = plot_cache(i);
    nexttile;
    [xw, vw] = downsample_cloud_local(P.x_wide, P.v_wide, 20000);
    scatter(xw, vw, 3, [0.80 0.80 0.80], 'filled'); hold on;
    scatter(P.x_selected, P.v_selected, 6, P.weight_selected, 'filled');
    xline(P.x_domain(1), 'k-.', 'LineWidth', 1.0);
    xline(P.x_domain(2), 'k-.', 'LineWidth', 1.0);
    xline(P.dynamic_query(1), '-', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.1);
    xline(P.dynamic_query(2), '-', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.1);
    yline(P.threshold, '--', 'Color', [0.35 0.35 0.35]);
    xlabel('x relative to candidate center (mm)');
    ylabel('Voltage (V)');
    title(sprintf('CH%d gradient-gated static cloud', P.sensor_id));
    cb = colorbar;
    cb.Label.String = 'Gradient weight';
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 8.5, 'TickDir', 'in', 'Box', 'on');
end
exportgraphics(fig, fullfile(figure_dir, sprintf('Step01_00_SelectedClouds_B%d_%s_20250527.png', ...
    target_blade, sensor_tag)), 'Resolution', 300);
exportgraphics(fig, fullfile(figure_dir, sprintf('Step01_00_SelectedClouds_B%d_%s_20250527.pdf', ...
    target_blade, sensor_tag)), 'ContentType', 'vector');
end

function [x_ds, v_ds] = downsample_cloud_local(x, v, max_points)
keep = isfinite(x) & isfinite(v);
x = x(keep);
v = v(keep);
if numel(x) > max_points
    idx = round(linspace(1, numel(x), max_points));
else
    idx = 1:numel(x);
end
x_ds = x(idx);
v_ds = v(idx);
end
