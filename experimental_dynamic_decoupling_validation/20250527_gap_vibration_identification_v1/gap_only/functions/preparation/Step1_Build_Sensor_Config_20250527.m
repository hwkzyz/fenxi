function Sensor_Config = Step1_Build_Sensor_Config_20250527(show_plots)
%STEP1_BUILD_SENSOR_CONFIG_20250527 Build Step1-style fingerprint config for 20250527.

if nargin < 1
    show_plots = true;
end

cfg = Get_20250527_BTT_Config();
case_data = Extract_BTT_Features_20250527(cfg.low_speed_case, cfg);

if ~exist(cfg.reference_output_dir, 'dir')
    mkdir(cfg.reference_output_dir);
end

anchor_sid = cfg.reference_sensor_id;
anchor_feat = case_data.channels(anchor_sid).peak_features(:);
if numel(anchor_feat) < cfg.blades_num
    error('Reference sensor CH%d does not have enough pulses.', anchor_sid);
end

limit = min(numel(anchor_feat), cfg.blades_num * cfg.reference_search_revs);
[~, max_idx] = max(anchor_feat(1:limit));
if max_idx < cfg.blades_num
    max_idx = max_idx + cfg.blades_num;
end
anchor_start_idx = min(max(max_idx, 1), numel(anchor_feat) - cfg.blades_num + 1);
ref_fp = anchor_feat(anchor_start_idx:(anchor_start_idx + cfg.blades_num - 1));

fingerprints = containers.Map('KeyType', 'double', 'ValueType', 'any');
target_indices = containers.Map('KeyType', 'double', 'ValueType', 'double');
summary_rows = repmat(struct( ...
    'sensor_id', NaN, ...
    'start_index', NaN, ...
    'best_corr', NaN, ...
    'pulse_count', NaN), numel(cfg.sensor_ids), 1);

fingerprints(anchor_sid) = ref_fp(:).';
target_indices(anchor_sid) = anchor_start_idx;

for i = 1:numel(cfg.sensor_ids)
    sid = cfg.sensor_ids(i);
    feats = case_data.channels(sid).peak_features(:);
    if numel(feats) < cfg.blades_num
        error('CH%d does not have enough low-speed pulses for fingerprint alignment.', sid);
    end

    [best_idx, best_corr] = match_reference_fingerprint_local( ...
        feats, ref_fp, cfg.blades_num, cfg.max_alignment_search_starts);
    aligned_fp = feats(best_idx:(best_idx + cfg.blades_num - 1)).';

    fingerprints(sid) = aligned_fp;
    target_indices(sid) = best_idx;

    summary_rows(i).sensor_id = sid;
    summary_rows(i).start_index = best_idx;
    summary_rows(i).best_corr = best_corr;
    summary_rows(i).pulse_count = numel(feats);
end

std_angles_start_edge = zeros(max([cfg.sensor_ids, cfg.opr_id]), cfg.blades_num);
std_angles_start_edge(:) = NaN;
if numel(case_data.opr_times) > cfg.blades_num
    spd_t = case_data.opr_times(1:(end - cfg.blades_num));
    spd_v = 360 ./ (case_data.opr_times((cfg.blades_num + 1):end) - case_data.opr_times(1:(end - cfg.blades_num)));
    F_omega_deg = griddedInterpolant(spd_t, spd_v, 'linear', 'nearest');

    for sid = cfg.sensor_ids
        meas_times = case_data.channels(sid).arrival_times(:);
        start_idx = target_indices(sid);

        for blade_id = 1:cfg.blades_num
            rel_angles = [];
            for lap = 0:(cfg.max_laps_process - 1)
                pulse_idx = start_idx + lap * cfg.blades_num + (blade_id - 1);
                if pulse_idx > numel(meas_times)
                    break;
                end

                t_meas = meas_times(pulse_idx);
                idx_prev_opr = find(case_data.opr_times < t_meas, 1, 'last');
                if isempty(idx_prev_opr)
                    continue;
                end

                t_ref = case_data.opr_times(idx_prev_opr);
                t_grid = linspace(t_ref, t_meas, 10);
                rel_angles(end + 1) = trapz(t_grid, F_omega_deg(t_grid)); %#ok<AGROW>
            end

            if ~isempty(rel_angles)
                std_angles_start_edge(sid, blade_id) = mean(rel_angles);
            end
        end
    end
end

opr_reference = build_opr_center_reference_for_case_local(cfg, cfg.low_speed_case);
std_angles_opr_center = std_angles_start_edge - opr_reference.phase_shift_deg;

Sensor_Config = struct();
Sensor_Config.ReferenceCase = cfg.low_speed_case;
Sensor_Config.Sensor_IDs = cfg.sensor_ids;
Sensor_Config.Capacitance_IDs = cfg.capacitance_ids;
Sensor_Config.Eddy_Current_IDs = cfg.eddy_current_ids;
Sensor_Config.OPR_ID = cfg.opr_id;
Sensor_Config.Blades_Num = cfg.blades_num;
Sensor_Config.Pinlv = cfg.pinlv;
Sensor_Config.Gap_Points = cfg.gap_points;
Sensor_Config.OPR_Threshold = cfg.opr_threshold;
Sensor_Config.Sensor_Thresholds = cfg.sensor_thresholds;
Sensor_Config.Fingerprints = fingerprints;
Sensor_Config.Target_Indices = target_indices;
Sensor_Config.Standard_Relative_Angles = std_angles_opr_center;
Sensor_Config.Standard_Relative_Angles_OPRCenter = std_angles_opr_center;
Sensor_Config.Standard_Relative_Angles_StartEdge = std_angles_start_edge;
Sensor_Config.Standard_Relative_Angles_Reference = 'opr_pulse_center';
Sensor_Config.Standard_Relative_Angles_StartEdge_Reference = 'opr_threshold_start_edge';
Sensor_Config.OPRReference = opr_reference;

summary_table = struct2table(summary_rows);
disp(summary_table);

save(fullfile(cfg.reference_output_dir, 'Sensor_Config_20250527.mat'), 'Sensor_Config');
save(fullfile(cfg.reference_output_dir, 'Sensor_Config.mat'), 'Sensor_Config');
writetable(summary_table, fullfile(cfg.reference_output_dir, 'Sensor_Config_Summary_20250527.csv'));

fprintf('Saved Sensor_Config to:\n  %s\n', fullfile(cfg.reference_output_dir, 'Sensor_Config_20250527.mat'));
fprintf('Standard_Relative_Angles now use OPR pulse-center reference; start-edge values are kept in Standard_Relative_Angles_StartEdge.\n');
fprintf('OPR center shift: %.6f deg, %.6f mm, %.3f us.\n', ...
    opr_reference.phase_shift_deg, opr_reference.phase_shift_mm, ...
    1e6 * opr_reference.median_center_minus_start_s);

if show_plots
    plot_step1_results_local(case_data, Sensor_Config, cfg, summary_table);
end
end


function [best_idx, best_corr] = match_reference_fingerprint_local(feats, ref_fp, blades_num, max_starts)
max_start = min(max_starts, numel(feats) - blades_num + 1);
best_idx = 1;
best_corr = -inf;

for idx = 1:max_start
    seg = feats(idx:(idx + blades_num - 1));
    c = safe_corr_local(ref_fp(:), seg(:));
    if c > best_corr
        best_corr = c;
        best_idx = idx;
    end
end

if ~isfinite(best_corr)
    best_corr = NaN;
end
end


function c = safe_corr_local(a, b)
if numel(a) ~= numel(b) || numel(a) < 2
    c = -inf;
    return;
end

if all(abs(a - a(1)) < eps) || all(abs(b - b(1)) < eps)
    c = -inf;
    return;
end

r = corrcoef(a, b);
if numel(r) < 4 || ~isfinite(r(1, 2))
    c = -inf;
else
    c = r(1, 2);
end
end


function plot_step1_results_local(case_data, Sensor_Config, cfg, summary_table)
sensor_ids = Sensor_Config.Sensor_IDs;
blades_num = Sensor_Config.Blades_Num;
n_show = min(4 * blades_num, min_pulse_count_local(case_data, sensor_ids));

fig1 = figure('Name', '20250527 Step1 Fingerprint Alignment', ...
    'Color', 'w', 'Position', [80, 80, 1400, 900], 'NumberTitle', 'off');
tiledlayout(fig1, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

anchor_sid = cfg.reference_sensor_id;
ref_fp = Sensor_Config.Fingerprints(anchor_sid);
for i = 1:numel(sensor_ids)
    sid = sensor_ids(i);
    feats = case_data.channels(sid).peak_features(:);
    start_idx = Sensor_Config.Target_Indices(sid);
    pulse_idx = 1:n_show;

    nexttile;
    hold on;
    grid on;
    box on;
    plot(pulse_idx, feats(pulse_idx), '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.0, ...
        'DisplayName', 'Raw pulse feature');
    sel_idx = start_idx:(start_idx + blades_num - 1);
    plot(sel_idx, feats(sel_idx), 'o-', 'Color', [0 0.45 0.74], 'LineWidth', 1.5, ...
        'MarkerFaceColor', [0 0.45 0.74], 'DisplayName', 'Selected fingerprint');
    if sid ~= anchor_sid
        plot(sel_idx, ref_fp(:), 's--', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.2, ...
            'MarkerFaceColor', [0.85 0.33 0.10], 'DisplayName', 'Reference CH1');
    end
    xlabel('Pulse index');
    ylabel('Peak feature');
    title(sprintf('CH%d fingerprint alignment', sid));
    legend('Location', 'best');
end

fig2 = figure('Name', '20250527 Step1 Standard Relative Angles', ...
    'Color', 'w', 'Position', [120, 120, 1200, 520], 'NumberTitle', 'off');
tiledlayout(fig2, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
angle_mat = Sensor_Config.Standard_Relative_Angles(sensor_ids, :);
imagesc(angle_mat);
axis tight;
grid on;
colorbar;
set(gca, 'XTick', 1:blades_num, 'YTick', 1:numel(sensor_ids), ...
    'YTickLabel', compose('CH%d', sensor_ids));
xlabel('Blade ID');
ylabel('Sensor');
    title('Standard relative angles (OPR center, deg)');

nexttile;
hold on;
grid on;
box on;
colors = lines(numel(sensor_ids));
for i = 1:numel(sensor_ids)
    sid = sensor_ids(i);
    plot(1:blades_num, angle_mat(i, :), 'o-', 'LineWidth', 1.5, ...
        'Color', colors(i, :), 'MarkerFaceColor', colors(i, :), ...
        'DisplayName', sprintf('CH%d', sid));
end
xlabel('Blade ID');
ylabel('Angle relative to previous OPR center (deg)');
title('Per-sensor blade angle curves, OPR-center reference');
legend('Location', 'best');

figure('Name', '20250527 Step1 Low-Speed RPM', ...
    'Color', 'w', 'Position', [160, 160, 1200, 420], 'NumberTitle', 'off');
hold on;
grid on;
box on;
plot(case_data.omega_time_s, case_data.omega_rpm, '-', 'LineWidth', 1.2, ...
    'Color', [0.2 0.2 0.2]);
xlabel('Time (s)');
ylabel('RPM');
title(sprintf('Low-speed RPM trend: %s', cfg.low_speed_case), 'Interpreter', 'none');

figure('Name', '20250527 Step1 Summary Table', ...
    'Color', 'w', 'Position', [200, 200, 800, 220], 'NumberTitle', 'off');
uitable('Data', table2cell(summary_table), ...
    'ColumnName', summary_table.Properties.VariableNames, ...
    'Units', 'normalized', ...
    'Position', [0 0 1 1]);
end


function n = min_pulse_count_local(case_data, sensor_ids)
counts = zeros(numel(sensor_ids), 1);
for i = 1:numel(sensor_ids)
    counts(i) = numel(case_data.channels(sensor_ids(i)).peak_features);
end
n = min(counts);
end


function opr_reference = build_opr_center_reference_for_case_local(cfg, case_name)
case_dir = fullfile(cfg.dataset_root, case_name);
file_ids = list_case_file_ids_local(case_dir, cfg.opr_id);
center_times = [];
start_times = [];
last_end_time = [];
tail = struct('t', [], 'v', []);

for iFile = 1:numel(file_ids)
    file_id = file_ids(iFile);
    [t_local, v_local] = load_raw_channel_local(case_dir, cfg.opr_id, file_id, cfg.pinlv);
    if isempty(last_end_time)
        offset = 0;
    else
        offset = last_end_time + 1 / cfg.pinlv - t_local(1);
    end
    t_global = t_local + offset;
    last_end_time = t_global(end);

    if iFile == 1 && cfg.initial_trim_points > 0
        keep_idx = (cfg.initial_trim_points + 1):numel(t_global);
        t_global = t_global(keep_idx);
        v_local = v_local(keep_idx);
    end

    [seg_center, seg_start, tail] = extract_opr_centers_local( ...
        t_global, v_local, cfg.opr_threshold, cfg.gap_points, tail);
    center_times = [center_times; seg_center(:)]; %#ok<AGROW>
    start_times = [start_times; seg_start(:)]; %#ok<AGROW>
end

opr_reference = build_opr_reference_from_center_start_local( ...
    center_times, start_times, cfg.blades_num, cfg.r_tip_mm);
opr_reference.mode = 'multi_threshold_center';
opr_reference.standard_angle_reference = 'opr_pulse_center';
opr_reference.source_case = case_name;
opr_reference.source_note = 'Standard_Relative_Angles is converted from OPR threshold start-edge to multi-threshold OPR pulse center.';
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


function [t_sec, v] = load_raw_channel_local(case_dir, channel_id, file_id, pinlv)
filepath = fullfile(case_dir, sprintf('4-%d-%d.mat', channel_id, file_id));
loaded = load(filepath);
fn = fieldnames(loaded);
raw = loaded.(fn{1});
raw(raw(:, 1) == 0, :) = [];
t_sec = raw(:, 1) / pinlv;
v = raw(:, 2);
end


function [center_times, start_times, next_tail] = extract_opr_centers_local(t, v, threshold, gap_points, tail)
if nargin < 5 || isempty(tail)
    tail = struct('t', [], 'v', []);
end
if ~isempty(tail.t)
    t = [tail.t; t];
    v = [tail.v; v];
end

idx_above = find(v > threshold);
center_times = [];
start_times = [];
next_tail = struct('t', [], 'v', []);
if isempty(idx_above)
    return;
end

jumps = find(diff(idx_above) > gap_points);
seg_start = [idx_above(1); idx_above(jumps + 1)];
seg_end = [idx_above(jumps); idx_above(end)];

if numel(v) - idx_above(end) < gap_points
    carry_start = max(1, seg_start(end) - gap_points);
    next_tail.t = t(carry_start:end);
    next_tail.v = v(carry_start:end);
    seg_start(end) = [];
    seg_end(end) = [];
end

for iSeg = 1:numel(seg_start)
    a = seg_start(iSeg);
    b = seg_end(iSeg);
    start_times(end + 1, 1) = interpolate_threshold_time_local(t, v, a, threshold); %#ok<AGROW>
    center_times(end + 1, 1) = compute_multithreshold_center_local(t, v, a, b); %#ok<AGROW>
end
end


function t0 = interpolate_threshold_time_local(t, v, idx, threshold)
if idx <= 1
    t0 = t(idx);
    return;
end
t1 = t(idx - 1);
t2 = t(idx);
v1 = v(idx - 1);
v2 = v(idx);
if abs(v2 - v1) < eps
    t0 = t2;
else
    t0 = t1 + (threshold - v1) * (t2 - t1) / (v2 - v1);
end
end


function t_center = compute_multithreshold_center_local(t, v, a, b)
pad = max(8, round(0.25 * (b - a + 1)));
a0 = max(1, a - pad);
b0 = min(numel(v), b + pad);
t_seg = t(a0:b0);
v_seg = v(a0:b0);
try
    v_smooth = smooth(v_seg, 16);
catch
    v_smooth = smoothdata(v_seg, 'movmean', 16);
end
[peak_val, i_peak] = max(v_smooth);
base_val = median(v_smooth(v_smooth <= prctile(v_smooth, 30)), 'omitnan');
if ~isfinite(base_val)
    base_val = min(v_smooth);
end
levels = base_val + [0.30 0.40 0.50 0.60 0.70] * (peak_val - base_val);
centers = nan(numel(levels), 1);
for i = 1:numel(levels)
    level = levels(i);
    i_rise = find(v_smooth(1:i_peak) >= level, 1, 'first');
    i_fall_rel = find(v_smooth(i_peak:end) <= level, 1, 'first');
    if isempty(i_rise) || isempty(i_fall_rel)
        continue;
    end
    i_fall = i_peak + i_fall_rel - 1;
    t_rise = interpolate_level_crossing_local(t_seg, v_smooth, i_rise, level);
    t_fall = interpolate_level_crossing_local(t_seg, v_smooth, i_fall, level);
    centers(i) = 0.5 * (t_rise + t_fall);
end
t_center = median(centers, 'omitnan');
if ~isfinite(t_center)
    t_center = t_seg(i_peak);
end
end


function tx = interpolate_level_crossing_local(t, v, idx, level)
i1 = max(1, idx - 1);
i2 = idx;
if i1 == i2 || abs(v(i2) - v(i1)) < eps
    tx = t(idx);
else
    tx = t(i1) + (level - v(i1)) * (t(i2) - t(i1)) / (v(i2) - v(i1));
end
end


function opr_reference = build_opr_reference_from_center_start_local(center_times, start_times, pulses_per_rev, r_tip_mm)
opr_reference = struct('phase_shift_deg', 0, 'phase_shift_mm', 0, ...
    'median_center_minus_start_s', 0);
n = min(numel(center_times) - pulses_per_rev, numel(start_times));
if n < 1
    return;
end
dt_center = center_times(1:n) - start_times(1:n);
dt_rev = center_times((1:n) + pulses_per_rev) - center_times(1:n);
valid = isfinite(dt_center) & isfinite(dt_rev) & dt_rev > eps;
if ~any(valid)
    return;
end
shift_deg = 360 * dt_center(valid) ./ dt_rev(valid);
opr_reference.phase_shift_deg = median(shift_deg, 'omitnan');
opr_reference.phase_shift_mm = opr_reference.phase_shift_deg * (pi / 180) * r_tip_mm;
opr_reference.median_center_minus_start_s = median(dt_center(valid), 'omitnan');
end
