function summary_table = Step2_Extract_JiluBlade_20250527(target_cases, show_plots)
%STEP2_EXTRACT_JILUBLADE_20250527 Original-like Step2 for 20250527 dynamic cases.
% This version follows the 20251222 Step2 idea more closely:
%   1. Read absolute sample-index waveforms directly from exported MAT files
%   2. Segment pulses with threshold + gap logic and keep the last segment as tail
%   3. Compute probe arrival time by half-area integration inside each pulse window
%   4. Determine blade order from fitted-peak fingerprints using Step1 Sensor_Config
%   5. Save jilublade_probe*.mat / jiluOPR.mat / omega.mat and visualize the result

cfg = Get_20250527_BTT_Config();
if nargin < 1 || isempty(target_cases)
    target_cases = cfg.dynamic_cases;
elseif ischar(target_cases) || isstring(target_cases)
    target_cases = cellstr(target_cases);
end
if nargin < 2
    show_plots = true;
end

config_path = fullfile(cfg.reference_output_dir, 'Sensor_Config_20250527.mat');
if ~isfile(config_path)
    Step1_Build_Sensor_Config_20250527();
end

loaded = load(config_path, 'Sensor_Config');
Sensor_Config = loaded.Sensor_Config;

all_rows = [];

for iCase = 1:numel(target_cases)
    case_name = target_cases{iCase};
    try
        case_data = process_case_step2_local(case_name, cfg, Sensor_Config);
    catch ME
        warning('Skipping case %s: %s', case_name, ME.message);
        continue;
    end

    output_dir = fullfile(cfg.output_root, case_name);
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    save_opr_products_local(output_dir, case_data);

    case_rows = repmat(struct( ...
        'case_name', '', ...
        'sensor_id', NaN, ...
        'pulse_count', NaN, ...
        'first_blade_id', NaN, ...
        'match_revolution_id', NaN, ...
        'valid_revolution_count', NaN, ...
        'best_shift', NaN, ...
        'best_corr', NaN), numel(cfg.sensor_ids), 1);

    for iSensor = 1:numel(cfg.sensor_ids)
        sid = cfg.sensor_ids(iSensor);
        ch = case_data.channels(sid);
        if isempty(ch.jilublade)
            continue;
        end

        jilublade = ch.jilublade;
        save(fullfile(output_dir, sprintf('jilublade_probe%d.mat', sid)), 'jilublade');

        case_rows(iSensor).case_name = char(case_name);
        case_rows(iSensor).sensor_id = sid;
        case_rows(iSensor).pulse_count = size(jilublade, 1);
        case_rows(iSensor).first_blade_id = jilublade(1, 4);
        case_rows(iSensor).match_revolution_id = ch.match_revolution_id;
        case_rows(iSensor).valid_revolution_count = ch.valid_revolution_count;
        case_rows(iSensor).best_shift = ch.best_shift;
        case_rows(iSensor).best_corr = ch.best_corr;
    end

    case_summary = struct2table(case_rows);
    writetable(case_summary, fullfile(output_dir, 'Step2_Summary_20250527.csv'));
    save(fullfile(output_dir, 'Step2_Summary_20250527.mat'), 'case_summary');
    all_rows = [all_rows; case_rows]; %#ok<AGROW>

    fprintf('\nCase %s completed. Output dir:\n  %s\n', case_name, output_dir);

    if show_plots
        plot_step2_case_local(case_data, cfg, case_name, case_summary);
    end
end

summary_table = struct2table(all_rows);
disp(summary_table);
end


function case_data = process_case_step2_local(case_name, cfg, Sensor_Config)
case_dir = fullfile(cfg.dataset_root, case_name);
if ~isfolder(case_dir)
    error('Case directory not found: %s', case_dir);
end

fprintf('>>> [Step2][%s] Loading case directory: %s\n', case_name, case_dir);

file_ids = list_case_file_ids_local(case_dir, cfg.opr_id);
if isempty(file_ids)
    error('No OPR files found in %s.', case_dir);
end

max_channel_id = max([cfg.sensor_ids, cfg.opr_id]);
channels(max_channel_id) = struct( ...
    'start_times', [], ...
    'end_times', [], ...
    'arrival_times', [], ...
    'peak_features', [], ...
    'raw_peak_features', [], ...
    'fitted_peak_features', [], ...
    'jilublade', [], ...
    'match_revolution_id', NaN, ...
    'relative_angles_deg', [], ...
    'prev_opr_index', [], ...
    'revolution_index', [], ...
    'revolution_shift_path', [], ...
    'revolution_score_path', [], ...
    'valid_revolution_count', NaN, ...
    'best_shift', NaN, ...
    'best_corr', NaN, ...
    'raw_first_lap_t', [], ...
    'raw_first_lap_v', [], ...
    'fingerprint_peak_times', [], ...
    'fingerprint_peak_features', [], ...
    'aligned_fingerprint', [], ...
    'current_fingerprint', [], ...
    'ref_fingerprint', [], ...
    'selected_peak_idx', []);

opr_result = process_opr_local(case_dir, file_ids, cfg);
fprintf('>>> [Step2][%s] OPR pulses=%d, revolutions≈%d\n', ...
    case_name, numel(opr_result.arrival_times), floor(numel(opr_result.arrival_times) / cfg.blades_num));
channels(cfg.opr_id).start_times = opr_result.start_times;
channels(cfg.opr_id).end_times = opr_result.end_times;
channels(cfg.opr_id).arrival_times = opr_result.arrival_times;

for iSensor = 1:numel(cfg.sensor_ids)
    sid = cfg.sensor_ids(iSensor);
    sensor_result = process_probe_local(case_dir, file_ids, sid, cfg, Sensor_Config, opr_result.arrival_times, case_name);
    channels(sid) = sensor_result;
end

[omega_time_s, omega_rad_s, omega_rpm] = compute_speed_local(opr_result.arrival_times, cfg.blades_num);

case_data = struct();
case_data.case_name = case_name;
case_data.case_dir = case_dir;
case_data.file_ids = file_ids;
case_data.channels = channels;
case_data.jiluOPR = [opr_result.arrival_times, opr_result.start_times, opr_result.end_times];
case_data.opr_times = opr_result.arrival_times;
case_data.omega_time_s = omega_time_s;
case_data.omega_rad_s = omega_rad_s;
case_data.omega_rpm = omega_rpm;
end


function opr_result = process_opr_local(case_dir, file_ids, cfg)
tail = [];
gap_points = 1e4;
start_times = [];
end_times = [];
arrival_times = [];

for iFile = 1:numel(file_ids)
    file_id = file_ids(iFile);
    raw = load_raw_case_channel_local(case_dir, cfg.opr_id, file_id);
    if isempty(raw)
        continue;
    end

    if iFile == 1 && cfg.initial_trim_points > 0 && size(raw, 1) > cfg.initial_trim_points
        raw(1:cfg.initial_trim_points, :) = [];
    end

    if ~isempty(tail)
        raw = [tail; raw]; %#ok<AGROW>
    end

    [segments, tail, gap_points] = segment_signal_original_like_local(raw, cfg.opr_threshold, gap_points);
    if isempty(segments)
        continue;
    end

    start_times = [start_times; segments.start_sample / cfg.pinlv]; %#ok<AGROW>
    end_times = [end_times; segments.end_sample / cfg.pinlv]; %#ok<AGROW>
    center_times = nan(numel(segments.start_idx), 1);
    for iSeg = 1:numel(segments.start_idx)
        center_times(iSeg) = compute_opr_multithreshold_center_local( ...
            raw, segments.start_idx(iSeg), segments.end_idx(iSeg), cfg);
    end
    arrival_times = [arrival_times; center_times]; %#ok<AGROW>
end

opr_result = struct();
opr_result.start_times = start_times(:);
opr_result.end_times = end_times(:);
opr_result.arrival_times = arrival_times(:);
end


function sensor_result = process_probe_local(case_dir, file_ids, sid, cfg, Sensor_Config, opr_times, case_name)
tail = [];
gap_points = 1e4;
rows = zeros(1e4, 4);
peak_features = zeros(1e4, 1);
raw_peak_features = zeros(1e4, 1);
fitted_peak_features = zeros(1e4, 1);
row_count = 0;
raw_first_lap_data = [];

fprintf('>>> [Step2][%s] CH%d: reading waveform files...\n', case_name, sid);

for iFile = 1:numel(file_ids)
    file_id = file_ids(iFile);
    if ~isfile(fullfile(case_dir, sprintf('4-%d-%d.mat', sid, file_id)))
        continue;
    end

    raw = load_raw_case_channel_local(case_dir, sid, file_id);
    if isempty(raw)
        continue;
    end

    if iFile == 1 && cfg.initial_trim_points > 0 && size(raw, 1) > cfg.initial_trim_points
        raw(1:cfg.initial_trim_points, :) = [];
        raw_first_lap_data = raw;
    elseif iFile == 1
        raw_first_lap_data = raw;
    end

    if ~isempty(tail)
        raw = [tail; raw]; %#ok<AGROW>
    end

    [segments, tail, gap_points] = segment_signal_original_like_local(raw, cfg.sensor_thresholds(sid), gap_points);
    if isempty(segments)
        continue;
    end

    for iSeg = 1:numel(segments.start_idx)
        a = segments.start_idx(iSeg);
        b = segments.end_idx(iSeg);
        [t_start, t_end, t_arrival, feat, feat_raw, feat_fitted] = ...
            compute_half_area_arrival_local(raw, a, b, cfg.pinlv, cfg.sensor_thresholds(sid));

        row_count = row_count + 1;
        if row_count > size(rows, 1)
            rows = [rows; zeros(size(rows, 1), 4)]; %#ok<AGROW>
            peak_features = [peak_features; zeros(size(peak_features, 1), 1)]; %#ok<AGROW>
            raw_peak_features = [raw_peak_features; zeros(size(raw_peak_features, 1), 1)]; %#ok<AGROW>
            fitted_peak_features = [fitted_peak_features; zeros(size(fitted_peak_features, 1), 1)]; %#ok<AGROW>
        end

        rows(row_count, 1:3) = [t_start, t_end, t_arrival];
        peak_features(row_count) = feat;
        raw_peak_features(row_count) = feat_raw;
        fitted_peak_features(row_count) = feat_fitted;
    end
end

rows = rows(1:row_count, :);
peak_features = peak_features(1:row_count);
raw_peak_features = raw_peak_features(1:row_count);
fitted_peak_features = fitted_peak_features(1:row_count);

fprintf('>>> [Step2][%s] CH%d: extracted pulses=%d\n', case_name, sid, row_count);

[blade_ids, match_info, fp_diag] = assign_blade_ids_by_revolution_local( ...
    raw_first_lap_data, rows(:, 3), peak_features, sid, cfg, Sensor_Config, opr_times, case_name);

if ~isempty(rows)
    rows(:, 4) = blade_ids;
end

sensor_result = struct();
sensor_result.start_times = rows(:, 1);
sensor_result.end_times = rows(:, 2);
sensor_result.arrival_times = rows(:, 3);
sensor_result.peak_features = peak_features(:);
sensor_result.raw_peak_features = raw_peak_features(:);
sensor_result.fitted_peak_features = fitted_peak_features(:);
sensor_result.jilublade = rows;
sensor_result.match_revolution_id = match_info.match_revolution_id;
sensor_result.relative_angles_deg = match_info.relative_angles_deg;
sensor_result.prev_opr_index = match_info.prev_opr_index;
sensor_result.revolution_index = match_info.revolution_index;
sensor_result.revolution_shift_path = match_info.revolution_shift_path;
sensor_result.revolution_score_path = match_info.revolution_score_path;
sensor_result.valid_revolution_count = match_info.valid_revolution_count;
sensor_result.best_shift = match_info.best_shift;
sensor_result.best_corr = match_info.best_corr;
sensor_result.raw_first_lap_t = fp_diag.raw_t;
sensor_result.raw_first_lap_v = fp_diag.raw_v;
sensor_result.fingerprint_peak_times = fp_diag.peak_times;
sensor_result.fingerprint_peak_features = fp_diag.peak_features;
sensor_result.aligned_fingerprint = fp_diag.aligned_fp;
sensor_result.current_fingerprint = fp_diag.current_fp;
sensor_result.ref_fingerprint = fp_diag.ref_fp;
sensor_result.selected_peak_idx = fp_diag.selected_peak_idx;

fprintf('>>> [Step2][%s] CH%d: valid revolutions=%d, matched revolution=%d, shift=%d, corr=%.4f\n', ...
    case_name, sid, sensor_result.valid_revolution_count, sensor_result.match_revolution_id, ...
    sensor_result.best_shift, sensor_result.best_corr);
end


function raw = load_raw_case_channel_local(case_dir, channel_id, file_id)
filepath = fullfile(case_dir, sprintf('4-%d-%d.mat', channel_id, file_id));
if ~isfile(filepath)
    raw = [];
    return;
end

loaded = load(filepath);
fn = fieldnames(loaded);
raw = loaded.(fn{1});
raw(raw(:, 1) == 0, :) = [];
end


function t_center = compute_opr_multithreshold_center_local(raw, a, b, cfg)
pad = max(8, round(0.25 * (b - a + 1)));
a0 = max(1, a - pad);
b0 = min(size(raw, 1), b + pad);
t = raw(a0:b0, 1) / cfg.pinlv;
v = raw(a0:b0, 2);
try
    vs = smooth(v, 16);
catch
    vs = smoothdata(v, 'movmean', 16);
end
[peakVal, iPeak] = max(vs);
baseVal = median(vs(vs <= prctile(vs, 30)), 'omitnan');
if ~isfinite(baseVal)
    baseVal = min(vs);
end
levels = baseVal + [0.30 0.40 0.50 0.60 0.70] .* max(peakVal - baseVal, eps);
centers = nan(numel(levels), 1);
for k = 1:numel(levels)
    tr = crossing_time_local(t(1:iPeak), vs(1:iPeak), levels(k), 'rising');
    tf = crossing_time_local(t(iPeak:end), vs(iPeak:end), levels(k), 'falling');
    if isfinite(tr) && isfinite(tf)
        centers(k) = 0.5 * (tr + tf);
    end
end
t_center = median(centers, 'omitnan');
if ~isfinite(t_center)
    t_center = 0.5 * (raw(a, 1) + raw(b, 1)) / cfg.pinlv;
end
end


function tc = crossing_time_local(t, v, level, direction)
tc = NaN;
t = t(:);
v = v(:);
if numel(t) < 2
    return;
end
if strcmpi(direction, 'rising')
    idx = find(v(1:end-1) < level & v(2:end) >= level, 1, 'first');
else
    idx = find(v(1:end-1) >= level & v(2:end) < level, 1, 'last');
end
if isempty(idx)
    return;
end
dv = v(idx + 1) - v(idx);
if abs(dv) < eps
    tc = t(idx);
else
    alpha = (level - v(idx)) / dv;
    tc = t(idx) + alpha * (t(idx + 1) - t(idx));
end
end


function [segments, tail, next_gap] = segment_signal_original_like_local(raw, threshold, gap_points)
segments = struct('start_idx', [], 'end_idx', [], 'start_sample', [], 'end_sample', []);
tail = [];
next_gap = gap_points;

if isempty(raw)
    return;
end

sig = raw(:, 2);
try
    sig_smooth = smooth(sig, 16);
catch
    sig_smooth = smoothdata(sig, 'movmean', 16);
end

chase = find(sig_smooth > threshold);
if isempty(chase)
    return;
end

sample_order = raw(chase, 1);
seg_starts = 1;
seg_ends = [];

for ii = 1:(numel(sample_order) - 1)
    c = sample_order(ii + 1) - sample_order(ii);
    if c > next_gap
        seg_ends(end + 1, 1) = ii; %#ok<AGROW>
        seg_starts(end + 1, 1) = ii + 1; %#ok<AGROW>
        next_gap = 0.6 * c;
    end
end
seg_ends(end + 1, 1) = numel(sample_order);

tail_point = chase(seg_starts(end)) - floor(next_gap / 2);
if tail_point > 0 && tail_point < size(raw, 1)
    tail = raw(tail_point:end, :);
end

if numel(seg_starts) < 2
    return;
end

comp_starts = seg_starts(1:end-1);
comp_ends = seg_ends(1:end-1);

segments.start_idx = chase(comp_starts);
segments.end_idx = chase(comp_ends);
segments.start_sample = raw(segments.start_idx, 1);
segments.end_sample = raw(segments.end_idx, 1);
end


function [t_start, t_end, t_arrival, feat, feat_raw, feat_fitted] = compute_half_area_arrival_local(raw, a, b, pinlv, threshold)
expand_pts = floor((b - a) * 0.3);
a = max(1, a - expand_pts);
b = min(size(raw, 1), b + expand_pts);

time_s = raw(a:b, 1) / pinlv;
data = raw(a:b, 2);

t_start = time_s(1);
t_end = time_s(end);
feat_raw = max(data);
feat_fitted = fit_peak_feature_local(time_s, data, threshold);
if isfinite(feat_fitted)
    feat = feat_fitted;
else
    feat = feat_raw;
end

if numel(time_s) < 2
    t_arrival = t_start;
    return;
end

dt = time_s(2) - time_s(1);
leiji = zeros(numel(data), 1);
sumj = 0;
for i = 1:numel(data)
    leiji(i) = sumj;
    sumj = sumj + data(i) * dt;
end

maxu = 0.5 * max(sumj);
z1 = find(leiji < maxu, 1, 'last');
if ~isempty(z1) && z1 < numel(time_s)
    t_arrival = 0.5 * (time_s(z1) + time_s(z1 + 1));
else
    t_arrival = t_start;
end
end


function peak_feature = fit_peak_feature_local(t_seg, v_seg, threshold)
peak_feature = NaN;
if isempty(t_seg) || isempty(v_seg)
    return;
end

smooth_seg = smooth_signal_for_peak_local(v_seg);
valid_mask = smooth_seg > threshold;
if sum(valid_mask) < 4
    peak_feature = max(smooth_seg);
    return;
end

t_fit = t_seg(valid_mask);
v_fit = smooth_seg(valid_mask);
mu = mean(t_fit);
unique_count = numel(unique(t_fit));
order = min(3, unique_count - 1);
if order < 1
    peak_feature = max(v_fit);
    return;
end

try
    warn_state_1 = warning('off', 'MATLAB:polyfit:RepeatedPointsOrRescale');
    warn_state_2 = warning('off', 'MATLAB:polyfit:PolyNotUnique');
    warn_state_3 = warning('off', 'MATLAB:singularMatrix');
    warn_state_4 = warning('off', 'MATLAB:nearlySingularMatrix');
    cleanup_obj = onCleanup(@() restore_polyfit_warnings_local( ...
        warn_state_1, warn_state_2, warn_state_3, warn_state_4)); %#ok<NASGU>

    p = polyfit(t_fit - mu, v_fit, order);
    dt = linspace(min(t_fit) - mu, max(t_fit) - mu, 100);
    vf = polyval(p, dt);
    peak_feature = max(vf);
catch
    peak_feature = max(smooth_seg);
end
end


function smooth_v = smooth_signal_for_peak_local(v)
if numel(v) >= 21
    smooth_v = sgolayfilt(v, 3, 21);
elseif numel(v) >= 5
    span = max(3, 2 * floor(numel(v) / 4) + 1);
    smooth_v = smoothdata(v, 'movmean', span);
else
    smooth_v = v;
end
end


function [blade_ids, match_info, fp_diag] = assign_blade_ids_by_revolution_local( ...
    raw_first_lap_data, arrival_times, peak_features, sid, cfg, Sensor_Config, opr_times, case_name)
total_pulses = numel(arrival_times);
blade_ids = nan(total_pulses, 1);

match_info = struct( ...
    'match_revolution_id', NaN, ...
    'relative_angles_deg', nan(total_pulses, 1), ...
    'prev_opr_index', nan(total_pulses, 1), ...
    'revolution_index', nan(total_pulses, 1), ...
    'revolution_shift_path', [], ...
    'revolution_score_path', [], ...
    'valid_revolution_count', 0, ...
    'best_shift', NaN, ...
    'best_corr', NaN);

fp_diag = struct('raw_t', [], 'raw_v', [], 'peak_times', [], 'peak_features', [], ...
    'current_fp', [], 'aligned_fp', [], 'ref_fp', [], 'selected_peak_idx', []);

if isempty(raw_first_lap_data) || isempty(arrival_times) || isempty(peak_features) || numel(opr_times) <= cfg.blades_num
    return;
end

t_raw = raw_first_lap_data(:, 1) / cfg.pinlv;
v_raw = raw_first_lap_data(:, 2);
fp_diag.raw_t = t_raw;
fp_diag.raw_v = v_raw;

try
    v_smooth = sgolayfilt(v_raw, 3, 21);
catch
    v_smooth = smooth(v_raw, 16);
end

[peak_times, fitted_peaks] = extract_poly_centroid_v2_local( ...
    t_raw, v_raw, v_smooth, cfg.sensor_thresholds(sid), cfg.gap_points, 5);
fp_diag.peak_times = peak_times(:);
fp_diag.peak_features = fitted_peaks(:);

if isKey(Sensor_Config.Fingerprints, sid)
    ref_fp = Sensor_Config.Fingerprints(sid);
else
    ref_fp = Sensor_Config.Fingerprints(cfg.reference_sensor_id);
end
fp_diag.ref_fp = ref_fp(:).';

[rel_angles_deg, prev_opr_index, revolution_index] = compute_pulse_reference_geometry_local(arrival_times, opr_times, cfg);
match_info.relative_angles_deg = rel_angles_deg;
match_info.prev_opr_index = prev_opr_index;
match_info.revolution_index = revolution_index;

valid_mask = ~isnan(revolution_index);
valid_rev_ids = unique(revolution_index(valid_mask));
if isempty(valid_rev_ids)
    return;
end

n_rev = numel(valid_rev_ids);
rev_best_shift = nan(n_rev, 1);
rev_best_corr = nan(n_rev, 1);
rev_pulse_indices = cell(n_rev, 1);
rev_is_valid = false(n_rev, 1);

fprintf('>>> [Step2][%s] CH%d: matching same-revolution 6-peak vectors...\n', case_name, sid);

for iRev = 1:n_rev
    rev_id = valid_rev_ids(iRev);
    idx = find(revolution_index == rev_id);
    rev_pulse_indices{iRev} = idx(:);
    if numel(idx) ~= cfg.blades_num
        continue;
    end

    rev_is_valid(iRev) = true;
    rev_feats = peak_features(idx(:)).';
    best_corr = -inf;
    best_shift = NaN;

    for shift = 0:(cfg.blades_num - 1)
        shifted_fp = circshift(rev_feats(:), -shift);
        feature_corr = safe_corr_local(ref_fp(:), shifted_fp(:));
        if feature_corr > best_corr
            best_corr = feature_corr;
            best_shift = shift;
        end
    end

    rev_best_shift(iRev) = best_shift;
    rev_best_corr(iRev) = best_corr;

    if mod(iRev, cfg.step2_progress_rev_interval) == 0 || iRev == n_rev
        fprintf('>>> [Step2][%s] CH%d: scored revolutions %d/%d\n', case_name, sid, iRev, n_rev);
    end
end

match_info.valid_revolution_count = sum(rev_is_valid);
match_info.revolution_shift_path = [valid_rev_ids(:), rev_best_shift(:)];
match_info.revolution_score_path = [valid_rev_ids(:), rev_best_corr(:)];

candidate_corr = rev_best_corr;
candidate_corr(~rev_is_valid) = -inf;
[best_corr, best_rev_pos] = max(candidate_corr);
if ~isfinite(best_corr)
    return;
end

best_shift = rev_best_shift(best_rev_pos);
best_rev_id = valid_rev_ids(best_rev_pos);
match_info.match_revolution_id = best_rev_id;
match_info.best_shift = best_shift;
match_info.best_corr = best_corr;

best_blade_seq = mod((0:(cfg.blades_num - 1)) - best_shift, cfg.blades_num) + 1;

for iRev = 1:n_rev
    idx = rev_pulse_indices{iRev};
    if isempty(idx)
        continue;
    end

    if numel(idx) == cfg.blades_num
        blade_ids(idx) = best_blade_seq(:);
    end
end

blade_ids = fill_missing_blade_ids_local(blade_ids, cfg.blades_num);

sel_rev_idx = best_rev_pos;
if isempty(sel_rev_idx) || ~isfinite(sel_rev_idx)
    fp_diag.current_fp = peak_features(1:min(cfg.blades_num, numel(peak_features))).';
    fp_diag.aligned_fp = fp_diag.current_fp;
    fp_diag.selected_peak_idx = 1:min(cfg.blades_num, numel(peak_features));
else
    sel_idx = rev_pulse_indices{sel_rev_idx};
    fp_diag.current_fp = peak_features(sel_idx(:)).';
    fp_diag.aligned_fp = circshift(fp_diag.current_fp(:), -best_shift).';
    fp_diag.selected_peak_idx = sel_idx(:).';
end
end


function [rel_angles_deg, prev_opr_index, revolution_index] = compute_pulse_reference_geometry_local(arrival_times, opr_times, cfg)
rel_angles_deg = nan(size(arrival_times));
prev_opr_index = nan(size(arrival_times));
revolution_index = nan(size(arrival_times));

if numel(opr_times) <= cfg.blades_num
    return;
end

spd_t = opr_times(1:(end - cfg.blades_num));
spd_v = 360 ./ (opr_times((cfg.blades_num + 1):end) - opr_times(1:(end - cfg.blades_num)));
F_omega_deg = griddedInterpolant(spd_t, spd_v, 'linear', 'nearest');

for i = 1:numel(arrival_times)
    t_meas = arrival_times(i);
    idx_prev_opr = find(opr_times < t_meas, 1, 'last');
    if isempty(idx_prev_opr)
        continue;
    end

    t_ref = opr_times(idx_prev_opr);
    t_grid = linspace(t_ref, t_meas, 10);
    rel_angles_deg(i) = trapz(t_grid, F_omega_deg(t_grid));
    prev_opr_index(i) = idx_prev_opr;
    revolution_index(i) = floor((idx_prev_opr - 1) / cfg.blades_num) + 1;
end
end


function blade_ids = fill_missing_blade_ids_local(blade_ids, blades_num)
valid_idx = find(isfinite(blade_ids));
if isempty(valid_idx)
    return;
end

first_valid = valid_idx(1);
for i = first_valid - 1:-1:1
    blade_ids(i) = mod(blade_ids(i + 1) - 2, blades_num) + 1;
end

for i = first_valid + 1:numel(blade_ids)
    if ~isfinite(blade_ids(i))
        blade_ids(i) = mod(blade_ids(i - 1), blades_num) + 1;
    end
end
end


function [cen_t, max_v] = extract_poly_centroid_v2_local(t_raw, v_raw, v_smooth, thres, gap_pts, poly_order)
cen_t = [];
max_v = [];
chase = find(v_smooth > thres);
if isempty(chase)
    return;
end

jumps = find(diff(chase) > gap_pts);
starts = [chase(1); chase(jumps + 1)];
ends = [chase(jumps); chase(end)];

for i = 1:length(starts)
    idx_s = starts(i);
    idx_e = ends(i);

    t_seg = t_raw(idx_s:idx_e);
    v_seg = v_raw(idx_s:idx_e);
    valid_mask = v_seg > thres;

    if sum(valid_mask) < poly_order + 1
        [mx, mi] = max(v_seg);
        if ~isempty(mi)
            cen_t(end + 1, 1) = t_seg(mi); %#ok<AGROW>
            max_v(end + 1, 1) = mx; %#ok<AGROW>
        end
        continue;
    end

    t_fit = t_seg(valid_mask);
    v_fit = v_seg(valid_mask);
    mu = mean(t_fit);

    try
        poly_deg = min([poly_order, 3, sum(valid_mask) - 2, numel(unique(t_fit)) - 1]);
        if poly_deg < 1
            [mx, mi] = max(v_seg);
            cen_t(end + 1, 1) = t_seg(mi); %#ok<AGROW>
            max_v(end + 1, 1) = mx; %#ok<AGROW>
            continue;
        end

        warn_state_1 = warning('off', 'MATLAB:polyfit:RepeatedPointsOrRescale');
        warn_state_2 = warning('off', 'MATLAB:polyfit:PolyNotUnique');
        warn_state_3 = warning('off', 'MATLAB:singularMatrix');
        warn_state_4 = warning('off', 'MATLAB:nearlySingularMatrix');
        cleanup_obj = onCleanup(@() restore_polyfit_warnings_local( ...
            warn_state_1, warn_state_2, warn_state_3, warn_state_4));

        p = polyfit(t_fit - mu, v_fit, poly_deg);
        dt = linspace(min(t_fit) - mu, max(t_fit) - mu, 100);
        vf = polyval(p, dt);
        tf = dt + mu;
        weights = vf - thres;
        weights(weights < 0) = 0;

        if sum(weights) > 1e-9
            cen_t(end + 1, 1) = sum(tf .* weights) / sum(weights); %#ok<AGROW>
            max_v(end + 1, 1) = max(vf); %#ok<AGROW>
        else
            [mx, mi] = max(vf);
            cen_t(end + 1, 1) = tf(mi); %#ok<AGROW>
            max_v(end + 1, 1) = mx; %#ok<AGROW>
        end
    catch
        [mx, mi] = max(v_seg);
        cen_t(end + 1, 1) = t_seg(mi); %#ok<AGROW>
        max_v(end + 1, 1) = mx; %#ok<AGROW>
    end
end
end


function file_ids = list_case_file_ids_local(case_dir, channel_id)
d = dir(fullfile(case_dir, sprintf('4-%d-*.mat', channel_id)));
file_ids = nan(numel(d), 1);
for i = 1:numel(d)
    tok = regexp(d(i).name, ['4-' num2str(channel_id) '-(\d+)\.mat'], 'tokens', 'once');
    if isempty(tok)
        continue;
    end
    file_ids(i) = str2double(tok{1});
end
file_ids = sort(unique(file_ids(~isnan(file_ids))));
end


function [omega_time_s, omega_rad_s, omega_rpm] = compute_speed_local(opr_times, blades_num)
if numel(opr_times) <= blades_num
    omega_time_s = [];
    omega_rad_s = [];
    omega_rpm = [];
    return;
end

rev_period = opr_times((blades_num + 1):end) - opr_times(1:(end - blades_num));
omega_time_s = opr_times(1:(end - blades_num));
omega_rad_s = 2 * pi ./ rev_period;
omega_rpm = 60 ./ rev_period;
end


function save_opr_products_local(output_dir, case_data)
jiluOPR = case_data.jiluOPR; %#ok<NASGU>
omega_time_s = case_data.omega_time_s; %#ok<NASGU>
omega_rad_s = case_data.omega_rad_s; %#ok<NASGU>
omega_rpm = case_data.omega_rpm; %#ok<NASGU>
omega_zong = omega_rad_s; %#ok<NASGU>

save(fullfile(output_dir, 'jiluOPR.mat'), 'jiluOPR');
save(fullfile(output_dir, 'omega.mat'), 'omega_zong', 'omega_time_s', 'omega_rpm', 'omega_rad_s');
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


function plot_step2_case_local(case_data, cfg, case_name, case_summary)
sensor_ids = cfg.sensor_ids;
blades_num = cfg.blades_num;
colors = lines(blades_num);

figure('Name', sprintf('20250527 Step2 RPM - %s', case_name), ...
    'Color', 'w', 'Position', [80, 80, 1200, 380], 'NumberTitle', 'off');
hold on;
grid on;
box on;
plot(case_data.omega_time_s, case_data.omega_rpm, 'k-', 'LineWidth', 1.2);
xlabel('Time (s)');
ylabel('RPM');
title(sprintf('OPR-derived RPM trend: %s', case_name), 'Interpreter', 'none');

fig2 = figure('Name', sprintf('20250527 Step2 Blade ID Sequence - %s', case_name), ...
    'Color', 'w', 'Position', [100, 100, 1450, 900], 'NumberTitle', 'off');
tiledlayout(fig2, numel(sensor_ids), 1, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:numel(sensor_ids)
    sid = sensor_ids(i);
    ch = case_data.channels(sid);
    if isempty(ch.jilublade)
        continue;
    end

    n_show = min(size(ch.jilublade, 1), 120);
    nexttile;
    hold on;
    grid on;
    box on;
    for blade_id = 1:blades_num
        mask = ch.jilublade(1:n_show, 4) == blade_id;
        scatter(ch.jilublade(mask, 3), ch.jilublade(mask, 4), 24, colors(blade_id, :), ...
            'filled', 'DisplayName', sprintf('B%d', blade_id));
    end
    xlabel('Arrival time (s)');
    ylabel('Blade ID');
    ylim([0.5, blades_num + 0.5]);
    title(sprintf('CH%d first %d labeled arrivals | match rev=%d | valid rev=%d | shift=%d | corr=%.4f', ...
        sid, n_show, ch.match_revolution_id, ch.valid_revolution_count, ch.best_shift, ch.best_corr));
    legend('Location', 'eastoutside');
end

fig3 = figure('Name', sprintf('20250527 Step2 First-Lap Fingerprint - %s', case_name), ...
    'Color', 'w', 'Position', [130, 130, 1500, 950], 'NumberTitle', 'off');
tiledlayout(fig3, numel(sensor_ids), 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:numel(sensor_ids)
    sid = sensor_ids(i);
    ch = case_data.channels(sid);

    nexttile;
    hold on;
    grid on;
    box on;
    sel_idx = ch.selected_peak_idx;
    if ~isempty(sel_idx)
        [outline_t, outline_v] = build_selected_pulse_outline_local( ...
            ch.start_times(sel_idx), ch.arrival_times(sel_idx), ch.end_times(sel_idx), ch.peak_features(sel_idx));
        plot(outline_t, outline_v, 'Color', [0.75 0.75 0.75], 'LineWidth', 1.0);
        first_ids = ch.jilublade(sel_idx, 4);
        for k = 1:blades_num
            plot(ch.arrival_times(sel_idx(k)), ch.peak_features(sel_idx(k)), 'o', ...
                'MarkerSize', 7, 'LineWidth', 1.6, 'Color', colors(first_ids(k), :));
            text(ch.arrival_times(sel_idx(k)), ch.peak_features(sel_idx(k)) * 1.01, ...
                sprintf('B%d', first_ids(k)), 'Color', colors(first_ids(k), :), ...
                'FontWeight', 'bold', 'HorizontalAlignment', 'center');
        end
        xlim([min(ch.start_times(sel_idx)), max(ch.end_times(sel_idx))]);
    end
    xlabel('Time (s)');
    ylabel('Voltage / feature');
    title(sprintf('CH%d selected revolution waveform and labeling', sid));

    nexttile;
    hold on;
    grid on;
    box on;
    if ~isempty(ch.current_fingerprint)
        plot(1:blades_num, normalize_vec_local(ch.current_fingerprint), 'o-', ...
            'LineWidth', 1.3, 'Color', [0.50 0.50 0.50], 'MarkerFaceColor', [0.50 0.50 0.50], ...
            'DisplayName', 'Raw current fingerprint');
    end
    if ~isempty(ch.aligned_fingerprint)
        plot(1:blades_num, normalize_vec_local(ch.aligned_fingerprint), 'o-', ...
            'LineWidth', 1.6, 'Color', [0 0.45 0.74], 'MarkerFaceColor', [0 0.45 0.74], ...
            'DisplayName', 'Aligned current fingerprint');
    end
    if ~isempty(ch.ref_fingerprint)
        plot(1:blades_num, normalize_vec_local(ch.ref_fingerprint), 's--', ...
            'LineWidth', 1.3, 'Color', [0.85 0.33 0.10], 'MarkerFaceColor', [0.85 0.33 0.10], ...
            'DisplayName', 'Reference fingerprint');
    end
    xlabel('Blade position in one lap');
    ylabel('Normalized feature');
    title(sprintf('CH%d fingerprint match | match rev=%d | shift=%d | corr=%.4f', ...
        sid, ch.match_revolution_id, ch.best_shift, ch.best_corr));
    legend('Location', 'best');
end

fig4 = figure('Name', sprintf('20250527 Step2 Revolution Shift Path - %s', case_name), ...
    'Color', 'w', 'Position', [160, 160, 1450, 900], 'NumberTitle', 'off');
tiledlayout(fig4, numel(sensor_ids), 1, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:numel(sensor_ids)
    sid = sensor_ids(i);
    ch = case_data.channels(sid);
    nexttile;
    hold on;
    grid on;
    box on;
    if ~isempty(ch.revolution_shift_path)
        plot(ch.revolution_shift_path(:, 1), ch.revolution_shift_path(:, 2), 'o-', ...
            'LineWidth', 1.2, 'MarkerSize', 4, 'Color', [0 0.45 0.74], ...
            'DisplayName', 'Per-revolution local best shift');
    end
    if ~isnan(ch.best_shift)
        yline(ch.best_shift, '--', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.4, ...
            'DisplayName', sprintf('Global best shift = %d', ch.best_shift));
    end
    xlabel('Revolution index');
    ylabel('Local best shift');
    ylim([-0.5, blades_num - 0.5]);
    title(sprintf('CH%d per-revolution local best shift', sid));
    legend('Location', 'eastoutside');
end

figure('Name', sprintf('20250527 Step2 Summary Table - %s', case_name), ...
    'Color', 'w', 'Position', [180, 180, 900, 220], 'NumberTitle', 'off');
uitable('Data', sanitize_table_cells_local(table2cell(case_summary)), ...
    'ColumnName', case_summary.Properties.VariableNames, ...
    'Units', 'normalized', ...
    'Position', [0 0 1 1]);
end


function y = normalize_vec_local(x)
x = x(:);
if isempty(x)
    y = x;
    return;
end
x = x - mean(x);
s = std(x);
if s < eps
    s = 1;
end
y = x / s;
end


function [outline_t, outline_v] = build_selected_pulse_outline_local(start_times, arrival_times, end_times, peak_features)
outline_t = [];
outline_v = [];
for i = 1:numel(arrival_times)
    outline_t = [outline_t; start_times(i); arrival_times(i); end_times(i); nan]; %#ok<AGROW>
    outline_v = [outline_v; 0; peak_features(i); 0; nan]; %#ok<AGROW>
end
end


function cells_out = sanitize_table_cells_local(cells_in)
cells_out = cells_in;
for i = 1:numel(cells_out)
    value = cells_out{i};
    if isstring(value)
        if isscalar(value)
            cells_out{i} = char(value);
        else
            cells_out{i} = char(join(value, ", "));
        end
    elseif ismissing(value)
        cells_out{i} = '';
    end
end
end


function restore_polyfit_warnings_local(varargin)
for i = 1:nargin
    warning(varargin{i}.state, varargin{i}.identifier);
end
end
