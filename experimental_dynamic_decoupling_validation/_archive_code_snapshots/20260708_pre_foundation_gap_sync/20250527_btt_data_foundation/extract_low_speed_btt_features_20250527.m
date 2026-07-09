function case_data = extract_low_speed_btt_features_20250527(cfg, sensor_ids)
%EXTRACT_LOW_SPEED_BTT_FEATURES_20250527 Extract low-speed OPR/probe features.
%
% Output is method-neutral: OPR times, per-sensor pulse times, pulse feature
% values, start/end times, and speed trace. It does not assign blade IDs.
%
% 2026-07 revision:
%   - OPR timing supports rising-edge or multi-threshold pulse center.
%   - Probe timing supports half-area, centroid, peak, and legacy polynomial
%     centroid definitions.
%   - MAT time-index mode supports local, absolute, and auto detection.

if nargin < 1 || isempty(cfg)
    cfg = BTTDataConfig_20250527();
end
if nargin < 2 || isempty(sensor_ids)
    sensor_ids = cfg.sensor_ids;
end
sensor_ids = unique(sensor_ids(:).', 'stable');

cfg = fill_missing_extraction_defaults_local(cfg);

case_name = cfg.low_speed_case;
case_dir = fullfile(cfg.dataset_root, case_name);
if ~isfolder(case_dir)
    error('Low-speed case directory not found: %s', case_dir);
end

file_ids = list_case_file_ids_local(case_dir, cfg.opr_id);
if isempty(file_ids)
    error('No OPR files found in %s.', case_dir);
end

max_channel_id = max([sensor_ids, cfg.opr_id]);
channels(max_channel_id) = struct( ...
    'arrival_times', [], ...
    'peak_features', [], ...
    'start_times', [], ...
    'end_times', [], ...
    'arrival_method', '');
raw_preview(max_channel_id) = struct('t', [], 'v', []);
tail_data(max_channel_id) = struct('t', [], 'v', []);

file_ranges = repmat(struct('file_id', NaN, ...
    'raw_t_start', NaN, 'raw_t_end', NaN, ...
    'offset', NaN, 't_start', NaN, 't_end', NaN, ...
    'time_index_mode_detected', ''), numel(file_ids), 1);

opr_times = [];
last_global_end_time = [];
last_raw_end_time = [];
time_mode_detected_all = cell(numel(file_ids), 1);

for iFile = 1:numel(file_ids)
    file_id = file_ids(iFile);
    [t_opr_local, v_opr_local] = load_raw_channel_local(case_dir, cfg.opr_id, file_id, cfg.sample_rate_hz);
    if isempty(t_opr_local)
        continue;
    end

    [offset, mode_detected] = compute_file_time_offset_local( ...
        t_opr_local, last_global_end_time, last_raw_end_time, cfg);
    t_opr_global = t_opr_local + offset;
    last_global_end_time = t_opr_global(end);
    last_raw_end_time = t_opr_local(end);
    time_mode_detected_all{iFile} = mode_detected;

    if iFile == 1
        raw_preview(cfg.opr_id).t = t_opr_global(:);
        raw_preview(cfg.opr_id).v = v_opr_local(:);
    end

    file_ranges(iFile).file_id = file_id;
    file_ranges(iFile).raw_t_start = t_opr_local(1);
    file_ranges(iFile).raw_t_end = t_opr_local(end);
    file_ranges(iFile).offset = offset;
    file_ranges(iFile).t_start = t_opr_global(1);
    file_ranges(iFile).t_end = t_opr_global(end);
    file_ranges(iFile).time_index_mode_detected = mode_detected;

    if iFile == 1 && cfg.initial_trim_points > 0 && numel(t_opr_global) > cfg.initial_trim_points
        keep_idx = (cfg.initial_trim_points + 1):numel(t_opr_global);
        t_opr_global = t_opr_global(keep_idx);
        v_opr_local = v_opr_local(keep_idx);
    end

    [opr_arrivals, tail_data(cfg.opr_id)] = extract_opr_pulses_local( ...
        t_opr_global, v_opr_local, cfg, tail_data(cfg.opr_id));
    opr_times = [opr_times; opr_arrivals]; %#ok<AGROW>

    for sid = sensor_ids
        if ~channel_file_exists_local(case_dir, sid, file_id)
            continue;
        end

        [t_local, v_local] = load_raw_channel_local(case_dir, sid, file_id, cfg.sample_rate_hz);
        t_global = t_local + offset;
        if iFile == 1
            raw_preview(sid).t = t_global(:);
            raw_preview(sid).v = v_local(:);
        end

        if iFile == 1 && cfg.initial_trim_points > 0 && numel(t_global) > cfg.initial_trim_points
            keep_idx = (cfg.initial_trim_points + 1):numel(t_global);
            t_global = t_global(keep_idx);
            v_local = v_local(keep_idx);
        end

        [pulse_info, tail_data(sid)] = extract_sensor_pulses_local( ...
            t_global, v_local, cfg.sensor_thresholds(sid), cfg, tail_data(sid));

        channels(sid).arrival_times = [channels(sid).arrival_times; pulse_info.arrival_times]; %#ok<AGROW>
        channels(sid).peak_features = [channels(sid).peak_features; pulse_info.peak_features]; %#ok<AGROW>
        channels(sid).start_times = [channels(sid).start_times; pulse_info.start_times]; %#ok<AGROW>
        channels(sid).end_times = [channels(sid).end_times; pulse_info.end_times]; %#ok<AGROW>
        channels(sid).arrival_method = cfg.probe_arrival_method;
    end
end

opr_times = sort(opr_times(:));
[omega_time_s, omega_rad_s, omega_rpm] = compute_speed_local(opr_times, cfg.blades_num);

case_data = struct();
case_data.case_name = case_name;
case_data.case_dir = case_dir;
case_data.file_ids = file_ids;
case_data.sensor_ids = sensor_ids;
case_data.channels = channels;
case_data.opr_times = opr_times(:);
case_data.opr_timing_method = cfg.opr_timing_method;
case_data.probe_arrival_method = cfg.probe_arrival_method;
case_data.time_index_mode = cfg.time_index_mode;
case_data.time_index_mode_detected = unique_nonempty_strings_local(time_mode_detected_all);
case_data.omega_time_s = omega_time_s(:);
case_data.omega_rad_s = omega_rad_s(:);
case_data.omega_rpm = omega_rpm(:);
case_data.file_ranges = file_ranges;
case_data.raw_preview = raw_preview;
end

function cfg = fill_missing_extraction_defaults_local(cfg)
if ~isfield(cfg, 'opr_timing_method') || isempty(cfg.opr_timing_method)
    cfg.opr_timing_method = 'multi_threshold_center';
end
if ~isfield(cfg, 'opr_center_level_ratios') || isempty(cfg.opr_center_level_ratios)
    cfg.opr_center_level_ratios = [0.30 0.40 0.50 0.60 0.70];
end
if ~isfield(cfg, 'opr_center_smooth_span_points') || isempty(cfg.opr_center_smooth_span_points)
    cfg.opr_center_smooth_span_points = 16;
end
if ~isfield(cfg, 'opr_center_pad_fraction') || isempty(cfg.opr_center_pad_fraction)
    cfg.opr_center_pad_fraction = 0.25;
end
if ~isfield(cfg, 'opr_center_min_pad_points') || isempty(cfg.opr_center_min_pad_points)
    cfg.opr_center_min_pad_points = 8;
end
if ~isfield(cfg, 'probe_arrival_method') || isempty(cfg.probe_arrival_method)
    cfg.probe_arrival_method = 'half_area';
end
if ~isfield(cfg, 'probe_smooth_span_points') || isempty(cfg.probe_smooth_span_points)
    cfg.probe_smooth_span_points = 21;
end
if ~isfield(cfg, 'time_index_mode') || isempty(cfg.time_index_mode)
    cfg.time_index_mode = 'auto';
end
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

function [t_sec, v] = load_raw_channel_local(case_dir, channel_id, file_id, sample_rate_hz)
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
t_sec = raw(:, 1) / sample_rate_hz;
v = raw(:, 2);
end

function tf = channel_file_exists_local(case_dir, channel_id, file_id)
tf = isfile(fullfile(case_dir, sprintf('4-%d-%d.mat', channel_id, file_id)));
end

function [offset, mode_detected] = compute_file_time_offset_local(t_local, last_global_end_time, last_raw_end_time, cfg)
mode = lower(strtrim(cfg.time_index_mode));
if isempty(last_global_end_time)
    offset = 0;
    if strcmp(mode, 'auto')
        mode_detected = 'first_file';
    else
        mode_detected = mode;
    end
    return;
end

sample_dt = 1 / cfg.sample_rate_hz;
switch mode
    case 'absolute'
        offset = 0;
        mode_detected = 'absolute';
    case 'local'
        offset = last_global_end_time + sample_dt - t_local(1);
        mode_detected = 'local_concatenated';
    case 'auto'
        if t_local(1) > last_raw_end_time
            offset = 0;
            mode_detected = 'absolute_detected';
        else
            offset = last_global_end_time + sample_dt - t_local(1);
            mode_detected = 'local_detected_concatenated';
        end
    otherwise
        error('cfg.time_index_mode must be auto, local, or absolute. Current: %s', cfg.time_index_mode);
end
end

function [opr_arrivals, next_tail] = extract_opr_pulses_local(t, v, cfg, tail)
if nargin < 4 || isempty(tail)
    tail = struct('t', [], 'v', []);
end
if ~isempty(tail.t)
    t = [tail.t; t];
    v = [tail.v; v];
end

idx_above = find(v > cfg.opr_threshold);
opr_arrivals = [];
next_tail = struct('t', [], 'v', []);
if isempty(idx_above)
    return;
end
[seg_start, seg_end] = split_segments_local(idx_above, cfg.gap_points);
if numel(v) - idx_above(end) < cfg.gap_points
    carry_start = max(1, seg_start(end) - cfg.gap_points);
    next_tail.t = t(carry_start:end);
    next_tail.v = v(carry_start:end);
    seg_start(end) = [];
    seg_end(end) = [];
end
if isempty(seg_start)
    return;
end

opr_arrivals = nan(numel(seg_start), 1);
for i = 1:numel(seg_start)
    switch lower(strtrim(cfg.opr_timing_method))
        case 'rising_edge'
            opr_arrivals(i) = compute_rising_edge_time_local( ...
                t, v, seg_start(i), cfg.opr_threshold);
        case {'multi_threshold_center', 'center', 'opr_center'}
            opr_arrivals(i) = compute_opr_multithreshold_center_local( ...
                t, v, seg_start(i), seg_end(i), cfg);
        otherwise
            error('Unsupported cfg.opr_timing_method: %s', cfg.opr_timing_method);
    end
end
end

function t_cross = compute_rising_edge_time_local(t, v, idx, threshold)
if idx <= 1
    t_cross = t(idx);
    return;
end
t1 = t(idx - 1);
t2 = t(idx);
v1 = v(idx - 1);
v2 = v(idx);
if abs(v2 - v1) < eps
    t_cross = t2;
else
    t_cross = t1 + (threshold - v1) * (t2 - t1) / (v2 - v1);
end
end

function tCenter = compute_opr_multithreshold_center_local(t, v, a, b, cfg)
pad = max(cfg.opr_center_min_pad_points, round(cfg.opr_center_pad_fraction * (b - a + 1)));
a0 = max(1, a - pad);
b0 = min(numel(v), b + pad);
t_seg = t(a0:b0);
v_seg = v(a0:b0);
v_s = smooth_signal_local(v_seg, cfg.opr_center_smooth_span_points);

[peak_val, i_peak] = max(v_s);
base_candidates = v_s(v_s <= prctile(v_s, 30));
if isempty(base_candidates)
    base_val = min(v_s);
else
    base_val = median(base_candidates, 'omitnan');
end
height = max(peak_val - base_val, eps);
levels = base_val + cfg.opr_center_level_ratios(:) * height;
centers = nan(numel(levels), 1);
for k = 1:numel(levels)
    tr = crossing_time_level_local(t_seg(1:i_peak), v_s(1:i_peak), levels(k), 'rising');
    tf = crossing_time_level_local(t_seg(i_peak:end), v_s(i_peak:end), levels(k), 'falling');
    if isfinite(tr) && isfinite(tf)
        centers(k) = 0.5 * (tr + tf);
    end
end
if any(isfinite(centers))
    tCenter = median(centers, 'omitnan');
else
    tCenter = t_seg(i_peak);
end
end

function t_cross = crossing_time_level_local(t, v, level, direction)
t_cross = NaN;
if numel(t) < 2
    return;
end
switch lower(direction)
    case 'rising'
        idx = find(v(1:end-1) <= level & v(2:end) >= level, 1, 'last');
    case 'falling'
        idx = find(v(1:end-1) >= level & v(2:end) <= level, 1, 'first');
    otherwise
        idx = [];
end
if isempty(idx)
    return;
end
v1 = v(idx);
v2 = v(idx + 1);
t1 = t(idx);
t2 = t(idx + 1);
if abs(v2 - v1) < eps
    t_cross = 0.5 * (t1 + t2);
else
    t_cross = t1 + (level - v1) * (t2 - t1) / (v2 - v1);
end
end

function [pulse_info, next_tail] = extract_sensor_pulses_local(t, v, threshold, cfg, tail)
if nargin < 5 || isempty(tail)
    tail = struct('t', [], 'v', []);
end
if ~isempty(tail.t)
    t = [tail.t; t];
    v = [tail.v; v];
end

smooth_v = smooth_signal_local(v, cfg.probe_smooth_span_points);
idx_above = find(smooth_v > threshold);
pulse_info = struct('arrival_times', [], 'peak_features', [], 'start_times', [], 'end_times', []);
next_tail = struct('t', [], 'v', []);
if isempty(idx_above)
    return;
end

[seg_start, seg_end] = split_segments_local(idx_above, cfg.gap_points);
if numel(v) - idx_above(end) < cfg.gap_points
    carry_start = max(1, seg_start(end) - cfg.gap_points);
    next_tail.t = t(carry_start:end);
    next_tail.v = v(carry_start:end);
    seg_start(end) = [];
    seg_end(end) = [];
end
if isempty(seg_start)
    return;
end

nSeg = numel(seg_start);
arrival_times = nan(nSeg, 1);
peak_features = nan(nSeg, 1);
start_times = nan(nSeg, 1);
end_times = nan(nSeg, 1);
for i = 1:nSeg
    idx_range = seg_start(i):seg_end(i);
    t_seg = t(idx_range);
    v_seg = v(idx_range);
    smooth_seg = smooth_v(idx_range);
    [arrival_times(i), peak_features(i)] = compute_probe_arrival_local( ...
        t_seg, v_seg, smooth_seg, threshold, cfg.probe_arrival_method);
    start_times(i) = t_seg(1);
    end_times(i) = t_seg(end);
end

pulse_info.arrival_times = arrival_times;
pulse_info.peak_features = peak_features;
pulse_info.start_times = start_times;
pulse_info.end_times = end_times;
end

function [arrival_time, peak_feature] = compute_probe_arrival_local(t_seg, v_seg, smooth_seg, threshold, method)
[peak_feature, idx_max] = max(v_seg);
method = lower(strtrim(method));
y = smooth_seg - threshold;
y(y < 0) = 0;
valid = y > 0;
if sum(valid) < 2
    arrival_time = t_seg(idx_max);
    return;
end

switch method
    case 'half_area'
        t_valid = t_seg(valid);
        y_valid = y(valid);
        area_cum = cumtrapz(t_valid, y_valid);
        total_area = area_cum(end);
        if total_area > 0
            arrival_time = interp1(area_cum, t_valid, 0.5 * total_area, 'linear', 'extrap');
        else
            arrival_time = t_seg(idx_max);
        end
    case 'threshold_centroid'
        arrival_time = sum(t_seg(:) .* y(:)) / sum(y(:));
    case 'peak'
        arrival_time = t_seg(idx_max);
    case 'polynomial_centroid'
        arrival_time = polynomial_centroid_arrival_local(t_seg, v_seg, smooth_seg, threshold);
    otherwise
        error('Unsupported cfg.probe_arrival_method: %s', method);
end
end

function arrival_time = polynomial_centroid_arrival_local(t_seg, v_seg, smooth_seg, threshold)
valid_mask = smooth_seg > threshold;
if sum(valid_mask) < 4
    [~, idx_max] = max(v_seg);
    arrival_time = t_seg(idx_max);
    return;
end

t_fit = t_seg(valid_mask);
v_fit = v_seg(valid_mask);
mu = mean(t_fit);
order = min(3, numel(unique(t_fit)) - 1);
if order < 1
    [~, idx_max] = max(v_seg);
    arrival_time = t_seg(idx_max);
    return;
end

warn_state = warning('off', 'MATLAB:polyfit:RepeatedPointsOrRescale');
try
    p = polyfit(t_fit - mu, v_fit, order);
    dt = linspace(min(t_fit) - mu, max(t_fit) - mu, 100);
    vf = polyval(p, dt);
    tf = dt + mu;
    weights = vf - threshold;
    weights(weights < 0) = 0;
    if sum(weights) > 0
        arrival_time = sum(tf .* weights) / sum(weights);
    else
        [~, idx_max] = max(vf);
        arrival_time = tf(idx_max);
    end
catch
    [~, idx_max] = max(v_seg);
    arrival_time = t_seg(idx_max);
end
warning(warn_state);
end

function [seg_start, seg_end] = split_segments_local(idx_above, gap_points)
jumps = find(diff(idx_above) > gap_points);
seg_start = [idx_above(1); idx_above(jumps + 1)];
seg_end = [idx_above(jumps); idx_above(end)];
end

function smooth_v = smooth_signal_local(v, span)
if nargin < 2 || isempty(span)
    span = 21;
end
span = max(3, round(span));
if mod(span, 2) == 0
    span = span + 1;
end
span = min(span, numel(v));
if mod(span, 2) == 0
    span = span - 1;
end
if span >= 5 && exist('sgolayfilt', 'file') == 2
    smooth_v = sgolayfilt(v, 3, span);
elseif span >= 3 && exist('smoothdata', 'file') == 2
    smooth_v = smoothdata(v, 'movmean', span);
else
    smooth_v = v;
end
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

function out = unique_nonempty_strings_local(c)
if isempty(c)
    out = {};
    return;
end
mask = ~cellfun(@isempty, c);
out = unique(c(mask), 'stable');
end
