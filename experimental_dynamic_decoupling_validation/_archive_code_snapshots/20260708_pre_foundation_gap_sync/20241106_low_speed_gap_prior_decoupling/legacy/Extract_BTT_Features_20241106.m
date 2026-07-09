function case_data = Extract_BTT_Features_20241106(case_name, cfg, opts)
%EXTRACT_BTT_FEATURES_20241106 Extract OPR edges and probe arrivals for one case.

if nargin < 2 || isempty(cfg)
    cfg = Get_20241106_BTT_Config();
end
if nargin < 3 || isempty(opts)
    opts = struct();
end

case_dir = fullfile(cfg.dataset_root, case_name);
if ~isfolder(case_dir)
    error('Case directory not found: %s', case_dir);
end

file_ids = list_case_file_ids(case_dir, cfg.opr_id);
if isempty(file_ids)
    error('No OPR files found in %s.', case_dir);
end

if isfield(opts, 'sensor_ids') && ~isempty(opts.sensor_ids)
    sensor_ids = unique(opts.sensor_ids(:).', 'stable');
else
    sensor_ids = cfg.sensor_ids;
end

max_channel_id = max([sensor_ids, cfg.opr_id]);
channels(max_channel_id) = struct( ...
    'arrival_times', [], ...
    'peak_features', [], ...
    'start_times', [], ...
    'end_times', []);

tail_data(max_channel_id) = struct('t', [], 'v', []);
file_ranges = repmat(struct('file_id', NaN, 'offset', NaN, 't_start', NaN, 't_end', NaN), numel(file_ids), 1);

opr_times = [];
last_end_time = [];

for iFile = 1:numel(file_ids)
    file_id = file_ids(iFile);
    [t_opr_local, v_opr_local] = load_raw_channel(case_dir, cfg.opr_id, file_id, cfg.pinlv);
    if isempty(last_end_time)
        offset = 0;
    else
        offset = last_end_time + 1 / cfg.pinlv - t_opr_local(1);
    end

    t_opr_global = t_opr_local + offset;
    last_end_time = t_opr_global(end);
    file_ranges(iFile).file_id = file_id;
    file_ranges(iFile).offset = offset;
    file_ranges(iFile).t_start = t_opr_global(1);
    file_ranges(iFile).t_end = t_opr_global(end);

    if iFile == 1 && cfg.initial_trim_points > 0
        keep_idx = (cfg.initial_trim_points + 1):numel(t_opr_global);
        t_opr_global = t_opr_global(keep_idx);
        v_opr_local = v_opr_local(keep_idx);
    end

    [opr_arrivals, tail_data(cfg.opr_id)] = extract_opr_edges_local( ...
        t_opr_global, v_opr_local, cfg.opr_threshold, cfg.gap_points, tail_data(cfg.opr_id));
    opr_times = [opr_times; opr_arrivals]; %#ok<AGROW>

    for sid = sensor_ids
        if ~channel_file_exists_local(case_dir, sid, file_id)
            continue;
        end

        [t_local, v_local] = load_raw_channel(case_dir, sid, file_id, cfg.pinlv);
        t_global = t_local + offset;

        if iFile == 1 && cfg.initial_trim_points > 0
            keep_idx = (cfg.initial_trim_points + 1):numel(t_global);
            t_global = t_global(keep_idx);
            v_local = v_local(keep_idx);
        end

        [pulse_info, next_tail] = extract_sensor_pulses_local( ...
            t_global, v_local, cfg.sensor_thresholds(sid), cfg.gap_points, tail_data(sid));
        tail_data(sid) = next_tail;

        if isempty(pulse_info.arrival_times)
            continue;
        end

        channels(sid).arrival_times = [channels(sid).arrival_times; pulse_info.arrival_times]; %#ok<AGROW>
        channels(sid).peak_features = [channels(sid).peak_features; pulse_info.peak_features]; %#ok<AGROW>
        channels(sid).start_times = [channels(sid).start_times; pulse_info.start_times]; %#ok<AGROW>
        channels(sid).end_times = [channels(sid).end_times; pulse_info.end_times]; %#ok<AGROW>
    end
end

[omega_time_s, omega_rad_s, omega_rpm] = compute_speed_local(opr_times, cfg.opr_pulses_per_rev);

case_data = struct();
case_data.case_name = case_name;
case_data.case_dir = case_dir;
case_data.file_ids = file_ids;
case_data.sensor_ids = sensor_ids(:).';
case_data.channels = channels;
case_data.opr_times = opr_times(:);
case_data.omega_time_s = omega_time_s(:);
case_data.omega_rad_s = omega_rad_s(:);
case_data.omega_rpm = omega_rpm(:);
case_data.file_ranges = file_ranges;
end


function file_ids = list_case_file_ids(case_dir, channel_id)
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


function [t_sec, v] = load_raw_channel(case_dir, channel_id, file_id, pinlv)
filepath = fullfile(case_dir, sprintf('4-%d-%d.mat', channel_id, file_id));
if ~isfile(filepath)
    error('Missing channel file: %s', filepath);
end

loaded = load(filepath);
fn = fieldnames(loaded);
raw = loaded.(fn{1});
raw(raw(:, 1) == 0, :) = [];

t_sec = raw(:, 1) / pinlv;
v = raw(:, 2);
end


function tf = channel_file_exists_local(case_dir, channel_id, file_id)
tf = isfile(fullfile(case_dir, sprintf('4-%d-%d.mat', channel_id, file_id)));
end


function [pulse_info, next_tail] = extract_sensor_pulses_local(t, v, threshold, gap_points, tail)
if nargin < 5 || isempty(tail)
    tail.t = [];
    tail.v = [];
end

if ~isempty(tail.t)
    t = [tail.t; t];
    v = [tail.v; v];
end

smooth_v = smooth_signal_local(v);
idx_above = find(smooth_v > threshold);

pulse_info = struct('arrival_times', [], 'peak_features', [], 'start_times', [], 'end_times', []);
next_tail = struct('t', [], 'v', []);
if isempty(idx_above)
    return;
end

[seg_start, seg_end] = split_segments_local(idx_above, gap_points);

if numel(v) - idx_above(end) < gap_points
    carry_start = max(1, seg_start(end) - gap_points);
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

    [arrival_times(i), peak_features(i)] = fit_arrival_time_local(t_seg, v_seg, smooth_seg, threshold);
    start_times(i) = t_seg(1);
    end_times(i) = t_seg(end);
end

pulse_info.arrival_times = arrival_times;
pulse_info.peak_features = peak_features;
pulse_info.start_times = start_times;
pulse_info.end_times = end_times;
end


function [opr_arrivals, next_tail] = extract_opr_edges_local(t, v, threshold, gap_points, tail)
if nargin < 5 || isempty(tail)
    tail.t = [];
    tail.v = [];
end

if ~isempty(tail.t)
    t = [tail.t; t];
    v = [tail.v; v];
end

idx_above = find(v > threshold);
opr_arrivals = [];
next_tail = struct('t', [], 'v', []);
if isempty(idx_above)
    return;
end

[seg_start, ~] = split_segments_local(idx_above, gap_points);

if numel(v) - idx_above(end) < gap_points
    carry_start = max(1, seg_start(end) - gap_points);
    next_tail.t = t(carry_start:end);
    next_tail.v = v(carry_start:end);
    seg_start(end) = [];
end

if isempty(seg_start)
    return;
end

opr_arrivals = nan(numel(seg_start), 1);
for i = 1:numel(seg_start)
    idx = seg_start(i);
    if idx <= 1
        opr_arrivals(i) = t(idx);
    else
        t1 = t(idx - 1);
        t2 = t(idx);
        v1 = v(idx - 1);
        v2 = v(idx);
        if abs(v2 - v1) < eps
            opr_arrivals(i) = t2;
        else
            opr_arrivals(i) = t1 + (threshold - v1) * (t2 - t1) / (v2 - v1);
        end
    end
end
end


function [seg_start, seg_end] = split_segments_local(idx_above, gap_points)
jumps = find(diff(idx_above) > gap_points);
seg_start = [idx_above(1); idx_above(jumps + 1)];
seg_end = [idx_above(jumps); idx_above(end)];
end


function smooth_v = smooth_signal_local(v)
if numel(v) >= 21
    smooth_v = sgolayfilt(v, 3, 21);
elseif numel(v) >= 5
    span = max(3, 2 * floor(numel(v) / 4) + 1);
    smooth_v = smoothdata(v, 'movmean', span);
else
    smooth_v = v;
end
end


function [arrival_time, peak_feature] = fit_arrival_time_local(t_seg, v_seg, smooth_seg, threshold)
valid_mask = smooth_seg > threshold;

if sum(valid_mask) < 4
    [peak_feature, idx_max] = max(v_seg);
    arrival_time = t_seg(idx_max);
    return;
end

t_fit = t_seg(valid_mask);
v_fit = v_seg(valid_mask);
mu = mean(t_fit);
unique_count = numel(unique(t_fit));
order = min(3, unique_count - 1);

if order < 1
    [peak_feature, idx_max] = max(v_seg);
    arrival_time = t_seg(idx_max);
    return;
end

try
    warning('off', 'MATLAB:polyfit:RepeatedPointsOrRescale');
    p = polyfit(t_fit - mu, v_fit, order);
    warning('on', 'MATLAB:polyfit:RepeatedPointsOrRescale');
    dt = linspace(min(t_fit) - mu, max(t_fit) - mu, 100);
    vf = polyval(p, dt);
    tf = dt + mu;
    weights = vf - threshold;
    weights(weights < 0) = 0;

    if sum(weights) > 0
        arrival_time = sum(tf .* weights) / sum(weights);
        peak_feature = max(vf);
    else
        [peak_feature, idx_max] = max(vf);
        arrival_time = tf(idx_max);
    end
catch
    warning('on', 'MATLAB:polyfit:RepeatedPointsOrRescale');
    [peak_feature, idx_max] = max(v_seg);
    arrival_time = t_seg(idx_max);
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


