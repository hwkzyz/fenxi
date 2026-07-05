clc; clear; close all;

%STEP02_EXTRACT_DYNAMIC_BTT_20241106 Method-neutral dynamic BTT extraction.
%
% Improved version for the 20241106 data-foundation route.
%
% Main improvements over the previous draft:
%   1. Dynamic OPR timing uses the same explicit definition as Step01 by
%      default: threshold start edge for revolution segmentation.
%   2. jiluOPR is saved as [opr_arrival_or_center, opr_start, opr_end]
%      with metadata, so downstream scripts do not accidentally use the
%      threshold start edge as the OPR reference.
%   3. Probe arrival time uses polynomial-centroid timing by default, matching
%      the current 20241106 low-speed reference extractor.
%   4. Gap, OPR timing, probe-arrival method, and time-index interpretation
%      are controlled by BTTDataConfig_20241106 rather than hard-coded.
%   5. The final tail segment of the last file is flushed instead of being
%      silently dropped.
%   6. Dynamic blade labeling uses same-revolution fingerprint matching,
%      reports second-best/gap/shift-consistency quality indicators, and
%      does not blindly fill missing blade IDs unless explicitly enabled.
%   7. Full case_data and metadata are saved for diagnosis in addition to
%      legacy-compatible jiluOPR/omega/jilublade files.
%
% Usage:
%   Run this script directly after Step01.

cfg = BTTDataConfig_20241106();
cfg = apply_step02_defaults_local(cfg);

%% Step02 run settings
show_plots = true;
target_cases = cfg.dynamic_cases;

%% Step02 extraction and labeling settings
step02_gap_points = cfg.step02_gap_points;
step02_segment_pad_fraction = cfg.step02_segment_pad_fraction;
step02_segment_pad_min_points = cfg.step02_segment_pad_min_points;
fingerprint_min_corr = cfg.fingerprint_min_corr;
fingerprint_min_corr_gap = cfg.fingerprint_min_corr_gap;
fingerprint_min_shift_consistency = cfg.fingerprint_min_shift_consistency;
fingerprint_quality_policy = cfg.fingerprint_quality_policy;
step02_fill_missing_blade_ids = cfg.step02_fill_missing_blade_ids;
step02_revolution_vote_min_corr = cfg.step02_revolution_vote_min_corr;
step02_progress_rev_interval = cfg.step02_progress_rev_interval;

cfg.step02_gap_points = step02_gap_points;
cfg.step02_segment_pad_fraction = step02_segment_pad_fraction;
cfg.step02_segment_pad_min_points = step02_segment_pad_min_points;
cfg.fingerprint_min_corr = fingerprint_min_corr;
cfg.fingerprint_min_corr_gap = fingerprint_min_corr_gap;
cfg.fingerprint_min_shift_consistency = fingerprint_min_shift_consistency;
cfg.fingerprint_quality_policy = fingerprint_quality_policy;
cfg.step02_fill_missing_blade_ids = step02_fill_missing_blade_ids;
cfg.step02_revolution_vote_min_corr = step02_revolution_vote_min_corr;
cfg.step02_progress_rev_interval = step02_progress_rev_interval;

if ischar(target_cases) || isstring(target_cases)
    target_cases = cellstr(target_cases);
end

config_path = fullfile(cfg.step01_output_dir, 'Sensor_Config_20241106.mat');
if ~isfile(config_path)
    error('Step01 reference not found. Please run Step01_Build_LowSpeed_Reference_20241106.m first:\n  %s', config_path);
end

loaded = load(config_path, 'Sensor_Config');
Sensor_Config = loaded.Sensor_Config;
check_step01_step02_consistency_local(cfg, Sensor_Config);

all_rows = [];

for iCase = 1:numel(target_cases)
    case_name = char(target_cases{iCase});
    try
        case_data = process_case_step02_local(case_name, cfg, Sensor_Config);
    catch ME
        warning('Skipping case %s: %s', case_name, ME.message);
        continue;
    end

    output_dir = fullfile(cfg.step02_output_dir, case_name);
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    metadata = build_step02_metadata_local(cfg, Sensor_Config, case_name);
    save_opr_products_local(output_dir, case_data, metadata);

    case_rows = repmat(struct( ...
        'case_name', '', ...
        'sensor_id', NaN, ...
        'pulse_count', NaN, ...
        'labeled_pulse_count', NaN, ...
        'unlabeled_pulse_count', NaN, ...
        'first_blade_id', NaN, ...
        'match_revolution_id', NaN, ...
        'valid_revolution_count', NaN, ...
        'best_shift', NaN, ...
        'best_corr', NaN, ...
        'second_best_corr', NaN, ...
        'corr_gap', NaN, ...
        'shift_consistency', NaN, ...
        'best_shift_support', NaN, ...
        'second_best_shift_support', NaN, ...
        'shift_support_gap', NaN, ...
        'quality_status', ''), numel(cfg.sensor_ids), 1);

    for iSensor = 1:numel(cfg.sensor_ids)
        sid = cfg.sensor_ids(iSensor);
        ch = case_data.channels(sid);
        if isempty(ch.jilublade)
            continue;
        end

        jilublade = ch.jilublade; %#ok<NASGU>
        save(fullfile(output_dir, sprintf('jilublade_probe%d.mat', sid)), ...
            'jilublade', 'metadata');

        case_rows(iSensor).case_name = case_name;
        case_rows(iSensor).sensor_id = sid;
        case_rows(iSensor).pulse_count = size(ch.jilublade, 1);
        case_rows(iSensor).labeled_pulse_count = sum(isfinite(ch.jilublade(:, 4)));
        case_rows(iSensor).unlabeled_pulse_count = sum(~isfinite(ch.jilublade(:, 4)));
        first_labeled_idx = find(isfinite(ch.jilublade(:, 4)), 1, 'first');
        if ~isempty(first_labeled_idx)
            case_rows(iSensor).first_blade_id = ch.jilublade(first_labeled_idx, 4);
        end
        case_rows(iSensor).match_revolution_id = ch.match_revolution_id;
        case_rows(iSensor).valid_revolution_count = ch.valid_revolution_count;
        case_rows(iSensor).best_shift = ch.best_shift;
        case_rows(iSensor).best_corr = ch.best_corr;
        case_rows(iSensor).second_best_corr = ch.second_best_corr;
        case_rows(iSensor).corr_gap = ch.corr_gap;
        case_rows(iSensor).shift_consistency = ch.shift_consistency;
        case_rows(iSensor).best_shift_support = ch.best_shift_support;
        case_rows(iSensor).second_best_shift_support = ch.second_best_shift_support;
        case_rows(iSensor).shift_support_gap = ch.shift_support_gap;
        case_rows(iSensor).quality_status = ch.quality_status;
    end

    case_summary = struct2table(case_rows);
    writetable(case_summary, fullfile(output_dir, 'Step02_Summary_20241106.csv'));
    save(fullfile(output_dir, 'Step02_Summary_20241106.mat'), 'case_summary', 'metadata');

    cfg_snapshot = cfg; %#ok<NASGU>
    save(fullfile(output_dir, 'Step02_Dynamic_BTT_Extraction_20241106.mat'), ...
        'case_data', 'Sensor_Config', 'cfg_snapshot', 'metadata', '-v7.3');

    all_rows = [all_rows; case_rows]; %#ok<AGROW>

    fprintf('\nCase %s completed. Output dir:\n  %s\n', case_name, output_dir);

    if show_plots || cfg.save_figures
        plot_step02_case_local(case_data, cfg, case_name, case_summary, show_plots);
    end
end

summary_table = struct2table(all_rows);
if ~isempty(summary_table)
    disp(summary_table);
end


function cfg = apply_step02_defaults_local(cfg)
% Defaults are added here so this Step02 remains compatible with older
% versions of BTTDataConfig_20241106.

if ~isfield(cfg, 'opr_timing_method') || isempty(cfg.opr_timing_method)
    cfg.opr_timing_method = 'threshold_start_edge';
end
if ~isfield(cfg, 'opr_center_level_ratios') || isempty(cfg.opr_center_level_ratios)
    cfg.opr_center_level_ratios = [0.30 0.40 0.50 0.60 0.70];
end
if ~isfield(cfg, 'probe_arrival_method') || isempty(cfg.probe_arrival_method)
    cfg.probe_arrival_method = 'polynomial_peak_fit';
end
if ~isfield(cfg, 'probe_arrival_baseline_mode') || isempty(cfg.probe_arrival_baseline_mode)
    cfg.probe_arrival_baseline_mode = 'edge_median';
end
if ~isfield(cfg, 'probe_arrival_integral_signal') || isempty(cfg.probe_arrival_integral_signal)
    cfg.probe_arrival_integral_signal = 'baseline_corrected';
end
if ~isfield(cfg, 'time_index_mode') || isempty(cfg.time_index_mode)
    cfg.time_index_mode = 'auto';
end
if ~isfield(cfg, 'step02_gap_points') || isempty(cfg.step02_gap_points)
    cfg.step02_gap_points = cfg.gap_points;
end
if ~isfield(cfg, 'step02_segment_pad_fraction') || isempty(cfg.step02_segment_pad_fraction)
    cfg.step02_segment_pad_fraction = 0.30;
end
if ~isfield(cfg, 'step02_segment_pad_min_points') || isempty(cfg.step02_segment_pad_min_points)
    cfg.step02_segment_pad_min_points = 4;
end
if ~isfield(cfg, 'fingerprint_min_corr') || isempty(cfg.fingerprint_min_corr)
    cfg.fingerprint_min_corr = 0.85;
end
if ~isfield(cfg, 'fingerprint_min_corr_gap') || isempty(cfg.fingerprint_min_corr_gap)
    cfg.fingerprint_min_corr_gap = 0.05;
end
if ~isfield(cfg, 'fingerprint_min_shift_consistency') || isempty(cfg.fingerprint_min_shift_consistency)
    cfg.fingerprint_min_shift_consistency = 0.60;
end
if ~isfield(cfg, 'fingerprint_quality_policy') || isempty(cfg.fingerprint_quality_policy)
    cfg.fingerprint_quality_policy = 'warn';  % 'warn' or 'error'
end
if ~isfield(cfg, 'step02_fill_missing_blade_ids') || isempty(cfg.step02_fill_missing_blade_ids)
    cfg.step02_fill_missing_blade_ids = false;
end
if ~isfield(cfg, 'step02_revolution_vote_min_corr') || isempty(cfg.step02_revolution_vote_min_corr)
    cfg.step02_revolution_vote_min_corr = -inf;
end
if ~isfield(cfg, 'step02_progress_rev_interval') || isempty(cfg.step02_progress_rev_interval)
    cfg.step02_progress_rev_interval = 200;
end
end


function check_step01_step02_consistency_local(cfg, Sensor_Config)
if isfield(Sensor_Config, 'OPR_Timing_Method')
    if ~strcmpi(string(Sensor_Config.OPR_Timing_Method), string(cfg.opr_timing_method))
        warning(['Step01/Step02 OPR timing mismatch: Step01=%s, Step02=%s. ' ...
            'Coordinate bias may appear in downstream displacement/observation layers.'], ...
            string(Sensor_Config.OPR_Timing_Method), string(cfg.opr_timing_method));
    end
elseif isfield(Sensor_Config, 'Standard_Relative_Angles_Reference')
    ref = lower(string(Sensor_Config.Standard_Relative_Angles_Reference));
    if contains(ref, 'rising') && strcmpi(cfg.opr_timing_method, 'multi_threshold_center')
        warning(['Step01 Sensor_Config was built with %s, but Step02 uses %s. ' ...
            'Rebuild Step01 with OPR center before using OPRCenterStd downstream.'], ...
            Sensor_Config.Standard_Relative_Angles_Reference, cfg.opr_timing_method);
    end
end

if isfield(Sensor_Config, 'Probe_Arrival_Method')
    if ~strcmpi(string(Sensor_Config.Probe_Arrival_Method), string(cfg.probe_arrival_method))
        warning('Step01/Step02 probe arrival method mismatch: Step01=%s, Step02=%s.', ...
            string(Sensor_Config.Probe_Arrival_Method), string(cfg.probe_arrival_method));
    end
end
end


function metadata = build_step02_metadata_local(cfg, Sensor_Config, case_name)
metadata = struct();
metadata.dataset = cfg.dataset;
metadata.case_name = case_name;
metadata.created_by = mfilename;
metadata.created_on = datestr(now, 31);
metadata.opr_timing_method = cfg.opr_timing_method;
metadata.opr_center_level_ratios = cfg.opr_center_level_ratios;
metadata.probe_arrival_method = cfg.probe_arrival_method;
metadata.probe_arrival_baseline_mode = cfg.probe_arrival_baseline_mode;
metadata.probe_arrival_integral_signal = cfg.probe_arrival_integral_signal;
metadata.time_index_mode = cfg.time_index_mode;
metadata.step02_gap_points = cfg.step02_gap_points;
metadata.jiluOPR_columns = {'opr_arrival_or_center_s', 'opr_start_s', 'opr_end_s'};
metadata.jilublade_columns = {'pulse_start_s', 'pulse_end_s', 'probe_arrival_s', 'blade_id'};
metadata.Sensor_Config_CreatedBy = get_optional_field_local(Sensor_Config, 'CreatedBy', '');
metadata.Sensor_Config_CreatedOn = get_optional_field_local(Sensor_Config, 'CreatedOn', '');
end


function value = get_optional_field_local(s, name, default_value)
if isstruct(s) && isfield(s, name)
    value = s.(name);
else
    value = default_value;
end
end


function case_data = process_case_step02_local(case_name, cfg, Sensor_Config)
case_dir = fullfile(cfg.dataset_root, case_name);
if ~isfolder(case_dir)
    error('Case directory not found: %s', case_dir);
end

fprintf('>>> [Step02][%s] Loading case directory: %s\n', case_name, case_dir);

file_ids = list_case_file_ids_local(case_dir, cfg.opr_id);
if isempty(file_ids)
    error('No OPR files found in %s.', case_dir);
end

file_index_map = build_file_index_map_local(case_dir, file_ids, cfg);

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
    'second_best_corr', NaN, ...
    'corr_gap', NaN, ...
    'shift_consistency', NaN, ...
    'best_shift_support', NaN, ...
    'second_best_shift_support', NaN, ...
    'shift_support_gap', NaN, ...
    'quality_status', '', ...
    'raw_first_lap_t', [], ...
    'raw_first_lap_v', [], ...
    'fingerprint_peak_times', [], ...
    'fingerprint_peak_features', [], ...
    'aligned_fingerprint', [], ...
    'current_fingerprint', [], ...
    'ref_fingerprint', [], ...
    'selected_peak_idx', [], ...
    'missing_fill_count', NaN, ...
    'unlabeled_count', NaN);

opr_result = process_opr_local(case_dir, file_ids, file_index_map, cfg);
fprintf('>>> [Step02][%s] OPR pulses=%d, revolutions~=%d\n', ...
    case_name, numel(opr_result.arrival_times), floor(numel(opr_result.arrival_times) / cfg.opr_pulses_per_rev));
channels(cfg.opr_id).start_times = opr_result.start_times;
channels(cfg.opr_id).end_times = opr_result.end_times;
channels(cfg.opr_id).arrival_times = opr_result.arrival_times;

for iSensor = 1:numel(cfg.sensor_ids)
    sid = cfg.sensor_ids(iSensor);
    sensor_result = process_probe_local(case_dir, file_ids, file_index_map, sid, ...
        cfg, Sensor_Config, opr_result.arrival_times, case_name);
    channels(sid) = sensor_result;
end

[omega_time_s, omega_rad_s, omega_rpm] = compute_speed_local(opr_result.arrival_times, cfg.opr_pulses_per_rev);

case_data = struct();
case_data.case_name = case_name;
case_data.case_dir = case_dir;
case_data.file_ids = file_ids;
case_data.file_index_map = file_index_map;
case_data.channels = channels;
case_data.jiluOPR = [opr_result.arrival_times, opr_result.start_times, opr_result.end_times];
case_data.opr_times = opr_result.arrival_times;
case_data.opr_start_times = opr_result.start_times;
case_data.opr_end_times = opr_result.end_times;
case_data.omega_time_s = omega_time_s;
case_data.omega_rad_s = omega_rad_s;
case_data.omega_rpm = omega_rpm;
end


function file_index_map = build_file_index_map_local(case_dir, file_ids, cfg)
file_index_map = repmat(struct( ...
    'file_id', NaN, ...
    'raw_first_sample', NaN, ...
    'raw_last_sample', NaN, ...
    'offset_samples', NaN, ...
    'global_first_sample', NaN, ...
    'global_last_sample', NaN, ...
    'mode_used', ''), numel(file_ids), 1);

last_global_sample = NaN;
last_raw_sample = NaN;
for iFile = 1:numel(file_ids)
    file_id = file_ids(iFile);
    raw0 = load_raw_case_channel_no_offset_local(case_dir, cfg.opr_id, file_id);
    if isempty(raw0)
        continue;
    end
    raw_first = raw0(1, 1);
    raw_last = raw0(end, 1);

    mode_used = lower(strtrim(cfg.time_index_mode));
    if strcmpi(mode_used, 'auto')
        if iFile == 1 || ~isfinite(last_raw_sample)
            mode_used = 'absolute';
        elseif raw_first > last_raw_sample
            mode_used = 'absolute';
        else
            mode_used = 'local';
        end
    end

    if strcmpi(mode_used, 'absolute')
        offset = 0;
    elseif strcmpi(mode_used, 'local')
        if iFile == 1 || ~isfinite(last_global_sample)
            offset = 0;
        else
            offset = last_global_sample + 1 - raw_first;
        end
    else
        error('Unsupported cfg.time_index_mode: %s', cfg.time_index_mode);
    end

    global_first = raw_first + offset;
    global_last = raw_last + offset;

    file_index_map(iFile).file_id = file_id;
    file_index_map(iFile).raw_first_sample = raw_first;
    file_index_map(iFile).raw_last_sample = raw_last;
    file_index_map(iFile).offset_samples = offset;
    file_index_map(iFile).global_first_sample = global_first;
    file_index_map(iFile).global_last_sample = global_last;
    file_index_map(iFile).mode_used = mode_used;

    last_global_sample = global_last;
    last_raw_sample = raw_last;
end
end


function opr_result = process_opr_local(case_dir, file_ids, file_index_map, cfg)
tail = [];
start_times = [];
end_times = [];
arrival_times = [];

for iFile = 1:numel(file_ids)
    file_id = file_ids(iFile);
    raw = load_raw_case_channel_local(case_dir, cfg.opr_id, file_id, file_index_map);
    if isempty(raw)
        continue;
    end

    if iFile == 1 && cfg.initial_trim_points > 0 && size(raw, 1) > cfg.initial_trim_points
        raw(1:cfg.initial_trim_points, :) = [];
    end

    if ~isempty(tail)
        raw = [tail; raw]; %#ok<AGROW>
        raw = unique_raw_samples_local(raw);
    end

    is_last_file = (iFile == numel(file_ids));
    [segments, tail] = segment_signal_local(raw, cfg.opr_threshold, cfg.step02_gap_points, is_last_file);
    if isempty(segments.start_idx)
        continue;
    end

    for iSeg = 1:numel(segments.start_idx)
        a = segments.start_idx(iSeg);
        b = segments.end_idx(iSeg);
        [t_start, t_end, t_arrival] = compute_opr_arrival_local(raw, a, b, cfg);
        start_times(end + 1, 1) = t_start; %#ok<AGROW>
        end_times(end + 1, 1) = t_end; %#ok<AGROW>
        arrival_times(end + 1, 1) = t_arrival; %#ok<AGROW>
    end
end

[start_times, end_times, arrival_times] = clean_monotonic_triplets_local(start_times, end_times, arrival_times);

opr_result = struct();
opr_result.start_times = start_times(:);
opr_result.end_times = end_times(:);
opr_result.arrival_times = arrival_times(:);
end


function sensor_result = process_probe_local(case_dir, file_ids, file_index_map, sid, cfg, Sensor_Config, opr_times, case_name)
tail = [];
rows = zeros(1e4, 4);
peak_features = zeros(1e4, 1);
raw_peak_features = zeros(1e4, 1);
fitted_peak_features = zeros(1e4, 1);
row_count = 0;
raw_first_lap_data = [];

fprintf('>>> [Step02][%s] CH%d: reading waveform files...\n', case_name, sid);

for iFile = 1:numel(file_ids)
    file_id = file_ids(iFile);
    if ~isfile(fullfile(case_dir, sprintf('4-%d-%d.mat', sid, file_id)))
        continue;
    end

    raw = load_raw_case_channel_local(case_dir, sid, file_id, file_index_map);
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
        raw = unique_raw_samples_local(raw);
    end

    is_last_file = (iFile == numel(file_ids));
    [segments, tail] = segment_signal_local(raw, cfg.sensor_thresholds(sid), cfg.step02_gap_points, is_last_file);
    if isempty(segments.start_idx)
        continue;
    end

    for iSeg = 1:numel(segments.start_idx)
        a = segments.start_idx(iSeg);
        b = segments.end_idx(iSeg);
        [t_start, t_end, t_arrival, feat, feat_raw, feat_fitted] = ...
            compute_probe_arrival_local(raw, a, b, cfg, cfg.sensor_thresholds(sid));

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

fprintf('>>> [Step02][%s] CH%d: extracted pulses=%d\n', case_name, sid, row_count);

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
sensor_result.second_best_corr = match_info.second_best_corr;
sensor_result.corr_gap = match_info.corr_gap;
sensor_result.shift_consistency = match_info.shift_consistency;
sensor_result.best_shift_support = match_info.best_shift_support;
sensor_result.second_best_shift_support = match_info.second_best_shift_support;
sensor_result.shift_support_gap = match_info.shift_support_gap;
sensor_result.quality_status = match_info.quality_status;
sensor_result.raw_first_lap_t = fp_diag.raw_t;
sensor_result.raw_first_lap_v = fp_diag.raw_v;
sensor_result.fingerprint_peak_times = fp_diag.peak_times;
sensor_result.fingerprint_peak_features = fp_diag.peak_features;
sensor_result.aligned_fingerprint = fp_diag.aligned_fp;
sensor_result.current_fingerprint = fp_diag.current_fp;
sensor_result.ref_fingerprint = fp_diag.ref_fp;
sensor_result.selected_peak_idx = fp_diag.selected_peak_idx;
sensor_result.missing_fill_count = match_info.missing_fill_count;
sensor_result.unlabeled_count = sum(~isfinite(rows(:, 4)));

fprintf(['>>> [Step02][%s] CH%d: valid revolutions=%d, match revolution=%d, ' ...
    'shift=%d, corr=%.4f, gap=%.4f, shift consistency=%.3f, quality=%s\n'], ...
    case_name, sid, sensor_result.valid_revolution_count, sensor_result.match_revolution_id, ...
    sensor_result.best_shift, sensor_result.best_corr, sensor_result.corr_gap, ...
    sensor_result.shift_consistency, sensor_result.quality_status);
end


function raw = load_raw_case_channel_no_offset_local(case_dir, channel_id, file_id)
filepath = fullfile(case_dir, sprintf('4-%d-%d.mat', channel_id, file_id));
if ~isfile(filepath)
    raw = [];
    return;
end
loaded = load(filepath);
fn = fieldnames(loaded);
raw = loaded.(fn{1});
raw(raw(:, 1) == 0, :) = [];
raw = sortrows(raw, 1);
end


function raw = load_raw_case_channel_local(case_dir, channel_id, file_id, file_index_map)
raw = load_raw_case_channel_no_offset_local(case_dir, channel_id, file_id);
if isempty(raw)
    return;
end
idx = find([file_index_map.file_id] == file_id, 1, 'first');
if ~isempty(idx) && isfinite(file_index_map(idx).offset_samples)
    raw(:, 1) = raw(:, 1) + file_index_map(idx).offset_samples;
end
raw = sortrows(raw, 1);
end


function raw = unique_raw_samples_local(raw)
if isempty(raw)
    return;
end
raw = sortrows(raw, 1);
[~, ia] = unique(raw(:, 1), 'stable');
raw = raw(ia, :);
end


function [segments, tail] = segment_signal_local(raw, threshold, gap_points, is_last_file)
segments = struct('start_idx', [], 'end_idx', [], 'start_sample', [], 'end_sample', []);
tail = [];

if isempty(raw)
    return;
end

sig = raw(:, 2);
try
    sig_smooth = smooth(sig, 16);
catch
    sig_smooth = smoothdata(sig, 'movmean', 16);
end

idx_above = find(sig_smooth > threshold);
if isempty(idx_above)
    return;
end

sample_order = raw(idx_above, 1);
jumps = find(diff(sample_order) > gap_points);
seg_start_pos = [1; jumps(:) + 1];
seg_end_pos = [jumps(:); numel(idx_above)];

if isempty(seg_start_pos)
    return;
end

if is_last_file
    complete_pos = 1:numel(seg_start_pos);
else
    complete_pos = 1:(numel(seg_start_pos) - 1);
    last_start_idx = idx_above(seg_start_pos(end));
    tail_start = max(1, last_start_idx - floor(gap_points / 2));
    tail = raw(tail_start:end, :);
end

if isempty(complete_pos)
    return;
end

segments.start_idx = idx_above(seg_start_pos(complete_pos));
segments.end_idx = idx_above(seg_end_pos(complete_pos));
segments.start_sample = raw(segments.start_idx, 1);
segments.end_sample = raw(segments.end_idx, 1);
end


function [t_start, t_end, t_arrival] = compute_opr_arrival_local(raw, a, b, cfg)
fs = cfg.sample_rate_hz;
t_start = raw(a, 1) / fs;
t_end = raw(b, 1) / fs;
method = lower(strtrim(cfg.opr_timing_method));

switch method
    case {'rising_edge', 'threshold_rising_edge', 'start_edge', 'threshold_start_edge'}
        t_arrival = threshold_crossing_start_local(raw, a, cfg.opr_threshold, fs);
    case {'multi_threshold_center', 'opr_center', 'center'}
        t_arrival = opr_multithreshold_center_local(raw, a, b, cfg);
    otherwise
        warning('Unknown cfg.opr_timing_method=%s. Falling back to threshold_start_edge.', cfg.opr_timing_method);
        t_arrival = threshold_crossing_start_local(raw, a, cfg.opr_threshold, fs);
end

if ~isfinite(t_arrival)
    t_arrival = t_start;
end
end


function t_cross = threshold_crossing_start_local(raw, a, threshold, fs)
if a <= 1
    t_cross = raw(a, 1) / fs;
    return;
end
s1 = raw(a - 1, 1) / fs;
s2 = raw(a, 1) / fs;
v1 = raw(a - 1, 2);
v2 = raw(a, 2);
if abs(v2 - v1) < eps
    t_cross = s2;
else
    t_cross = s1 + (threshold - v1) * (s2 - s1) / (v2 - v1);
end
end


function t_center = opr_multithreshold_center_local(raw, a, b, cfg)
fs = cfg.sample_rate_hz;
pad = max(8, round(0.25 * (b - a + 1)));
a0 = max(1, a - pad);
b0 = min(size(raw, 1), b + pad);
t = raw(a0:b0, 1) / fs;
v = raw(a0:b0, 2);

try
    vs = smooth(v, 16);
catch
    vs = smoothdata(v, 'movmean', 16);
end

[peak_val, i_peak] = max(vs);
base_val = median(vs(vs <= prctile(vs, 30)), 'omitnan');
if ~isfinite(base_val)
    base_val = min(vs);
end
amp = peak_val - base_val;
if amp <= eps
    t_center = raw(a, 1) / fs;
    return;
end
levels = base_val + cfg.opr_center_level_ratios(:).' .* amp;
centers = nan(numel(levels), 1);
for k = 1:numel(levels)
    tr = crossing_time_local(t(1:i_peak), vs(1:i_peak), levels(k), 'rising');
    tf = crossing_time_local(t(i_peak:end), vs(i_peak:end), levels(k), 'falling');
    if isfinite(tr) && isfinite(tf)
        centers(k) = 0.5 * (tr + tf);
    end
end

t_center = median(centers, 'omitnan');
if ~isfinite(t_center)
    t_center = threshold_crossing_start_local(raw, a, cfg.opr_threshold, fs);
end
end


function tc = crossing_time_local(t, v, level, direction)
tc = NaN;
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
    tc = 0.5 * (t1 + t2);
else
    tc = t1 + (level - v1) * (t2 - t1) / (v2 - v1);
end
end


function [t_start, t_end, t_arrival, feat, feat_raw, feat_fitted] = compute_probe_arrival_local(raw, a, b, cfg, threshold)
fs = cfg.sample_rate_hz;
t_start = raw(a, 1) / fs;
t_end = raw(b, 1) / fs;

pad_pts = max(cfg.step02_segment_pad_min_points, floor((b - a + 1) * cfg.step02_segment_pad_fraction));
a_pad = max(1, a - pad_pts);
b_pad = min(size(raw, 1), b + pad_pts);
time_s = raw(a_pad:b_pad, 1) / fs;
data = raw(a_pad:b_pad, 2);

feat_raw = max(data);
feat_fitted = fit_peak_feature_local(time_s, data, threshold);
if isfinite(feat_fitted)
    feat = feat_fitted;
else
    feat = feat_raw;
end

method = lower(strtrim(cfg.probe_arrival_method));
switch method
    case {'half_area', 'baseline_half_area'}
        t_arrival = baseline_corrected_half_area_local(time_s, data, threshold, cfg);
    case {'threshold_centroid', 'centroid'}
        t_arrival = threshold_centroid_arrival_local(time_s, data, threshold);
    case {'peak', 'peak_time'}
        [~, i_max] = max(data);
        t_arrival = time_s(i_max);
    case {'polynomial_centroid', 'poly_centroid', 'polynomial_peak_fit'}
        [t_arrival, ~] = polynomial_centroid_arrival_local(time_s, data, threshold);
    otherwise
        warning('Unknown cfg.probe_arrival_method=%s. Falling back to polynomial_peak_fit.', cfg.probe_arrival_method);
        [t_arrival, ~] = polynomial_centroid_arrival_local(time_s, data, threshold);
end

if ~isfinite(t_arrival)
    [~, i_max] = max(data);
    t_arrival = time_s(i_max);
end
end


function t_arrival = baseline_corrected_half_area_local(time_s, data, threshold, cfg)
t_arrival = NaN;
if numel(time_s) < 2
    return;
end
baseline = estimate_baseline_local(data, cfg.probe_arrival_baseline_mode);
mode = lower(strtrim(cfg.probe_arrival_integral_signal));
switch mode
    case {'baseline_corrected', 'baseline'}
        pulse = data - baseline;
    case {'threshold_excess', 'threshold'}
        pulse = data - threshold;
    otherwise
        pulse = data - baseline;
end
pulse(~isfinite(pulse)) = 0;
pulse(pulse < 0) = 0;
if sum(pulse) <= eps
    pulse = data - min(data);
    pulse(pulse < 0) = 0;
end
if sum(pulse) <= eps
    return;
end
cum_area = cumtrapz(time_s, pulse);
total_area = cum_area(end);
if total_area <= eps
    return;
end
target = 0.5 * total_area;
idx = find(cum_area >= target, 1, 'first');
if isempty(idx)
    return;
elseif idx == 1
    t_arrival = time_s(1);
else
    a0 = cum_area(idx - 1);
    a1 = cum_area(idx);
    if abs(a1 - a0) < eps
        t_arrival = 0.5 * (time_s(idx - 1) + time_s(idx));
    else
        ratio = (target - a0) / (a1 - a0);
        t_arrival = time_s(idx - 1) + ratio * (time_s(idx) - time_s(idx - 1));
    end
end
end


function baseline = estimate_baseline_local(data, mode)
mode = lower(strtrim(mode));
n = numel(data);
if n == 0
    baseline = 0;
    return;
end
switch mode
    case {'edge_median', 'edge'}
        m = max(1, round(0.15 * n));
        edge_data = [data(1:m); data(max(1, n - m + 1):n)];
        baseline = median(edge_data, 'omitnan');
    case {'low_percentile', 'p20'}
        baseline = prctile(data, 20);
    case {'threshold_zero', 'zero'}
        baseline = 0;
    otherwise
        baseline = median(data(data <= prctile(data, 30)), 'omitnan');
end
if ~isfinite(baseline)
    baseline = min(data);
end
end


function t_arrival = threshold_centroid_arrival_local(time_s, data, threshold)
weights = data - threshold;
weights(~isfinite(weights)) = 0;
weights(weights < 0) = 0;
if sum(weights) <= eps
    t_arrival = NaN;
else
    t_arrival = sum(time_s(:) .* weights(:)) / sum(weights(:));
end
end


function [arrival_time, peak_feature] = polynomial_centroid_arrival_local(time_s, data, threshold)
arrival_time = NaN;
peak_feature = NaN;
if isempty(time_s) || isempty(data)
    return;
end
try
    smooth_seg = smooth_signal_for_peak_local(data);
catch
    smooth_seg = data;
end
valid_mask = smooth_seg > threshold;
if sum(valid_mask) < 4
    [peak_feature, idx_max] = max(data);
    arrival_time = time_s(idx_max);
    return;
end

t_fit = time_s(valid_mask);
v_fit = smooth_seg(valid_mask);
mu = mean(t_fit);
order = min(3, numel(unique(t_fit)) - 1);
if order < 1
    [peak_feature, idx_max] = max(data);
    arrival_time = time_s(idx_max);
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
    tf = dt + mu;
    weights = vf - threshold;
    weights(weights < 0) = 0;
    peak_feature = max(vf);
    if sum(weights) > 0
        arrival_time = sum(tf .* weights) / sum(weights);
    else
        [~, idx_max] = max(vf);
        arrival_time = tf(idx_max);
    end
catch
    [peak_feature, idx_max] = max(data);
    arrival_time = time_s(idx_max);
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
if numel(v) >= 21 && exist('sgolayfilt', 'file') == 2
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
    'best_corr', NaN, ...
    'second_best_corr', NaN, ...
    'corr_gap', NaN, ...
    'shift_consistency', NaN, ...
    'best_shift_support', NaN, ...
    'second_best_shift_support', NaN, ...
    'shift_support_gap', NaN, ...
    'quality_status', 'unmatched', ...
    'missing_fill_count', 0);

fp_diag = struct('raw_t', [], 'raw_v', [], 'peak_times', [], 'peak_features', [], ...
    'current_fp', [], 'aligned_fp', [], 'ref_fp', [], 'selected_peak_idx', []);

if isempty(arrival_times) || isempty(peak_features) || numel(opr_times) <= cfg.blades_num
    return;
end

if ~isempty(raw_first_lap_data)
    t_raw = raw_first_lap_data(:, 1) / cfg.sample_rate_hz;
    v_raw = raw_first_lap_data(:, 2);
    fp_diag.raw_t = t_raw;
    fp_diag.raw_v = v_raw;
    try
        v_smooth = smooth_signal_for_peak_local(v_raw);
        [peak_times, fitted_peaks] = extract_poly_centroid_v2_local( ...
            t_raw, v_raw, v_smooth, cfg.sensor_thresholds(sid), cfg.step02_gap_points, 5);
        fp_diag.peak_times = peak_times(:);
        fp_diag.peak_features = fitted_peaks(:);
    catch
        fp_diag.peak_times = [];
        fp_diag.peak_features = [];
    end
end

if isKey(Sensor_Config.Fingerprints, sid)
    ref_fp = Sensor_Config.Fingerprints(sid);
else
    ref_fp = Sensor_Config.Fingerprints(cfg.reference_sensor_id);
end
ref_fp = ref_fp(:).';
fp_diag.ref_fp = ref_fp;

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
rev_second_corr = nan(n_rev, 1);
rev_corr_gap = nan(n_rev, 1);
rev_pulse_indices = cell(n_rev, 1);
rev_is_valid = false(n_rev, 1);
rev_all_corr = nan(n_rev, cfg.blades_num);

fprintf('>>> [Step02][%s] CH%d: matching same-revolution 6-peak vectors...\n', case_name, sid);

for iRev = 1:n_rev
    rev_id = valid_rev_ids(iRev);
    idx = find(revolution_index == rev_id);
    idx = idx(:);
    rev_pulse_indices{iRev} = idx;
    if numel(idx) ~= cfg.blades_num
        continue;
    end

    rev_feats = peak_features(idx).';
    corr_by_shift = nan(1, cfg.blades_num);
    for shift = 0:(cfg.blades_num - 1)
        shifted_fp = circshift(rev_feats(:), -shift);
        corr_by_shift(shift + 1) = safe_corr_local(ref_fp(:), shifted_fp(:));
    end

    [sorted_corr, sorted_idx] = sort(corr_by_shift, 'descend');
    best_corr = sorted_corr(1);
    best_shift = sorted_idx(1) - 1;
    second_corr = NaN;
    if numel(sorted_corr) >= 2
        second_corr = sorted_corr(2);
    end

    rev_is_valid(iRev) = isfinite(best_corr);
    rev_best_shift(iRev) = best_shift;
    rev_best_corr(iRev) = best_corr;
    rev_second_corr(iRev) = second_corr;
    rev_corr_gap(iRev) = best_corr - second_corr;
    rev_all_corr(iRev, :) = corr_by_shift;

    if mod(iRev, cfg.step02_progress_rev_interval) == 0 || iRev == n_rev
        fprintf('>>> [Step02][%s] CH%d: scored revolutions %d/%d\n', case_name, sid, iRev, n_rev);
    end
end

match_info.valid_revolution_count = sum(rev_is_valid);
match_info.revolution_shift_path = [valid_rev_ids(:), rev_best_shift(:)];
match_info.revolution_score_path = [valid_rev_ids(:), rev_best_corr(:), rev_second_corr(:), rev_corr_gap(:)];

if ~any(rev_is_valid)
    return;
end

[best_shift, shift_consistency, shift_score] = choose_global_shift_local( ...
    rev_best_shift, rev_best_corr, rev_is_valid, cfg);
match_info.best_shift = best_shift;
match_info.shift_consistency = shift_consistency;

shift_mask = rev_is_valid & rev_best_shift == best_shift;
if any(shift_mask)
    [best_corr, local_pos] = max(rev_best_corr .* double(shift_mask) + (-inf) .* double(~shift_mask)); %#ok<NASGU>
    candidates = find(shift_mask);
    [best_corr, j] = max(rev_best_corr(candidates));
    best_rev_pos = candidates(j);
else
    [best_corr, best_rev_pos] = max(rev_best_corr);
end
match_info.match_revolution_id = valid_rev_ids(best_rev_pos);
match_info.best_corr = best_corr;

% Per-revolution ambiguity at the selected reference revolution.
% These values stay in the correlation scale [-1, 1] and are suitable for
% the fingerprint_min_corr_gap threshold.
match_info.second_best_corr = rev_second_corr(best_rev_pos);
match_info.corr_gap = rev_corr_gap(best_rev_pos);

% Global shift support is a separate vote-style diagnostic across all valid
% revolutions. It is normalized to [0, 1], so it is not confused with a
% correlation coefficient.
shift_support = shift_score;
shift_support(~isfinite(shift_support)) = 0;
total_support = sum(shift_support);
if total_support > eps
    shift_support = shift_support ./ total_support;
end
if isfinite(best_shift) && best_shift >= 0 && best_shift < cfg.blades_num
    match_info.best_shift_support = shift_support(best_shift + 1);
end
[sorted_shift_support, ~] = sort(shift_support, 'descend');
if numel(sorted_shift_support) >= 2
    match_info.second_best_shift_support = sorted_shift_support(2);
    match_info.shift_support_gap = sorted_shift_support(1) - sorted_shift_support(2);
end

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

if cfg.step02_fill_missing_blade_ids
    before_missing = sum(~isfinite(blade_ids));
    blade_ids = fill_missing_blade_ids_conservative_local(blade_ids, cfg.blades_num);
    after_missing = sum(~isfinite(blade_ids));
    match_info.missing_fill_count = before_missing - after_missing;
end

match_info.quality_status = evaluate_dynamic_match_quality_local(match_info, cfg);
if strcmpi(match_info.quality_status, 'weak') || strcmpi(match_info.quality_status, 'ambiguous')
    msg = sprintf('CH%d dynamic fingerprint match is %s: corr=%.4f, gap=%.4f, shiftConsistency=%.3f.', ...
        sid, match_info.quality_status, match_info.best_corr, match_info.corr_gap, match_info.shift_consistency);
    if strcmpi(cfg.fingerprint_quality_policy, 'error')
        error('%s', msg);
    else
        warning('%s', msg);
    end
end

sel_idx = rev_pulse_indices{best_rev_pos};
if ~isempty(sel_idx) && numel(sel_idx) == cfg.blades_num
    fp_diag.current_fp = peak_features(sel_idx(:)).';
    fp_diag.aligned_fp = circshift(fp_diag.current_fp(:), -best_shift).';
    fp_diag.selected_peak_idx = sel_idx(:).';
else
    n = min(cfg.blades_num, numel(peak_features));
    fp_diag.current_fp = peak_features(1:n).';
    fp_diag.aligned_fp = fp_diag.current_fp;
    fp_diag.selected_peak_idx = 1:n;
end
end


function [best_shift, shift_consistency, shift_score] = choose_global_shift_local(rev_best_shift, rev_best_corr, rev_is_valid, cfg)
shift_score = zeros(1, cfg.blades_num);
valid_idx = find(rev_is_valid & isfinite(rev_best_shift) & isfinite(rev_best_corr));
if isempty(valid_idx)
    best_shift = NaN;
    shift_consistency = NaN;
    return;
end

corr_use = rev_best_corr(valid_idx);
shift_use = rev_best_shift(valid_idx);
qualified = corr_use >= cfg.step02_revolution_vote_min_corr;
if ~any(qualified)
    qualified = true(size(corr_use));
end
corr_use = corr_use(qualified);
shift_use = shift_use(qualified);

weights = corr_use;
weights(~isfinite(weights)) = 0;
weights = max(weights, 0);
if all(weights <= eps)
    weights = ones(size(weights));
end
for k = 1:numel(shift_use)
    s = shift_use(k);
    if isfinite(s) && s >= 0 && s < cfg.blades_num
        shift_score(s + 1) = shift_score(s + 1) + weights(k);
    end
end
[~, best_idx] = max(shift_score);
best_shift = best_idx - 1;
shift_consistency = sum(shift_use == best_shift) / numel(shift_use);
end


function quality = evaluate_dynamic_match_quality_local(match_info, cfg)
quality = 'good';
if ~isfinite(match_info.best_corr)
    quality = 'unmatched';
    return;
end
if match_info.best_corr < cfg.fingerprint_min_corr
    quality = 'weak';
    return;
end
if isfinite(match_info.corr_gap) && match_info.corr_gap < cfg.fingerprint_min_corr_gap
    quality = 'ambiguous';
    return;
end
if isfinite(match_info.shift_consistency) && match_info.shift_consistency < cfg.fingerprint_min_shift_consistency
    quality = 'ambiguous';
    return;
end
end


function [rel_angles_deg, prev_opr_index, revolution_index] = compute_pulse_reference_geometry_local(arrival_times, opr_times, cfg)
rel_angles_deg = nan(size(arrival_times));
prev_opr_index = nan(size(arrival_times));
revolution_index = nan(size(arrival_times));

opr_step = cfg.opr_pulses_per_rev;
if numel(opr_times) <= opr_step
    return;
end

rev_period = opr_times((opr_step + 1):end) - opr_times(1:(end - opr_step));
valid_rev = isfinite(rev_period) & rev_period > eps;
spd_t = opr_times(1:(end - opr_step));
spd_t = spd_t(valid_rev);
spd_v = 360 ./ rev_period(valid_rev);
if isempty(spd_t)
    return;
end
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
    revolution_index(i) = floor((idx_prev_opr - 1) / opr_step) + 1;
end
end


function blade_ids = fill_missing_blade_ids_conservative_local(blade_ids, blades_num)
% Only fills isolated single missing labels between two labels that are
% consistent with the expected sequence. Long missing blocks are left as NaN.
for i = 2:(numel(blade_ids) - 1)
    if isfinite(blade_ids(i))
        continue;
    end
    if isfinite(blade_ids(i - 1)) && isfinite(blade_ids(i + 1))
        expected_now = mod(blade_ids(i - 1), blades_num) + 1;
        expected_next = mod(expected_now, blades_num) + 1;
        if blade_ids(i + 1) == expected_next
            blade_ids(i) = expected_now;
        end
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

sample_order = chase(:);
jumps = find(diff(sample_order) > gap_pts);
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
            warn_state_1, warn_state_2, warn_state_3, warn_state_4)); %#ok<NASGU>

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
opr_times = opr_times(:);
valid = isfinite(opr_times);
opr_times = opr_times(valid);
opr_times = unique(opr_times, 'stable');
if numel(opr_times) <= blades_num
    omega_time_s = [];
    omega_rad_s = [];
    omega_rpm = [];
    return;
end

rev_period = opr_times((blades_num + 1):end) - opr_times(1:(end - blades_num));
valid_period = isfinite(rev_period) & rev_period > eps;
omega_time_s = opr_times(1:(end - blades_num));
omega_time_s = omega_time_s(valid_period);
omega_rad_s = 2 * pi ./ rev_period(valid_period);
omega_rpm = 60 ./ rev_period(valid_period);
end


function save_opr_products_local(output_dir, case_data, metadata)
jiluOPR = case_data.jiluOPR; %#ok<NASGU>
omega_time_s = case_data.omega_time_s; %#ok<NASGU>
omega_rad_s = case_data.omega_rad_s; %#ok<NASGU>
omega_rpm = case_data.omega_rpm; %#ok<NASGU>
omega_zong = omega_rad_s; %#ok<NASGU>

save(fullfile(output_dir, 'jiluOPR.mat'), 'jiluOPR', 'metadata');
save(fullfile(output_dir, 'omega.mat'), 'omega_zong', 'omega_time_s', 'omega_rpm', 'omega_rad_s', 'metadata');
end


function [start_times, end_times, arrival_times] = clean_monotonic_triplets_local(start_times, end_times, arrival_times)
if isempty(arrival_times)
    return;
end
valid = isfinite(start_times) & isfinite(end_times) & isfinite(arrival_times);
start_times = start_times(valid);
end_times = end_times(valid);
arrival_times = arrival_times(valid);
[arrival_times, order] = sort(arrival_times);
start_times = start_times(order);
end_times = end_times(order);
keep = [true; diff(arrival_times) > eps];
start_times = start_times(keep);
end_times = end_times(keep);
arrival_times = arrival_times(keep);
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


function plot_step02_case_local(case_data, cfg, case_name, case_summary, show_plots)
sensor_ids = cfg.sensor_ids;
blades_num = cfg.blades_num;
colors = lines(blades_num);
case_fig_dir = fullfile(cfg.step02_figure_dir, case_name);

fig1 = figure('Name', sprintf('20241106 Step02 RPM - %s', case_name), ...
    'Color', 'w', 'Position', [80, 80, 1200, 380], 'NumberTitle', 'off', ...
    'Visible', visibility_state_local(show_plots));
hold on; grid on; box on;
plot(case_data.omega_time_s, case_data.omega_rpm, 'k-', 'LineWidth', 1.2);
xlabel('Time (s)'); ylabel('RPM');
title(sprintf('OPR-derived RPM trend: %s', case_name), 'Interpreter', 'none');
save_step02_figure_local(fig1, case_fig_dir, cfg, 'Step02_RPM_Trend');

fig2 = figure('Name', sprintf('20241106 Step02 Blade ID Sequence - %s', case_name), ...
    'Color', 'w', 'Position', [100, 100, 1450, 900], 'NumberTitle', 'off', ...
    'Visible', visibility_state_local(show_plots));
tiledlayout(fig2, numel(sensor_ids), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:numel(sensor_ids)
    sid = sensor_ids(i);
    ch = case_data.channels(sid);
    nexttile; hold on; grid on; box on;
    if isempty(ch.jilublade)
        title(sprintf('CH%d no extracted pulses', sid));
        continue;
    end
    n_show = min(size(ch.jilublade, 1), 120);
    for blade_id = 1:blades_num
        mask = ch.jilublade(1:n_show, 4) == blade_id;
        scatter(ch.jilublade(mask, 3), ch.jilublade(mask, 4), 24, colors(blade_id, :), ...
            'filled', 'DisplayName', sprintf('B%d', blade_id));
    end
    unlabeled = isnan(ch.jilublade(1:n_show, 4));
    if any(unlabeled)
        scatter(ch.jilublade(unlabeled, 3), zeros(sum(unlabeled), 1), 18, [0.4 0.4 0.4], ...
            'filled', 'DisplayName', 'unlabeled');
    end
    xlabel('Arrival time (s)'); ylabel('Blade ID');
    ylim([0, blades_num + 0.5]);
    title(sprintf('CH%d first %d arrivals | match rev=%d | valid rev=%d | shift=%d | corr=%.4f | quality=%s', ...
        sid, n_show, ch.match_revolution_id, ch.valid_revolution_count, ch.best_shift, ch.best_corr, ch.quality_status));
    legend('Location', 'eastoutside');
end
save_step02_figure_local(fig2, case_fig_dir, cfg, 'Step02_BladeID_Sequence');

fig3 = figure('Name', sprintf('20241106 Step02 Fingerprint - %s', case_name), ...
    'Color', 'w', 'Position', [130, 130, 1500, 950], 'NumberTitle', 'off', ...
    'Visible', visibility_state_local(show_plots));
tiledlayout(fig3, numel(sensor_ids), 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:numel(sensor_ids)
    sid = sensor_ids(i);
    ch = case_data.channels(sid);

    nexttile; hold on; grid on; box on;
    sel_idx = ch.selected_peak_idx;
    if ~isempty(sel_idx) && numel(sel_idx) >= blades_num
        sel_idx = sel_idx(1:blades_num);
        [outline_t, outline_v] = build_selected_pulse_outline_local( ...
            ch.start_times(sel_idx), ch.arrival_times(sel_idx), ch.end_times(sel_idx), ch.peak_features(sel_idx));
        plot(outline_t, outline_v, 'Color', [0.75 0.75 0.75], 'LineWidth', 1.0);
        first_ids = ch.jilublade(sel_idx, 4);
        for k = 1:numel(sel_idx)
            if isfinite(first_ids(k))
                plot(ch.arrival_times(sel_idx(k)), ch.peak_features(sel_idx(k)), 'o', ...
                    'MarkerSize', 7, 'LineWidth', 1.6, 'Color', colors(first_ids(k), :));
                text(ch.arrival_times(sel_idx(k)), ch.peak_features(sel_idx(k)) * 1.01, ...
                    sprintf('B%d', first_ids(k)), 'Color', colors(first_ids(k), :), ...
                    'FontWeight', 'bold', 'HorizontalAlignment', 'center');
            end
        end
        xlim([min(ch.start_times(sel_idx)), max(ch.end_times(sel_idx))]);
    end
    xlabel('Time (s)'); ylabel('Voltage / feature');
    title(sprintf('CH%d selected revolution waveform and labeling', sid));

    nexttile; hold on; grid on; box on;
    if ~isempty(ch.current_fingerprint)
        plot(1:numel(ch.current_fingerprint), normalize_vec_local(ch.current_fingerprint), 'o-', ...
            'LineWidth', 1.3, 'Color', [0.50 0.50 0.50], 'MarkerFaceColor', [0.50 0.50 0.50], ...
            'DisplayName', 'Raw current fingerprint');
    end
    if ~isempty(ch.aligned_fingerprint)
        plot(1:numel(ch.aligned_fingerprint), normalize_vec_local(ch.aligned_fingerprint), 'o-', ...
            'LineWidth', 1.6, 'Color', [0 0.45 0.74], 'MarkerFaceColor', [0 0.45 0.74], ...
            'DisplayName', 'Aligned current fingerprint');
    end
    if ~isempty(ch.ref_fingerprint)
        plot(1:numel(ch.ref_fingerprint), normalize_vec_local(ch.ref_fingerprint), 's--', ...
            'LineWidth', 1.3, 'Color', [0.85 0.33 0.10], 'MarkerFaceColor', [0.85 0.33 0.10], ...
            'DisplayName', 'Reference fingerprint');
    end
    xlabel('Blade position in one lap'); ylabel('Normalized feature');
    title(sprintf('CH%d fingerprint | shift=%d | corr=%.4f | gap=%.4f', ...
        sid, ch.best_shift, ch.best_corr, ch.corr_gap));
    legend('Location', 'best');
end
save_step02_figure_local(fig3, case_fig_dir, cfg, 'Step02_Fingerprint_Match');

fig4 = figure('Name', sprintf('20241106 Step02 Revolution Shift Path - %s', case_name), ...
    'Color', 'w', 'Position', [160, 160, 1450, 900], 'NumberTitle', 'off', ...
    'Visible', visibility_state_local(show_plots));
tiledlayout(fig4, numel(sensor_ids), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:numel(sensor_ids)
    sid = sensor_ids(i);
    ch = case_data.channels(sid);
    nexttile; hold on; grid on; box on;
    if ~isempty(ch.revolution_shift_path)
        plot(ch.revolution_shift_path(:, 1), ch.revolution_shift_path(:, 2), 'o-', ...
            'LineWidth', 1.2, 'MarkerSize', 4, 'Color', [0 0.45 0.74], ...
            'DisplayName', 'Per-revolution local best shift');
    end
    if ~isnan(ch.best_shift)
        yline(ch.best_shift, '--', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.4, ...
            'DisplayName', sprintf('Global best shift = %d', ch.best_shift));
    end
    xlabel('Revolution index'); ylabel('Local best shift');
    ylim([-0.5, blades_num - 0.5]);
    title(sprintf('CH%d per-revolution shift | consistency=%.3f | quality=%s', ...
        sid, ch.shift_consistency, ch.quality_status));
    legend('Location', 'eastoutside');
end
save_step02_figure_local(fig4, case_fig_dir, cfg, 'Step02_Revolution_Shift_Path');

fig5 = figure('Name', sprintf('20241106 Step02 Summary Table - %s', case_name), ...
    'Color', 'w', 'Position', [180, 180, 1100, 260], 'NumberTitle', 'off', ...
    'Visible', visibility_state_local(show_plots));
uitable('Data', sanitize_table_cells_local(table2cell(case_summary)), ...
    'ColumnName', case_summary.Properties.VariableNames, ...
    'Units', 'normalized', 'Position', [0 0 1 1]);
save_step02_figure_local(fig5, case_fig_dir, cfg, 'Step02_Summary_Table');

if ~show_plots
    close([fig1 fig2 fig3 fig4 fig5]);
end
end


function save_step02_figure_local(fig, figure_dir, cfg, tag)
if ~cfg.save_figures
    return;
end
if exist(figure_dir, 'dir') ~= 7
    mkdir(figure_dir);
end
png_file = fullfile(figure_dir, [tag, '.png']);
pdf_file = fullfile(figure_dir, [tag, '.pdf']);
try
    exportgraphics(fig, png_file, 'Resolution', 300);
    exportgraphics(fig, pdf_file, 'ContentType', 'vector');
catch
    saveas(fig, png_file);
    saveas(fig, pdf_file);
end
end


function state = visibility_state_local(show_plots)
if show_plots
    state = 'on';
else
    state = 'off';
end
end


function y = normalize_vec_local(x)
x = x(:);
if isempty(x)
    y = x;
    return;
end
x = x - mean(x, 'omitnan');
s = std(x, 'omitnan');
if s < eps || ~isfinite(s)
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
            cells_out{i} = char(join(value, ', '));
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

