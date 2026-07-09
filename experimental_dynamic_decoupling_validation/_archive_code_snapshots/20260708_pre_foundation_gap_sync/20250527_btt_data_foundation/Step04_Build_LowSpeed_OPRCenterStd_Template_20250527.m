clc; clear; close all;

%STEP04_BUILD_LOWSPEED_OPRCENTERSTD_TEMPLATE_20250527 Build low-speed waveform templates.
%
% This Step04 is the bridge between the method-neutral data foundation
% (Step01-Step03) and later direct waveform inversion.
%
% Purpose:
%   1. Load the Step01 low-speed reference: blade fingerprints, internal
%      blade numbering, OPR reference, and standard relative angles.
%   2. Re-read the low-speed raw BTT waveforms from cfg.low_speed_case.
%   3. Map each low-speed waveform sample to the OPRCenterStd spatial
%      coordinate x = (theta_actual - theta_standard) * R.
%   4. Aggregate all low-speed pulse point clouds into non-parametric
%      sensor-blade templates T_{s,b}(x), including template gradients,
%      coverage weights, and trusted spatial domains.
%   5. Save a reusable Template structure for later dynamic waveform map
%      construction and direct-template synchronous identification.
%
% This step intentionally does NOT identify vibration parameters. It only
% builds the low-speed forward model V ~= T_{s,b}(x).
%
% Usage:
%   Run this script directly after Step01.

%% Step04 run settings
show_plots = true;

cfg = BTTDataConfig_20250527();
cfg = apply_step04_defaults_local(cfg);

%% Step04 template-building settings
step04_standard_angle_source = cfg.step04_standard_angle_source;
step04_pulse_pad_fraction = cfg.step04_pulse_pad_fraction;
step04_pulse_min_pad_points = cfg.step04_pulse_min_pad_points;
step04_max_points_per_pulse = cfg.step04_max_points_per_pulse;
step04_max_pulses_per_sensor = cfg.step04_max_pulses_per_sensor;
step04_max_laps_per_sensor = cfg.step04_max_laps_per_sensor;
step04_x_collect_abs_limit_mm = cfg.step04_x_collect_abs_limit_mm;

step04_grid_dx_mm = cfg.step04_grid_dx_mm;
step04_min_bin_count = cfg.step04_min_bin_count;
step04_smooth_span_bins = cfg.step04_smooth_span_bins;

step04_gradient_min_ratio = cfg.step04_gradient_min_ratio;
step04_amplitude_min_ratio = cfg.step04_amplitude_min_ratio;
step04_min_half_width_mm = cfg.step04_min_half_width_mm;
step04_max_half_width_mm = cfg.step04_max_half_width_mm;
step04_domain_center_mode = cfg.step04_domain_center_mode;

step04_min_template_points = cfg.step04_min_template_points;
step04_min_template_pulses = cfg.step04_min_template_pulses;
step04_min_domain_width_mm = cfg.step04_min_domain_width_mm;
step04_plot_max_templates_per_figure = cfg.step04_plot_max_templates_per_figure;

cfg.step04_standard_angle_source = step04_standard_angle_source;
cfg.step04_pulse_pad_fraction = step04_pulse_pad_fraction;
cfg.step04_pulse_min_pad_points = step04_pulse_min_pad_points;
cfg.step04_max_points_per_pulse = step04_max_points_per_pulse;
cfg.step04_max_pulses_per_sensor = step04_max_pulses_per_sensor;
cfg.step04_max_laps_per_sensor = step04_max_laps_per_sensor;
cfg.step04_x_collect_abs_limit_mm = step04_x_collect_abs_limit_mm;
cfg.step04_grid_dx_mm = step04_grid_dx_mm;
cfg.step04_min_bin_count = step04_min_bin_count;
cfg.step04_smooth_span_bins = step04_smooth_span_bins;
cfg.step04_gradient_min_ratio = step04_gradient_min_ratio;
cfg.step04_amplitude_min_ratio = step04_amplitude_min_ratio;
cfg.step04_min_half_width_mm = step04_min_half_width_mm;
cfg.step04_max_half_width_mm = step04_max_half_width_mm;
cfg.step04_domain_center_mode = step04_domain_center_mode;
cfg.step04_min_template_points = step04_min_template_points;
cfg.step04_min_template_pulses = step04_min_template_pulses;
cfg.step04_min_domain_width_mm = step04_min_domain_width_mm;
cfg.step04_plot_max_templates_per_figure = step04_plot_max_templates_per_figure;

sensor_config_path = fullfile(cfg.step01_output_dir, 'Sensor_Config_20250527.mat');
low_speed_feature_path = fullfile(cfg.step01_output_dir, 'LowSpeed_Features_20250527.mat');
if ~isfile(sensor_config_path) || ~isfile(low_speed_feature_path)
    error('Missing Step01 low-speed reference. Please run Step01_Build_LowSpeed_Reference_20250527.m first.');
end

loaded_cfg = load(sensor_config_path, 'Sensor_Config');
Sensor_Config = loaded_cfg.Sensor_Config;
loaded_features = load(low_speed_feature_path, 'case_data');
case_data = loaded_features.case_data;

check_step04_inputs_local(cfg, Sensor_Config, case_data);
std_angles = resolve_standard_angles_local(Sensor_Config, cfg.step04_standard_angle_source);

fprintf('\n=== Step04: low-speed OPRCenterStd non-parametric template ===\n');
fprintf('Dataset: %s\n', cfg.dataset);
fprintf('Low-speed case: %s\n', cfg.low_speed_case);
fprintf('Sensors: %s\n', mat2str(cfg.sensor_ids));
fprintf('Blades: 1..%d\n', cfg.blades_num);
fprintf('OPR timing reference: %s\n', get_string_field_local(Sensor_Config, 'OPR_Timing_Method', 'unknown'));
fprintf('Standard-angle source: %s\n', cfg.step04_standard_angle_source);
fprintf('OPR events/rev: %d (%s)\n', ...
    cfg.step04_opr_events_per_revolution, ...
    get_string_field_local(cfg, 'opr_events_per_revolution_policy', 'unspecified'));

[F_omega_deg_s, speed_diagnostic] = build_low_speed_speed_interpolant_local( ...
    case_data.opr_times(:), cfg.step04_opr_events_per_revolution);
speed_diagnostic.opr_events_per_revolution = cfg.step04_opr_events_per_revolution;
speed_diagnostic.opr_events_per_revolution_source = 'cfg.step04_opr_events_per_revolution';
point_cloud = build_low_speed_point_cloud_local(cfg, Sensor_Config, case_data, std_angles, F_omega_deg_s);
Template = aggregate_template_from_point_cloud_local(cfg, Sensor_Config, case_data, std_angles, point_cloud, speed_diagnostic, sensor_config_path, low_speed_feature_path);

out_dir = cfg.step04_output_dir;
fig_dir = cfg.step04_figure_dir;
if exist(out_dir, 'dir') ~= 7
    mkdir(out_dir);
end
if cfg.save_figures && exist(fig_dir, 'dir') ~= 7
    mkdir(fig_dir);
end

template_file = fullfile(out_dir, sprintf('Template_OPRCenterStd_LowSpeed_AllBlades_S%s_%s.mat', ...
    sprintf('%d', cfg.sensor_ids), cfg.dataset));
low_speed_feature_path_out = fullfile(out_dir, sprintf('LowSpeed_Template_SourceData_%s.mat', cfg.dataset));
LowSpeedTemplateSourceData = struct();
LowSpeedTemplateSourceData.dataset = cfg.dataset;
LowSpeedTemplateSourceData.low_speed_case = cfg.low_speed_case;
LowSpeedTemplateSourceData.created_by = mfilename;
LowSpeedTemplateSourceData.created_on = datestr(now, 31);
LowSpeedTemplateSourceData.sensor_config_path = sensor_config_path;
LowSpeedTemplateSourceData.low_speed_feature_path = low_speed_feature_path;
LowSpeedTemplateSourceData.point_cloud_note = ...
    'Saved for calibration-method comparison; Step04 template algorithm is unchanged.';
summary_table = Template.Summary_Table; %#ok<NASGU>
metadata = Template.Metadata; %#ok<NASGU>
save(template_file, 'Template', 'metadata', '-v7.3');
save(low_speed_feature_path_out, 'LowSpeedTemplateSourceData', 'point_cloud', 'speed_diagnostic', '-v7.3');
writetable(summary_table, fullfile(out_dir, sprintf('Template_OPRCenterStd_Summary_%s.csv', cfg.dataset)));

fprintf('>>> [Step04] Saved template:\n  %s\n', template_file);
fprintf('>>> [Step04] Saved low-speed template source data:\n  %s\n', low_speed_feature_path_out);
fprintf('>>> [Step04] Template entries: %d, good/usable entries: %d\n', ...
    numel(Template.SensorBlade), nnz(ismember(string(Template.Summary_Table.quality_status), ["good", "usable"])));

if show_plots || cfg.save_figures
    plot_step04_template_local(Template, cfg, show_plots);
end


function cfg = apply_step04_defaults_local(cfg)
if ~isfield(cfg, 'step04_opr_events_per_revolution') || ...
        isempty(cfg.step04_opr_events_per_revolution)
    if isfield(cfg, 'opr_events_per_revolution') && ...
            ~isempty(cfg.opr_events_per_revolution)
        cfg.step04_opr_events_per_revolution = cfg.opr_events_per_revolution;
    else
        cfg.step04_opr_events_per_revolution = cfg.blades_num;
    end
end

% Output paths.
if ~isfield(cfg, 'step04_output_dir') || isempty(cfg.step04_output_dir)
    cfg.step04_output_dir = fullfile(cfg.output_root, 'step04_low_speed_template');
end
if ~isfield(cfg, 'step04_figure_dir') || isempty(cfg.step04_figure_dir)
    cfg.step04_figure_dir = fullfile(cfg.figure_root, 'step04_low_speed_template');
end

% Standard-angle field used as template coordinate reference.
if ~isfield(cfg, 'step04_standard_angle_source') || isempty(cfg.step04_standard_angle_source)
    cfg.step04_standard_angle_source = 'selected';  % selected, mean, median
end

% Pulse-window and point-cloud controls.
if ~isfield(cfg, 'step04_pulse_pad_fraction') || isempty(cfg.step04_pulse_pad_fraction)
    cfg.step04_pulse_pad_fraction = 0.30;
end
if ~isfield(cfg, 'step04_pulse_min_pad_points') || isempty(cfg.step04_pulse_min_pad_points)
    cfg.step04_pulse_min_pad_points = 40;
end
if ~isfield(cfg, 'step04_max_points_per_pulse') || isempty(cfg.step04_max_points_per_pulse)
    cfg.step04_max_points_per_pulse = 450;
end
if ~isfield(cfg, 'step04_max_pulses_per_sensor') || isempty(cfg.step04_max_pulses_per_sensor)
    cfg.step04_max_pulses_per_sensor = inf;
end
if ~isfield(cfg, 'step04_max_laps_per_sensor') || isempty(cfg.step04_max_laps_per_sensor)
    cfg.step04_max_laps_per_sensor = cfg.max_laps_process;
end
if ~isfield(cfg, 'step04_x_collect_abs_limit_mm') || isempty(cfg.step04_x_collect_abs_limit_mm)
    cfg.step04_x_collect_abs_limit_mm = 12.0;
end

% Template grid and smoothing.
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

% Trusted-domain selection. These values follow the spirit of the earlier
% GradientXRange030 OPRCenterStd route, but are kept here as explicit Step04
% tuning parameters.
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
if ~isfield(cfg, 'step04_domain_center_mode') || isempty(cfg.step04_domain_center_mode)
    cfg.step04_domain_center_mode = 'zero';  % zero or peak
end

% Quality gates.
if ~isfield(cfg, 'step04_min_template_points') || isempty(cfg.step04_min_template_points)
    cfg.step04_min_template_points = 1500;
end
if ~isfield(cfg, 'step04_min_template_pulses') || isempty(cfg.step04_min_template_pulses)
    cfg.step04_min_template_pulses = 20;
end
if ~isfield(cfg, 'step04_min_domain_width_mm') || isempty(cfg.step04_min_domain_width_mm)
    cfg.step04_min_domain_width_mm = 2.0 * cfg.step04_min_half_width_mm;
end
if ~isfield(cfg, 'step04_plot_max_templates_per_figure') || isempty(cfg.step04_plot_max_templates_per_figure)
    cfg.step04_plot_max_templates_per_figure = 6;
end
end


function check_step04_inputs_local(cfg, Sensor_Config, case_data)
if ~isfield(case_data, 'opr_times') || ...
        numel(case_data.opr_times) <= cfg.step04_opr_events_per_revolution
    error('Step04 requires low-speed OPR times from Step01 LowSpeed_Features_20250527.mat.');
end
if ~isfield(Sensor_Config, 'Target_Indices')
    error('Sensor_Config.Target_Indices is missing. Re-run improved Step01.');
end
if ~isfield(Sensor_Config, 'Standard_Relative_Angles') && ...
        ~isfield(Sensor_Config, 'Standard_Relative_Angles_Mean') && ...
        ~isfield(Sensor_Config, 'Standard_Relative_Angles_Median')
    error('Sensor_Config does not contain standard relative angles. Re-run Step01.');
end
if isfield(Sensor_Config, 'OPR_Timing_Method') && isfield(cfg, 'opr_timing_method')
    if ~strcmpi(strtrim(Sensor_Config.OPR_Timing_Method), strtrim(cfg.opr_timing_method))
        warning('Step04:OPRReferenceMismatch', ['cfg.opr_timing_method=%s, but Sensor_Config.OPR_Timing_Method=%s. ' ...
            'The low-speed template coordinate may not match later dynamic coordinates.'], ...
            cfg.opr_timing_method, Sensor_Config.OPR_Timing_Method);
    end
end
if isfield(case_data, 'opr_timing_method') && isfield(Sensor_Config, 'OPR_Timing_Method')
    if ~strcmpi(strtrim(case_data.opr_timing_method), strtrim(Sensor_Config.OPR_Timing_Method))
        warning('Step04:LowSpeedReferenceMismatch', ...
            'LowSpeed_Features OPR method (%s) differs from Sensor_Config (%s).', ...
            case_data.opr_timing_method, Sensor_Config.OPR_Timing_Method);
    end
end
for sid = cfg.sensor_ids
    if sid > numel(case_data.channels) || isempty(case_data.channels(sid).arrival_times)
        warning('Step04:MissingLowSpeedChannel', 'Low-speed features for CH%d are missing or empty.', sid);
    end
end
end


function std_angles = resolve_standard_angles_local(Sensor_Config, source_name)
source_name = lower(strtrim(source_name));
switch source_name
    case 'mean'
        if isfield(Sensor_Config, 'Standard_Relative_Angles_Mean')
            std_angles = Sensor_Config.Standard_Relative_Angles_Mean;
        else
            std_angles = Sensor_Config.Standard_Relative_Angles;
        end
    case 'median'
        if isfield(Sensor_Config, 'Standard_Relative_Angles_Median')
            std_angles = Sensor_Config.Standard_Relative_Angles_Median;
        else
            std_angles = Sensor_Config.Standard_Relative_Angles;
        end
    otherwise
        std_angles = Sensor_Config.Standard_Relative_Angles;
end
end


function [F_omega_deg_s, diagnostic] = build_low_speed_speed_interpolant_local( ...
    opr_times, opr_events_per_revolution)
opr_times = opr_times(:);
if numel(opr_times) <= opr_events_per_revolution
    error('Not enough OPR events (%d) for events/rev=%d.', ...
        numel(opr_times), opr_events_per_revolution);
end
rev_period = opr_times((opr_events_per_revolution + 1):end) - ...
    opr_times(1:(end - opr_events_per_revolution));
speed_time = opr_times(1:(end - opr_events_per_revolution));
speed_deg_s = 360 ./ max(rev_period, eps);
valid = isfinite(speed_time) & isfinite(speed_deg_s) & speed_deg_s > 0;
if nnz(valid) < 2
    error('Not enough valid OPR periods to build low-speed speed interpolant.');
end
F_omega_deg_s = griddedInterpolant(speed_time(valid), speed_deg_s(valid), 'linear', 'nearest');
diagnostic = struct();
diagnostic.speed_time_s = speed_time(valid);
diagnostic.speed_deg_s = speed_deg_s(valid);
diagnostic.rpm = speed_deg_s(valid) / 360 * 60;
end


function point_cloud = build_low_speed_point_cloud_local(cfg, Sensor_Config, case_data, std_angles, F_omega_deg_s)
case_dir = fullfile(cfg.dataset_root, cfg.low_speed_case);
if ~isfolder(case_dir)
    error('Low-speed raw data folder not found: %s', case_dir);
end

file_ranges = case_data.file_ranges;
file_ids = case_data.file_ids(:).';
opr_times = case_data.opr_times(:);

max_entries = numel(cfg.sensor_ids) * cfg.blades_num;
entry_template = struct( ...
    'sensor_id', NaN, 'blade_id', NaN, ...
    'x_mm', [], 'v', [], 't_s', [], 'pulse_index', [], 'lap_index', [], ...
    'baseline_by_pulse', [], 'feature_by_pulse', [], 'pulse_count', 0, 'point_count', 0);
point_cloud = repmat(entry_template, max_entries, 1);
idx_entry = 0;
for sid = cfg.sensor_ids
    for blade_id = 1:cfg.blades_num
        idx_entry = idx_entry + 1;
        point_cloud(idx_entry).sensor_id = sid;
        point_cloud(idx_entry).blade_id = blade_id;
    end
end

fprintf('>>> [Step04] Building low-speed x-V point cloud from raw waveforms...\n');
for iSensor = 1:numel(cfg.sensor_ids)
    sid = cfg.sensor_ids(iSensor);
    if sid > numel(case_data.channels) || isempty(case_data.channels(sid).arrival_times)
        continue;
    end
    ch = case_data.channels(sid);
    pulse_times = [ch.start_times(:), ch.end_times(:), ch.arrival_times(:)];
    pulse_count_total = size(pulse_times, 1);
    if pulse_count_total == 0
        continue;
    end

    start_idx = get_target_index_local(Sensor_Config, sid);
    if ~isfinite(start_idx)
        warning('Step04:MissingTargetIndex', 'CH%d has no target index in Sensor_Config; skipping.', sid);
        continue;
    end

    max_pulses = min([pulse_count_total, cfg.step04_max_pulses_per_sensor, ...
        start_idx + cfg.step04_max_laps_per_sensor * cfg.blades_num - 1]);
    if ~isfinite(max_pulses)
        max_pulses = pulse_count_total;
    end
    max_pulses = min(pulse_count_total, floor(max_pulses));

    fprintf('>>> [Step04] CH%d: processing %d/%d low-speed pulses.\n', sid, max_pulses, pulse_count_total);

    % Read files one by one. This avoids loading the whole raw channel into
    % memory and uses the same time-index offsets detected by Step01.
    for iFile = 1:numel(file_ids)
        file_id = file_ids(iFile);
        filepath = fullfile(case_dir, sprintf('4-%d-%d.mat', sid, file_id));
        if ~isfile(filepath)
            continue;
        end
        raw = load_raw_case_channel_local(filepath, cfg.sample_rate_hz);
        if isempty(raw.t)
            continue;
        end
        offset = find_file_offset_local(file_ranges, file_id);
        t_global = raw.t + offset;
        v = raw.v;

        if iFile == 1 && cfg.initial_trim_points > 0 && numel(t_global) > cfg.initial_trim_points
            keep = (cfg.initial_trim_points + 1):numel(t_global);
            t_global = t_global(keep);
            v = v(keep);
        end

        if isempty(t_global)
            continue;
        end
        t_file_start = t_global(1);
        t_file_end = t_global(end);

        pulse_idx_in_file = find(pulse_times(1:max_pulses, 2) >= t_file_start & ...
            pulse_times(1:max_pulses, 1) <= t_file_end);
        if isempty(pulse_idx_in_file)
            continue;
        end

        for k = 1:numel(pulse_idx_in_file)
            pulse_idx = pulse_idx_in_file(k);
            start_t = pulse_times(pulse_idx, 1);
            end_t = pulse_times(pulse_idx, 2);
            arrival_t = pulse_times(pulse_idx, 3); %#ok<NASGU>
            if ~isfinite(start_t) || ~isfinite(end_t) || end_t <= start_t
                continue;
            end

            blade_id = mod((pulse_idx - start_idx), cfg.blades_num) + 1;
            lap_index = floor((pulse_idx - start_idx) / cfg.blades_num) + 1;
            if blade_id < 1 || blade_id > cfg.blades_num || lap_index < 1
                continue;
            end
            if sid > size(std_angles, 1) || ~isfinite(std_angles(sid, blade_id))
                continue;
            end

            pulse_width = end_t - start_t;
            pad_s = max(cfg.step04_pulse_min_pad_points / cfg.sample_rate_hz, ...
                cfg.step04_pulse_pad_fraction * pulse_width);
            win_start = start_t - pad_s;
            win_end = end_t + pad_s;
            idx = find(t_global >= win_start & t_global <= win_end);
            if numel(idx) < 5
                continue;
            end
            if numel(idx) > cfg.step04_max_points_per_pulse
                idx = idx(unique(round(linspace(1, numel(idx), cfg.step04_max_points_per_pulse))));
            end

            t_seg = t_global(idx);
            v_seg = v(idx);
            [x_mm, keep_mask] = map_time_to_x_local(t_seg, opr_times, F_omega_deg_s, std_angles(sid, blade_id), cfg);
            if ~any(keep_mask)
                continue;
            end
            t_seg = t_seg(keep_mask);
            v_seg = v_seg(keep_mask);
            x_mm = x_mm(keep_mask);

            within = isfinite(x_mm) & isfinite(v_seg) & abs(x_mm) <= cfg.step04_x_collect_abs_limit_mm;
            if ~any(within)
                continue;
            end
            x_mm = x_mm(within);
            v_seg = v_seg(within);
            t_seg = t_seg(within);

            baseline = estimate_pulse_baseline_local(v_seg);
            feature = max(v_seg) - baseline;
            entry_idx = sensor_blade_entry_index_local(cfg, sid, blade_id);
            point_cloud(entry_idx).x_mm = [point_cloud(entry_idx).x_mm; x_mm(:)]; %#ok<AGROW>
            point_cloud(entry_idx).v = [point_cloud(entry_idx).v; v_seg(:)]; %#ok<AGROW>
            point_cloud(entry_idx).t_s = [point_cloud(entry_idx).t_s; t_seg(:)]; %#ok<AGROW>
            point_cloud(entry_idx).pulse_index = [point_cloud(entry_idx).pulse_index; repmat(pulse_idx, numel(x_mm), 1)]; %#ok<AGROW>
            point_cloud(entry_idx).lap_index = [point_cloud(entry_idx).lap_index; repmat(lap_index, numel(x_mm), 1)]; %#ok<AGROW>
            point_cloud(entry_idx).baseline_by_pulse = [point_cloud(entry_idx).baseline_by_pulse; baseline]; %#ok<AGROW>
            point_cloud(entry_idx).feature_by_pulse = [point_cloud(entry_idx).feature_by_pulse; feature]; %#ok<AGROW>
        end
    end
end

for i = 1:numel(point_cloud)
    point_cloud(i).point_count = numel(point_cloud(i).x_mm);
    point_cloud(i).pulse_count = numel(unique(point_cloud(i).pulse_index));
end
end


function raw = load_raw_case_channel_local(filepath, sample_rate_hz)
loaded = load(filepath);
fn = fieldnames(loaded);
if isempty(fn)
    raw = struct('t', [], 'v', []);
    return;
end
arr = loaded.(fn{1});
if isempty(arr) || size(arr, 2) < 2
    raw = struct('t', [], 'v', []);
    return;
end
arr(arr(:, 1) == 0, :) = [];
raw = struct();
raw.t = arr(:, 1) / sample_rate_hz;
raw.v = arr(:, 2);
end


function offset = find_file_offset_local(file_ranges, file_id)
offset = 0;
if isempty(file_ranges)
    return;
end
ids = [file_ranges.file_id];
idx = find(ids == file_id, 1, 'first');
if ~isempty(idx) && isfield(file_ranges, 'offset') && isfinite(file_ranges(idx).offset)
    offset = file_ranges(idx).offset;
end
end


function idx = sensor_blade_entry_index_local(cfg, sid, blade_id)
sensor_pos = find(cfg.sensor_ids == sid, 1, 'first');
idx = (sensor_pos - 1) * cfg.blades_num + blade_id;
end


function start_idx = get_target_index_local(Sensor_Config, sid)
start_idx = NaN;
try
    if isKey(Sensor_Config.Target_Indices, sid)
        start_idx = Sensor_Config.Target_Indices(sid);
    end
catch
    start_idx = NaN;
end
end


function [x_mm, keep_mask] = map_time_to_x_local(t_seg, opr_times, F_omega_deg_s, theta_std_deg, cfg)
t_seg = t_seg(:);
prev_idx = find_previous_opr_indices_local(opr_times, t_seg);
keep_mask = isfinite(prev_idx) & prev_idx >= 1 & prev_idx <= numel(opr_times);
x_mm = nan(size(t_seg));
if ~any(keep_mask)
    return;
end
idx_valid = find(keep_mask);
for ii = 1:numel(idx_valid)
    i = idx_valid(ii);
    t_ref = opr_times(prev_idx(i));
    t_now = t_seg(i);
    if ~isfinite(t_now) || t_now < t_ref
        keep_mask(i) = false;
        continue;
    end
    theta_actual = integrate_angle_local(t_ref, t_now, F_omega_deg_s);
    theta_diff = wrap_to_180_local(theta_actual - theta_std_deg);
    x_mm(i) = theta_diff * (pi / 180) * cfg.r_tip_mm;
end
keep_mask = keep_mask & isfinite(x_mm);
end


function prev_idx = find_previous_opr_indices_local(opr_times, t)
prev_idx = nan(size(t));
edges = [-inf; opr_times(:); inf];
bin = discretize(t(:), edges);
idx = bin - 1;
idx(idx < 1 | idx > numel(opr_times)) = NaN;
prev_idx(:) = idx;
end


function theta_deg = integrate_angle_local(t0, t1, F_omega_deg_s)
if t1 <= t0
    theta_deg = 0;
    return;
end
% Five nodes are sufficient for the low-speed quasi-steady reference and are
% much faster than dense numerical integration for every raw waveform point.
t_grid = linspace(t0, t1, 5);
theta_deg = trapz(t_grid, F_omega_deg_s(t_grid));
end


function baseline = estimate_pulse_baseline_local(v)
v = v(:);
if isempty(v)
    baseline = NaN;
    return;
end
n = numel(v);
n_edge = max(3, min(floor(0.15 * n), 50));
edge_values = [v(1:n_edge); v((n-n_edge+1):n)];
baseline = median(edge_values, 'omitnan');
if ~isfinite(baseline)
    baseline = prctile(v, 10);
end
end


function Template = aggregate_template_from_point_cloud_local(cfg, Sensor_Config, case_data, std_angles, point_cloud, speed_diagnostic, sensor_config_path, low_speed_feature_path)
entry_count = numel(point_cloud);
SensorBlade = repmat(make_empty_template_entry_local(), 1, entry_count);
summary_rows = repmat(struct( ...
    'sensor_id', NaN, 'blade_id', NaN, ...
    'pulse_count', NaN, 'point_count', NaN, ...
    'x_min_mm', NaN, 'x_max_mm', NaN, ...
    'xc_mm', NaN, 'xc_reference_source', '', ...
    'domain_left_mm', NaN, 'domain_right_mm', NaN, 'domain_width_mm', NaN, ...
    'baseline_v', NaN, 'amplitude_v', NaN, 'max_abs_gradient_v_per_mm', NaN, ...
    'valid_bin_count', NaN, 'coverage_fraction', NaN, ...
    'quality_status', ''), entry_count, 1);

for i = 1:entry_count
    pc = point_cloud(i);
    sid = pc.sensor_id;
    blade_id = pc.blade_id;
    [entry, row] = aggregate_one_sensor_blade_local(cfg, pc, sid, blade_id);
    SensorBlade(i) = entry;
    summary_rows(i) = row;
end

summary_table = struct2table(summary_rows);

Template = struct();
Template.Dataset = cfg.dataset;
Template.Case_Name = cfg.low_speed_case;
Template.CreatedBy = mfilename;
Template.CreatedOn = datestr(now, 31);
Template.Route_Dir = cfg.route_dir;
Template.Step01_Sensor_Config_Path = sensor_config_path;
Template.LowSpeed_Feature_Path = low_speed_feature_path;
Template.Sensor_IDs = cfg.sensor_ids;
Template.Blades_Num = cfg.blades_num;
Template.OPR_Events_Per_Revolution = cfg.step04_opr_events_per_revolution;
Template.OPR_Events_Per_Revolution_Policy = ...
    get_string_field_local(cfg, 'opr_events_per_revolution_policy', 'unspecified');
Template.R_Tip_mm = cfg.r_tip_mm;
Template.OPR_Timing_Method = get_string_field_local(Sensor_Config, 'OPR_Timing_Method', get_string_field_local(case_data, 'opr_timing_method', 'unknown'));
Template.Probe_Arrival_Method = get_string_field_local(Sensor_Config, 'Probe_Arrival_Method', get_string_field_local(case_data, 'probe_arrival_method', 'unknown'));
Template.Time_Index_Mode = get_string_field_local(Sensor_Config, 'Time_Index_Mode', get_string_field_local(case_data, 'time_index_mode', 'unknown'));
Template.Time_Index_Mode_Detected = get_string_field_local(Sensor_Config, 'Time_Index_Mode_Detected', get_string_field_local(case_data, 'time_index_mode_detected', 'unknown'));
Template.Blade_ID_Definition = get_string_field_local(Sensor_Config, 'Blade_ID_Definition', 'internal Blade 1 from low-speed reference');
Template.Standard_Relative_Angles = std_angles;
if isfield(Sensor_Config, 'Standard_Relative_Angles_Std')
    Template.Standard_Relative_Angles_Std = Sensor_Config.Standard_Relative_Angles_Std;
end
if isfield(Sensor_Config, 'Standard_Relative_Angles_Count')
    Template.Standard_Relative_Angles_Count = Sensor_Config.Standard_Relative_Angles_Count;
end
Template.Speed_Diagnostic = speed_diagnostic;
Template.PointCloud_Settings = struct( ...
    'pulse_pad_fraction', cfg.step04_pulse_pad_fraction, ...
    'max_points_per_pulse', cfg.step04_max_points_per_pulse, ...
    'x_collect_abs_limit_mm', cfg.step04_x_collect_abs_limit_mm);
Template.Template_Settings = struct( ...
    'grid_dx_mm', cfg.step04_grid_dx_mm, ...
    'min_bin_count', cfg.step04_min_bin_count, ...
    'smooth_span_bins', cfg.step04_smooth_span_bins, ...
    'template_xc_policy', 'top85_weighted_centroid_on_effective_template', ...
    'gradient_min_ratio', cfg.step04_gradient_min_ratio, ...
    'amplitude_min_ratio', cfg.step04_amplitude_min_ratio, ...
    'min_half_width_mm', cfg.step04_min_half_width_mm, ...
    'max_half_width_mm', cfg.step04_max_half_width_mm, ...
    'domain_center_mode', cfg.step04_domain_center_mode);
Template.SensorBlade = SensorBlade;
Template.Summary_Table = summary_table;
Template.Metadata = build_template_metadata_local(Template, cfg);
end


function empty = make_empty_template_entry_local()
empty = struct( ...
    'sensor_id', NaN, 'blade_id', NaN, ...
    'x_grid', [], 'v_grid', [], 'v_grid_raw', [], 'v_grid_baseline_removed', [], ...
    'dv_dx', [], 'weight_grid', [], 'count_grid', [], 'valid_grid_mask', [], ...
    'domain_effective_mask', [], 'domain_mask', [], ...
    'x_domain', [NaN NaN], 'x_coordinate_frame', 'OPRCenterStd_abs', ...
    'x_center_subtracted', false, 'xc', NaN, 'xc_mm', NaN, ...
    'xc_abs', NaN, 'xc_reference_source', '', ...
    'x_grid_abs', [], 'x_domain_abs', [NaN NaN], ...
    'baseline', NaN, 'amplitude', NaN, ...
    'pulse_count', 0, 'point_count', 0, 'quality_status', '', ...
    'x_points_preview', [], 'x_abs_points_preview', [], ...
    'v_points_preview', []);
end


function [entry, row] = aggregate_one_sensor_blade_local(cfg, pc, sid, blade_id)
entry = make_empty_template_entry_local();
entry.sensor_id = sid;
entry.blade_id = blade_id;

row = struct( ...
    'sensor_id', sid, 'blade_id', blade_id, ...
    'pulse_count', pc.pulse_count, 'point_count', pc.point_count, ...
    'x_min_mm', NaN, 'x_max_mm', NaN, ...
    'xc_mm', NaN, 'xc_reference_source', '', ...
    'domain_left_mm', NaN, 'domain_right_mm', NaN, 'domain_width_mm', NaN, ...
    'baseline_v', NaN, 'amplitude_v', NaN, 'max_abs_gradient_v_per_mm', NaN, ...
    'valid_bin_count', NaN, 'coverage_fraction', NaN, ...
    'quality_status', 'empty');

x = pc.x_mm(:);
v = pc.v(:);
valid = isfinite(x) & isfinite(v);
x = x(valid);
v = v(valid);
entry.point_count = numel(x);
entry.pulse_count = pc.pulse_count;
if isempty(x)
    entry.quality_status = 'empty';
    return;
end

row.x_min_mm = min(x);
row.x_max_mm = max(x);

% Robust x range for grid construction. Avoid letting a few extreme points
% create a very wide sparse template.
x_lo = prctile(x, 0.5);
x_hi = prctile(x, 99.5);
x_lo = max(x_lo, -cfg.step04_x_collect_abs_limit_mm);
x_hi = min(x_hi, cfg.step04_x_collect_abs_limit_mm);
if x_hi <= x_lo
    x_lo = min(x);
    x_hi = max(x);
end
if x_hi <= x_lo
    entry.quality_status = 'degenerate_x_range';
    row.quality_status = entry.quality_status;
    return;
end

edges = x_lo:cfg.step04_grid_dx_mm:x_hi;
if numel(edges) < 5
    edges = linspace(x_lo, x_hi, 20);
end
x_grid = 0.5 * (edges(1:end-1) + edges(2:end));
x_grid = x_grid(:);
bin = discretize(x, edges);
n_bins = numel(x_grid);
v_grid_raw = nan(n_bins, 1);
count_grid = zeros(n_bins, 1);
for ib = 1:n_bins
    mask = bin == ib;
    count_grid(ib) = nnz(mask);
    if count_grid(ib) > 0
        v_grid_raw(ib) = median(v(mask), 'omitnan');
    end
end

valid_bin_mask = count_grid >= cfg.step04_min_bin_count & isfinite(v_grid_raw);
if nnz(valid_bin_mask) < 5
    entry.quality_status = 'insufficient_coverage';
    row.quality_status = entry.quality_status;
    row.valid_bin_count = nnz(valid_bin_mask);
    row.coverage_fraction = mean(valid_bin_mask);
    return;
end

v_interp = v_grid_raw;
missing = ~isfinite(v_interp);
if any(missing)
    v_interp(missing) = interp1(x_grid(valid_bin_mask), v_grid_raw(valid_bin_mask), ...
        x_grid(missing), 'linear', 'extrap');
end
v_smooth = smooth_template_local(v_interp, cfg.step04_smooth_span_bins);

baseline = estimate_template_baseline_local(x, v, v_smooth);
amplitude = max(v_smooth, [], 'omitnan') - baseline;
if ~isfinite(amplitude) || amplitude <= 0
    amplitude = max(v_smooth, [], 'omitnan') - min(v_smooth, [], 'omitnan');
end
if ~isfinite(amplitude)
    amplitude = 0;
end

dv_dx = gradient(v_smooth, x_grid);
max_abs_grad = max(abs(dv_dx(valid_bin_mask)), [], 'omitnan');
if ~isfinite(max_abs_grad)
    max_abs_grad = 0;
end

[x_domain, domain_mask, effective_mask] = determine_trusted_domain_local(cfg, x_grid(:), v_smooth(:), dv_dx(:), count_grid(:), valid_bin_mask(:), baseline, amplitude, max_abs_grad);
weight_grid = build_template_weight_grid_local(count_grid(:), valid_bin_mask(:), domain_mask(:), dv_dx(:), max_abs_grad);
[xc_mm, xc_source] = estimate_template_xc_from_grid_local( ...
    x_grid(:), v_smooth(:), baseline, valid_bin_mask(:), effective_mask(:));
if ~isfinite(xc_mm)
    xc_mm = 0;
    xc_source = 'zero_fallback_invalid_template_xc';
end

x_model = x - xc_mm;
x_grid_model = x_grid(:) - xc_mm;
x_domain_model = x_domain - xc_mm;

entry.x_grid = x_grid_model(:);
entry.x_coordinate_frame = 'OPRCenterStd_centered';
entry.x_center_subtracted = true;
entry.xc = xc_mm;
entry.xc_mm = xc_mm;
entry.xc_abs = xc_mm;
entry.xc_reference_source = xc_source;
entry.x_grid_abs = x_grid(:);
entry.v_grid = v_smooth(:);
entry.v_grid_raw = v_grid_raw(:);
entry.v_grid_baseline_removed = v_smooth(:) - baseline;
entry.dv_dx = dv_dx(:);
entry.weight_grid = weight_grid(:);
entry.count_grid = count_grid(:);
entry.valid_grid_mask = valid_bin_mask(:);
entry.domain_effective_mask = effective_mask(:);
entry.domain_mask = domain_mask(:);
entry.x_domain = x_domain_model;
entry.x_domain_abs = x_domain;
entry.baseline = baseline;
entry.amplitude = amplitude;
entry.quality_status = evaluate_template_quality_local(cfg, pc.pulse_count, numel(x), x_domain, valid_bin_mask);

preview_n = min(numel(x), 4000);
if numel(x) > preview_n
    idx = unique(round(linspace(1, numel(x), preview_n)));
else
    idx = 1:numel(x);
end
entry.x_points_preview = x(idx);
entry.x_abs_points_preview = x(idx);
entry.x_points_preview = x_model(idx);
entry.v_points_preview = v(idx);

row.x_min_mm = min(x_model);
row.x_max_mm = max(x_model);
row.domain_left_mm = x_domain_model(1);
row.domain_right_mm = x_domain_model(2);
row.domain_width_mm = x_domain_model(2) - x_domain_model(1);
row.xc_mm = xc_mm;
row.xc_reference_source = xc_source;
row.baseline_v = baseline;
row.amplitude_v = amplitude;
row.max_abs_gradient_v_per_mm = max_abs_grad;
row.valid_bin_count = nnz(valid_bin_mask);
row.coverage_fraction = mean(valid_bin_mask);
row.quality_status = entry.quality_status;
end


function [xc, source] = estimate_template_xc_from_grid_local( ...
    x_grid, v_grid, baseline, valid_bin_mask, effective_mask)
% Static sensor-frame center used only as metadata/eta prior.  The template
% itself remains in the absolute OPRCenterStd x frame.
xc = NaN;
source = 'invalid_template_center';

x = x_grid(:);
v = v_grid(:);
valid = isfinite(x) & isfinite(v);

if nargin >= 4 && numel(valid_bin_mask) == numel(x)
    valid = valid & valid_bin_mask(:);
end
if nargin >= 5 && numel(effective_mask) == numel(x)
    valid = valid & effective_mask(:);
end

if nnz(valid) < 5
    return;
end

xv = x(valid);
vv = v(valid);

if ~isfinite(baseline)
    baseline = median(vv, 'omitnan');
end

w = max(vv - baseline, 0);
if nnz(isfinite(w) & w > 0) < 5
    w = max(vv - min(vv, [], 'omitnan'), 0);
end

wp = w(isfinite(w) & w > 0);
if numel(wp) < 5
    [~, imax] = max(vv);
    xc = xv(imax);
    source = 'template_peak_fallback_xc';
    return;
end

thr = prctile(wp, 85);
high = isfinite(w) & w >= thr;

if nnz(high) < 3
    [~, imax] = max(w);
    xc = xv(imax);
    source = 'template_peak_fallback_xc';
    return;
end

den = sum(w(high), 'omitnan');
if den <= eps
    xc = median(xv(high), 'omitnan');
    source = 'template_top85_median_xc';
else
    xc = sum(xv(high) .* w(high), 'omitnan') ./ den;
    source = 'template_top85_weighted_centroid_xc';
end
end


function v_smooth = smooth_template_local(v, span)
v = v(:);
if numel(v) < 3
    v_smooth = v;
    return;
end
span = min(span, numel(v));
if mod(span, 2) == 0
    span = span - 1;
end
span = max(span, 3);
try
    v_smooth = smoothdata(v, 'movmedian', span, 'omitnan');
    v_smooth = smoothdata(v_smooth, 'movmean', span, 'omitnan');
catch
    v_smooth = movmedian(v, span, 'omitnan');
    v_smooth = movmean(v_smooth, span, 'omitnan');
end
end


function baseline = estimate_template_baseline_local(x, v, v_grid)
% Prefer the low-voltage percentile because not every template has clean
% far-field edge coverage after the OPRCenterStd coordinate transform.
if isempty(v)
    baseline = NaN;
    return;
end
baseline = prctile(v, 10);
if ~isfinite(baseline)
    baseline = median(v, 'omitnan');
end
if ~isfinite(baseline) && ~isempty(v_grid)
    baseline = min(v_grid, [], 'omitnan');
end
end


function [x_domain, domain_mask, effective] = determine_trusted_domain_local(cfg, x_grid, v_grid, dv_dx, count_grid, valid_bin_mask, baseline, amplitude, max_abs_grad)
%DETERMINE_TRUSTED_DOMAIN_LOCAL Select a practical template trust interval.
%
% Important: this is a coverage-and-gradient trust interval, not a plot-only
% decoration. The previous version used (amplitude OR gradient) and then found
% a contiguous region around x=0. For a bell-shaped pulse this makes almost
% the whole high-amplitude plateau "trusted", so the domain often expands to
% the configured max half-width and becomes visually/physically uninformative.
%
% The corrected logic follows the earlier GradientXRange030 route more closely:
%   1. Find bins where the template gradient is informative.
%   2. Optionally require the signal amplitude to be above a low threshold.
%   3. Use the outermost informative bins on the left and right of the chosen
%      center to define the span.
%   4. Enforce practical min/max half-widths and clip to available coverage.

x_grid = x_grid(:);
v_grid = v_grid(:);
dv_dx = dv_dx(:);
valid_bin_mask = valid_bin_mask(:);
count_grid = count_grid(:); %#ok<NASGU>

domain_mask = false(size(x_grid));
effective = false(size(x_grid));
x_domain = [NaN NaN];

if isempty(x_grid) || nnz(valid_bin_mask) < 5
    return;
end

amp_norm = zeros(size(x_grid));
if isfinite(amplitude) && amplitude > 0 && isfinite(baseline)
    amp_norm = max(v_grid - baseline, 0) ./ amplitude;
end

grad_mask = false(size(x_grid));
if isfinite(max_abs_grad) && max_abs_grad > 0
    grad_mask = abs(dv_dx) >= cfg.step04_gradient_min_ratio * max_abs_grad;
end

% Amplitude threshold is an auxiliary guard, not a substitute for gradient.
effective = valid_bin_mask & grad_mask;
if isfinite(cfg.step04_amplitude_min_ratio) && cfg.step04_amplitude_min_ratio > 0
    effective = effective & amp_norm >= cfg.step04_amplitude_min_ratio;
end

switch lower(strtrim(cfg.step04_domain_center_mode))
    case 'peak'
        [~, i_center] = max(v_grid);
        if isempty(i_center) || ~isfinite(i_center)
            [~, i_center] = min(abs(x_grid));
        end
    otherwise
        [~, i_center] = min(abs(x_grid));
end
if isempty(i_center) || ~isfinite(i_center)
    return;
end
center_x = 0;
if strcmpi(cfg.step04_domain_center_mode, 'peak')
    center_x = x_grid(i_center);
end

left_candidates = find(x_grid < center_x & effective);
right_candidates = find(x_grid > center_x & effective);

% Fallback: if one side lacks gradient bins, use amplitude/coverage support.
if isempty(left_candidates) || isempty(right_candidates)
    amp_support = valid_bin_mask & amp_norm >= max(cfg.step04_amplitude_min_ratio, 0.05);
    left_candidates = find(x_grid < center_x & amp_support);
    right_candidates = find(x_grid > center_x & amp_support);
end
if isempty(left_candidates) || isempty(right_candidates)
    return;
end

left_L = max(center_x - x_grid(left_candidates(1)), 0);
right_L = max(x_grid(right_candidates(end)) - center_x, 0);

% Practical constraints inherited from the OPRCenterStd GradientXRange route.
left_L = max(left_L, cfg.step04_min_half_width_mm);
right_L = max(right_L, cfg.step04_min_half_width_mm);
if isfinite(cfg.step04_max_half_width_mm) && cfg.step04_max_half_width_mm > 0
    left_L = min(left_L, cfg.step04_max_half_width_mm);
    right_L = min(right_L, cfg.step04_max_half_width_mm);
end

coverage_left = min(x_grid(valid_bin_mask));
coverage_right = max(x_grid(valid_bin_mask));
left_x = max(center_x - left_L, coverage_left);
right_x = min(center_x + right_L, coverage_right);

if right_x <= left_x
    return;
end

x_domain = [left_x, right_x];
domain_mask = x_grid >= left_x & x_grid <= right_x & valid_bin_mask;
end

function weight_grid = build_template_weight_grid_local(count_grid, valid_bin_mask, domain_mask, dv_dx, max_abs_grad)
count_grid = count_grid(:);
valid_bin_mask = valid_bin_mask(:);
domain_mask = domain_mask(:);
dv_dx = dv_dx(:);
weight_grid = zeros(size(count_grid));
if isempty(count_grid)
    return;
end
max_count = max(count_grid(valid_bin_mask), [], 'omitnan');
if ~isfinite(max_count) || max_count <= 0
    max_count = 1;
end
coverage_w = sqrt(max(count_grid, 0) ./ max_count);
if isfinite(max_abs_grad) && max_abs_grad > 0
    gradient_w = 0.25 + 0.75 * min(abs(dv_dx) ./ max_abs_grad, 1);
else
    gradient_w = ones(size(count_grid));
end
weight_grid(valid_bin_mask) = coverage_w(valid_bin_mask) .* gradient_w(valid_bin_mask);
weight_grid(~domain_mask) = 0.25 * weight_grid(~domain_mask);
weight_grid(~isfinite(weight_grid)) = 0;
end


function quality = evaluate_template_quality_local(cfg, pulse_count, point_count, x_domain, valid_bin_mask)
if point_count <= 0
    quality = 'empty';
    return;
end
if pulse_count < cfg.step04_min_template_pulses || point_count < cfg.step04_min_template_points
    quality = 'weak_points';
    return;
end
if any(~isfinite(x_domain)) || x_domain(2) <= x_domain(1)
    quality = 'no_trusted_domain';
    return;
end
if (x_domain(2) - x_domain(1)) < cfg.step04_min_domain_width_mm
    quality = 'narrow_domain';
    return;
end
coverage_fraction = mean(valid_bin_mask);
if coverage_fraction < 0.20
    quality = 'sparse_coverage';
elseif coverage_fraction < 0.35
    quality = 'usable';
else
    quality = 'good';
end
end


function metadata = build_template_metadata_local(Template, cfg)
metadata = struct();
metadata.dataset = Template.Dataset;
metadata.case_name = Template.Case_Name;
metadata.created_by = Template.CreatedBy;
metadata.created_on = Template.CreatedOn;
metadata.description = ['Low-speed non-parametric OPRCenterStd sensor-blade templates. ' ...
    'Each template maps spatial coordinate x_mm to measured sensor voltage.'];
metadata.forward_model = 'V ~= T_{sensor,blade}(x_mm)';
metadata.intended_next_steps = {'Step05 dynamic waveform map', 'Step06 direct-template synchronous waveform inversion'};
metadata.cfg_step04 = Template.Template_Settings;
metadata.point_cloud_settings = Template.PointCloud_Settings;
metadata.opr_timing_method = Template.OPR_Timing_Method;
metadata.probe_arrival_method = Template.Probe_Arrival_Method;
metadata.blade_id_definition = Template.Blade_ID_Definition;
metadata.sensor_ids = cfg.sensor_ids;
metadata.blades_num = cfg.blades_num;
metadata.opr_events_per_revolution = Template.OPR_Events_Per_Revolution;
metadata.opr_events_per_revolution_policy = Template.OPR_Events_Per_Revolution_Policy;
end


function plot_step04_template_local(Template, cfg, show_plots)
fig_dir = cfg.step04_figure_dir;
if exist(fig_dir, 'dir') ~= 7
    mkdir(fig_dir);
end
sensor_ids = Template.Sensor_IDs;
blades_num = Template.Blades_Num;
colors = lines(blades_num);

% Figure 1: clean template and gradient curves by sensor.  The raw fitting
% points are intentionally shown in a separate figure to avoid confusing the
% trusted-domain markers with real waveform data.
for iSensor = 1:numel(sensor_ids)
    sid = sensor_ids(iSensor);
    fig = figure('Name', sprintf('Step04 clean low-speed templates CH%d', sid), ...
        'Color', 'w', 'Position', [80, 80, 1400, 850], 'NumberTitle', 'off', ...
        'Visible', visibility_state_local(show_plots));
    tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

    ax1 = nexttile;
    hold(ax1, 'on'); grid(ax1, 'on'); box(ax1, 'on');
    for blade_id = 1:blades_num
        entry = get_template_entry_local(Template, sid, blade_id);
        if isempty(entry) || isempty(entry.x_grid)
            continue;
        end
        plot(ax1, entry.x_grid, entry.v_grid, '-', 'LineWidth', 1.6, ...
            'Color', colors(blade_id, :), ...
            'DisplayName', sprintf('B%d (%s)', blade_id, entry.quality_status));
    end
    yl = ylim(ax1);
    for blade_id = 1:blades_num
        entry = get_template_entry_local(Template, sid, blade_id);
        if isempty(entry) || any(~isfinite(entry.x_domain))
            continue;
        end
        draw_domain_boundaries_local(ax1, entry.x_domain, yl, colors(blade_id, :));
    end
    xlabel(ax1, 'x in OPRCenterStd frame (mm)');
    ylabel(ax1, 'Voltage (V)');
    title(ax1, sprintf('CH%d low-speed non-parametric templates', sid));
    legend(ax1, 'Location', 'eastoutside');

    ax2 = nexttile;
    hold(ax2, 'on'); grid(ax2, 'on'); box(ax2, 'on');
    for blade_id = 1:blades_num
        entry = get_template_entry_local(Template, sid, blade_id);
        if isempty(entry) || isempty(entry.x_grid)
            continue;
        end
        plot(ax2, entry.x_grid, entry.dv_dx, '-', 'LineWidth', 1.3, ...
            'Color', colors(blade_id, :), 'DisplayName', sprintf('B%d', blade_id));
    end
    yline(ax2, 0, 'k:', 'HandleVisibility', 'off');
    yl = ylim(ax2);
    for blade_id = 1:blades_num
        entry = get_template_entry_local(Template, sid, blade_id);
        if isempty(entry) || any(~isfinite(entry.x_domain))
            continue;
        end
        draw_domain_boundaries_local(ax2, entry.x_domain, yl, colors(blade_id, :));
        if isfield(entry, 'domain_effective_mask') && ~isempty(entry.domain_effective_mask)
            eff = entry.domain_effective_mask(:) & entry.valid_grid_mask(:);
            plot(ax2, entry.x_grid(eff), zeros(nnz(eff), 1), '.', ...
                'Color', colors(blade_id, :), 'MarkerSize', 5, 'HandleVisibility', 'off');
        end
    end
    xlabel(ax2, 'x in OPRCenterStd frame (mm)');
    ylabel(ax2, 'dV/dx (V/mm)');
    title(ax2, sprintf('CH%d template gradients and trusted-domain boundaries', sid));
    legend(ax2, 'Location', 'eastoutside');

    save_step04_figure_local(fig, fig_dir, cfg, sprintf('Step04_Template_Curves_Clean_CH%d', sid));
    if ~show_plots
        close(fig);
    end
end

% Figure 2: raw x-V point clouds used for template aggregation.  Each tile
% shows the actual low-speed samples retained for one sensor-blade template,
% with the median/smoothed template overlaid.  This is the key diagnostic for
% whether the non-parametric template is supported by the original data.
for iSensor = 1:numel(sensor_ids)
    sid = sensor_ids(iSensor);
    fig_raw = figure('Name', sprintf('Step04 raw template points CH%d', sid), ...
        'Color', 'w', 'Position', [90, 70, 1500, 880], 'NumberTitle', 'off', ...
        'Visible', visibility_state_local(show_plots));
    tiledlayout(fig_raw, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    for blade_id = 1:blades_num
        ax = nexttile;
        hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
        entry = get_template_entry_local(Template, sid, blade_id);
        if isempty(entry) || isempty(entry.x_grid)
            title(ax, sprintf('CH%d B%d: empty', sid, blade_id));
            continue;
        end

        if ~isempty(entry.x_points_preview)
            plot(ax, entry.x_points_preview, entry.v_points_preview, '.', ...
                'Color', [0.72 0.72 0.72], 'MarkerSize', 2, ...
                'DisplayName', 'raw retained points');
        end

        yl_now = ylim(ax);
        if all(isfinite(entry.x_domain))
            shade_domain_local(ax, entry.x_domain, yl_now, colors(blade_id, :));
        end
        if ~isempty(entry.x_points_preview)
            plot(ax, entry.x_points_preview, entry.v_points_preview, '.', ...
                'Color', [0.72 0.72 0.72], 'MarkerSize', 2, ...
                'HandleVisibility', 'off');
        end
        plot(ax, entry.x_grid, entry.v_grid_raw, '-', 'LineWidth', 0.8, ...
            'Color', [0.35 0.35 0.35], 'DisplayName', 'binned median');
        plot(ax, entry.x_grid, entry.v_grid, '-', 'LineWidth', 1.8, ...
            'Color', colors(blade_id, :), 'DisplayName', 'smoothed template');
        if isfinite(entry.baseline)
            yline(ax, entry.baseline, ':', 'Color', [0.25 0.25 0.25], ...
                'HandleVisibility', 'off');
        end
        if all(isfinite(entry.x_domain))
            yl_now = ylim(ax);
            draw_domain_boundaries_local(ax, entry.x_domain, yl_now, colors(blade_id, :));
        end
        xlabel(ax, 'x in OPRCenterStd frame (mm)');
        ylabel(ax, 'Voltage (V)');
        title(ax, sprintf('CH%d B%d | pulses=%d, points=%d, domain=%.2f mm', ...
            sid, blade_id, entry.pulse_count, entry.point_count, diff(entry.x_domain)));
        legend(ax, 'Location', 'best');
    end

    save_step04_figure_local(fig_raw, fig_dir, cfg, sprintf('Step04_RawPointCloud_TemplateFit_CH%d', sid));
    if ~show_plots
        close(fig_raw);
    end
end

% Figure 3: quality summary heatmaps.
S = Template.Summary_Table;
fig2 = figure('Name', 'Step04 template quality summary', 'Color', 'w', ...
    'Position', [100, 100, 1300, 650], 'NumberTitle', 'off', ...
    'Visible', visibility_state_local(show_plots));
tiledlayout(fig2, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
[pulse_mat, width_mat, amp_mat] = summary_to_matrices_local(S, sensor_ids, blades_num);

nexttile;
imagesc(1:blades_num, 1:numel(sensor_ids), pulse_mat);
set(gca, 'YTick', 1:numel(sensor_ids), 'YTickLabel', compose('CH%d', sensor_ids), 'XTick', 1:blades_num);
xlabel('Blade ID'); ylabel('Sensor'); title('Pulse count'); colorbar;

nexttile;
imagesc(1:blades_num, 1:numel(sensor_ids), width_mat);
set(gca, 'YTick', 1:numel(sensor_ids), 'YTickLabel', compose('CH%d', sensor_ids), 'XTick', 1:blades_num);
xlabel('Blade ID'); ylabel('Sensor'); title('Trusted domain width (mm)'); colorbar;

nexttile;
imagesc(1:blades_num, 1:numel(sensor_ids), amp_mat);
set(gca, 'YTick', 1:numel(sensor_ids), 'YTickLabel', compose('CH%d', sensor_ids), 'XTick', 1:blades_num);
xlabel('Blade ID'); ylabel('Sensor'); title('Template amplitude (V)'); colorbar;
save_step04_figure_local(fig2, fig_dir, cfg, 'Step04_Template_Quality_Heatmap');
if ~show_plots
    close(fig2);
end

% Figure 4: summary table.
fig3 = figure('Name', 'Step04 template summary table', 'Color', 'w', ...
    'Position', [120, 120, 1500, 520], 'NumberTitle', 'off', ...
    'Visible', visibility_state_local(show_plots));
show_cols = {'sensor_id', 'blade_id', 'pulse_count', 'point_count', 'domain_width_mm', ...
    'baseline_v', 'amplitude_v', 'max_abs_gradient_v_per_mm', 'quality_status'};
show_cols = show_cols(ismember(show_cols, S.Properties.VariableNames));
uitable('Data', sanitize_table_cells_local(table2cell(S(:, show_cols))), ...
    'ColumnName', show_cols, 'Units', 'normalized', 'Position', [0 0 1 1]);
save_step04_figure_local(fig3, fig_dir, cfg, 'Step04_Template_Summary_Table');
if ~show_plots
    close(fig3);
end
end

function shade_domain_local(ax, x_domain, yl, color)
if any(~isfinite(x_domain)) || diff(x_domain) <= 0 || any(~isfinite(yl)) || diff(yl) <= 0
    return;
end
patch(ax, [x_domain(1), x_domain(2), x_domain(2), x_domain(1)], ...
    [yl(1), yl(1), yl(2), yl(2)], color, ...
    'FaceAlpha', 0.06, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end

function draw_domain_boundaries_local(ax, x_domain, yl, color)
if any(~isfinite(x_domain)) || diff(x_domain) <= 0 || any(~isfinite(yl)) || diff(yl) <= 0
    return;
end
plot(ax, [x_domain(1), x_domain(1)], yl, ':', 'Color', color, ...
    'LineWidth', 0.8, 'HandleVisibility', 'off');
plot(ax, [x_domain(2), x_domain(2)], yl, ':', 'Color', color, ...
    'LineWidth', 0.8, 'HandleVisibility', 'off');
end

function entry = get_template_entry_local(Template, sid, blade_id)
entry = [];
for i = 1:numel(Template.SensorBlade)
    if Template.SensorBlade(i).sensor_id == sid && Template.SensorBlade(i).blade_id == blade_id
        entry = Template.SensorBlade(i);
        return;
    end
end
end


function [pulse_mat, width_mat, amp_mat] = summary_to_matrices_local(S, sensor_ids, blades_num)
pulse_mat = nan(numel(sensor_ids), blades_num);
width_mat = nan(numel(sensor_ids), blades_num);
amp_mat = nan(numel(sensor_ids), blades_num);
for iSensor = 1:numel(sensor_ids)
    sid = sensor_ids(iSensor);
    for blade_id = 1:blades_num
        idx = find(S.sensor_id == sid & S.blade_id == blade_id, 1, 'first');
        if isempty(idx)
            continue;
        end
        pulse_mat(iSensor, blade_id) = S.pulse_count(idx);
        width_mat(iSensor, blade_id) = S.domain_width_mm(idx);
        amp_mat(iSensor, blade_id) = S.amplitude_v(idx);
    end
end
end


function save_step04_figure_local(fig, figure_dir, cfg, tag)
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
    elseif isnumeric(value) && isscalar(value)
        if isfinite(value)
            cells_out{i} = sprintf('%.6g', value);
        else
            cells_out{i} = '';
        end
    end
end
end


function s = get_string_field_local(S, name, default_value)
s = default_value;
if isstruct(S) && isfield(S, name)
    value = S.(name);
    if ischar(value)
        s = value;
    elseif isstring(value) && isscalar(value)
        s = char(value);
    elseif isnumeric(value)
        s = mat2str(value);
    end
end
end


function ang = wrap_to_180_local(ang)
ang = mod(ang + 180, 360) - 180;
end
