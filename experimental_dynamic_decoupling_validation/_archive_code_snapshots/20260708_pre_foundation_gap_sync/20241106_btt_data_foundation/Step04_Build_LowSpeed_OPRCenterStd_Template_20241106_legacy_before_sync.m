clc; clear; close all;

%STEP04_BUILD_LOWSPEED_OPRCENTERSTD_TEMPLATE_20241106 Build low-speed waveform templates.
%
% This Step04 is the bridge between the method-neutral data foundation
% (Step01-Step03) and later direct waveform inversion.
%
% Purpose:
%   1. Load the Step01 low-speed reference: blade fingerprints, internal
%      blade numbering, OPR reference, and standard relative angles.
%   2. Extract OPR centers and probe pulse windows from the low-speed raw
%      BTT waveforms in cfg.low_speed_case.
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

cfg = BTTDataConfig_20241106();
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
step04_query_safe_margin_mm = cfg.step04_query_safe_margin_mm;

step04_min_template_points = cfg.step04_min_template_points;
step04_min_template_pulses = cfg.step04_min_template_pulses;
step04_min_domain_width_mm = cfg.step04_min_domain_width_mm;
step04_plot_max_templates_per_figure = cfg.step04_plot_max_templates_per_figure;
step04_initial_gap_points = cfg.step04_initial_gap_points;
step04_opr_center_level_ratios = cfg.step04_opr_center_level_ratios;

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
cfg.step04_query_safe_margin_mm = step04_query_safe_margin_mm;
cfg.step04_min_template_points = step04_min_template_points;
cfg.step04_min_template_pulses = step04_min_template_pulses;
cfg.step04_min_domain_width_mm = step04_min_domain_width_mm;
cfg.step04_plot_max_templates_per_figure = step04_plot_max_templates_per_figure;
cfg.step04_initial_gap_points = step04_initial_gap_points;
cfg.step04_opr_center_level_ratios = step04_opr_center_level_ratios;

sensor_config_path = fullfile(cfg.step01_output_dir, 'Sensor_Config_20241106.mat');
if ~isfile(sensor_config_path)
    error('Missing Step01 low-speed reference. Please run Step01_Build_LowSpeed_Reference_20241106.m first.');
end

loaded_cfg = load(sensor_config_path, 'Sensor_Config');
Sensor_Config = loaded_cfg.Sensor_Config;

check_step04_reference_local(cfg, Sensor_Config);
std_angles = resolve_standard_angles_local(Sensor_Config, cfg.step04_standard_angle_source);
case_data = extract_low_speed_case_data_from_raw_local(cfg, Sensor_Config);
check_step04_inputs_local(cfg, Sensor_Config, case_data);

fprintf('\n=== Step04: low-speed OPRCenterStd non-parametric template ===\n');
fprintf('Dataset: %s\n', cfg.dataset);
fprintf('Low-speed case: %s\n', cfg.low_speed_case);
fprintf('Sensors: %s\n', mat2str(cfg.sensor_ids));
fprintf('Blades: 1..%d\n', cfg.blades_num);
fprintf('OPR timing reference: %s\n', get_string_field_local(Sensor_Config, 'OPR_Timing_Method', 'unknown'));
fprintf('Standard-angle source: %s\n', cfg.step04_standard_angle_source);

[F_omega_deg_s, speed_diagnostic] = build_low_speed_speed_interpolant_local( ...
    case_data.opr_times(:), cfg.opr_pulses_per_rev);
point_cloud = build_low_speed_point_cloud_local(cfg, Sensor_Config, case_data, std_angles, F_omega_deg_s);
Template = aggregate_template_from_point_cloud_local(cfg, Sensor_Config, case_data, std_angles, point_cloud, speed_diagnostic, sensor_config_path);
point_cloud = annotate_point_cloud_with_template_domain_local(point_cloud, Template);

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
low_speed_feature_path = fullfile(out_dir, sprintf('LowSpeed_Template_SourceData_%s.mat', cfg.dataset));
summary_table = Template.Summary_Table; %#ok<NASGU>
metadata = Template.Metadata; %#ok<NASGU>
LowSpeedTemplateSourceData = case_data; %#ok<NASGU>
save(template_file, 'Template', 'metadata', '-v7.3');
save(low_speed_feature_path, 'LowSpeedTemplateSourceData', 'point_cloud', 'speed_diagnostic', '-v7.3');
writetable(summary_table, fullfile(out_dir, sprintf('Template_OPRCenterStd_Summary_%s.csv', cfg.dataset)));

fprintf('>>> [Step04] Saved template:\n  %s\n', template_file);
fprintf('>>> [Step04] Saved low-speed template source data:\n  %s\n', low_speed_feature_path);
fprintf('>>> [Step04] Template entries: %d, good/usable entries: %d\n', ...
    numel(Template.SensorBlade), nnz(ismember(string(Template.Summary_Table.quality_status), ["good", "usable"])));

if show_plots || cfg.save_figures
    plot_step04_template_local(Template, cfg, show_plots);
end


function cfg = apply_step04_defaults_local(cfg)
% Output paths.
if ~isfield(cfg, 'step04_output_dir') || isempty(cfg.step04_output_dir)
    cfg.step04_output_dir = fullfile(cfg.output_root, 'step04_low_speed_template');
end
if ~isfield(cfg, 'step04_figure_dir') || isempty(cfg.step04_figure_dir)
    cfg.step04_figure_dir = fullfile(cfg.figure_root, 'step04_low_speed_template');
end

% Standard-angle field used as template coordinate reference.
if ~isfield(cfg, 'step04_standard_angle_source') || isempty(cfg.step04_standard_angle_source)
    cfg.step04_standard_angle_source = 'opr_center';  % opr_center, selected, mean, median
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
if ~isfield(cfg, 'step04_query_safe_margin_mm') || isempty(cfg.step04_query_safe_margin_mm)
    cfg.step04_query_safe_margin_mm = 0.20;
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
if ~isfield(cfg, 'step04_initial_gap_points') || isempty(cfg.step04_initial_gap_points)
    cfg.step04_initial_gap_points = 1e4;
end
if ~isfield(cfg, 'step04_opr_center_level_ratios') || isempty(cfg.step04_opr_center_level_ratios)
    cfg.step04_opr_center_level_ratios = [0.30 0.40 0.50 0.60 0.70];
end
end


function check_step04_reference_local(cfg, Sensor_Config)
if ~isfield(Sensor_Config, 'Target_Indices') || ~isa(Sensor_Config.Target_Indices, 'containers.Map')
    error('Step04 requires Sensor_Config.Target_Indices from Step01.');
end
if ~isfield(Sensor_Config, 'Fingerprints') || ~isa(Sensor_Config.Fingerprints, 'containers.Map')
    error('Step04 requires Sensor_Config.Fingerprints from Step01.');
end
if ~isfield(Sensor_Config, 'Standard_Relative_Angles') && ...
        ~isfield(Sensor_Config, 'Standard_Relative_Angles_OPRCenter') && ...
        ~isfield(Sensor_Config, 'Standard_Relative_Angles_Mean') && ...
        ~isfield(Sensor_Config, 'Standard_Relative_Angles_Median')
    error('Sensor_Config does not contain standard relative angles. Re-run Step01.');
end
for sid = cfg.sensor_ids(:).'
    if ~isKey(Sensor_Config.Target_Indices, sid)
        error('Sensor_Config.Target_Indices does not contain CH%d.', sid);
    end
    if ~isKey(Sensor_Config.Fingerprints, sid)
        error('Sensor_Config.Fingerprints does not contain CH%d.', sid);
    end
end
end


function check_step04_inputs_local(cfg, Sensor_Config, case_data)
if ~isfield(case_data, 'opr_times') || numel(case_data.opr_times) <= cfg.opr_pulses_per_rev
    error('Step04 requires low-speed OPR center times extracted from raw low-speed data.');
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
            'Step04 raw low-speed OPR method (%s) differs from Sensor_Config (%s).', ...
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
    case {'opr_center', 'oprcenter'}
        if isfield(Sensor_Config, 'Standard_Relative_Angles_OPRCenter')
            std_angles = Sensor_Config.Standard_Relative_Angles_OPRCenter;
        else
            std_angles = Sensor_Config.Standard_Relative_Angles;
        end
    case 'selected'
        std_angles = Sensor_Config.Standard_Relative_Angles;
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


function case_data = extract_low_speed_case_data_from_raw_local(cfg, Sensor_Config)
case_dir = fullfile(cfg.dataset_root, cfg.low_speed_case);
if ~isfolder(case_dir)
    error('Low-speed raw data folder not found:\n  %s', case_dir);
end

file_ids = list_case_file_ids_local(case_dir, cfg.opr_id);
if isempty(file_ids)
    error('No low-speed OPR raw files found under:\n  %s', case_dir);
end

fprintf('>>> [Step04] Extracting low-speed OPR centers from raw data...\n');
OPRRaw = extract_low_speed_opr_from_raw_local(case_dir, file_ids, cfg);

max_channel_id = max([cfg.sensor_ids(:); cfg.opr_id]);
empty_ch = struct( ...
    'sensor_id', NaN, ...
    'start_times', [], ...
    'end_times', [], ...
    'arrival_times', [], ...
    'pulse_peak_value', [], ...
    'pulse_feature_value', [], ...
    'raw_peak_value', [], ...
    'source_file_id', []);
channels = repmat(empty_ch, max_channel_id, 1);

for sid = cfg.sensor_ids(:).'
    fprintf('>>> [Step04] Extracting low-speed probe pulses CH%d...\n', sid);
    ProbeRaw = extract_low_speed_probe_from_raw_local(case_dir, file_ids, sid, cfg);
    channels(sid).sensor_id = sid;
    channels(sid).start_times = ProbeRaw.start_time_s;
    channels(sid).end_times = ProbeRaw.end_time_s;
    channels(sid).arrival_times = ProbeRaw.arrival_time_s;
    channels(sid).pulse_peak_value = ProbeRaw.pulse_peak_value;
    channels(sid).pulse_feature_value = ProbeRaw.pulse_feature_value;
    channels(sid).raw_peak_value = ProbeRaw.raw_peak_value;
    channels(sid).source_file_id = ProbeRaw.source_file_id;
end

case_data = struct();
case_data.dataset = cfg.dataset;
case_data.case_name = cfg.low_speed_case;
case_data.case_dir = case_dir;
case_data.file_ids = file_ids(:).';
case_data.file_ranges = struct('file_id', {}, 'offset', {});
case_data.opr_times_raw_extracted = OPRRaw.center_time_s(:);
if isfield(Sensor_Config, 'OPR_Center_Times') && ~isempty(Sensor_Config.OPR_Center_Times)
    case_data.opr_times = Sensor_Config.OPR_Center_Times(:);
    case_data.opr_time_source = 'Sensor_Config.OPR_Center_Times';
else
    case_data.opr_times = OPRRaw.center_time_s(:);
    case_data.opr_time_source = 'Step04 raw multi-threshold extraction';
end
case_data.opr_start_times = OPRRaw.start_time_s(:);
case_data.opr_end_times = OPRRaw.end_time_s(:);
case_data.opr_width_s = OPRRaw.width_s(:);
case_data.opr_center_quality = OPRRaw.center_quality(:);
case_data.channels = channels;
case_data.opr_timing_method = 'multi_threshold_center';
case_data.opr_raw_extraction_method = 'multi_threshold_width_center';
case_data.probe_arrival_method = 'half_area_arrival_from_raw_low_speed';
case_data.time_index_mode = get_string_field_local(Sensor_Config, 'Time_Index_Mode', 'raw_sample_index_divided_by_sample_rate');
case_data.time_index_mode_detected = 'raw_sample_index_divided_by_sample_rate';
case_data.created_by = mfilename;
case_data.created_on = datestr(now, 31);
case_data.note = ['Step04 extracted these low-speed timing features directly from raw ', ...
    'files. They are not legacy/NewFlow template outputs.'];
case_data.SourceInventoryTable = inspect_low_speed_source_inventory_local(case_dir, cfg);
end


function T = inspect_low_speed_source_inventory_local(case_dir, cfg)
rows = [];
channels = [cfg.sensor_ids(:); cfg.opr_id];
for channel_id = channels(:).'
    files = dir(fullfile(case_dir, sprintf('4-%d-*.mat', channel_id)));
    rows = [rows; struct( ... %#ok<AGROW>
        'channel_id', channel_id, ...
        'file_count', numel(files), ...
        'first_file', string(first_file_name_local(files)), ...
        'last_file', string(last_file_name_local(files)), ...
        'case_dir', string(case_dir))];
end
T = struct2table(rows);
end


function name = first_file_name_local(files)
if isempty(files)
    name = '';
else
    [~, order] = sort({files.name});
    name = files(order(1)).name;
end
end


function name = last_file_name_local(files)
if isempty(files)
    name = '';
else
    [~, order] = sort({files.name});
    name = files(order(end)).name;
end
end


function file_ids = list_case_file_ids_local(case_dir, channel_id)
files = dir(fullfile(case_dir, sprintf('4-%d-*.mat', channel_id)));
file_ids = nan(numel(files), 1);
for i = 1:numel(files)
    tok = regexp(files(i).name, sprintf('^4-%d-(\\d+)\\.mat$', channel_id), 'tokens', 'once');
    if isempty(tok)
        continue;
    end
    file_ids(i) = str2double(tok{1});
end
file_ids = sort(unique(file_ids(isfinite(file_ids))));
end


function OPRRaw = extract_low_speed_opr_from_raw_local(case_dir, file_ids, cfg)
tail = [];
gap_points = cfg.step04_initial_gap_points;
rows = [];
for iFile = 1:numel(file_ids)
    file_id = file_ids(iFile);
    raw = load_raw_channel_by_id_local(case_dir, cfg.opr_id, file_id);
    if isempty(raw)
        continue;
    end
    if iFile == 1 && cfg.initial_trim_points > 0 && size(raw, 1) > cfg.initial_trim_points
        raw(1:cfg.initial_trim_points, :) = [];
    end
    if ~isempty(tail)
        raw = [tail; raw]; %#ok<AGROW>
    end
    [segments, tail, gap_points] = segment_signal_local(raw, cfg.opr_threshold, gap_points, true);
    rows = [rows; build_low_speed_opr_rows_local(raw, segments, file_id, cfg)]; %#ok<AGROW>
end
if ~isempty(tail)
    [segments, ~, ~] = segment_signal_local(tail, cfg.opr_threshold, gap_points, false);
    rows = [rows; build_low_speed_opr_rows_local(tail, segments, file_ids(end), cfg)]; %#ok<AGROW>
end
if isempty(rows)
    error('No OPR pulses extracted from raw low-speed data.');
end
[~, order] = sort([rows.center_time_s]);
rows = rows(order);
OPRRaw = struct();
OPRRaw.center_time_s = [rows.center_time_s].';
OPRRaw.start_time_s = [rows.start_time_s].';
OPRRaw.end_time_s = [rows.end_time_s].';
OPRRaw.width_s = [rows.width_s].';
OPRRaw.peak_value = [rows.peak_value].';
OPRRaw.file_id = [rows.file_id].';
OPRRaw.center_quality = [rows.center_quality].';
fprintf('>>> [Step04] Low-speed OPR pulses: %d\n', numel(OPRRaw.center_time_s));
end


function rows = build_low_speed_opr_rows_local(raw, segments, file_id, cfg)
rows = [];
for iSeg = 1:numel(segments.start_idx)
    a = segments.start_idx(iSeg);
    b = segments.end_idx(iSeg);
    [center_t, peak_v, q] = compute_low_speed_opr_center_local(raw, a, b, cfg);
    start_t = raw(a, 1) / cfg.pinlv;
    end_t = raw(b, 1) / cfg.pinlv;
    rows = [rows; struct( ... %#ok<AGROW>
        'center_time_s', center_t, ...
        'start_time_s', start_t, ...
        'end_time_s', end_t, ...
        'width_s', end_t - start_t, ...
        'peak_value', peak_v, ...
        'file_id', file_id, ...
        'center_quality', q)];
end
end


function [center_t, peak_v, quality] = compute_low_speed_opr_center_local(raw, a, b, cfg)
expand_pts = max(3, floor((b - a + 1) * 0.4));
a2 = max(1, a - expand_pts);
b2 = min(size(raw, 1), b + expand_pts);
t = raw(a2:b2, 1) / cfg.pinlv;
v = raw(a2:b2, 2);
try
    v_smooth = smoothdata(v, 'movmean', min(21, max(3, 2 * floor(numel(v) / 8) + 1)));
catch
    v_smooth = v;
end

peak_v = max(v_smooth);
baseline = min(v_smooth);
amp = peak_v - baseline;
centers = nan(numel(cfg.step04_opr_center_level_ratios), 1);
if amp > 0 && numel(t) >= 3
    for i = 1:numel(cfg.step04_opr_center_level_ratios)
        level = baseline + cfg.step04_opr_center_level_ratios(i) * amp;
        centers(i) = pulse_width_center_at_level_local(t, v_smooth, level);
    end
end
centers = centers(isfinite(centers));
if isempty(centers)
    center_t = 0.5 * (raw(a, 1) + raw(b, 1)) / cfg.pinlv;
    quality = 0;
else
    center_t = median(centers);
    quality = numel(centers) / numel(cfg.step04_opr_center_level_ratios);
end
end


function center_t = pulse_width_center_at_level_local(t, v, level)
above = find(v >= level);
if isempty(above)
    center_t = NaN;
    return;
end
i1 = above(1);
i2 = above(end);
t_rise = interp_crossing_local(t, v, i1 - 1, i1, level);
t_fall = interp_crossing_local(t, v, i2, i2 + 1, level);
if ~isfinite(t_rise) || ~isfinite(t_fall) || t_fall < t_rise
    center_t = NaN;
else
    center_t = 0.5 * (t_rise + t_fall);
end
end


function tc = interp_crossing_local(t, v, i_left, i_right, level)
if i_left < 1 || i_right > numel(t)
    tc = NaN;
    return;
end
v1 = v(i_left);
v2 = v(i_right);
if abs(v2 - v1) < eps
    tc = 0.5 * (t(i_left) + t(i_right));
else
    alpha = (level - v1) / (v2 - v1);
    alpha = min(max(alpha, 0), 1);
    tc = t(i_left) + alpha * (t(i_right) - t(i_left));
end
end


function ProbeRaw = extract_low_speed_probe_from_raw_local(case_dir, file_ids, sid, cfg)
tail = [];
gap_points = cfg.step04_initial_gap_points;
rows = [];
threshold = cfg.sensor_thresholds(sid);
for iFile = 1:numel(file_ids)
    file_id = file_ids(iFile);
    raw = load_raw_channel_by_id_local(case_dir, sid, file_id);
    if isempty(raw)
        continue;
    end
    if iFile == 1 && cfg.initial_trim_points > 0 && size(raw, 1) > cfg.initial_trim_points
        raw(1:cfg.initial_trim_points, :) = [];
    end
    if ~isempty(tail)
        raw = [tail; raw]; %#ok<AGROW>
    end
    [segments, tail, gap_points] = segment_signal_local(raw, threshold, gap_points, true);
    rows = [rows; build_low_speed_probe_rows_local(raw, segments, file_id, sid, cfg, threshold)]; %#ok<AGROW>
end
if ~isempty(tail)
    [segments, ~, ~] = segment_signal_local(tail, threshold, gap_points, false);
    rows = [rows; build_low_speed_probe_rows_local(tail, segments, file_ids(end), sid, cfg, threshold)]; %#ok<AGROW>
end

if isempty(rows)
    warning('Step04:NoLowSpeedProbePulses', 'No low-speed probe pulses extracted for CH%d.', sid);
    ProbeRaw = struct('sensor_id', sid, 'start_time_s', [], 'end_time_s', [], ...
        'arrival_time_s', [], 'pulse_peak_value', [], 'pulse_feature_value', [], ...
        'raw_peak_value', [], 'source_file_id', []);
    return;
end

[~, order] = sort([rows.arrival_time_s]);
rows = rows(order);
ProbeRaw = struct();
ProbeRaw.sensor_id = sid;
ProbeRaw.start_time_s = [rows.start_time_s].';
ProbeRaw.end_time_s = [rows.end_time_s].';
ProbeRaw.arrival_time_s = [rows.arrival_time_s].';
ProbeRaw.pulse_peak_value = [rows.pulse_peak_value].';
ProbeRaw.pulse_feature_value = [rows.pulse_feature_value].';
ProbeRaw.raw_peak_value = [rows.raw_peak_value].';
ProbeRaw.source_file_id = [rows.source_file_id].';
fprintf('>>> [Step04] CH%d low-speed probe pulses: %d\n', sid, numel(ProbeRaw.arrival_time_s));
end


function rows = build_low_speed_probe_rows_local(raw, segments, file_id, sid, cfg, threshold)
rows = [];
for iSeg = 1:numel(segments.start_idx)
    a = segments.start_idx(iSeg);
    b = segments.end_idx(iSeg);
    [t_start, t_end, t_arrival, feat, feat_raw, feat_fitted] = ...
        compute_half_area_arrival_local(raw, a, b, cfg.pinlv, threshold);
    rows = [rows; struct( ... %#ok<AGROW>
        'sensor_id', sid, ...
        'start_time_s', t_start, ...
        'end_time_s', t_end, ...
        'arrival_time_s', t_arrival, ...
        'pulse_peak_value', feat_raw, ...
        'pulse_feature_value', feat, ...
        'raw_peak_value', feat_fitted, ...
        'source_file_id', file_id)];
end
end


function raw = load_raw_channel_by_id_local(case_dir, channel_id, file_id)
filepath = fullfile(case_dir, sprintf('4-%d-%d.mat', channel_id, file_id));
if ~isfile(filepath)
    raw = [];
    return;
end
loaded = load(filepath);
fn = fieldnames(loaded);
raw = [];
for i = 1:numel(fn)
    value = loaded.(fn{i});
    if isnumeric(value) && ismatrix(value) && size(value, 2) >= 2
        raw = value(:, 1:2);
        break;
    end
end
if isempty(raw)
    return;
end
raw(raw(:, 1) == 0, :) = [];
raw = raw(all(isfinite(raw), 2), :);
end


function [segments, tail, next_gap] = segment_signal_local(raw, threshold, gap_points, keep_last_as_tail)
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

if keep_last_as_tail
    tail_point = chase(seg_starts(end)) - floor(next_gap / 2);
    if tail_point > 0 && tail_point < size(raw, 1)
        tail = raw(tail_point:end, :);
    end
end

if keep_last_as_tail
    if numel(seg_starts) < 2
        return;
    end
    comp_starts = seg_starts(1:end-1);
    comp_ends = seg_ends(1:end-1);
else
    comp_starts = seg_starts;
    comp_ends = seg_ends;
end

segments.start_idx = chase(comp_starts);
segments.end_idx = chase(comp_ends);
segments.start_sample = raw(segments.start_idx, 1);
segments.end_sample = raw(segments.end_idx, 1);
end


function [t_start, t_end, t_arrival, feat, feat_raw, feat_fitted] = compute_half_area_arrival_local(raw, a, b, pinlv, threshold)
expand_pts = floor((b - a + 1) * 0.3);
a = max(1, a - expand_pts);
b = min(size(raw, 1), b + expand_pts);
t = raw(a:b, 1) / pinlv;
v = raw(a:b, 2);
t_start = t(1);
t_end = t(end);
feat_raw = max(v);
feat_fitted = fit_peak_feature_local(t, v, threshold);
if isfinite(feat_fitted)
    feat = feat_fitted;
else
    feat = feat_raw;
end
if numel(t) < 2
    t_arrival = t_start;
    return;
end
dt = median(diff(t));
area = cumsum(max(v, 0) * dt);
target = 0.5 * area(end);
idx = find(area >= target, 1, 'first');
if isempty(idx)
    t_arrival = t_start;
elseif idx == 1
    t_arrival = t(1);
else
    t_arrival = 0.5 * (t(idx - 1) + t(idx));
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
order = min(3, numel(unique(t_fit)) - 1);
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
    tt = linspace(min(t_fit) - mu, max(t_fit) - mu, 100);
    peak_feature = max(polyval(p, tt));
catch
    peak_feature = max(smooth_seg);
end
end


function smooth_v = smooth_signal_for_peak_local(v)
try
    if numel(v) >= 21
        smooth_v = sgolayfilt(v, 3, 21);
    elseif numel(v) >= 5
        span = max(3, 2 * floor(numel(v) / 4) + 1);
        smooth_v = smoothdata(v, 'movmean', span);
    else
        smooth_v = v;
    end
catch
    smooth_v = v;
end
end


function restore_polyfit_warnings_local(w1, w2, w3, w4)
warning(w1);
warning(w2);
warning(w3);
warning(w4);
end


function [F_omega_deg_s, diagnostic] = build_low_speed_speed_interpolant_local(opr_times, opr_pulses_per_rev)
opr_times = opr_times(:);
rev_period = opr_times((opr_pulses_per_rev + 1):end) - opr_times(1:(end - opr_pulses_per_rev));
speed_time = opr_times(1:(end - opr_pulses_per_rev));
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


function xc = estimate_low_speed_template_center_xc_local(x, v)
% Match the legacy NewFlow convention: detect a pulse center in absolute
% OPRCenterStd x, then build the formal template in x - xc coordinates.
x = x(:);
v = v(:);
valid = isfinite(x) & isfinite(v);
x = x(valid);
v = v(valid);
xc = NaN;
if numel(x) < 10
    return;
end
baseline = prctile(v, 20);
v_zero = max(v - baseline, 0);
if ~any(v_zero > 0)
    [~, imax] = max(v);
    xc = x(imax);
    return;
end
high_level = prctile(v_zero, 85);
high_mask = v_zero >= high_level & v_zero > 0;
if nnz(high_mask) >= 5 && sum(v_zero(high_mask), 'omitnan') > eps
    xc = sum(x(high_mask) .* v_zero(high_mask), 'omitnan') ./ ...
        sum(v_zero(high_mask), 'omitnan');
else
    [~, imax] = max(v_zero);
    xc = x(imax);
end
end


function selection = select_legacy_like_template_domain_local(cfg, x_abs, v, threshold)
% Legacy NewFlow Step05 selection:
%   wide cloud -> high-level centroid -> gradient-threshold span -> sgfit xc.
selection = struct( ...
    'ok', false, 'xc', NaN, 'xc_initial', NaN, ...
    'x_domain_rel', [NaN NaN], 'point_count', 0);
x_abs = x_abs(:);
v = v(:);
valid = isfinite(x_abs) & isfinite(v);
x_abs = x_abs(valid);
v = v(valid);
if numel(x_abs) < 50
    return;
end

baseline = estimate_template_baseline_local(x_abs, v, []);
v_zero = max(v - baseline, 0);
high_level = prctile(v_zero, 85);
high_mask = v_zero >= high_level & v_zero > 0;
if nnz(high_mask) >= 5 && sum(v_zero(high_mask), 'omitnan') > eps
    xc0 = sum(x_abs(high_mask) .* v_zero(high_mask), 'omitnan') ./ ...
        sum(v_zero(high_mask), 'omitnan');
else
    [~, imax] = max(v_zero);
    xc0 = x_abs(imax);
end
selection.xc_initial = xc0;

[left_l, right_l, ok_span] = estimate_legacy_gradient_span_local( ...
    cfg, x_abs, v, baseline, xc0);
if ~ok_span
    if isfinite(threshold)
        strict_mask = v >= threshold;
    else
        strict_mask = v_zero >= prctile(v_zero, 60);
    end
    dx_strict = x_abs(strict_mask) - xc0;
    left_dist = -dx_strict(dx_strict < 0);
    right_dist = dx_strict(dx_strict > 0);
    if isempty(left_dist) || isempty(right_dist)
        half_width = prctile(abs(dx_strict), 95);
        left_l = half_width;
        right_l = half_width;
    else
        left_l = prctile(left_dist, 99.5);
        right_l = prctile(right_dist, 99.5);
    end
    left_l = max(left_l - 0.02, cfg.step04_min_half_width_mm);
    right_l = max(right_l - 0.02, cfg.step04_min_half_width_mm);
    left_l = min(left_l, cfg.step04_max_half_width_mm);
    right_l = min(right_l, cfg.step04_max_half_width_mm);
end

trust_mask = x_abs >= xc0 - left_l & x_abs <= xc0 + right_l;
if isfinite(threshold) && nnz(trust_mask) < 30
    trust_mask = v >= threshold;
end
if nnz(trust_mask) < 30
    return;
end

xc = refine_center_by_sg_fit_local(x_abs(trust_mask), v(trust_mask), baseline, xc0);
if ~isfinite(xc)
    xc = xc0;
end
x_rel_selected = x_abs(trust_mask) - xc;
selection.ok = nnz(isfinite(x_rel_selected)) >= 30;
selection.xc = xc;
selection.point_count = nnz(trust_mask);
selection.x_domain_rel = [min(x_rel_selected, [], 'omitnan'), ...
    max(x_rel_selected, [], 'omitnan')];
if any(~isfinite(selection.x_domain_rel)) || selection.x_domain_rel(2) <= selection.x_domain_rel(1)
    selection.ok = false;
end
end


function [left_l, right_l, ok] = estimate_legacy_gradient_span_local(cfg, x, v, baseline, xc)
left_l = NaN;
right_l = NaN;
ok = false;
[x_grid, v_smooth] = build_denoised_static_profile_local(x, v, baseline);
if numel(x_grid) < 20
    return;
end
v_above = max(v_smooth(:), 0);
peak_v = max(v_above, [], 'omitnan');
if ~isfinite(peak_v) || peak_v <= 0
    return;
end
g_abs = abs(gradient(v_smooth(:), x_grid(:)));
g_threshold = estimate_noise_aware_gradient_threshold_local( ...
    g_abs, v_above, cfg.step04_gradient_min_ratio);
effective = g_abs >= g_threshold;
if isfield(cfg, 'step04_amplitude_min_ratio') && cfg.step04_amplitude_min_ratio > 0
    effective = effective & v_above >= cfg.step04_amplitude_min_ratio * peak_v;
end
left_candidates = find(x_grid(:) < xc & effective(:));
right_candidates = find(x_grid(:) > xc & effective(:));
if isempty(left_candidates) || isempty(right_candidates)
    return;
end
left_l = max(xc - x_grid(left_candidates(1)), 0);
right_l = max(x_grid(right_candidates(end)) - xc, 0);
left_l = max(left_l, cfg.step04_min_half_width_mm);
right_l = max(right_l, cfg.step04_min_half_width_mm);
if isfinite(cfg.step04_max_half_width_mm) && cfg.step04_max_half_width_mm > 0
    left_l = min(left_l, cfg.step04_max_half_width_mm);
    right_l = min(right_l, cfg.step04_max_half_width_mm);
end
ok = isfinite(left_l) && isfinite(right_l) && left_l > 0 && right_l > 0;
end


function [x_grid, v_smooth] = build_denoised_static_profile_local(x, v, baseline)
x = x(:);
v = v(:);
valid = isfinite(x) & isfinite(v);
x = x(valid);
v = v(valid);
if numel(x) < 20
    x_grid = [];
    v_smooth = [];
    return;
end
[x_sort, order] = sort(x);
v_sort = v(order) - baseline;
grid_n = 600;
edges = linspace(min(x_sort), max(x_sort), grid_n + 1);
x_grid = 0.5 * (edges(1:end-1) + edges(2:end));
bin_id = discretize(x_sort, edges);
valid_bin = isfinite(bin_id);
v_med = accumarray(bin_id(valid_bin), v_sort(valid_bin), [grid_n, 1], @median, NaN);
valid_grid = isfinite(v_med);
if nnz(valid_grid) < 10
    x_grid = [];
    v_smooth = [];
    return;
end
v_fill = v_med;
missing = ~valid_grid;
if any(missing)
    v_fill(missing) = interp1(x_grid(valid_grid), v_med(valid_grid), ...
        x_grid(missing), 'linear', 'extrap');
end
span = min(41, 2 * floor(numel(v_fill) / 8) + 1);
span = max(span, 5);
try
    v_smooth = smoothdata(v_fill, 'sgolay', span);
catch
    v_smooth = smoothdata(v_fill, 'movmean', span);
end
x_grid = x_grid(:);
v_smooth = v_smooth(:);
end


function g_threshold = estimate_noise_aware_gradient_threshold_local(g_abs, v_above, ratio)
g_abs = g_abs(:);
v_above = v_above(:);
finite = isfinite(g_abs) & isfinite(v_above);
if nnz(finite) < 5
    g_threshold = 0;
    return;
end
noise_mask = v_above <= prctile(v_above(finite), 35);
noise_level = median(g_abs(noise_mask & finite), 'omitnan') + ...
    3 * mad(g_abs(noise_mask & finite), 1);
if ~isfinite(noise_level)
    noise_level = 0;
end
g_max = max(g_abs(finite), [], 'omitnan');
if ~isfinite(g_max)
    g_max = 0;
end
g_threshold = max(noise_level, ratio * g_max);
end


function xc = refine_center_by_sg_fit_local(x, v, baseline, xc_initial)
x = x(:);
v = v(:);
finite = isfinite(x) & isfinite(v);
x = x(finite);
v = v(finite);
xc = xc_initial;
if numel(x) < 50
    return;
end
[peak_val, idx_peak] = max(v);
b0 = max(peak_val - baseline, 0.1);
w0 = max(std(x), 0.2);
n0 = 3.0;
xc0 = xc_initial;
if ~isfinite(xc0)
    xc0 = x(idx_peak);
end
sg_model = @(p, xx) p(1) .* exp(-abs((xx - p(4)) ./ ...
    max(p(2), 1e-6)).^max(p(3), 1e-6)) + baseline;
obj_fun = @(p) sum((v - sg_model(p, x)).^2);
p0 = [b0, w0, n0, xc0];
lb = [0.05, 0.05, 1.2, min(x) - 0.5];
ub = [max(10, 2 * b0 + 0.5), 6.0, 10.0, max(x) + 0.5];
try
    if exist('fmincon', 'file') == 2
        opts = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp');
        p_opt = fmincon(obj_fun, p0, [], [], [], [], lb, ub, [], opts);
    else
        penalty_obj = @(p) obj_fun(min(max(p, lb), ub));
        opts = optimset('Display', 'off', 'MaxIter', 400, 'MaxFunEvals', 1200);
        p_opt = fminsearch(penalty_obj, p0, opts);
        p_opt = min(max(p_opt, lb), ub);
    end
    if isfinite(p_opt(4))
        xc = p_opt(4);
    end
catch
    xc = xc_initial;
end
end


function Template = aggregate_template_from_point_cloud_local(cfg, Sensor_Config, case_data, std_angles, point_cloud, speed_diagnostic, sensor_config_path)
entry_count = numel(point_cloud);
SensorBlade = repmat(make_empty_template_entry_local(), 1, entry_count);
summary_rows = repmat(struct( ...
    'sensor_id', NaN, 'blade_id', NaN, ...
    'pulse_count', NaN, 'point_count', NaN, 'trusted_point_count', NaN, ...
    'xc_mm', NaN, 'xc_detected_mm', NaN, ...
    'x_min_mm', NaN, 'x_max_mm', NaN, ...
    'domain_left_mm', NaN, 'domain_right_mm', NaN, 'domain_width_mm', NaN, ...
    'query_safe_left_mm', NaN, 'query_safe_right_mm', NaN, 'query_safe_width_mm', NaN, ...
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
Template.LowSpeed_Template_Source_Mode = 'raw_low_speed_extraction_inside_step04';
Template.LowSpeed_Raw_Case_Dir = case_data.case_dir;
Template.Sensor_IDs = cfg.sensor_ids;
Template.Blades_Num = cfg.blades_num;
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
    'gradient_min_ratio', cfg.step04_gradient_min_ratio, ...
    'amplitude_min_ratio', cfg.step04_amplitude_min_ratio, ...
    'min_half_width_mm', cfg.step04_min_half_width_mm, ...
    'max_half_width_mm', cfg.step04_max_half_width_mm, ...
    'domain_center_mode', cfg.step04_domain_center_mode, ...
    'query_safe_margin_mm', cfg.step04_query_safe_margin_mm, ...
    'final_template_point_policy', 'trusted_domain_recalibration');
Template.Final_Template_Point_Policy = 'trusted_domain_recalibration';
Template.Source_Point_Policy = ['wide pulse-window point cloud is retained for ', ...
    'shape/domain diagnostics; final x_grid/v_grid are re-aggregated from ', ...
    'raw low-speed points inside x_domain.'];
Template.SensorBlade = SensorBlade;
Template.Summary_Table = summary_table;
Template.Metadata = build_template_metadata_local(Template, cfg);
end


function empty = make_empty_template_entry_local()
empty = struct( ...
    'sensor_id', NaN, 'blade_id', NaN, ...
    'xc', NaN, 'xc_detected', NaN, 'xc_reference_source', '', ...
    'x_coordinate_mode', '', 'x_wide_domain_abs', [NaN NaN], ...
    'x_grid', [], 'v_grid', [], 'v_grid_raw', [], 'v_grid_baseline_removed', [], ...
    'dv_dx', [], 'weight_grid', [], 'count_grid', [], 'valid_grid_mask', [], ...
    'domain_effective_mask', [], 'domain_mask', [], ...
    'wide_x_grid', [], 'wide_v_grid', [], 'wide_v_grid_raw', [], ...
    'wide_dv_dx', [], 'wide_weight_grid', [], 'wide_count_grid', [], ...
    'wide_valid_grid_mask', [], 'wide_domain_effective_mask', [], ...
    'wide_domain_mask', [], ...
    'x_wide_domain', [NaN NaN], 'x_domain', [NaN NaN], ...
    'x_fit_domain', [NaN NaN], 'x_query_safe_domain', [NaN NaN], ...
    'baseline', NaN, 'amplitude', NaN, ...
    'pulse_count', 0, 'point_count', 0, 'trusted_point_count', 0, ...
    'source_point_policy', '', 'final_template_point_policy', '', ...
    'quality_status', '', ...
    'x_points_preview', [], 'v_points_preview', []);
end


function [entry, row] = aggregate_one_sensor_blade_local(cfg, pc, sid, blade_id)
entry = make_empty_template_entry_local();
entry.sensor_id = sid;
entry.blade_id = blade_id;

row = struct( ...
    'sensor_id', sid, 'blade_id', blade_id, ...
    'pulse_count', pc.pulse_count, 'point_count', pc.point_count, ...
    'trusted_point_count', NaN, ...
    'xc_mm', NaN, 'xc_detected_mm', NaN, ...
    'x_min_mm', NaN, 'x_max_mm', NaN, ...
    'domain_left_mm', NaN, 'domain_right_mm', NaN, 'domain_width_mm', NaN, ...
    'query_safe_left_mm', NaN, 'query_safe_right_mm', NaN, 'query_safe_width_mm', NaN, ...
    'baseline_v', NaN, 'amplitude_v', NaN, 'max_abs_gradient_v_per_mm', NaN, ...
    'valid_bin_count', NaN, 'coverage_fraction', NaN, ...
    'quality_status', 'empty');

x_abs = pc.x_mm(:);
v = pc.v(:);
valid = isfinite(x_abs) & isfinite(v);
x_abs = x_abs(valid);
v = v(valid);
entry.point_count = numel(x_abs);
entry.pulse_count = pc.pulse_count;
if isempty(x_abs)
    entry.quality_status = 'empty';
    return;
end

selection_threshold = NaN;
if isfield(cfg, 'sensor_thresholds') && numel(cfg.sensor_thresholds) >= sid
    selection_threshold = cfg.sensor_thresholds(sid);
end
selection = select_legacy_like_template_domain_local( ...
    cfg, x_abs, v, selection_threshold);
xc = selection.xc;
if ~isfinite(xc)
    xc = 0;
end
x = x_abs - xc;

row.xc_mm = xc;
row.xc_detected_mm = xc;
row.x_min_mm = min(x);
row.x_max_mm = max(x);
entry.xc = xc;
entry.xc_detected = selection.xc_initial;
entry.xc_reference_source = 'detected_from_current_template_points';
entry.x_coordinate_mode = 'x_rel = x_opr_center_std_abs - xc';
entry.x_wide_domain_abs = [min(x_abs), max(x_abs)];
entry.x_wide_domain = [row.x_min_mm, row.x_max_mm];
entry.source_point_policy = 'wide pulse-window point cloud retained for domain selection and diagnostics';
entry.final_template_point_policy = 'trusted_domain_recalibration';

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

if selection.ok
    x_domain = selection.x_domain_rel;
    domain_mask = x_grid(:) >= x_domain(1) & x_grid(:) <= x_domain(2) & valid_bin_mask(:);
    effective_mask = domain_mask;
else
    [x_domain, domain_mask, effective_mask] = determine_trusted_domain_local(cfg, x_grid(:), v_smooth(:), dv_dx(:), count_grid(:), valid_bin_mask(:), baseline, amplitude, max_abs_grad);
end
weight_grid = build_template_weight_grid_local(count_grid(:), valid_bin_mask(:), domain_mask(:), dv_dx(:), max_abs_grad);

% Keep the full wide-window template for diagnostics, then publish the
% formal calibration template from trusted-domain raw points only.
entry.wide_x_grid = x_grid(:);
entry.wide_v_grid = v_smooth(:);
entry.wide_v_grid_raw = v_grid_raw(:);
entry.wide_dv_dx = dv_dx(:);
entry.wide_weight_grid = weight_grid(:);
entry.wide_count_grid = count_grid(:);
entry.wide_valid_grid_mask = valid_bin_mask(:);
entry.wide_domain_effective_mask = effective_mask(:);
entry.wide_domain_mask = domain_mask(:);

trusted_point_mask = isfinite(x) & isfinite(v) & ...
    x >= x_domain(1) & x <= x_domain(2);
entry.trusted_point_count = nnz(trusted_point_mask);

[formal_x_grid, formal_v_raw, formal_v_smooth, formal_count_grid, ...
    formal_valid_bin_mask, formal_dv_dx, formal_weight_grid] = ...
    build_formal_template_from_trusted_points_local( ...
    cfg, x(trusted_point_mask), v(trusted_point_mask), x_domain);

if nnz(formal_valid_bin_mask) < 5
    entry.quality_status = 'insufficient_trusted_coverage';
    entry.x_domain = x_domain;
    entry.x_fit_domain = x_domain;
    entry.x_query_safe_domain = shrink_domain_local(x_domain, cfg.step04_query_safe_margin_mm);
else
    entry.x_grid = formal_x_grid(:);
    entry.v_grid = formal_v_smooth(:);
    entry.v_grid_raw = formal_v_raw(:);
    entry.v_grid_baseline_removed = formal_v_smooth(:) - baseline;
    entry.dv_dx = formal_dv_dx(:);
    entry.weight_grid = formal_weight_grid(:);
    entry.count_grid = formal_count_grid(:);
    entry.valid_grid_mask = formal_valid_bin_mask(:);
    entry.domain_effective_mask = true(size(entry.x_grid));
    entry.domain_mask = true(size(entry.x_grid));
    entry.x_domain = x_domain;
    entry.x_fit_domain = x_domain;
    entry.x_query_safe_domain = shrink_domain_local(entry.x_domain, cfg.step04_query_safe_margin_mm);
    entry.quality_status = evaluate_template_quality_local( ...
        cfg, pc.pulse_count, entry.trusted_point_count, ...
        entry.x_domain, formal_valid_bin_mask);
end
entry.baseline = baseline;
entry.amplitude = amplitude;

preview_n = min(numel(x), 4000);
if numel(x) > preview_n
    idx = unique(round(linspace(1, numel(x), preview_n)));
else
    idx = 1:numel(x);
end
entry.x_points_preview = x(idx);
entry.v_points_preview = v(idx);

row.trusted_point_count = entry.trusted_point_count;
row.domain_left_mm = entry.x_domain(1);
row.domain_right_mm = entry.x_domain(2);
row.domain_width_mm = entry.x_domain(2) - entry.x_domain(1);
row.query_safe_left_mm = entry.x_query_safe_domain(1);
row.query_safe_right_mm = entry.x_query_safe_domain(2);
row.query_safe_width_mm = entry.x_query_safe_domain(2) - entry.x_query_safe_domain(1);
row.baseline_v = baseline;
row.amplitude_v = amplitude;
row.max_abs_gradient_v_per_mm = max_abs_grad;
row.valid_bin_count = nnz(entry.valid_grid_mask);
row.coverage_fraction = mean(entry.valid_grid_mask);
row.quality_status = entry.quality_status;
end


function [x_grid, v_grid_raw, v_smooth, count_grid, valid_bin_mask, dv_dx, weight_grid] = ...
    build_formal_template_from_trusted_points_local(cfg, x, v, x_domain)
x = x(:);
v = v(:);
valid = isfinite(x) & isfinite(v);
x = x(valid);
v = v(valid);

x_grid = [];
v_grid_raw = [];
v_smooth = [];
count_grid = [];
valid_bin_mask = false(0, 1);
dv_dx = [];
weight_grid = [];

if numel(x) < cfg.step04_min_template_points || any(~isfinite(x_domain)) || x_domain(2) <= x_domain(1)
    return;
end

edges = x_domain(1):cfg.step04_grid_dx_mm:x_domain(2);
if numel(edges) < 5
    edges = linspace(x_domain(1), x_domain(2), 20);
end
if edges(end) < x_domain(2)
    edges(end + 1) = x_domain(2); %#ok<AGROW>
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
    return;
end

v_interp = v_grid_raw;
missing = ~isfinite(v_interp);
if any(missing)
    v_interp(missing) = interp1(x_grid(valid_bin_mask), v_grid_raw(valid_bin_mask), ...
        x_grid(missing), 'linear', 'extrap');
end

v_smooth = smooth_template_local(v_interp, cfg.step04_smooth_span_bins);
dv_dx = gradient(v_smooth, x_grid);
max_abs_grad = max(abs(dv_dx(valid_bin_mask)), [], 'omitnan');
if ~isfinite(max_abs_grad)
    max_abs_grad = 0;
end

domain_mask = true(size(valid_bin_mask));
weight_grid = build_template_weight_grid_local( ...
    count_grid(:), valid_bin_mask(:), domain_mask(:), dv_dx(:), max_abs_grad);
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


function safe_domain = shrink_domain_local(x_domain, margin_mm)
safe_domain = [NaN NaN];
if any(~isfinite(x_domain)) || x_domain(2) <= x_domain(1)
    return;
end
if ~isfinite(margin_mm) || margin_mm < 0
    margin_mm = 0;
end
safe_domain = [x_domain(1) + margin_mm, x_domain(2) - margin_mm];
if safe_domain(2) <= safe_domain(1)
    center_x = mean(x_domain);
    safe_domain = [center_x, center_x];
end
end


function point_cloud = annotate_point_cloud_with_template_domain_local(point_cloud, Template)
for i = 1:numel(point_cloud)
    sid = point_cloud(i).sensor_id;
    blade_id = point_cloud(i).blade_id;
    tpl = get_template_entry_for_annotation_local(Template, sid, blade_id);
    n = numel(point_cloud(i).x_mm);
    point_cloud(i).is_inside_x_domain = false(n, 1);
    point_cloud(i).is_used_for_final_template = false(n, 1);
    point_cloud(i).x_domain = [NaN NaN];
    point_cloud(i).x_query_safe_domain = [NaN NaN];
    point_cloud(i).xc = NaN;
    point_cloud(i).x_coordinate_mode = '';
    point_cloud(i).source_point_policy = 'wide pulse-window point cloud retained for diagnostics';
    if isempty(tpl) || any(~isfinite(tpl.x_domain))
        continue;
    end
    x = point_cloud(i).x_mm(:);
    xc = 0;
    if isfield(tpl, 'xc') && isscalar(tpl.xc) && isfinite(tpl.xc)
        xc = tpl.xc;
    end
    x_rel = x - xc;
    inside = isfinite(x_rel) & x_rel >= tpl.x_domain(1) & x_rel <= tpl.x_domain(2);
    point_cloud(i).is_inside_x_domain = inside(:);
    point_cloud(i).is_used_for_final_template = inside(:);
    point_cloud(i).x_domain = tpl.x_domain;
    point_cloud(i).xc = xc;
    point_cloud(i).x_coordinate_mode = 'x_rel = x_mm - xc';
    if isfield(tpl, 'x_query_safe_domain')
        point_cloud(i).x_query_safe_domain = tpl.x_query_safe_domain;
    end
end
end


function entry = get_template_entry_for_annotation_local(Template, sid, blade_id)
entry = [];
if ~isfield(Template, 'SensorBlade')
    return;
end
for i = 1:numel(Template.SensorBlade)
    if Template.SensorBlade(i).sensor_id == sid && Template.SensorBlade(i).blade_id == blade_id
        entry = Template.SensorBlade(i);
        return;
    end
end
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
metadata.intended_next_steps = {'Step05 direct-template identification', 'method-specific gap/VP/VARPRO variants'};
metadata.cfg_step04 = Template.Template_Settings;
metadata.point_cloud_settings = Template.PointCloud_Settings;
metadata.source_point_policy = Template.Source_Point_Policy;
metadata.final_template_point_policy = Template.Final_Template_Point_Policy;
metadata.opr_timing_method = Template.OPR_Timing_Method;
metadata.probe_arrival_method = Template.Probe_Arrival_Method;
metadata.blade_id_definition = Template.Blade_ID_Definition;
metadata.sensor_ids = cfg.sensor_ids;
metadata.blades_num = cfg.blades_num;
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
    xlabel(ax1, 'x - xc in OPRCenterStd frame (mm)');
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
    xlabel(ax2, 'x - xc in OPRCenterStd frame (mm)');
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
        xlabel(ax, 'x - xc in OPRCenterStd frame (mm)');
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
