clc; clear; close all;

%STEP04_BUILD_LOWSPEED_OPRCENTERSTD_TEMPLATE_20251222 Build low-speed waveform templates.
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
show_plots = parse_bool_env_local('STEP04_SHOW_PLOTS', true);

cfg = BTTDataConfig_20251222();
cfg = apply_step04_defaults_local(cfg);
cfg.save_figures = parse_bool_env_local('STEP04_SAVE_FIGURES', cfg.save_figures);
cfg.step04_use_cached_source = parse_bool_env_local('STEP04_USE_CACHED_SOURCE', false);
stable_window_plan_env = strtrim(getenv('STEP04_STABLE_WINDOW_PLAN_FILE'));
if ~isempty(stable_window_plan_env)
    cfg.step04_stable_window_plan_file = stable_window_plan_env;
else
    canonical_plan = fullfile(cfg.route_dir, 'input', 'step04_stable_window_plan', ...
        sprintf('Step04_LowSpeed_StaticStableWindowPlan_%s.csv', cfg.dataset));
    if isfile(canonical_plan)
        cfg.step04_stable_window_plan_file = canonical_plan;
    end
end

sensor_config_path = fullfile(cfg.step01_output_dir, 'Sensor_Config_20251222.mat');
if ~isfile(sensor_config_path)
    error('Missing Step01 low-speed reference. Please run Step01_Build_LowSpeed_Reference_20251222.m first.');
end

loaded_cfg = load(sensor_config_path, 'Sensor_Config');
Sensor_Config = loaded_cfg.Sensor_Config;

check_step04_reference_local(cfg, Sensor_Config);
[std_angles, standard_angle_field] = resolve_standard_angles_local(Sensor_Config, cfg.step04_standard_angle_source);
cfg.step04_standard_angle_field = standard_angle_field;
source_cache_path = fullfile(cfg.step04_output_dir, sprintf('LowSpeed_Template_SourceData_%s.mat', cfg.dataset));
[case_data, point_cloud, speed_diagnostic, used_source_cache] = prepare_step04_source_data_local( ...
    cfg, Sensor_Config, std_angles, source_cache_path);
check_step04_inputs_local(cfg, Sensor_Config, case_data);

fprintf('\n=== Step04: low-speed OPRCenterStd non-parametric template ===\n');
fprintf('Dataset: %s\n', cfg.dataset);
fprintf('Low-speed case: %s\n', cfg.low_speed_case);
fprintf('Sensors: %s\n', mat2str(cfg.sensor_ids));
fprintf('Blades: 1..%d\n', cfg.blades_num);
fprintf('OPR timing reference: %s\n', get_string_field_local(Sensor_Config, 'OPR_Timing_Method', 'unknown'));
fprintf('Standard-angle source: %s\n', cfg.step04_standard_angle_source);
fprintf('Low-speed source cache used: %d\n', used_source_cache);
Template = aggregate_template_from_point_cloud_local(cfg, Sensor_Config, case_data, std_angles, point_cloud, speed_diagnostic, sensor_config_path);
validate_low_speed_template_contract_local(Template, cfg);

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
low_speed_feature_path = source_cache_path;
summary_table = Template.Summary_Table; %#ok<NASGU>
metadata = Template.Metadata; %#ok<NASGU>
LowSpeedTemplateSourceData = case_data; %#ok<NASGU>
save(template_file, 'Template', 'metadata', '-v7.3');
save(low_speed_feature_path, 'LowSpeedTemplateSourceData', 'point_cloud', 'speed_diagnostic', '-v7.3');
writetable(summary_table, fullfile(out_dir, sprintf('Template_OPRCenterStd_Summary_%s.csv', cfg.dataset)));
write_step04_diagnostic_tables_local(out_dir, cfg, Template, case_data, speed_diagnostic);

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
    cfg.step04_standard_angle_source = 'oprcenter';  % oprcenter, selected, mean, median
end
if ~isfield(cfg, 'opr_timing_reference') || isempty(cfg.opr_timing_reference)
    cfg.opr_timing_reference = 'rising_edge';
end
if ~isfield(cfg, 'opr_timing_method') || isempty(cfg.opr_timing_method)
    cfg.opr_timing_method = 'threshold_rising_edge';
end

% Pulse-window and point-cloud controls.
if ~isfield(cfg, 'step04_pulse_pad_fraction') || isempty(cfg.step04_pulse_pad_fraction)
    cfg.step04_pulse_pad_fraction = 0.30;
end
if ~isfield(cfg, 'step04_pulse_min_pad_points') || isempty(cfg.step04_pulse_min_pad_points)
    cfg.step04_pulse_min_pad_points = 40;
end
if ~isfield(cfg, 'step04_max_points_per_pulse') || isempty(cfg.step04_max_points_per_pulse)
    cfg.step04_max_points_per_pulse = inf;
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
if ~isfield(cfg, 'step04_stable_window_enable') || isempty(cfg.step04_stable_window_enable)
    cfg.step04_stable_window_enable = true;
end
if ~isfield(cfg, 'step04_stable_window_lap_counts') || isempty(cfg.step04_stable_window_lap_counts)
    cfg.step04_stable_window_lap_counts = [20 30 40 50];
end
if ~isfield(cfg, 'step04_stable_window_scan_step_laps') || isempty(cfg.step04_stable_window_scan_step_laps)
    cfg.step04_stable_window_scan_step_laps = 2;
end
if ~isfield(cfg, 'step04_stable_window_neighbor_span_laps') || isempty(cfg.step04_stable_window_neighbor_span_laps)
    cfg.step04_stable_window_neighbor_span_laps = 1;
end
if ~isfield(cfg, 'step04_stable_window_min_points') || isempty(cfg.step04_stable_window_min_points)
    cfg.step04_stable_window_min_points = 3e4;
end
if ~isfield(cfg, 'step04_stable_window_min_pulses') || isempty(cfg.step04_stable_window_min_pulses)
    cfg.step04_stable_window_min_pulses = 20;
end
if ~isfield(cfg, 'step04_stable_window_selection_policy') || isempty(cfg.step04_stable_window_selection_policy)
    cfg.step04_stable_window_selection_policy = 'earliest_short_stable_static_fit';
end
if ~isfield(cfg, 'step04_stable_window_width_deficit_limit') || isempty(cfg.step04_stable_window_width_deficit_limit)
    cfg.step04_stable_window_width_deficit_limit = 0.20;
end
if ~isfield(cfg, 'step04_stable_window_prefer_shorter') || isempty(cfg.step04_stable_window_prefer_shorter)
    cfg.step04_stable_window_prefer_shorter = true;
end
if ~isfield(cfg, 'step04_stable_window_max_start_lap') || isempty(cfg.step04_stable_window_max_start_lap)
    cfg.step04_stable_window_max_start_lap = 25;
end
if ~isfield(cfg, 'step04_stable_window_start_penalty_weight') || isempty(cfg.step04_stable_window_start_penalty_weight)
    cfg.step04_stable_window_start_penalty_weight = 2.0;
end
if ~isfield(cfg, 'step04_stable_window_early_score_tolerance') || isempty(cfg.step04_stable_window_early_score_tolerance)
    cfg.step04_stable_window_early_score_tolerance = 0.35;
end
if ~isfield(cfg, 'step04_stable_window_static_fit_max_points') || isempty(cfg.step04_stable_window_static_fit_max_points)
    cfg.step04_stable_window_static_fit_max_points = 80000;
end
if ~isfield(cfg, 'step04_stable_window_static_fit_center_max_points') || isempty(cfg.step04_stable_window_static_fit_center_max_points)
    cfg.step04_stable_window_static_fit_center_max_points = 5000;
end
if ~isfield(cfg, 'step04_stable_window_center_gap_limit_mm') || isempty(cfg.step04_stable_window_center_gap_limit_mm)
    cfg.step04_stable_window_center_gap_limit_mm = 0.030;
end
if ~isfield(cfg, 'step04_stable_window_local_xc_mad_limit_mm') || isempty(cfg.step04_stable_window_local_xc_mad_limit_mm)
    cfg.step04_stable_window_local_xc_mad_limit_mm = 0.012;
end
if ~isfield(cfg, 'step04_stable_window_plan_file') || isempty(cfg.step04_stable_window_plan_file)
    cfg.step04_stable_window_plan_file = '';
end
if ~isfield(cfg, 'step04_baseline_mode') || isempty(cfg.step04_baseline_mode)
    cfg.step04_baseline_mode = 'sensor_stream_low_state_histogram';
end
if ~isfield(cfg, 'step04_baseline_pad_points') || isempty(cfg.step04_baseline_pad_points)
    cfg.step04_baseline_pad_points = 4000;
end
if ~isfield(cfg, 'step04_baseline_trim_fraction') || isempty(cfg.step04_baseline_trim_fraction)
    cfg.step04_baseline_trim_fraction = 0.05;
end
if ~isfield(cfg, 'step04_baseline_max_samples_per_sensor') || isempty(cfg.step04_baseline_max_samples_per_sensor)
    cfg.step04_baseline_max_samples_per_sensor = 5e5;
end
if ~isfield(cfg, 'step04_baseline_low_state_threshold_v') || isempty(cfg.step04_baseline_low_state_threshold_v)
    cfg.step04_baseline_low_state_threshold_v = cfg.sensor_threshold_default;
end
if ~isfield(cfg, 'step04_wide_threshold_offset_v') || isempty(cfg.step04_wide_threshold_offset_v)
    cfg.step04_wide_threshold_offset_v = 0.05;
end
if ~isfield(cfg, 'step04_center_seed_peak_fraction') || isempty(cfg.step04_center_seed_peak_fraction)
    cfg.step04_center_seed_peak_fraction = 0.60;
end
if ~isfield(cfg, 'step04_center_seed_min_points') || isempty(cfg.step04_center_seed_min_points)
    cfg.step04_center_seed_min_points = 8;
end
if ~isfield(cfg, 'step04_baseline_histogram_bins') || isempty(cfg.step04_baseline_histogram_bins)
    cfg.step04_baseline_histogram_bins = 240;
end
if ~isfield(cfg, 'step04_baseline_histogram_half_width_v') || isempty(cfg.step04_baseline_histogram_half_width_v)
    cfg.step04_baseline_histogram_half_width_v = 0.03;
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
if ~isfield(cfg, 'step04_stable_window_target_width_mm') || isempty(cfg.step04_stable_window_target_width_mm)
    cfg.step04_stable_window_target_width_mm = max(2 * cfg.step04_min_half_width_mm, 6.0);
end
if ~isfield(cfg, 'step04_domain_center_mode') || isempty(cfg.step04_domain_center_mode)
    cfg.step04_domain_center_mode = 'zero';  % zero or peak
end
if ~isfield(cfg, 'step04_center_mode') || isempty(cfg.step04_center_mode)
    cfg.step04_center_mode = 'sgfit';  % centroid or sgfit
end
if ~isfield(cfg, 'step04_center_seed_source') || isempty(cfg.step04_center_seed_source)
    cfg.step04_center_seed_source = 'stable_window_xc';
end
if ~isfield(cfg, 'step04_center_refine_max_shift_mm') || isempty(cfg.step04_center_refine_max_shift_mm)
    cfg.step04_center_refine_max_shift_mm = 0.08;
end
if ~isfield(cfg, 'step04_xrange_mode') || isempty(cfg.step04_xrange_mode)
    cfg.step04_xrange_mode = 'threshold';  % threshold or energy
end
if ~isfield(cfg, 'step04_template_point_selection_mode') || isempty(cfg.step04_template_point_selection_mode)
    cfg.step04_template_point_selection_mode = 'wide_window';  % wide_window or gradient_gate
end
if ~isfield(cfg, 'step04_domain_min_points_per_bin') || isempty(cfg.step04_domain_min_points_per_bin)
    cfg.step04_domain_min_points_per_bin = 30;
end
if ~isfield(cfg, 'step04_domain_min_laps_per_bin') || isempty(cfg.step04_domain_min_laps_per_bin)
    cfg.step04_domain_min_laps_per_bin = 10;
end
if ~isfield(cfg, 'step04_domain_noise_sigma_factor') || isempty(cfg.step04_domain_noise_sigma_factor)
    cfg.step04_domain_noise_sigma_factor = 3.0;
end
if ~isfield(cfg, 'step04_domain_std_abs_limit_v') || isempty(cfg.step04_domain_std_abs_limit_v)
    cfg.step04_domain_std_abs_limit_v = 0.03;
end
if ~isfield(cfg, 'step04_domain_repeat_ratio') || isempty(cfg.step04_domain_repeat_ratio)
    cfg.step04_domain_repeat_ratio = 0.12;
end
if ~isfield(cfg, 'step04_domain_boundary_margin_mm') || isempty(cfg.step04_domain_boundary_margin_mm)
    cfg.step04_domain_boundary_margin_mm = 0.10;
end
if ~isfield(cfg, 'step04_query_guard_mm') || isempty(cfg.step04_query_guard_mm)
    cfg.step04_query_guard_mm = 0.50;
end
if ~isfield(cfg, 'step04_final_template_domain_buffer_mm') || isempty(cfg.step04_final_template_domain_buffer_mm)
    cfg.step04_final_template_domain_buffer_mm = 0.0;
end
if ~isfield(cfg, 'step04_trust_quantile') || isempty(cfg.step04_trust_quantile)
    cfg.step04_trust_quantile = 0.995;
end
if ~isfield(cfg, 'step04_trust_edge_margin_mm') || isempty(cfg.step04_trust_edge_margin_mm)
    cfg.step04_trust_edge_margin_mm = 0.02;
end
if ~isfield(cfg, 'step04_gradient_energy_quantile') || isempty(cfg.step04_gradient_energy_quantile)
    cfg.step04_gradient_energy_quantile = 0.995;
end
if ~isfield(cfg, 'step04_gradient_edge_margin_mm') || isempty(cfg.step04_gradient_edge_margin_mm)
    cfg.step04_gradient_edge_margin_mm = 0.02;
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
if ~isfield(cfg, 'step04_opr_events_per_revolution_policy') || isempty(cfg.step04_opr_events_per_revolution_policy)
    cfg.step04_opr_events_per_revolution_policy = 'auto';
end
if ~isfield(cfg, 'step04_opr_events_per_revolution') || isempty(cfg.step04_opr_events_per_revolution)
    cfg.step04_opr_events_per_revolution = 'auto';
end
if isnumeric(cfg.step04_opr_events_per_revolution)
    cfg.step04_opr_events_per_revolution_source = 'fixed_known_from_experiment';
else
    cfg.step04_opr_events_per_revolution_source = char(string(cfg.step04_opr_events_per_revolution_policy));
end
if ~isfield(cfg, 'step04_low_speed_nominal_rpm') || isempty(cfg.step04_low_speed_nominal_rpm)
    cfg.step04_low_speed_nominal_rpm = 1000;
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
if ~isfield(case_data, 'opr_times') || numel(case_data.opr_times) <= cfg.blades_num
    error('Step04 requires low-speed OPR center times extracted from raw low-speed data.');
end
if isfield(cfg, 'pinlv') && isfield(cfg, 'sample_rate_hz')
    if abs(double(cfg.pinlv) - double(cfg.sample_rate_hz)) > 1e-9
        error('Step04:SampleRateMismatch', 'cfg.pinlv and cfg.sample_rate_hz mismatch.');
    end
end
if ~isfield(Sensor_Config, 'Target_Indices')
    error('Sensor_Config.Target_Indices is missing. Re-run improved Step01.');
end
if ~isfield(Sensor_Config, 'Standard_Relative_Angles_OPRCenter') && ...
        ~isfield(Sensor_Config, 'Standard_Relative_Angles') && ...
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


function [case_data, point_cloud, speed_diagnostic, used_cache] = prepare_step04_source_data_local( ...
    cfg, Sensor_Config, std_angles, source_cache_path)
used_cache = false;
case_data = struct();

if cfg.step04_use_cached_source && isfile(source_cache_path)
    try
        loaded = load(source_cache_path, 'LowSpeedTemplateSourceData');
        if isfield(loaded, 'LowSpeedTemplateSourceData')
            case_data = loaded.LowSpeedTemplateSourceData;
            expected_signature = build_step04_source_extraction_signature_local(cfg, case_data.file_ids);
            if is_step04_source_signature_match_local(case_data, expected_signature)
                used_cache = true;
                fprintf('>>> [Step04] Loaded low-speed raw/timing source cache:\n  %s\n', source_cache_path);
            else
                case_data = struct();
                fprintf('>>> [Step04] Ignoring stale low-speed raw/timing source cache:\n  %s\n', source_cache_path);
            end
        end
    catch ME
        warning('Step04:SourceCacheLoadFailed', 'Ignoring unreadable Step04 source cache: %s', ME.message);
        used_cache = false;
    end
end

if ~used_cache
    case_data = extract_low_speed_case_data_from_raw_local(cfg, Sensor_Config);
else
    case_data = refresh_cached_case_diagnostics_step04_local(cfg, case_data);
    fprintf('>>> [Step04] Rebuilding x-V point cloud from cached raw/timing source.\n');
end
[F_omega_deg_s, speed_diagnostic] = build_low_speed_speed_interpolant_local(case_data.opr_times(:), cfg);
case_data.MappingSignature = build_step04_mapping_signature_local(cfg, speed_diagnostic);
point_cloud = build_low_speed_point_cloud_local(cfg, Sensor_Config, case_data, std_angles, F_omega_deg_s);
end


function signature = build_step04_source_extraction_signature_local(cfg, file_ids)
signature = struct();
signature.dataset = cfg.dataset;
signature.low_speed_case = cfg.low_speed_case;
signature.file_ids = double(file_ids(:).');
signature.opr_threshold = cfg.opr_threshold;
signature.sensor_thresholds = serialize_sensor_thresholds_local(cfg);
signature.opr_center_level_ratios = cfg.step04_opr_center_level_ratios;
signature.initial_gap_points = cfg.step04_initial_gap_points;
signature.initial_trim_points = cfg.initial_trim_points;
signature.sample_rate_hz = cfg.sample_rate_hz;
signature.sensor_ids = cfg.sensor_ids(:).';
signature.opr_id = cfg.opr_id;
signature.baseline_mode = cfg.step04_baseline_mode;
signature.baseline_pad_points = cfg.step04_baseline_pad_points;
signature.baseline_trim_fraction = cfg.step04_baseline_trim_fraction;
signature.baseline_max_samples_per_sensor = cfg.step04_baseline_max_samples_per_sensor;
signature.baseline_low_state_threshold_v = cfg.step04_baseline_low_state_threshold_v;
signature.baseline_histogram_bins = cfg.step04_baseline_histogram_bins;
signature.baseline_histogram_half_width_v = cfg.step04_baseline_histogram_half_width_v;
signature.opr_timing_reference = get_field_or_default_local( ...
    cfg, 'opr_timing_reference', 'rising_edge');
signature.step04_standard_angle_source = get_field_or_default_local( ...
    cfg, 'step04_standard_angle_source', 'oprcenter');
end


function tf = is_step04_source_signature_match_local(case_data, expected_signature)
tf = isfield(case_data, 'SourceExtractionSignature') && ...
    isequaln(case_data.SourceExtractionSignature, expected_signature);
end


function S = serialize_sensor_thresholds_local(cfg)
sensor_ids = cfg.sensor_ids(:).';
thresholds = nan(size(sensor_ids));
for i = 1:numel(sensor_ids)
    sid = sensor_ids(i);
    if isfield(cfg, 'sensor_thresholds') && isa(cfg.sensor_thresholds, 'containers.Map') && isKey(cfg.sensor_thresholds, sid)
        thresholds(i) = cfg.sensor_thresholds(sid);
    elseif isfield(cfg, 'sensor_threshold_default')
        thresholds(i) = cfg.sensor_threshold_default;
    end
end
S = struct('sensor_ids', sensor_ids, 'thresholds', thresholds);
end


function signature = build_step04_mapping_signature_local(cfg, speed_diagnostic)
signature = struct();
signature.opr_events_per_revolution = get_numeric_field_local( ...
    speed_diagnostic, 'opr_events_per_revolution', NaN);
signature.standard_angle_source = cfg.step04_standard_angle_source;
signature.standard_angle_field = cfg.step04_standard_angle_field;
signature.r_tip_mm = cfg.r_tip_mm;
signature.pulse_pad_fraction = cfg.step04_pulse_pad_fraction;
signature.pulse_min_pad_points = cfg.step04_pulse_min_pad_points;
signature.x_collect_abs_limit_mm = cfg.step04_x_collect_abs_limit_mm;
signature.stable_window_enable = cfg.step04_stable_window_enable;
signature.stable_window_selection_policy = cfg.step04_stable_window_selection_policy;
signature.stable_window_plan_file = cfg.step04_stable_window_plan_file;
signature.center_seed_source = cfg.step04_center_seed_source;
signature.center_mode = cfg.step04_center_mode;
signature.center_refine_max_shift_mm = cfg.step04_center_refine_max_shift_mm;
signature.created_on = datestr(now, 31);
end


function case_data = refresh_cached_case_diagnostics_step04_local(cfg, case_data)
case_dir = fullfile(cfg.dataset_root, cfg.low_speed_case);
if ~isfield(case_data, 'case_dir') || isempty(case_data.case_dir)
    case_data.case_dir = case_dir;
else
    case_dir = case_data.case_dir;
end
if ~isfield(case_data, 'FileTimeContinuityTable') || ~istable(case_data.FileTimeContinuityTable)
    if isfolder(case_dir)
        case_data.FileTimeContinuityTable = inspect_file_time_continuity_local(case_dir, cfg);
    end
end
if ~isfield(case_data, 'SourceInventoryTable') || ~istable(case_data.SourceInventoryTable)
    if isfolder(case_dir)
        case_data.SourceInventoryTable = inspect_low_speed_source_inventory_local(case_dir, cfg);
    end
end
if ~isfield(case_data, 'LowSpeedBaseline') || ~istable(case_data.LowSpeedBaseline) || ...
        ~is_step04_baseline_table_compatible_local(case_data.LowSpeedBaseline, cfg)
    if isfolder(case_dir) && isfield(case_data, 'file_ids')
        case_data.LowSpeedBaseline = estimate_low_speed_sensor_baselines_local( ...
            case_dir, case_data.file_ids, cfg, case_data.FileTimeContinuityTable);
    end
end
end


function [std_angles, source_field] = resolve_standard_angles_local(Sensor_Config, source_name)
source_name = lower(strtrim(source_name));
switch source_name
    case {'selected', 'oprcenter', 'opr_center', 'oprcenterstd', 'opr_center_std'}
        if isfield(Sensor_Config, 'Standard_Relative_Angles_OPRCenter')
            std_angles = Sensor_Config.Standard_Relative_Angles_OPRCenter;
            source_field = 'Standard_Relative_Angles_OPRCenter';
        elseif isfield(Sensor_Config, 'Standard_Relative_Angles')
            warning('Step04:OPRCenterAnglesMissing', ...
                'Standard_Relative_Angles_OPRCenter is missing; falling back to Standard_Relative_Angles.');
            std_angles = Sensor_Config.Standard_Relative_Angles;
            source_field = 'Standard_Relative_Angles';
        else
            error('No usable standard relative angle matrix found.');
        end
    case {'startedge', 'start_edge', 'rising', 'rising_edge', 'threshold_rising_edge'}
        if isfield(Sensor_Config, 'Standard_Relative_Angles_StartEdge')
            std_angles = Sensor_Config.Standard_Relative_Angles_StartEdge;
            source_field = 'Standard_Relative_Angles_StartEdge';
        elseif isfield(Sensor_Config, 'Standard_Relative_Angles')
            warning('Step04:StartEdgeAnglesMissing', ...
                'Standard_Relative_Angles_StartEdge is missing; falling back to Standard_Relative_Angles.');
            std_angles = Sensor_Config.Standard_Relative_Angles;
            source_field = 'Standard_Relative_Angles';
        else
            error('No usable rising-edge standard relative angle matrix found.');
        end
    case 'mean'
        if isfield(Sensor_Config, 'Standard_Relative_Angles_Mean')
            std_angles = Sensor_Config.Standard_Relative_Angles_Mean;
            source_field = 'Standard_Relative_Angles_Mean';
        else
            std_angles = Sensor_Config.Standard_Relative_Angles;
            source_field = 'Standard_Relative_Angles';
        end
    case 'median'
        if isfield(Sensor_Config, 'Standard_Relative_Angles_Median')
            std_angles = Sensor_Config.Standard_Relative_Angles_Median;
            source_field = 'Standard_Relative_Angles_Median';
        else
            std_angles = Sensor_Config.Standard_Relative_Angles;
            source_field = 'Standard_Relative_Angles';
        end
    case 'legacy'
        std_angles = Sensor_Config.Standard_Relative_Angles;
        source_field = 'Standard_Relative_Angles';
    otherwise
        error('Unknown Step04 standard angle source: %s', source_name);
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
file_time_table = inspect_file_time_continuity_local(case_dir, cfg);
OPRRaw = extract_low_speed_opr_from_raw_local(case_dir, file_ids, cfg, file_time_table);

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
    ProbeRaw = extract_low_speed_probe_from_raw_local(case_dir, file_ids, sid, cfg, file_time_table);
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
case_data.file_ranges = build_file_ranges_from_time_table_local(file_time_table, cfg.opr_id);
case_data.FileTimeContinuityTable = file_time_table;
case_data.opr_times = OPRRaw.reference_time_s(:);
case_data.opr_reference_times = OPRRaw.reference_time_s(:);
case_data.opr_start_times = OPRRaw.start_time_s(:);
case_data.opr_end_times = OPRRaw.end_time_s(:);
case_data.opr_center_times = OPRRaw.center_time_s(:);
case_data.opr_width_s = OPRRaw.width_s(:);
case_data.opr_center_quality = OPRRaw.center_quality(:);
case_data.channels = channels;
case_data.opr_timing_method = OPRRaw.reference_method;
case_data.opr_center_method = 'multi_threshold_width_center';
case_data.probe_arrival_method = 'half_area_arrival_from_raw_low_speed';
case_data.time_index_mode = get_string_field_local(Sensor_Config, 'Time_Index_Mode', 'raw_sample_index_divided_by_sample_rate');
case_data.time_index_mode_detected = 'raw_sample_index_divided_by_sample_rate';
case_data.LowSpeedBaseline = estimate_low_speed_sensor_baselines_local(case_dir, file_ids, cfg, file_time_table);
case_data.created_by = mfilename;
case_data.created_on = datestr(now, 31);
case_data.note = ['Step04 extracted these low-speed timing features directly from raw ', ...
    'files. They are not legacy/NewFlow template outputs.'];
case_data.SourceInventoryTable = inspect_low_speed_source_inventory_local(case_dir, cfg);
case_data.SourceExtractionSignature = build_step04_source_extraction_signature_local(cfg, file_ids);
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


function T = inspect_file_time_continuity_local(case_dir, cfg)
rows = [];
channels = [cfg.sensor_ids(:); cfg.opr_id];
for channel_id = channels(:).'
    file_ids = list_case_file_ids_local(case_dir, channel_id);
    prev_corrected_last = NaN;
    for i = 1:numel(file_ids)
        file_id = file_ids(i);
        raw = load_raw_channel_by_id_local(case_dir, channel_id, file_id);
        row = struct('channel_id', channel_id, 'file_id', file_id, ...
            'first_sample', NaN, 'last_sample', NaN, 'first_time_s', NaN, ...
            'last_time_s', NaN, 'n_samples', NaN, 'is_global_continuous', false, ...
            'needs_offset', false, 'computed_offset_s', NaN);
        if isempty(raw)
            rows = [rows; row]; %#ok<AGROW>
            continue;
        end
        first_sample = raw(1, 1);
        last_sample = raw(end, 1);
        row.first_sample = first_sample;
        row.last_sample = last_sample;
        row.first_time_s = first_sample / cfg.sample_rate_hz;
        row.last_time_s = last_sample / cfg.sample_rate_hz;
        row.n_samples = size(raw, 1);
        if i == 1
            row.is_global_continuous = true;
            row.needs_offset = false;
            row.computed_offset_s = 0;
            prev_corrected_last = last_sample;
        else
            row.is_global_continuous = first_sample > prev_corrected_last;
            row.needs_offset = ~row.is_global_continuous;
            if row.needs_offset
                corrected_first = prev_corrected_last + 1;
                offset_samples = corrected_first - first_sample;
                row.computed_offset_s = offset_samples / cfg.sample_rate_hz;
                prev_corrected_last = last_sample + offset_samples;
            else
                row.computed_offset_s = 0;
                prev_corrected_last = last_sample;
            end
        end
        rows = [rows; row]; %#ok<AGROW>
    end
end
if isempty(rows)
    T = table();
else
    T = struct2table(rows);
end
end


function file_ranges = build_file_ranges_from_time_table_local(T, channel_id)
file_ranges = struct('file_id', {}, 'offset', {});
if isempty(T) || ~istable(T)
    return;
end
mask = T.channel_id == channel_id;
Tc = T(mask, :);
file_ranges = repmat(struct('file_id', NaN, 'offset', NaN), height(Tc), 1);
for i = 1:height(Tc)
    file_ranges(i).file_id = Tc.file_id(i);
    file_ranges(i).offset = Tc.computed_offset_s(i);
end
end


function offset_s = get_channel_file_offset_local(T, channel_id, file_id)
offset_s = 0;
if isempty(T) || ~istable(T)
    return;
end
mask = T.channel_id == channel_id & T.file_id == file_id;
if any(mask)
    value = T.computed_offset_s(find(mask, 1, 'first'));
    if isfinite(value)
        offset_s = value;
    end
end
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


function tf = is_step04_baseline_table_compatible_local(T, cfg)
tf = false;
if ~istable(T) || ~all(ismember({'sensor_id','baseline_v','baseline_mode'}, T.Properties.VariableNames))
    return;
end
if height(T) < numel(cfg.sensor_ids)
    return;
end
mode = string(T.baseline_mode);
mode = mode(~ismissing(mode));
if isempty(mode) || any(~strcmpi(mode, string(cfg.step04_baseline_mode)))
    return;
end
tf = all(ismember(cfg.sensor_ids(:), T.sensor_id(:))) && all(isfinite(T.baseline_v));
end


function T = estimate_low_speed_sensor_baselines_local(case_dir, file_ids, cfg, file_time_table)
rows = repmat(struct( ...
    'sensor_id', NaN, ...
    'baseline_v', NaN, ...
    'baseline_mode', "", ...
    'base_ratio', NaN, ...
    'sample_count', NaN, ...
    'base_sample_count', NaN, ...
    'threshold_v', NaN, ...
    'pad_points', NaN), numel(cfg.sensor_ids), 1);

if ~any(strcmpi(string(cfg.step04_baseline_mode), ...
        ["sensor_stream_masked_global", "sensor_stream_low_state_histogram"]))
    for i = 1:numel(cfg.sensor_ids)
        sid = cfg.sensor_ids(i);
        rows(i).sensor_id = sid;
        rows(i).baseline_mode = string(cfg.step04_baseline_mode);
        rows(i).threshold_v = get_sensor_threshold_local(cfg, sid);
        rows(i).pad_points = cfg.step04_baseline_pad_points;
    end
    T = struct2table(rows);
    return;
end

    fprintf('>>> [Step04] Estimating low-speed baselines from raw sensor streams (mode=%s)...\n', ...
        cfg.step04_baseline_mode);
for i = 1:numel(cfg.sensor_ids)
    sid = cfg.sensor_ids(i);
    threshold = get_sensor_threshold_local(cfg, sid);
    diag = estimate_sensor_stream_baseline_from_files_local( ...
        case_dir, file_ids, sid, cfg, file_time_table, threshold);

    rows(i).sensor_id = sid;
    rows(i).baseline_v = diag.baseline;
    rows(i).baseline_mode = string(cfg.step04_baseline_mode);
    rows(i).base_ratio = diag.base_ratio;
    rows(i).sample_count = diag.sample_count;
    rows(i).base_sample_count = diag.base_sample_count;
    rows(i).threshold_v = threshold;
    rows(i).pad_points = cfg.step04_baseline_pad_points;

    fprintf('>>> [Step04] CH%d global baseline: %.6g V (base ratio %.2f%%, samples %d)\n', ...
        sid, rows(i).baseline_v, 100 * rows(i).base_ratio, rows(i).sample_count);
end
T = struct2table(rows);
end


function diag = estimate_sensor_stream_baseline_from_files_local(case_dir, file_ids, sid, cfg, file_time_table, threshold)
base_samples = [];
sample_count = 0;
base_sample_count = 0;
max_samples = cfg.step04_baseline_max_samples_per_sensor;
mode = lower(strtrim(string(cfg.step04_baseline_mode)));
for iFile = 1:numel(file_ids)
    file_id = file_ids(iFile);
    raw = load_raw_channel_by_id_local(case_dir, sid, file_id);
    if isempty(raw)
        continue;
    end
    offset_s = get_channel_file_offset_local(file_time_table, sid, file_id);
    raw(:, 1) = raw(:, 1) + offset_s * cfg.sample_rate_hz;
    if iFile == 1 && cfg.initial_trim_points > 0 && size(raw, 1) > cfg.initial_trim_points
        raw(1:cfg.initial_trim_points, :) = [];
    end
    v_here = raw(:, 2);
    v_here = v_here(isfinite(v_here));
    sample_count = sample_count + numel(v_here);
    if isempty(v_here)
        continue;
    end
    if mode == "sensor_stream_low_state_histogram"
        low_thr = cfg.step04_baseline_low_state_threshold_v;
        if ~isfinite(low_thr)
            low_thr = threshold;
        end
        base_here = v_here(v_here < low_thr);
    else
        base_mask = build_masked_baseline_mask_step04_local( ...
            v_here, threshold, cfg.step04_baseline_pad_points);
        base_here = v_here(base_mask);
    end
    base_sample_count = base_sample_count + numel(base_here);
    if isempty(base_here)
        continue;
    end
    if isfinite(max_samples) && max_samples > 0
        remaining = max(0, round(max_samples) - numel(base_samples));
        if remaining <= 0
            continue;
        end
        if numel(base_here) > remaining
            idx = unique(round(linspace(1, numel(base_here), remaining)));
            base_here = base_here(idx);
        end
    end
    base_samples = [base_samples; base_here(:)]; %#ok<AGROW>
end

if mode == "sensor_stream_low_state_histogram"
    baseline = estimate_low_state_histogram_baseline_step04_local(base_samples, cfg);
else
    trimmed_samples = trim_samples_step04_local(base_samples, cfg.step04_baseline_trim_fraction);
    if isempty(trimmed_samples)
        trimmed_samples = base_samples;
    end
    if isempty(trimmed_samples)
        baseline = NaN;
    else
        baseline = mean(trimmed_samples, 'omitnan');
    end
end
diag = struct();
diag.baseline = baseline;
diag.base_ratio = base_sample_count / max(sample_count, 1);
diag.sample_count = sample_count;
diag.base_sample_count = base_sample_count;
end


function baseline = estimate_low_state_histogram_baseline_step04_local(samples, cfg)
samples = samples(isfinite(samples));
if isempty(samples)
    baseline = NaN;
    return;
end
bin_count = max(40, round(cfg.step04_baseline_histogram_bins));
[counts, edges] = histcounts(samples, bin_count);
if isempty(counts) || ~any(counts > 0)
    baseline = median(samples, 'omitnan');
    return;
end
[~, imax] = max(counts);
centers = 0.5 * (edges(1:end-1) + edges(2:end));
mode_center = centers(imax);
half_width = cfg.step04_baseline_histogram_half_width_v;
near = samples(abs(samples - mode_center) <= half_width);
if numel(near) >= 20
    baseline = mean(trim_samples_step04_local(near, cfg.step04_baseline_trim_fraction), 'omitnan');
else
    baseline = mode_center;
end
end


function is_base = build_masked_baseline_mask_step04_local(v_full, threshold, pad_points)
v_full = v_full(:);
idx_above = find(v_full > threshold);
if isempty(idx_above)
    is_base = true(size(v_full));
else
    jumps = find(diff(idx_above) > 3);
    seg_start = [idx_above(1); idx_above(jumps + 1)];
    seg_end = [idx_above(jumps); idx_above(end)];
    is_base = true(size(v_full));
    pad_points = max(0, round(pad_points));
    for i = 1:numel(seg_start)
        lo = max(1, seg_start(i) - pad_points);
        hi = min(numel(v_full), seg_end(i) + pad_points);
        is_base(lo:hi) = false;
    end
end
end


function samples_trim = trim_samples_step04_local(samples, frac)
samples = samples(isfinite(samples));
if isempty(samples)
    samples_trim = samples;
    return;
end
samples = sort(samples(:));
trim_n = floor(max(0, min(frac, 0.45)) * numel(samples));
lo = 1 + trim_n;
hi = numel(samples) - trim_n;
if lo > hi
    samples_trim = samples;
else
    samples_trim = samples(lo:hi);
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


function OPRRaw = extract_low_speed_opr_from_raw_local(case_dir, file_ids, cfg, file_time_table)
tail = [];
gap_points = cfg.step04_initial_gap_points;
rows = [];
for iFile = 1:numel(file_ids)
    file_id = file_ids(iFile);
    raw = load_raw_channel_by_id_local(case_dir, cfg.opr_id, file_id);
    if isempty(raw)
        continue;
    end
    offset_s = get_channel_file_offset_local(file_time_table, cfg.opr_id, file_id);
    raw(:, 1) = raw(:, 1) + offset_s * cfg.sample_rate_hz;
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
OPRRaw.reference_time_s = select_opr_reference_times_step04_local(OPRRaw, cfg);
OPRRaw.reference_method = resolve_opr_reference_method_step04_local(cfg);
fprintf('>>> [Step04] Low-speed OPR pulses: %d\n', numel(OPRRaw.center_time_s));
fprintf('>>> [Step04] Low-speed OPR timing reference: %s\n', OPRRaw.reference_method);
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

function reference_time_s = select_opr_reference_times_step04_local(OPRRaw, cfg)
ref = lower(strtrim(string(get_field_or_default_local( ...
    cfg, 'opr_timing_reference', 'rising_edge'))));
switch ref
    case {"rising_edge", "start_edge", "start", "threshold_rising_edge"}
        reference_time_s = OPRRaw.start_time_s(:);
    case {"center", "opr_center", "pulse_center", "multi_threshold_width_center"}
        reference_time_s = OPRRaw.center_time_s(:);
    case {"falling_edge", "end_edge", "end"}
        reference_time_s = OPRRaw.end_time_s(:);
    otherwise
        error('Unsupported Step04 OPR timing reference: %s.', ref);
end
end

function method = resolve_opr_reference_method_step04_local(cfg)
ref = lower(strtrim(string(get_field_or_default_local( ...
    cfg, 'opr_timing_reference', 'rising_edge'))));
switch ref
    case {"rising_edge", "start_edge", "start", "threshold_rising_edge"}
        method = 'threshold_rising_edge';
    case {"center", "opr_center", "pulse_center", "multi_threshold_width_center"}
        method = 'multi_threshold_width_center';
    case {"falling_edge", "end_edge", "end"}
        method = 'threshold_falling_edge';
    otherwise
        method = char(ref);
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


function ProbeRaw = extract_low_speed_probe_from_raw_local(case_dir, file_ids, sid, cfg, file_time_table)
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
    offset_s = get_channel_file_offset_local(file_time_table, sid, file_id);
    raw(:, 1) = raw(:, 1) + offset_s * cfg.sample_rate_hz;
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
baseline = estimate_pulse_baseline_local(v);
v_eff = max(v - baseline, 0);
area = cumsum(v_eff * dt);
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


function [F_omega_deg_s, diagnostic] = build_low_speed_speed_interpolant_local(opr_times, cfg)
opr_times = opr_times(:);
[events_per_rev, semantics_table] = resolve_opr_events_per_revolution_step04_local(opr_times, cfg);
rev_period = opr_times((events_per_rev + 1):end) - opr_times(1:(end - events_per_rev));
speed_time = 0.5 .* (opr_times((events_per_rev + 1):end) + opr_times(1:(end - events_per_rev)));
speed_deg_s = 360 ./ max(rev_period, eps);
valid = isfinite(speed_time) & isfinite(speed_deg_s) & speed_deg_s > 0;
if nnz(valid) < 2
    error('Not enough valid OPR periods to build low-speed speed interpolant.');
end
F_omega_deg_s = griddedInterpolant(speed_time(valid), speed_deg_s(valid), 'linear', 'nearest');
diagnostic = struct();
diagnostic.opr_count = numel(opr_times);
diagnostic.opr_events_per_revolution = events_per_rev;
diagnostic.opr_semantics_table = semantics_table;
diagnostic.speed_time_s = speed_time(valid);
diagnostic.speed_deg_s = speed_deg_s(valid);
diagnostic.rpm = speed_deg_s(valid) / 360 * 60;
end


function [events_per_rev, T] = resolve_opr_events_per_revolution_step04_local(opr_times, cfg)
candidates = unique([1, cfg.blades_num], 'stable');
rows = repmat(struct('events_per_revolution', NaN, 'median_rpm', NaN, ...
    'mad_rpm', NaN, 'rpm_error_to_nominal', NaN, 'valid_period_count', NaN, ...
    'selected', false), numel(candidates), 1);
for i = 1:numel(candidates)
    epr = candidates(i);
    rows(i).events_per_revolution = epr;
    if numel(opr_times) <= epr
        rows(i).valid_period_count = 0;
        continue;
    end
    rev_period = opr_times((epr + 1):end) - opr_times(1:(end - epr));
    rpm = 60 ./ max(rev_period, eps);
    valid = isfinite(rpm) & rpm > 0;
    rows(i).valid_period_count = nnz(valid);
    if any(valid)
        rows(i).median_rpm = median(rpm(valid), 'omitnan');
        rows(i).mad_rpm = mad(rpm(valid), 1);
        rows(i).rpm_error_to_nominal = abs(rows(i).median_rpm - cfg.step04_low_speed_nominal_rpm);
    end
end
T = struct2table(rows);

setting = cfg.step04_opr_events_per_revolution;
if isnumeric(setting)
    events_per_rev = double(setting);
else
    setting = lower(strtrim(string(setting)));
    if setting == "auto"
        valid = isfinite(T.rpm_error_to_nominal) & T.valid_period_count >= 2;
        if ~any(valid)
            error('Cannot infer Step04 OPR events per revolution.');
        end
        valid_idx = find(valid);
        [~, best_local] = min(T.rpm_error_to_nominal(valid_idx));
        events_per_rev = T.events_per_revolution(valid_idx(best_local));
    else
        events_per_rev = str2double(setting);
        if ~isfinite(events_per_rev)
            error('cfg.step04_opr_events_per_revolution must be numeric or auto.');
        end
    end
end
T.selected = T.events_per_revolution == events_per_rev;
fprintf('>>> [Step04] OPR semantics: selected %g event/rev, median low-speed RPM %.3f\n', ...
    events_per_rev, T.median_rpm(T.selected));
end


function point_cloud = build_low_speed_point_cloud_local(cfg, Sensor_Config, case_data, std_angles, F_omega_deg_s)
case_dir = fullfile(cfg.dataset_root, cfg.low_speed_case);
if ~isfolder(case_dir)
    error('Low-speed raw data folder not found: %s', case_dir);
end

file_ids = case_data.file_ids(:).';
opr_times = case_data.opr_times(:);
if isfield(case_data, 'FileTimeContinuityTable') && istable(case_data.FileTimeContinuityTable)
    file_time_table = case_data.FileTimeContinuityTable;
else
    warning('Step04:MissingFileTimeContinuityTable', ...
        'FileTimeContinuityTable is missing; sensor raw offsets default to zero.');
    file_time_table = table();
end

max_entries = numel(cfg.sensor_ids) * cfg.blades_num;
entry_template = struct( ...
    'sensor_id', NaN, 'blade_id', NaN, ...
    'x_mm', [], 'v', [], 't_s', [], 'pulse_index', [], 'lap_index', [], ...
    'baseline_by_pulse', [], 'feature_by_pulse', [], ...
    'sensor_global_baseline', NaN, 'sensor_baseline_mode', '', ...
    'pulse_count', 0, 'point_count', 0, ...
    'stable_window_enabled', false, 'stable_window_start_lap', NaN, ...
    'stable_window_end_lap', NaN, 'stable_window_lap_count', NaN, ...
    'stable_window_score', NaN, 'stable_window_candidate_count', 0, ...
    'stable_window_reason', '', 'stable_window_candidates', table());
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
    [sensor_baseline, sensor_baseline_mode] = get_low_speed_sensor_baseline_local(case_data, sid);
    for blade_id = 1:cfg.blades_num
        entry_idx = sensor_blade_entry_index_local(cfg, sid, blade_id);
        point_cloud(entry_idx).sensor_global_baseline = sensor_baseline;
        point_cloud(entry_idx).sensor_baseline_mode = sensor_baseline_mode;
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
        offset = get_channel_file_offset_local(file_time_table, sid, file_id);
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
            arrival_t = pulse_times(pulse_idx, 3);
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
            if isfinite(cfg.step04_max_points_per_pulse) && numel(idx) > cfg.step04_max_points_per_pulse
                idx = idx(unique(round(linspace(1, numel(idx), cfg.step04_max_points_per_pulse))));
            end

            t_seg = t_global(idx);
            v_seg = v(idx);
            [x_mm, keep_mask] = map_time_to_x_local(t_seg, opr_times, F_omega_deg_s, ...
                std_angles(sid, blade_id), cfg, arrival_t);
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

            baseline = estimate_point_cloud_baseline_local(cfg, v_seg, sensor_baseline);
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

if cfg.step04_stable_window_enable
    point_cloud = apply_stable_window_selection_to_point_cloud_local(cfg, point_cloud);
end
end


function point_cloud = apply_stable_window_selection_to_point_cloud_local(cfg, point_cloud)
fprintf('>>> [Step04] Selecting stable low-speed template windows independently from extracted point clouds...\n');
external_plan = load_external_stable_window_plan_step04_local(cfg);
for i = 1:numel(point_cloud)
    pc = point_cloud(i);
    if pc.pulse_count < cfg.step04_stable_window_min_pulses || isempty(pc.x_mm)
        point_cloud(i).stable_window_reason = 'insufficient_pulses_for_window_scan';
        continue;
    end

    [external_start_lap, external_end_lap, external_reason] = ...
        resolve_external_stable_window_for_point_cloud_local(external_plan, pc);
    if isfinite(external_start_lap) && isfinite(external_end_lap)
        point_cloud(i) = subset_point_cloud_by_lap_window_step04_local( ...
            pc, external_start_lap, external_end_lap);
        point_cloud(i).stable_window_enabled = true;
        point_cloud(i).stable_window_start_lap = external_start_lap;
        point_cloud(i).stable_window_end_lap = external_end_lap;
        point_cloud(i).stable_window_lap_count = external_end_lap - external_start_lap + 1;
        point_cloud(i).stable_window_score = NaN;
        point_cloud(i).stable_window_candidate_count = 0;
        point_cloud(i).stable_window_reason = char(external_reason);
        point_cloud(i).stable_window_candidates = empty_stable_window_candidate_table_step04_local();

        fprintf('>>> [Step04] CH%d B%d stable window from external plan: laps %d-%d (%d pulses), reason=%s\n', ...
            pc.sensor_id, pc.blade_id, external_start_lap, external_end_lap, ...
            point_cloud(i).pulse_count, char(external_reason));
        continue;
    end

    candidates = build_stable_window_candidates_step04_local(cfg, pc);
    point_cloud(i).stable_window_candidate_count = height(candidates);
    point_cloud(i).stable_window_candidates = candidates;

    if isempty(candidates)
        point_cloud(i).stable_window_reason = 'no_candidate_window';
        continue;
    end

    candidates = score_stable_window_candidates_step04_local(cfg, candidates);
    point_cloud(i).stable_window_candidates = candidates;
    chosen = choose_stable_window_candidate_step04_local(candidates, cfg);
    if isempty(chosen) || ~isfinite(chosen.StartLap)
        point_cloud(i).stable_window_reason = 'no_finite_candidate_score';
        continue;
    end

    point_cloud(i) = subset_point_cloud_by_lap_window_step04_local(pc, chosen.StartLap, chosen.EndLap);
    point_cloud(i).stable_window_enabled = true;
    point_cloud(i).stable_window_start_lap = chosen.StartLap;
    point_cloud(i).stable_window_end_lap = chosen.EndLap;
    point_cloud(i).stable_window_lap_count = chosen.LapCount;
    point_cloud(i).stable_window_score = chosen.ChosenScore;
    point_cloud(i).stable_window_candidate_count = height(candidates);
    point_cloud(i).stable_window_reason = char(chosen.Reason);
    point_cloud(i).stable_window_candidates = candidates;

    fprintf('>>> [Step04] CH%d B%d stable window: laps %d-%d (%d pulses), score %.3g, reason=%s\n', ...
        pc.sensor_id, pc.blade_id, chosen.StartLap, chosen.EndLap, point_cloud(i).pulse_count, ...
        chosen.ChosenScore, char(chosen.Reason));
end
end


function plan = load_external_stable_window_plan_step04_local(cfg)
plan = table();
if ~isfield(cfg, 'step04_stable_window_plan_file') || ...
        isempty(cfg.step04_stable_window_plan_file)
    return;
end
file = char(cfg.step04_stable_window_plan_file);
if isempty(strtrim(file))
    return;
end
if ~isfile(file)
    error('Step04 external stable-window plan file does not exist: %s', file);
end
[~, ~, ext] = fileparts(file);
switch lower(ext)
    case '.mat'
        loaded = load(file);
        if isfield(loaded, 'window_plan')
            plan = loaded.window_plan;
        elseif isfield(loaded, 'FixedWindowPlan')
            plan = loaded.FixedWindowPlan;
        elseif isfield(loaded, 'plan')
            plan = loaded.plan;
        else
            names = fieldnames(loaded);
            for i = 1:numel(names)
                if istable(loaded.(names{i}))
                    plan = loaded.(names{i});
                    break;
                end
            end
        end
    otherwise
        plan = readtable(file);
end
if isempty(plan)
    error('Step04 external stable-window plan is empty or unreadable: %s', file);
end
required = {'BladeID','SensorID','LapCount'};
for i = 1:numel(required)
    if ~ismember(required{i}, plan.Properties.VariableNames)
        error('Step04 external stable-window plan missing required column: %s', required{i});
    end
end
if ~ismember('StartLap', plan.Properties.VariableNames) && ...
        ~ismember('StartOffsetLap', plan.Properties.VariableNames)
    error('Step04 external stable-window plan must contain StartLap or StartOffsetLap.');
end
fprintf('>>> [Step04] External stable-window plan enabled: %s\n', file);
end


function [start_lap, end_lap, reason] = resolve_external_stable_window_for_point_cloud_local(plan, pc)
start_lap = NaN;
end_lap = NaN;
reason = "external_stable_window_plan";
if isempty(plan)
    return;
end
mask = plan.BladeID == pc.blade_id & plan.SensorID == pc.sensor_id;
if ~any(mask)
    return;
end
row = plan(find(mask, 1, 'first'), :);
if ismember('StartLap', row.Properties.VariableNames)
    start_lap = row.StartLap;
else
    start_lap = row.StartOffsetLap + 1;
end
lap_count = row.LapCount;
end_lap = start_lap + lap_count - 1;
if ismember('SelectionRule', row.Properties.VariableNames)
    reason = "external_plan_" + string(row.SelectionRule);
end
if ~isfinite(start_lap) || ~isfinite(end_lap) || end_lap < start_lap
    start_lap = NaN;
    end_lap = NaN;
end
end


function candidates = build_stable_window_candidates_step04_local(cfg, pc)
candidates = empty_stable_window_candidate_table_step04_local();
laps = unique(pc.lap_index(isfinite(pc.lap_index)));
laps = laps(:);
if isempty(laps)
    return;
end
min_lap = min(laps);
max_lap = max(laps);
scan_step = max(1, round(cfg.step04_stable_window_scan_step_laps));
row_count = 0;
rows = repmat(make_empty_stable_window_candidate_row_step04_local(), 0, 1);
pulse_stats = build_pulse_stats_for_stable_window_step04_local(cfg, pc);
if isempty(pulse_stats)
    return;
end

for lap_count_raw = cfg.step04_stable_window_lap_counts(:).'
    lap_count = round(lap_count_raw);
    if lap_count < cfg.step04_stable_window_min_pulses || max_lap - min_lap + 1 < lap_count
        continue;
    end
    max_start = max_lap - lap_count + 1;
    start_laps = min_lap:scan_step:max_start;
    if isempty(start_laps) || start_laps(end) ~= max_start
        start_laps = [start_laps, max_start]; %#ok<AGROW>
    end

    for start_lap = start_laps
        end_lap = start_lap + lap_count - 1;
        stat_mask = pulse_stats.lap_index >= start_lap & pulse_stats.lap_index <= end_lap;
        row = make_empty_stable_window_candidate_row_step04_local();
        row.sensor_id = pc.sensor_id;
        row.blade_id = pc.blade_id;
        row.StartLap = start_lap;
        row.EndLap = end_lap;
        row.LapCount = lap_count;
        row.PulseCount = nnz(stat_mask);
        row.WidePointCount = sum(pulse_stats.point_count(stat_mask), 'omitnan');
        row.Reason = "candidate";

        if row.PulseCount < cfg.step04_stable_window_min_pulses || row.WidePointCount < cfg.step04_stable_window_min_points
            row.Reason = "insufficient_window_points";
            row_count = row_count + 1;
            rows(row_count, 1) = row; %#ok<AGROW>
            continue;
        end

        valid_center = stat_mask & isfinite(pulse_stats.xc_mm);
        valid_feature = stat_mask & isfinite(pulse_stats.feature_v);
        if nnz(valid_center) < cfg.step04_stable_window_min_pulses
            row.Reason = "insufficient_valid_pulse_centers";
            row_count = row_count + 1;
            rows(row_count, 1) = row; %#ok<AGROW>
            continue;
        end

        point_mask = pc.lap_index >= start_lap & pc.lap_index <= end_lap & ...
            isfinite(pc.x_mm) & isfinite(pc.v);
        x_win = pc.x_mm(point_mask);
        v_win = pc.v(point_mask);
        pulse_win = pc.pulse_index(point_mask);
        if numel(x_win) > cfg.step04_stable_window_static_fit_max_points
            decim = ceil(numel(x_win) / cfg.step04_stable_window_static_fit_max_points);
            x_quality = x_win(1:decim:end);
            v_quality = v_win(1:decim:end);
            pulse_quality = pulse_win(1:decim:end);
        else
            x_quality = x_win;
            v_quality = v_win;
            pulse_quality = pulse_win;
        end
        row.SelectedPointCount = numel(x_win);
        row.xc_seed_mm = median(pulse_stats.xc_mm(valid_center), 'omitnan');
        [row.xc_mm, row.center_gap_mm, row.RMSE, row.BoundRatio, ...
            row.LocalXcMAD, row.LeftGradientRatio, row.RightGradientRatio] = ...
            estimate_stable_window_static_quality_step04_local(cfg, pc, x_quality, v_quality, pulse_quality, row.xc_seed_mm);
        row.PulseXcMAD = row.LocalXcMAD;
        if ~isfinite(row.xc_mm)
            row.xc_mm = row.xc_seed_mm;
        end
        row.SelectedWidthMm = prctile(x_win, 95) - prctile(x_win, 5);
        if ~isfinite(row.SelectedWidthMm)
            row.SelectedWidthMm = median(pulse_stats.support_width_mm(valid_center), 'omitnan');
        end
        row.WidthDeficit = max(0, cfg.step04_stable_window_target_width_mm - row.SelectedWidthMm) ./ ...
            max(cfg.step04_stable_window_target_width_mm, eps);
        row.PointDeficit = max(0, cfg.step04_stable_window_min_points - row.SelectedPointCount) ./ ...
            max(cfg.step04_stable_window_min_points, eps);
        row.Reason = "static_fit_scored_candidate";
        row_count = row_count + 1;
        rows(row_count, 1) = row; %#ok<AGROW>
    end
end

if row_count > 0
    candidates = struct2table(rows(1:row_count));
    candidates.LocalXcMAD = compute_candidate_local_xc_mad_step04_local( ...
        candidates.xc_mm, candidates.StartLap, candidates.LapCount, ...
        cfg.step04_stable_window_neighbor_span_laps);
end
end


function pulse_stats = build_pulse_stats_for_stable_window_step04_local(cfg, pc)
pulse_stats = table();
ids = unique(pc.pulse_index(:), 'stable');
if isempty(ids)
    return;
end
threshold = get_sensor_threshold_local(cfg, pc.sensor_id);
sensor_baseline = NaN;
if isfield(pc, 'sensor_global_baseline')
    sensor_baseline = pc.sensor_global_baseline;
end
n = numel(ids);
pulse_id = nan(n, 1);
lap_index = nan(n, 1);
point_count = zeros(n, 1);
xc_mm = nan(n, 1);
support_width_mm = nan(n, 1);
baseline_v = nan(n, 1);
feature_v = nan(n, 1);

for i = 1:n
    pid = ids(i);
    mask = pc.pulse_index == pid;
    x = pc.x_mm(mask);
    v = pc.v(mask);
    finite = isfinite(x) & isfinite(v);
    x = x(finite);
    v = v(finite);
    pulse_id(i) = pid;
    point_count(i) = numel(x);
    lap_values = pc.lap_index(mask);
    lap_index(i) = median(lap_values(isfinite(lap_values)), 'omitnan');
    if numel(x) < 10
        continue;
    end
    baseline = estimate_point_cloud_baseline_local(cfg, v, sensor_baseline);
    baseline_v(i) = baseline;
    vz = max(v - baseline, 0);
    feature_v(i) = max(v, [], 'omitnan') - baseline;
    wide = isolate_main_pulse_step04_local(v, baseline + cfg.step04_wide_threshold_offset_v);
    if nnz(wide) < 5
        wide = isolate_main_pulse_step04_local(v, threshold);
    end
    if nnz(wide) < 5 && max(vz) > 0
        wide = vz >= prctile(vz, 75);
    end
    if nnz(wide) >= 3
        xc_mm(i) = estimate_center_seed_from_points_step04_local(x(wide), vz(wide), cfg);
        support_width_mm(i) = prctile(x(wide), 95) - prctile(x(wide), 5);
    end
end

pulse_stats = table(pulse_id, lap_index, point_count, xc_mm, support_width_mm, ...
    baseline_v, feature_v);
end


function [xc_fit, center_gap, rmse, bound_ratio, local_xc_mad, left_ratio, right_ratio] = ...
    estimate_stable_window_static_quality_step04_local(cfg, pc, x_win, v_win, pulse_win, xc_seed)
xc_fit = NaN;
center_gap = NaN;
rmse = NaN;
bound_ratio = NaN;
local_xc_mad = NaN;
left_ratio = NaN;
right_ratio = NaN;
if isempty(x_win) || isempty(v_win) || ~isfinite(xc_seed)
    return;
end
sensor_baseline = NaN;
if isfield(pc, 'sensor_global_baseline')
    sensor_baseline = pc.sensor_global_baseline;
end
baseline = estimate_point_cloud_baseline_local(cfg, v_win, sensor_baseline);
if ~isfinite(baseline)
    baseline = estimate_template_selection_baseline_local(v_win);
end
threshold = get_sensor_threshold_local(cfg, pc.sensor_id);
strict_mask = isolate_main_pulse_by_pulse_step04_local(v_win, pulse_win, threshold);
if nnz(strict_mask) < 30
    vz = max(v_win - baseline, 0);
    strict_mask = vz >= prctile(vz, 60);
end
center_x = x_win(strict_mask);
center_v = v_win(strict_mask);
if numel(center_x) < 30
    center_x = x_win;
    center_v = v_win;
end
max_center_points = max(1000, round(get_field_or_default_local( ...
    cfg, 'step04_stable_window_static_fit_center_max_points', 5000)));
if numel(center_x) > max_center_points
    [center_x, center_v] = downsample_static_fit_cloud_step04_local(center_x, center_v, max_center_points);
end
if nnz(strict_mask) >= 30
    xc_fit = refine_center_by_sg_fit_step04_local(center_x, center_v, baseline, xc_seed);
else
    xc_fit = refine_center_by_sg_fit_step04_local(center_x, center_v, baseline, xc_seed);
end
if ~isfinite(xc_fit)
    xc_fit = xc_seed;
end
center_gap = abs(xc_fit - xc_seed);
weight = build_gradient_weight_from_wide_cloud_step04_local(x_win, x_win, v_win, baseline, 0.05);
[rmse, bound_ratio] = estimate_static_fit_quality_step04_local(x_win - xc_fit, v_win, weight);
local_xc_mad = estimate_local_center_mad_step04_local(x_win, v_win, pulse_win, baseline, threshold);
[left_ratio, right_ratio] = estimate_gradient_side_balance_step04_local(x_win, weight, xc_fit);
end


function [x_ds, v_ds] = downsample_static_fit_cloud_step04_local(x, v, max_points)
x = x(:);
v = v(:);
finite = isfinite(x) & isfinite(v);
x = x(finite);
v = v(finite);
if numel(x) <= max_points
    x_ds = x;
    v_ds = v;
    return;
end
[x_sort, idx] = sort(x);
v_sort = v(idx);
bin_count = min(max_points, max(200, round(max_points / 2)));
edges = linspace(min(x_sort), max(x_sort), bin_count + 1).';
bin_id = discretize(x_sort, edges);
valid = ~isnan(bin_id);
x_grid = 0.5 * (edges(1:end-1) + edges(2:end));
v_med = accumarray(bin_id(valid), v_sort(valid), [bin_count, 1], @median, NaN);
ok = isfinite(v_med);
x_ds = x_grid(ok);
v_ds = v_med(ok);
if numel(x_ds) < 50
    step = ceil(numel(x_sort) / max_points);
    x_ds = x_sort(1:step:end);
    v_ds = v_sort(1:step:end);
end
end


function local_mad = compute_candidate_local_xc_mad_step04_local(xc_vals, start_laps, lap_counts, neighbor_span)
xc_vals = xc_vals(:);
start_laps = start_laps(:);
lap_counts = lap_counts(:);
n = numel(xc_vals);
neighbor_span = max(1, round(neighbor_span));
local_mad = nan(n, 1);
lap_unique = unique(lap_counts(isfinite(lap_counts)));
lap_tol = 0;
if numel(lap_unique) >= 2
    lap_diff = diff(sort(lap_unique));
    lap_tol = max(1, round(median(lap_diff, 'omitnan')));
end
for k = 1:n
    idx = isfinite(xc_vals) & abs(start_laps - start_laps(k)) <= neighbor_span & ...
        abs(lap_counts - lap_counts(k)) <= lap_tol;
    if nnz(idx) < 3
        idx = isfinite(xc_vals) & abs(start_laps - start_laps(k)) <= neighbor_span;
    end
    if nnz(idx) < 3
        idx = isfinite(xc_vals);
    end
    local_block = xc_vals(idx);
    local_center = median(local_block, 'omitnan');
    local_mad(k) = median(abs(local_block - local_center), 'omitnan');
end
end


function mask = isolate_main_pulse_by_pulse_step04_local(v, pulse_id, threshold)
mask = false(size(v));
finite = isfinite(v) & isfinite(pulse_id);
ids = unique(pulse_id(finite), 'stable');
for i = 1:numel(ids)
    id_mask = pulse_id == ids(i);
    local_mask = isolate_main_pulse_step04_local(v(id_mask), threshold);
    idx = find(id_mask);
    mask(idx(local_mask)) = true;
end
end


function mask = isolate_main_pulse_step04_local(v, threshold)
mask = false(size(v));
idx_above = find(v > threshold);
if isempty(idx_above)
    return;
end
jumps = find(diff(idx_above) > 3);
starts = [1; jumps + 1];
ends = [jumps; numel(idx_above)];
[~, idx_peak] = max(v);
best_seg = 1;
for i = 1:numel(starts)
    if idx_peak >= idx_above(starts(i)) && idx_peak <= idx_above(ends(i))
        best_seg = i;
        break;
    end
end
mask(idx_above(starts(best_seg)):idx_above(ends(best_seg))) = true;
end


function xc = estimate_center_seed_from_points_step04_local(x, v_zero, cfg)
x = x(:);
v_zero = max(v_zero(:), 0);
finite = isfinite(x) & isfinite(v_zero);
x = x(finite);
v_zero = v_zero(finite);
if isempty(x)
    xc = NaN;
    return;
end
[peak_val, idx_peak] = max(v_zero);
if ~isfinite(peak_val) || peak_val <= 0
    xc = x(idx_peak);
    return;
end
mask = v_zero >= cfg.step04_center_seed_peak_fraction * peak_val;
if nnz(mask) < cfg.step04_center_seed_min_points || sum(v_zero(mask), 'omitnan') <= eps
    xc = x(idx_peak);
else
    xc = sum(x(mask) .* v_zero(mask), 'omitnan') ./ sum(v_zero(mask), 'omitnan');
end
end


function candidates = score_stable_window_candidates_step04_local(cfg, candidates)
if isempty(candidates)
    return;
end
rmse_ref = max(median(candidates.RMSE, 'omitnan'), 0.02);
xc_ref = max(median(candidates.LocalXcMAD, 'omitnan'), 0.002);
gap_ref = max(median(candidates.center_gap_mm, 'omitnan'), 0.002);
bound_ref = max(median(candidates.BoundRatio, 'omitnan'), 0.02);
lap_span = max(candidates.LapCount) - min(candidates.LapCount);
if lap_span <= 0
    lap_reward = zeros(height(candidates), 1);
else
    lap_reward = (candidates.LapCount - min(candidates.LapCount)) ./ lap_span;
end
balance_penalty = abs(candidates.LeftGradientRatio - candidates.RightGradientRatio);
balance_penalty(~isfinite(balance_penalty)) = 1;

candidates.IsStableCandidate = ...
    candidates.PulseCount >= cfg.step04_stable_window_min_pulses & ...
    candidates.SelectedPointCount >= cfg.step04_stable_window_min_points & ...
    candidates.WidthDeficit <= cfg.step04_stable_window_width_deficit_limit & ...
    candidates.LocalXcMAD <= max(3 * xc_ref, 0.012) & ...
    candidates.center_gap_mm <= max(3 * gap_ref, 0.030) & ...
    candidates.RMSE <= max(2.5 * rmse_ref, 0.10) & ...
    candidates.BoundRatio <= max(2.5 * bound_ref, 0.20);

start_vals = candidates.StartLap;
if all(~isfinite(start_vals)) || max(start_vals) <= min(start_vals)
    start_penalty = zeros(height(candidates), 1);
else
    start_penalty = (start_vals - min(start_vals)) ./ max(max(start_vals) - min(start_vals), eps);
end
if cfg.step04_stable_window_prefer_shorter
    lap_penalty = lap_reward;
else
    lap_penalty = -lap_reward;
end

policy = lower(strtrim(string(cfg.step04_stable_window_selection_policy)));
switch policy
    case {"legacy_width_first", "coverage_first"}
        score = 80 * candidates.WidthDeficit + ...
            12 * candidates.PointDeficit + ...
            1.5 * candidates.RMSE ./ rmse_ref + ...
            1.2 * candidates.LocalXcMAD ./ xc_ref + ...
            0.8 * candidates.center_gap_mm ./ gap_ref + ...
            0.8 * candidates.BoundRatio ./ bound_ref + ...
            0.8 * balance_penalty - ...
            1.0 * lap_reward;
    case {"earliest_short_stable_static_fit", "sg_style_static_fit"}
        score = robust_positive_score_step04_local(candidates.LocalXcMAD) + ...
            0.60 * robust_positive_score_step04_local(candidates.RMSE) + ...
            0.45 * robust_positive_score_step04_local(candidates.center_gap_mm) + ...
            0.25 * robust_positive_score_step04_local(candidates.BoundRatio) + ...
            0.20 * candidates.WidthDeficit + ...
            0.10 * candidates.PointDeficit + ...
            0.08 * balance_penalty + ...
            0.05 * start_penalty + ...
            0.20 * lap_penalty;
    otherwise
        score = 2.5 * candidates.LocalXcMAD ./ xc_ref + ...
            1.8 * candidates.RMSE ./ rmse_ref + ...
            1.2 * candidates.center_gap_mm ./ gap_ref + ...
            1.0 * candidates.WidthDeficit + ...
            0.8 * candidates.PointDeficit + ...
            0.6 * candidates.BoundRatio ./ bound_ref + ...
            0.4 * balance_penalty + ...
            0.25 * lap_penalty + ...
            cfg.step04_stable_window_start_penalty_weight * start_penalty;
end
score(~isfinite(score)) = inf;
score(~candidates.IsStableCandidate) = score(~candidates.IsStableCandidate) + 25;
candidates.ChosenScore = score;
end


function chosen = choose_stable_window_candidate_step04_local(candidates, cfg)
if isempty(candidates)
    chosen = table();
    return;
end
valid = candidates(isfinite(candidates.ChosenScore), :);
if isempty(valid)
    chosen = candidates(1, :);
    return;
end
stable = valid(valid.IsStableCandidate, :);
if isempty(stable)
    stable = valid;
end
policy = lower(strtrim(string(cfg.step04_stable_window_selection_policy)));
switch policy
    case {"legacy_width_first", "coverage_first"}
        stable = sortrows(stable, {'WidthDeficit', 'ChosenScore', 'LapCount'}, {'ascend', 'ascend', 'descend'});
    case {"earliest_short_stable_static_fit", "sg_style_static_fit"}
        early_limit = get_field_or_default_local(cfg, 'step04_stable_window_max_start_lap', inf);
        static = stable;
        if isfinite(early_limit) && ismember('StartLap', static.Properties.VariableNames)
            early_static = static(static.StartLap <= early_limit, :);
            if ~isempty(early_static)
                static = early_static;
            end
        end
        if ismember('center_gap_mm', static.Properties.VariableNames)
            center_gap_limit = get_field_or_default_local(cfg, 'step04_stable_window_center_gap_limit_mm', 0.030);
            center_ok = static.center_gap_mm <= center_gap_limit | ~isfinite(static.center_gap_mm);
            if any(center_ok)
                static = static(center_ok, :);
            end
        end
        if ismember('LocalXcMAD', static.Properties.VariableNames)
            local_mad_limit = get_field_or_default_local(cfg, 'step04_stable_window_local_xc_mad_limit_mm', 0.012);
            local_ok = static.LocalXcMAD <= local_mad_limit | ~isfinite(static.LocalXcMAD);
            if any(local_ok)
                static = static(local_ok, :);
            end
        end
        if ~isempty(static)
            stable = sortrows(static, ...
                {'StartLap', 'LapCount', 'ChosenScore', 'RMSE', 'LocalXcMAD', 'center_gap_mm'}, ...
                {'ascend', 'ascend', 'ascend', 'ascend', 'ascend', 'ascend'});
        else
            stable = sortrows(stable, {'ChosenScore', 'StartLap', 'RMSE'}, ...
                {'ascend', 'ascend', 'ascend'});
        end
    otherwise
        early_limit = get_field_or_default_local(cfg, 'step04_stable_window_max_start_lap', inf);
        score_tol = get_field_or_default_local(cfg, 'step04_stable_window_early_score_tolerance', 0);
        if isfinite(early_limit) && ismember('StartLap', stable.Properties.VariableNames)
            early_mask = stable.StartLap <= early_limit;
            if any(early_mask)
                best_score = min(stable.ChosenScore, [], 'omitnan');
                early = stable(early_mask & ...
                    stable.ChosenScore <= best_score + score_tol, :);
                if ~isempty(early)
                    stable = early;
                end
            end
        end
        if cfg.step04_stable_window_prefer_shorter
            stable = sortrows(stable, {'ChosenScore', 'LapCount', 'StartLap', 'WidthDeficit'}, ...
                {'ascend', 'ascend', 'ascend', 'ascend'});
        else
            stable = sortrows(stable, {'ChosenScore', 'StartLap', 'WidthDeficit'}, ...
                {'ascend', 'ascend', 'ascend'});
        end
end
chosen = stable(1, :);
end


function pc_out = subset_point_cloud_by_lap_window_step04_local(pc, start_lap, end_lap)
point_mask = pc.lap_index >= start_lap & pc.lap_index <= end_lap;
pc_out = pc;
pc_out.x_mm = pc.x_mm(point_mask);
pc_out.v = pc.v(point_mask);
pc_out.t_s = pc.t_s(point_mask);
pc_out.pulse_index = pc.pulse_index(point_mask);
pc_out.lap_index = pc.lap_index(point_mask);

all_pulses = unique(pc.pulse_index(:), 'stable');
keep_pulses = unique(pc_out.pulse_index(:), 'stable');
pulse_keep = ismember(all_pulses, keep_pulses);
if numel(pc.baseline_by_pulse) == numel(all_pulses)
    pc_out.baseline_by_pulse = pc.baseline_by_pulse(pulse_keep);
else
    pc_out.baseline_by_pulse = pc.baseline_by_pulse;
end
if isfield(pc, 'sensor_global_baseline')
    pc_out.sensor_global_baseline = pc.sensor_global_baseline;
end
if isfield(pc, 'sensor_baseline_mode')
    pc_out.sensor_baseline_mode = pc.sensor_baseline_mode;
end
if numel(pc.feature_by_pulse) == numel(all_pulses)
    pc_out.feature_by_pulse = pc.feature_by_pulse(pulse_keep);
else
    pc_out.feature_by_pulse = pc.feature_by_pulse;
end
pc_out.point_count = numel(pc_out.x_mm);
pc_out.pulse_count = numel(unique(pc_out.pulse_index));
end


function T = collect_stable_window_candidates_local(point_cloud)
T = empty_stable_window_candidate_table_step04_local();
for i = 1:numel(point_cloud)
    if istable(point_cloud(i).stable_window_candidates) && ~isempty(point_cloud(i).stable_window_candidates)
        T = [T; point_cloud(i).stable_window_candidates]; %#ok<AGROW>
    end
end
end


function row = make_empty_stable_window_candidate_row_step04_local()
row = struct( ...
    'sensor_id', NaN, 'blade_id', NaN, ...
    'StartLap', NaN, 'EndLap', NaN, 'LapCount', NaN, ...
    'PulseCount', NaN, 'WidePointCount', NaN, 'SelectedPointCount', NaN, ...
    'xc_mm', NaN, 'xc_seed_mm', NaN, 'center_gap_mm', NaN, 'SelectedWidthMm', NaN, ...
    'RMSE', NaN, 'BoundRatio', NaN, 'LocalXcMAD', NaN, 'PulseXcMAD', NaN, ...
    'LeftGradientRatio', NaN, 'RightGradientRatio', NaN, ...
    'WidthDeficit', NaN, 'PointDeficit', NaN, ...
    'IsStableCandidate', false, 'ChosenScore', NaN, 'Reason', "");
end


function T = empty_stable_window_candidate_table_step04_local()
T = struct2table(repmat(make_empty_stable_window_candidate_row_step04_local(), 0, 1));
end


function [rmse, bound_ratio] = estimate_static_fit_quality_step04_local(x_rel, v, w)
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


function s = robust_positive_score_step04_local(x)
x = x(:);
x_med = median(x, 'omitnan');
x_mad = median(abs(x - x_med), 'omitnan');
if ~isfinite(x_mad) || x_mad <= eps
    x_mad = std(x, 'omitnan');
end
if ~isfinite(x_mad) || x_mad <= eps
    s = zeros(size(x));
else
    s = max(0, (x - x_med) ./ x_mad);
end
end


function local_xc_mad = estimate_local_center_mad_step04_local(x, v, pulse_id, baseline, threshold)
finite = isfinite(x) & isfinite(v) & isfinite(pulse_id);
x = x(finite);
v = v(finite);
pulse_id = pulse_id(finite);
ids = unique(pulse_id(:), 'stable');
xc = nan(numel(ids), 1);
for i = 1:numel(ids)
    mask = pulse_id == ids(i);
    xv = x(mask);
    vv = v(mask);
    vz = max(vv - baseline, 0);
    high = vv >= threshold;
    if nnz(high) < 5 && ~isempty(vz)
        high = vz >= prctile(vz, 75);
    end
    if nnz(high) >= 3 && sum(vz(high), 'omitnan') > eps
        xc(i) = sum(xv(high) .* vz(high), 'omitnan') / sum(vz(high), 'omitnan');
    end
end
local_xc_mad = mad(xc, 1);
end


function [left_ratio, right_ratio] = estimate_gradient_side_balance_step04_local(x_selected, weight_selected, xc)
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


function [baseline, mode] = get_low_speed_sensor_baseline_local(case_data, sid)
baseline = NaN;
mode = '';
if ~isfield(case_data, 'LowSpeedBaseline') || ~istable(case_data.LowSpeedBaseline)
    return;
end
T = case_data.LowSpeedBaseline;
if ~ismember('sensor_id', T.Properties.VariableNames) || ~ismember('baseline_v', T.Properties.VariableNames)
    return;
end
mask = T.sensor_id == sid;
if ~any(mask)
    return;
end
baseline = T.baseline_v(find(mask, 1, 'first'));
if ismember('baseline_mode', T.Properties.VariableNames)
    mode = char(string(T.baseline_mode(find(mask, 1, 'first'))));
else
    mode = 'sensor_stream_masked_global';
end
end


function baseline = estimate_point_cloud_baseline_local(cfg, v, sensor_baseline)
if any(strcmpi(string(cfg.step04_baseline_mode), ...
        ["sensor_stream_masked_global", "sensor_stream_low_state_histogram"])) && isfinite(sensor_baseline)
    baseline = sensor_baseline;
else
    baseline = estimate_pulse_baseline_local(v);
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


function [x_mm, keep_mask] = map_time_to_x_local(t_seg, opr_times, F_omega_deg_s, theta_std_deg, cfg, t_arrival)
t_seg = t_seg(:);
idx_ref = find(opr_times < t_arrival, 1, 'last');
keep_mask = false(size(t_seg));
x_mm = nan(size(t_seg));
if isempty(idx_ref) || ~isfinite(t_arrival)
    return;
end
% Match the original 20251222 OPRCenterStd template route: all samples from
% one blade-passage waveform use the OPR center immediately before that
% blade arrival as the angular reference.
t_ref = opr_times(idx_ref);
for i = 1:numel(t_seg)
    t_now = t_seg(i);
    if ~isfinite(t_now)
        continue;
    end
    theta_actual = integrate_angle_local(t_ref, t_now, F_omega_deg_s);
    theta_diff = wrap_to_180_local(theta_actual - theta_std_deg);
    x_mm(i) = theta_diff * (pi / 180) * cfg.r_tip_mm;
    keep_mask(i) = isfinite(x_mm(i));
end
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
% Keep the sign when t1 is before t0. This matches the original 20251222
% OPRCenterStd low-speed template route and preserves the leading side of
% blade passages close to the OPR reference.
t_grid = linspace(t0, t1, 10);
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


function Template = aggregate_template_from_point_cloud_local(cfg, Sensor_Config, case_data, std_angles, point_cloud, speed_diagnostic, sensor_config_path)
entry_count = numel(point_cloud);
SensorBlade = repmat(make_empty_template_entry_local(), 1, entry_count);
summary_rows = repmat(struct( ...
    'sensor_id', NaN, 'blade_id', NaN, ...
    'pulse_count', NaN, 'point_count', NaN, 'wide_point_count', NaN, ...
    'stable_window_enabled', false, 'stable_window_start_lap', NaN, ...
    'stable_window_end_lap', NaN, 'stable_window_lap_count', NaN, ...
    'stable_window_score', NaN, 'stable_window_candidate_count', NaN, ...
    'stable_window_reason', '', ...
    'x_abs_min_mm', NaN, 'x_abs_max_mm', NaN, ...
    'x_min_mm', NaN, 'x_max_mm', NaN, ...
    'xc_mm', NaN, 'xc_detected_mm', NaN, 'xc_seed_mm', NaN, ...
    'xc_refined_mm', NaN, 'xc_refine_shift_mm', NaN, ...
    'xc_refine_status', '', 'xc_reference_source', '', ...
    'selection_mode', '', ...
    'domain_left_mm', NaN, 'domain_right_mm', NaN, 'domain_width_mm', NaN, ...
    'support_left_mm', NaN, 'support_right_mm', NaN, ...
    'response_left_mm', NaN, 'response_right_mm', NaN, ...
    'query_safe_left_mm', NaN, 'query_safe_right_mm', NaN, ...
    'domain_status', '', 'domain_peak_amp_v', NaN, 'domain_noise_sigma_v', NaN, ...
    'source_point_policy', '', 'final_point_count', NaN, ...
    'baseline_v', NaN, 'baseline_source', '', 'threshold_v', NaN, ...
    'amplitude_v', NaN, 'max_abs_gradient_v_per_mm', NaN, ...
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
stable_window_table = collect_stable_window_candidates_local(point_cloud);

Template = struct();
Template.SchemaVersion = 'LowSpeedOPRCenterStdTemplate/v1';
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
Template.OPR_Timing_Method = get_string_field_local(case_data, 'opr_timing_method', get_string_field_local(Sensor_Config, 'OPR_Timing_Method', 'unknown'));
Template.Step01_OPR_Timing_Method = get_string_field_local(Sensor_Config, 'OPR_Timing_Method', 'unknown');
Template.Step04_OPR_Timing_Method = get_string_field_local(case_data, 'opr_timing_method', 'unknown');
Template.Step01_Probe_Arrival_Method = get_string_field_local(Sensor_Config, 'Probe_Arrival_Method', 'unknown');
Template.Step04_Probe_Arrival_Method = get_string_field_local(case_data, 'probe_arrival_method', 'unknown');
Template.Probe_Arrival_Method = Template.Step04_Probe_Arrival_Method;
Template.Time_Index_Mode = get_string_field_local(Sensor_Config, 'Time_Index_Mode', get_string_field_local(case_data, 'time_index_mode', 'unknown'));
Template.Time_Index_Mode_Detected = get_string_field_local(Sensor_Config, 'Time_Index_Mode_Detected', get_string_field_local(case_data, 'time_index_mode_detected', 'unknown'));
Template.OPR_Events_Per_Revolution = get_numeric_field_local(speed_diagnostic, 'opr_events_per_revolution', cfg.step04_opr_events_per_revolution);
Template.OPR_Events_Per_Revolution_Source = cfg.step04_opr_events_per_revolution_source;
Template.Blade_ID_Definition = get_string_field_local(Sensor_Config, 'Blade_ID_Definition', 'internal Blade 1 from low-speed reference');
Template.Standard_Angle_Source = cfg.step04_standard_angle_source;
Template.Standard_Angle_Field = cfg.step04_standard_angle_field;
Template.Standard_Relative_Angles = std_angles;
if isfield(Sensor_Config, 'Standard_Relative_Angles_Std')
    Template.Standard_Relative_Angles_Std = Sensor_Config.Standard_Relative_Angles_Std;
end
if isfield(Sensor_Config, 'Standard_Relative_Angles_Count')
    Template.Standard_Relative_Angles_Count = Sensor_Config.Standard_Relative_Angles_Count;
end
Template.Speed_Diagnostic = speed_diagnostic;
if isfield(case_data, 'LowSpeedBaseline')
    Template.LowSpeedBaseline = case_data.LowSpeedBaseline;
end
Template.PointCloud_Settings = struct( ...
    'pulse_pad_fraction', cfg.step04_pulse_pad_fraction, ...
    'max_points_per_pulse', cfg.step04_max_points_per_pulse, ...
    'x_collect_abs_limit_mm', cfg.step04_x_collect_abs_limit_mm, ...
    'stable_window_enable', cfg.step04_stable_window_enable, ...
    'stable_window_lap_counts', cfg.step04_stable_window_lap_counts, ...
    'stable_window_scan_step_laps', cfg.step04_stable_window_scan_step_laps, ...
    'stable_window_min_points', cfg.step04_stable_window_min_points, ...
    'stable_window_min_pulses', cfg.step04_stable_window_min_pulses, ...
    'stable_window_target_width_mm', cfg.step04_stable_window_target_width_mm, ...
    'stable_window_selection_policy', cfg.step04_stable_window_selection_policy, ...
    'stable_window_width_deficit_limit', cfg.step04_stable_window_width_deficit_limit, ...
    'stable_window_prefer_shorter', cfg.step04_stable_window_prefer_shorter, ...
    'stable_window_max_start_lap', cfg.step04_stable_window_max_start_lap, ...
    'stable_window_start_penalty_weight', cfg.step04_stable_window_start_penalty_weight, ...
    'stable_window_early_score_tolerance', cfg.step04_stable_window_early_score_tolerance, ...
    'baseline_mode', cfg.step04_baseline_mode, ...
    'baseline_pad_points', cfg.step04_baseline_pad_points, ...
    'baseline_trim_fraction', cfg.step04_baseline_trim_fraction, ...
    'baseline_low_state_threshold_v', cfg.step04_baseline_low_state_threshold_v, ...
    'baseline_histogram_bins', cfg.step04_baseline_histogram_bins, ...
    'baseline_histogram_half_width_v', cfg.step04_baseline_histogram_half_width_v);
Template.Template_Settings = struct( ...
    'grid_dx_mm', cfg.step04_grid_dx_mm, ...
    'min_bin_count', cfg.step04_min_bin_count, ...
    'smooth_span_bins', cfg.step04_smooth_span_bins, ...
    'gradient_min_ratio', cfg.step04_gradient_min_ratio, ...
    'amplitude_min_ratio', cfg.step04_amplitude_min_ratio, ...
    'min_half_width_mm', cfg.step04_min_half_width_mm, ...
    'max_half_width_mm', cfg.step04_max_half_width_mm, ...
    'domain_center_mode', cfg.step04_domain_center_mode, ...
    'center_mode', cfg.step04_center_mode, ...
    'center_seed_source', cfg.step04_center_seed_source, ...
    'center_refine_max_shift_mm', cfg.step04_center_refine_max_shift_mm, ...
    'xrange_mode', cfg.step04_xrange_mode, ...
    'template_point_selection_mode', cfg.step04_template_point_selection_mode, ...
    'domain_min_points_per_bin', cfg.step04_domain_min_points_per_bin, ...
    'domain_min_laps_per_bin', cfg.step04_domain_min_laps_per_bin, ...
    'domain_noise_sigma_factor', cfg.step04_domain_noise_sigma_factor, ...
    'domain_std_abs_limit_v', cfg.step04_domain_std_abs_limit_v, ...
    'domain_repeat_ratio', cfg.step04_domain_repeat_ratio, ...
    'domain_boundary_margin_mm', cfg.step04_domain_boundary_margin_mm, ...
    'query_guard_mm', cfg.step04_query_guard_mm, ...
    'final_template_domain_buffer_mm', cfg.step04_final_template_domain_buffer_mm, ...
    'trust_quantile', cfg.step04_trust_quantile, ...
    'trust_edge_margin_mm', cfg.step04_trust_edge_margin_mm, ...
    'gradient_energy_quantile', cfg.step04_gradient_energy_quantile, ...
    'gradient_edge_margin_mm', cfg.step04_gradient_edge_margin_mm);
Template.Final_Template_Point_Policy = 'trusted_domain_recalibration';
Template.Source_Point_Policy = 'raw_low_speed_waveform_to_signed_oprcenterstd_point_cloud';
Template.SensorBlade = SensorBlade;
Template.Summary_Table = summary_table;
Template.Stable_Window_Candidates = stable_window_table;
Template.Metadata = build_template_metadata_local(Template, cfg);
end


function empty = make_empty_template_entry_local()
empty = struct( ...
    'sensor_id', NaN, 'blade_id', NaN, ...
    'x_grid', [], 'v_grid', [], 'v_grid_raw', [], 'v_grid_baseline_removed', [], ...
    'dv_dx', [], 'weight_grid', [], 'count_grid', [], 'lap_count_grid', [], ...
    'v_std_grid', [], 'valid_grid_mask', [], ...
    'wide_x_grid', [], 'wide_v_grid', [], 'wide_dv_dx', [], 'wide_count_grid', [], ...
    'preliminary_x_grid', [], 'preliminary_v_grid', [], ...
    'preliminary_v_grid_raw', [], 'preliminary_count_grid', [], ...
    'preliminary_lap_count_grid', [], 'preliminary_v_std_grid', [], ...
    'domain_effective_mask', [], 'domain_mask', [], ...
    'x_collect_domain_mm', [NaN NaN], 'x_support_domain_mm', [NaN NaN], ...
    'x_response_domain_mm', [NaN NaN], 'x_domain', [NaN NaN], ...
    'x_collect_domain', [NaN NaN], 'x_support_domain', [NaN NaN], ...
    'x_response_domain', [NaN NaN], ...
    'x_query_safe_domain', [NaN NaN], 'domain_diagnostic', struct(), ...
    'xc', NaN, 'xc_mm', NaN, 'xc_detected', NaN, 'xc_detected_mm', NaN, ...
    'xc_seed', NaN, 'xc_refined', NaN, 'xc_refine_shift_mm', NaN, ...
    'xc_refine_status', '', ...
    'xc_reference_source', '', 'selection_mode', '', ...
    'baseline', NaN, 'baseline_source', '', 'threshold', NaN, 'amplitude', NaN, ...
    'pulse_count', 0, 'point_count', 0, 'wide_point_count', 0, ...
    'stable_window_enabled', false, 'stable_window_start_lap', NaN, ...
    'stable_window_end_lap', NaN, 'stable_window_lap_count', NaN, ...
    'stable_window_score', NaN, 'stable_window_candidate_count', 0, ...
    'stable_window_reason', '', ...
    'source_point_policy', '', 'final_template_point_policy', '', ...
    'final_point_count', 0, 'trusted_point_count', 0, ...
    'quality_status', '', 'x_points_preview', [], ...
    'x_used_points_preview', [], 'v_used_points_preview', [], ...
    'x_abs_points_preview', [], 'x_wide_points_preview', [], ...
    'v_points_preview', [], 'weight_points_preview', []);
end


function write_step04_diagnostic_tables_local(out_dir, cfg, Template, case_data, speed_diagnostic)
if isfield(speed_diagnostic, 'opr_semantics_table') && istable(speed_diagnostic.opr_semantics_table)
    writetable(speed_diagnostic.opr_semantics_table, fullfile(out_dir, ...
        sprintf('Step04_OPRSemanticsCheck_%s.csv', cfg.dataset)));
end
if isfield(case_data, 'FileTimeContinuityTable') && istable(case_data.FileTimeContinuityTable)
    writetable(case_data.FileTimeContinuityTable, fullfile(out_dir, ...
        sprintf('Step04_File_Time_Continuity_%s.csv', cfg.dataset)));
end
if isfield(case_data, 'SourceInventoryTable') && istable(case_data.SourceInventoryTable)
    writetable(case_data.SourceInventoryTable, fullfile(out_dir, ...
        sprintf('Step04_SourceInventory_%s.csv', cfg.dataset)));
end
if isfield(case_data, 'LowSpeedBaseline') && istable(case_data.LowSpeedBaseline)
    writetable(case_data.LowSpeedBaseline, fullfile(out_dir, ...
        sprintf('Step04_LowSpeedBaseline_%s.csv', cfg.dataset)));
end
if isfield(Template, 'LowSpeedBaseline') && istable(Template.LowSpeedBaseline)
    writetable(Template.LowSpeedBaseline, fullfile(out_dir, ...
        sprintf('Step04_LowSpeedBaseline_Final_%s.csv', cfg.dataset)));
end
if isfield(Template, 'Stable_Window_Candidates') && istable(Template.Stable_Window_Candidates)
    writetable(Template.Stable_Window_Candidates, fullfile(out_dir, ...
        sprintf('Step04_StableWindow_Candidates_%s.csv', cfg.dataset)));
end
S = Template.Summary_Table;
if istable(S)
    wanted = {'sensor_id','blade_id','stable_window_enabled','stable_window_start_lap', ...
        'stable_window_end_lap','stable_window_lap_count','stable_window_score', ...
        'stable_window_candidate_count','stable_window_reason','domain_left_mm', ...
        'domain_right_mm','domain_width_mm','query_safe_left_mm','query_safe_right_mm', ...
        'xc_mm','xc_detected_mm','xc_seed_mm','xc_refined_mm','xc_refine_shift_mm', ...
        'xc_refine_status','baseline_v','baseline_source','domain_status','quality_status'};
    keep = intersect(wanted, S.Properties.VariableNames, 'stable');
    writetable(S(:, keep), fullfile(out_dir, ...
        sprintf('Step04_StableWindow_Selected_%s.csv', cfg.dataset)));
end
end


function [entry, row] = aggregate_one_sensor_blade_local(cfg, pc, sid, blade_id)
entry = make_empty_template_entry_local();
entry.sensor_id = sid;
entry.blade_id = blade_id;

row = struct( ...
    'sensor_id', sid, 'blade_id', blade_id, ...
    'pulse_count', pc.pulse_count, 'point_count', pc.point_count, ...
    'wide_point_count', NaN, ...
    'stable_window_enabled', pc.stable_window_enabled, ...
    'stable_window_start_lap', pc.stable_window_start_lap, ...
    'stable_window_end_lap', pc.stable_window_end_lap, ...
    'stable_window_lap_count', pc.stable_window_lap_count, ...
    'stable_window_score', pc.stable_window_score, ...
    'stable_window_candidate_count', pc.stable_window_candidate_count, ...
    'stable_window_reason', pc.stable_window_reason, ...
    'x_abs_min_mm', NaN, 'x_abs_max_mm', NaN, ...
    'x_min_mm', NaN, 'x_max_mm', NaN, ...
    'xc_mm', NaN, 'xc_detected_mm', NaN, 'xc_reference_source', '', ...
    'selection_mode', '', ...
    'domain_left_mm', NaN, 'domain_right_mm', NaN, 'domain_width_mm', NaN, ...
    'support_left_mm', NaN, 'support_right_mm', NaN, ...
    'response_left_mm', NaN, 'response_right_mm', NaN, ...
    'query_safe_left_mm', NaN, 'query_safe_right_mm', NaN, ...
    'domain_status', '', 'domain_peak_amp_v', NaN, 'domain_noise_sigma_v', NaN, ...
    'source_point_policy', '', 'final_point_count', NaN, ...
    'baseline_v', NaN, 'threshold_v', NaN, ...
    'amplitude_v', NaN, 'max_abs_gradient_v_per_mm', NaN, ...
    'valid_bin_count', NaN, 'coverage_fraction', NaN, ...
    'quality_status', 'empty');

x_abs = pc.x_mm(:);
v = pc.v(:);
lap_abs = pc.lap_index(:);
pulse_abs = pc.pulse_index(:);
valid = isfinite(x_abs) & isfinite(v);
x_abs = x_abs(valid);
v = v(valid);
lap_abs = lap_abs(valid);
pulse_abs = pulse_abs(valid);
entry.point_count = numel(x_abs);
entry.wide_point_count = numel(x_abs);
entry.pulse_count = pc.pulse_count;
entry.stable_window_enabled = pc.stable_window_enabled;
entry.stable_window_start_lap = pc.stable_window_start_lap;
entry.stable_window_end_lap = pc.stable_window_end_lap;
entry.stable_window_lap_count = pc.stable_window_lap_count;
entry.stable_window_score = pc.stable_window_score;
entry.stable_window_candidate_count = pc.stable_window_candidate_count;
entry.stable_window_reason = pc.stable_window_reason;
if isempty(x_abs)
    entry.quality_status = 'empty';
    return;
end

row.x_abs_min_mm = min(x_abs);
row.x_abs_max_mm = max(x_abs);
row.wide_point_count = numel(x_abs);

[selection, center_info] = select_low_speed_template_points_local(cfg, sid, x_abs, v, pc, lap_abs, pulse_abs);
if ~isfinite(selection.xc) || isempty(selection.x_selected)
    entry.quality_status = 'invalid_xcenter';
    row.quality_status = entry.quality_status;
    return;
end

xc = selection.xc;
x_selected_abs = selection.x_selected(:);
v_selected = selection.v_selected(:);
point_weight = selection.weight_selected(:);
lap_selected = selection.lap_selected(:);
x = x_selected_abs - xc;
v = v_selected;

row.x_min_mm = min(x);
row.x_max_mm = max(x);
row.xc_mm = xc;
row.xc_detected_mm = selection.xc_detected;
row.xc_seed_mm = selection.xc_seed;
row.xc_refined_mm = selection.xc_refined;
row.xc_refine_shift_mm = selection.xc_refine_shift_mm;
row.xc_refine_status = selection.xc_refine_status;
row.xc_reference_source = center_info.reference_source;
row.selection_mode = selection.mode;
row.point_count = numel(x);
row.baseline_v = selection.baseline;
row.baseline_source = selection.baseline_source;
row.threshold_v = selection.threshold;
entry.point_count = numel(x);

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
weight_grid_raw = nan(n_bins, 1);
v_std_grid = nan(n_bins, 1);
lap_count_grid = zeros(n_bins, 1);
count_grid = zeros(n_bins, 1);
for ib = 1:n_bins
    mask = bin == ib;
    count_grid(ib) = nnz(mask);
    if count_grid(ib) > 0
        v_grid_raw(ib) = median(v(mask), 'omitnan');
        v_std_grid(ib) = std(v(mask), 'omitnan');
        weight_grid_raw(ib) = mean(point_weight(mask), 'omitnan');
        laps_here = lap_selected(mask);
        lap_count_grid(ib) = numel(unique(laps_here(isfinite(laps_here))));
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

baseline = selection.baseline;
if ~isfinite(baseline)
    baseline = estimate_template_baseline_local(x, v, v_smooth);
end
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

Domain = determine_trusted_domain_local(cfg, x_grid(:), v_smooth(:), dv_dx(:), ...
    count_grid(:), lap_count_grid(:), v_std_grid(:), valid_bin_mask(:), baseline, amplitude, max_abs_grad);
x_domain = Domain.x_domain;

preliminary_x_grid = x_grid(:);
preliminary_v_grid = v_smooth(:);
preliminary_v_grid_raw = v_grid_raw(:);
preliminary_count_grid = count_grid(:);
preliminary_lap_count_grid = lap_count_grid(:);
preliminary_v_std_grid = v_std_grid(:);

source_point_policy = 'wide_point_cloud_preliminary_only';
final_point_count = 0;
final_used_mask = false(size(x));
if strcmpi(string(Domain.status), "ok")
    fit_domain = [x_domain(1) - cfg.step04_final_template_domain_buffer_mm, ...
        x_domain(2) + cfg.step04_final_template_domain_buffer_mm];
    fit_mask = x >= fit_domain(1) & x <= fit_domain(2) & isfinite(x) & isfinite(v);
    if nnz(fit_mask) >= cfg.step04_min_template_points
        FinalGrid = build_template_grid_from_points_step04_local(cfg, x(fit_mask), v(fit_mask), ...
            point_weight(fit_mask), lap_selected(fit_mask), x_domain);
        if FinalGrid.ok
            x_grid = FinalGrid.x_grid;
            v_grid_raw = FinalGrid.v_grid_raw;
            v_smooth = FinalGrid.v_smooth;
            count_grid = FinalGrid.count_grid;
            lap_count_grid = FinalGrid.lap_count_grid;
            v_std_grid = FinalGrid.v_std_grid;
            weight_grid_raw = FinalGrid.weight_grid_raw;
            valid_bin_mask = FinalGrid.valid_bin_mask;
            dv_dx = gradient(v_smooth, x_grid);
            max_abs_grad = max(abs(dv_dx(valid_bin_mask)), [], 'omitnan');
            if ~isfinite(max_abs_grad)
                max_abs_grad = 0;
            end
            source_point_policy = 'trusted_domain_recalibration';
            final_point_count = nnz(fit_mask);
            final_used_mask = fit_mask;
        end
    end
end

domain_mask = x_grid >= x_domain(1) & x_grid <= x_domain(2) & valid_bin_mask;
effective_mask = domain_mask;
weight_grid = build_template_weight_grid_local(count_grid(:), valid_bin_mask(:), domain_mask(:), dv_dx(:), max_abs_grad, weight_grid_raw(:));

entry.x_grid = x_grid(:);
entry.v_grid = v_smooth(:);
entry.v_grid_raw = v_grid_raw(:);
entry.v_grid_baseline_removed = v_smooth(:) - baseline;
entry.dv_dx = dv_dx(:);
entry.weight_grid = weight_grid(:);
entry.count_grid = count_grid(:);
entry.lap_count_grid = lap_count_grid(:);
entry.v_std_grid = v_std_grid(:);
entry.preliminary_x_grid = preliminary_x_grid;
entry.preliminary_v_grid = preliminary_v_grid;
entry.preliminary_v_grid_raw = preliminary_v_grid_raw;
entry.preliminary_count_grid = preliminary_count_grid;
entry.preliminary_lap_count_grid = preliminary_lap_count_grid;
entry.preliminary_v_std_grid = preliminary_v_std_grid;
entry.wide_x_grid = preliminary_x_grid;
entry.wide_v_grid = preliminary_v_grid;
entry.wide_dv_dx = gradient(preliminary_v_grid, preliminary_x_grid);
entry.wide_count_grid = preliminary_count_grid;
entry.valid_grid_mask = valid_bin_mask(:);
entry.domain_effective_mask = effective_mask(:);
entry.domain_mask = domain_mask(:);
entry.x_domain = x_domain;
entry.x_collect_domain_mm = [-cfg.step04_x_collect_abs_limit_mm, cfg.step04_x_collect_abs_limit_mm];
entry.x_support_domain_mm = Domain.x_support_domain_mm;
entry.x_response_domain_mm = Domain.x_response_domain_mm;
entry.x_collect_domain = entry.x_collect_domain_mm;
entry.x_support_domain = entry.x_support_domain_mm;
entry.x_response_domain = entry.x_response_domain_mm;
entry.x_query_safe_domain = Domain.x_query_safe_domain;
entry.domain_diagnostic = Domain;
entry.xc = xc;
entry.xc_mm = xc;
entry.xc_detected = selection.xc_detected;
entry.xc_detected_mm = selection.xc_detected;
entry.xc_seed = selection.xc_seed;
entry.xc_refined = selection.xc_refined;
entry.xc_refine_shift_mm = selection.xc_refine_shift_mm;
entry.xc_refine_status = selection.xc_refine_status;
entry.xc_reference_source = center_info.reference_source;
entry.selection_mode = selection.mode;
entry.baseline = baseline;
entry.baseline_source = selection.baseline_source;
entry.threshold = selection.threshold;
entry.amplitude = amplitude;
entry.source_point_policy = source_point_policy;
entry.final_template_point_policy = source_point_policy;
entry.final_point_count = final_point_count;
entry.trusted_point_count = final_point_count;
if strcmpi(source_point_policy, 'trusted_domain_recalibration')
    quality_point_count = final_point_count;
else
    quality_point_count = 0;
end
entry.quality_status = evaluate_template_quality_local(cfg, pc.pulse_count, quality_point_count, x_domain, valid_bin_mask);

preview_n = min(numel(x), 4000);
if numel(x) > preview_n
    idx = unique(round(linspace(1, numel(x), preview_n)));
else
    idx = 1:numel(x);
end
entry.x_points_preview = x(idx);
entry.x_abs_points_preview = x_selected_abs(idx);
entry.v_points_preview = v(idx);
entry.weight_points_preview = point_weight(idx);
used_idx_all = find(final_used_mask);
used_preview_n = min(numel(used_idx_all), 4000);
if used_preview_n > 0
    if numel(used_idx_all) > used_preview_n
        used_idx = used_idx_all(unique(round(linspace(1, numel(used_idx_all), used_preview_n))));
    else
        used_idx = used_idx_all;
    end
    entry.x_used_points_preview = x(used_idx);
    entry.v_used_points_preview = v(used_idx);
end
wide_preview_n = min(numel(x_abs), 4000);
if numel(x_abs) > wide_preview_n
    wide_idx = unique(round(linspace(1, numel(x_abs), wide_preview_n)));
else
    wide_idx = 1:numel(x_abs);
end
entry.x_wide_points_preview = x_abs(wide_idx) - xc;

row.domain_left_mm = x_domain(1);
row.domain_right_mm = x_domain(2);
row.domain_width_mm = x_domain(2) - x_domain(1);
row.support_left_mm = Domain.x_support_domain_mm(1);
row.support_right_mm = Domain.x_support_domain_mm(2);
row.response_left_mm = Domain.x_response_domain_mm(1);
row.response_right_mm = Domain.x_response_domain_mm(2);
row.query_safe_left_mm = Domain.x_query_safe_domain(1);
row.query_safe_right_mm = Domain.x_query_safe_domain(2);
row.domain_status = char(Domain.status);
row.domain_peak_amp_v = Domain.peak_amp;
row.domain_noise_sigma_v = Domain.noise_sigma;
row.source_point_policy = source_point_policy;
row.final_point_count = final_point_count;
row.baseline_v = baseline;
row.amplitude_v = amplitude;
row.max_abs_gradient_v_per_mm = max_abs_grad;
row.valid_bin_count = nnz(valid_bin_mask);
row.coverage_fraction = mean(valid_bin_mask);
row.quality_status = entry.quality_status;
end


function [P, info] = select_low_speed_template_points_local(cfg, sid, x_wide, v_wide, pc, lap_wide, pulse_wide)
if nargin < 6 || isempty(lap_wide)
    lap_wide = nan(size(x_wide));
end
if nargin < 7 || isempty(pulse_wide)
    pulse_wide = nan(size(x_wide));
end
finite = isfinite(x_wide) & isfinite(v_wide);
x_wide = x_wide(finite);
v_wide = v_wide(finite);
lap_wide = lap_wide(finite);
pulse_wide = pulse_wide(finite);

baseline = NaN;
baseline_source = 'pulse_edge_median';
if nargin >= 5 && isstruct(pc) && isfield(pc, 'sensor_global_baseline') && isfinite(pc.sensor_global_baseline) && ...
        any(strcmpi(string(cfg.step04_baseline_mode), ...
        ["sensor_stream_masked_global", "sensor_stream_low_state_histogram"]))
    baseline = pc.sensor_global_baseline;
    baseline_source = char(string(cfg.step04_baseline_mode));
elseif nargin >= 5 && isstruct(pc) && isfield(pc, 'baseline_by_pulse') && ~isempty(pc.baseline_by_pulse)
    baseline = median(pc.baseline_by_pulse(isfinite(pc.baseline_by_pulse)), 'omitnan');
    baseline_source = 'median_pulse_baseline';
end
if ~isfinite(baseline)
    baseline = estimate_template_selection_baseline_local(v_wide);
    baseline_source = 'template_low_percentile_fallback';
end
threshold = get_sensor_threshold_local(cfg, sid);
v_zero = max(v_wide - baseline, 0);

info = struct( ...
    'reference_source', 'weighted_high_voltage_centroid_xcenter', ...
    'baseline_for_center_v', baseline, ...
    'threshold_v', threshold, ...
    'high_point_count', 0);

P = struct( ...
    'mode', 'empty', ...
    'xc', NaN, ...
    'xc_detected', NaN, ...
    'xc_seed', NaN, ...
    'xc_refined', NaN, ...
    'xc_refine_shift_mm', NaN, ...
    'xc_refine_status', 'not_run', ...
    'x_wide', x_wide(:), ...
    'v_wide', v_wide(:), ...
    'x_selected', [], ...
    'v_selected', [], ...
    'weight_selected', [], ...
    'lap_selected', [], ...
    'pulse_selected', [], ...
    'baseline', baseline, ...
    'baseline_source', baseline_source, ...
    'threshold', threshold);

if numel(x_wide) < 50 || ~isfinite(baseline)
    info.reference_source = 'insufficient_points_for_selection';
    return;
end

high_level = prctile(v_zero, 85);
high_mask = v_zero >= high_level & v_zero > 0;
info.high_point_count = nnz(high_mask);
if nnz(high_mask) >= 5 && sum(v_zero(high_mask), 'omitnan') > eps
    xc = sum(x_wide(high_mask) .* v_zero(high_mask), 'omitnan') ./ ...
        sum(v_zero(high_mask), 'omitnan');
else
    [~, imax] = max(v_zero);
    xc = x_wide(imax);
    info.reference_source = 'voltage_peak_fallback_xcenter';
end
xc_detected = xc;
xc_seed = choose_step04_center_seed_local(cfg, pc, xc_detected);
xc = xc_seed;

strict_mask = v_wide >= threshold;
dx_strict = x_wide(strict_mask) - xc;
if nnz(strict_mask) < 30
    strict_mask = v_zero >= prctile(v_zero, 60);
    dx_strict = x_wide(strict_mask) - xc;
end

if strcmpi(cfg.step04_template_point_selection_mode, 'wide_window')
    trust_mask = true(size(x_wide));
    xrange_mode_used = 'wide_window';
else
    [left_L, right_L, xrange_ok, xrange_mode_used] = estimate_xrange_span_step04_local( ...
        x_wide, v_wide, baseline, xc, cfg);
    if ~xrange_ok
        left_dist = -dx_strict(dx_strict < 0);
        right_dist = dx_strict(dx_strict > 0);
        if isempty(left_dist) || isempty(right_dist)
            L = prctile(abs(dx_strict), 95);
            left_L = L;
            right_L = L;
        else
            left_L = prctile(left_dist, 100 * cfg.step04_trust_quantile);
            right_L = prctile(right_dist, 100 * cfg.step04_trust_quantile);
        end
        left_L = max(left_L - cfg.step04_trust_edge_margin_mm, 0.1);
        right_L = max(right_L - cfg.step04_trust_edge_margin_mm, 0.1);
    end

    if xrange_ok
        trust_mask = x_wide >= xc - left_L & x_wide <= xc + right_L;
    else
        trust_mask = strict_mask & x_wide >= xc - left_L & x_wide <= xc + right_L;
    end
    if nnz(trust_mask) < 30
        trust_mask = strict_mask;
    end
end

x_selected = x_wide(trust_mask);
v_selected = v_wide(trust_mask);
lap_selected = lap_wide(trust_mask);
pulse_selected = pulse_wide(trust_mask);
xc_refined = NaN;
xc_refine_shift = NaN;
xc_refine_status = 'not_run';
if strcmpi(cfg.step04_center_mode, 'sgfit')
    if strcmpi(cfg.step04_template_point_selection_mode, 'wide_window') && nnz(strict_mask) >= 30
        xc_refined = refine_center_by_sg_fit_step04_local(x_wide(strict_mask), v_wide(strict_mask), baseline, xc_seed);
    else
        xc_refined = refine_center_by_sg_fit_step04_local(x_selected, v_selected, baseline, xc_seed);
    end
    xc_refine_shift = xc_refined - xc_seed;
    if isfinite(xc_refined) && abs(xc_refine_shift) <= cfg.step04_center_refine_max_shift_mm
        xc = xc_refined;
        xc_refine_status = 'accepted';
    elseif isfinite(xc_refined)
        xc = xc_seed;
        xc_refine_status = 'rejected_large_shift';
        warning('Step04:CenterRefineLargeShift', ...
            ['CH%d B? center refinement shifted %.4f mm from the stable-window seed; ', ...
            'using seed %.4f mm instead of refined %.4f mm.'], ...
            sid, xc_refine_shift, xc_seed, xc_refined);
    else
        xc = xc_seed;
        xc_refine_status = 'invalid_refine_fallback_seed';
    end
end

weight_selected = build_gradient_weight_from_wide_cloud_step04_local( ...
    x_selected, x_wide, v_wide, baseline, 0.05);

if strcmpi(cfg.step04_template_point_selection_mode, 'wide_window')
    mode = 'raw_low_speed_wide_window_full_point_cloud_gradient_weight';
elseif xrange_ok
    mode = ['raw_low_speed_', xrange_mode_used, '_gate_gradient_weight'];
else
    mode = 'raw_low_speed_sg_like_strict_gate_trust_window_gradient_weight';
end
if strcmpi(cfg.step04_center_mode, 'sgfit')
    mode = [mode, '_sgfit_center'];
end

P.mode = mode;
P.xc = xc;
P.xc_detected = xc_detected;
P.xc_seed = xc_seed;
P.xc_refined = xc_refined;
P.xc_refine_shift_mm = xc_refine_shift;
P.xc_refine_status = xc_refine_status;
P.x_selected = x_selected(:);
P.v_selected = v_selected(:);
P.weight_selected = weight_selected(:);
P.lap_selected = lap_selected(:);
P.pulse_selected = pulse_selected(:);
P.baseline = baseline;
P.baseline_source = baseline_source;
P.threshold = threshold;
info.reference_source = 'detected_from_current_template_points';
end


function xc_seed = choose_step04_center_seed_local(cfg, pc, xc_detected)
xc_seed = xc_detected;
source = lower(strtrim(string(cfg.step04_center_seed_source)));
if source == "stable_window_xc" && isstruct(pc) && isfield(pc, 'stable_window_score') && ...
        isfield(pc, 'stable_window_candidates') && istable(pc.stable_window_candidates)
    C = pc.stable_window_candidates;
    if isfield(pc, 'stable_window_start_lap') && isfinite(pc.stable_window_start_lap) && ...
            isfield(pc, 'stable_window_end_lap') && isfinite(pc.stable_window_end_lap) && ...
            all(ismember({'StartLap','EndLap','xc_mm'}, C.Properties.VariableNames))
        mask = C.StartLap == pc.stable_window_start_lap & C.EndLap == pc.stable_window_end_lap;
        if any(mask)
            candidate_xc = C.xc_mm(find(mask, 1, 'first'));
            if isfinite(candidate_xc)
                xc_seed = candidate_xc;
            end
        end
    end
end
end


function threshold = get_sensor_threshold_local(cfg, sid)
threshold = cfg.sensor_threshold_default;
if isfield(cfg, 'sensor_thresholds') && isa(cfg.sensor_thresholds, 'containers.Map')
    try
        if isKey(cfg.sensor_thresholds, sid)
            threshold = cfg.sensor_thresholds(sid);
        end
    catch
        threshold = cfg.sensor_threshold_default;
    end
end
if ~isfinite(threshold)
    threshold = cfg.sensor_threshold_default;
end
end


function baseline = estimate_template_selection_baseline_local(v)
v = v(:);
if isempty(v)
    baseline = NaN;
    return;
end
low_mask = v <= prctile(v, 30);
if nnz(low_mask) >= 10
    baseline = median(v(low_mask), 'omitnan');
else
    baseline = prctile(v, 10);
end
if ~isfinite(baseline)
    baseline = median(v, 'omitnan');
end
end


function [left_L, right_L, ok, mode_used] = estimate_xrange_span_step04_local(x, v, baseline, xc, cfg)
if strcmpi(cfg.step04_xrange_mode, 'threshold')
    [left_L, right_L, ok] = estimate_gradient_threshold_span_step04_local( ...
        x, v, baseline, xc, cfg.step04_gradient_min_ratio, ...
        cfg.step04_amplitude_min_ratio, cfg.step04_min_half_width_mm, ...
        cfg.step04_max_half_width_mm);
    mode_used = 'gradient_threshold';
else
    [left_L, right_L, ok] = estimate_gradient_energy_span_step04_local( ...
        x, v, baseline, xc, cfg.step04_gradient_energy_quantile, ...
        cfg.step04_gradient_edge_margin_mm);
    mode_used = 'gradient_energy';
end
end


function [left_L, right_L, ok] = estimate_gradient_threshold_span_step04_local( ...
    x, v, baseline, xc, gradient_min_ratio, amplitude_min_ratio, min_half_width, max_half_width)
left_L = NaN;
right_L = NaN;
ok = false;
[x_grid, v_smooth] = build_denoised_static_profile_step04_local(x, v, baseline);
if numel(x_grid) < 20
    return;
end
amp_norm = normalize01_step04_local(v_smooth(:));
g_abs = abs(gradient(v_smooth(:), x_grid(:)));
g_threshold = estimate_noise_aware_gradient_threshold_step04_local(g_abs, amp_norm, gradient_min_ratio);
effective = g_abs >= g_threshold;
if amplitude_min_ratio > 0
    effective = effective & amp_norm >= amplitude_min_ratio;
end
left_candidates = find(x_grid(:) < xc & effective(:));
right_candidates = find(x_grid(:) > xc & effective(:));
if isempty(left_candidates) || isempty(right_candidates)
    return;
end
left_L = max(xc - x_grid(left_candidates(1)), 0);
right_L = max(x_grid(right_candidates(end)) - xc, 0);
left_L = max(left_L, min_half_width);
right_L = max(right_L, min_half_width);
if isfinite(max_half_width) && max_half_width > 0
    left_L = min(left_L, max_half_width);
    right_L = min(right_L, max_half_width);
end
ok = isfinite(left_L) && isfinite(right_L) && left_L > 0 && right_L > 0;
end


function [left_L, right_L, ok] = estimate_gradient_energy_span_step04_local( ...
    x, v, baseline, xc, energy_quantile, edge_margin_mm)
left_L = NaN;
right_L = NaN;
ok = false;
[x_grid, v_smooth] = build_denoised_static_profile_step04_local(x, v, baseline);
if numel(x_grid) < 20 || max(v_smooth) <= 0
    return;
end
g = abs(gradient(v_smooth, x_grid));
signal_gate = v_smooth >= prctile(v_smooth, 40);
energy = g(:) .* signal_gate(:) + 0.05 * max(g(:)) * normalize01_step04_local(v_smooth(:));
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


function [x_grid, v_smooth] = build_denoised_static_profile_step04_local(x, v, baseline)
finite = isfinite(x) & isfinite(v);
x = x(finite);
v_zero = max(v(finite) - baseline, 0);
x_grid = [];
v_smooth = [];
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
    x_grid = [];
    v_smooth = [];
    return;
end
v_fill = fillmissing(v_med, 'linear', 'EndValues', 'nearest');
span = min(41, 2 * floor((numel(v_fill) - 1) / 2) + 1);
if span >= 5
    v_smooth = smoothdata(v_fill, 'sgolay', span);
else
    v_smooth = v_fill;
end
end


function g_threshold = estimate_noise_aware_gradient_threshold_step04_local(g_abs, amp_norm, ratio)
g_abs = g_abs(:);
amp_norm = amp_norm(:);
finite = isfinite(g_abs) & isfinite(amp_norm);
g_abs = g_abs(finite);
amp_norm = amp_norm(finite);
if isempty(g_abs)
    g_threshold = Inf;
    return;
end
noise_gate = amp_norm <= 0.10;
signal_gate = amp_norm >= 0.20;
if nnz(noise_gate) >= 10
    g_noise = prctile(g_abs(noise_gate), 95);
else
    g_noise = prctile(g_abs, 10);
end
if nnz(signal_gate) >= 10
    g_signal = prctile(g_abs(signal_gate), 95);
else
    g_signal = prctile(g_abs, 95);
end
if ~isfinite(g_noise)
    g_noise = 0;
end
if ~isfinite(g_signal) || g_signal <= g_noise
    g_signal = max(g_abs);
end
g_threshold = g_noise + ratio * max(g_signal - g_noise, 0);
end


function y = normalize01_step04_local(x)
x = x(:);
xmin = min(x, [], 'omitnan');
xmax = max(x, [], 'omitnan');
if ~isfinite(xmin) || ~isfinite(xmax) || xmax <= xmin
    y = zeros(size(x));
else
    y = (x - xmin) ./ (xmax - xmin);
end
end


function xc = refine_center_by_sg_fit_step04_local(x, v, baseline, xc_initial)
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
        bounded_obj = @(p) obj_fun(min(max(p, lb), ub));
        opts = optimset('Display', 'off', 'MaxIter', 400, 'MaxFunEvals', 1200);
        p_opt = fminsearch(bounded_obj, p0, opts);
        p_opt = min(max(p_opt, lb), ub);
    end
    if isfinite(p_opt(4))
        xc = p_opt(4);
    else
        xc = xc_initial;
    end
catch
    xc = xc_initial;
end
end


function W = build_gradient_weight_from_wide_cloud_step04_local(x_query, x_wide, v_wide, baseline, weight_floor)
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


function [xc, info] = estimate_template_xcenter_local(x_abs, v)
% Estimate the pulse center in absolute OPRCenterStd coordinates. This
% mirrors the 20251222 low-speed rotating calibration route before the
% foundation SensorBlade template is built in the relative x = x_abs - xc
% coordinate.
x_abs = x_abs(:);
v = v(:);
finite = isfinite(x_abs) & isfinite(v);
x_abs = x_abs(finite);
v = v(finite);

info = struct( ...
    'baseline_for_center_v', NaN, ...
    'high_level_v', NaN, ...
    'high_point_count', 0, ...
    'reference_source', 'estimated_from_low_speed_point_cloud');

xc = NaN;
if numel(x_abs) < 5
    info.reference_source = 'insufficient_points_for_xcenter';
    return;
end

low_mask = v <= prctile(v, 30);
if nnz(low_mask) >= 5
    baseline = median(v(low_mask), 'omitnan');
else
    baseline = prctile(v, 10);
end
if ~isfinite(baseline)
    baseline = median(v, 'omitnan');
end
if ~isfinite(baseline)
    info.reference_source = 'invalid_baseline_for_xcenter';
    return;
end

v_zero = max(v - baseline, 0);
info.baseline_for_center_v = baseline;

if all(~isfinite(v_zero)) || max(v_zero, [], 'omitnan') <= 0
    [~, imax_raw] = max(v);
    xc = x_abs(imax_raw);
    info.high_point_count = 1;
    info.reference_source = 'raw_peak_fallback_xcenter';
    return;
end

high_level = prctile(v_zero, 85);
high_mask = v_zero >= high_level & v_zero > 0 & isfinite(x_abs);
info.high_level_v = high_level;
info.high_point_count = nnz(high_mask);

if nnz(high_mask) >= 5 && sum(v_zero(high_mask), 'omitnan') > eps
    xc = sum(x_abs(high_mask) .* v_zero(high_mask), 'omitnan') ./ ...
        sum(v_zero(high_mask), 'omitnan');
    info.reference_source = 'weighted_high_voltage_centroid_xcenter';
else
    [~, imax] = max(v_zero);
    xc = x_abs(imax);
    info.reference_source = 'voltage_peak_fallback_xcenter';
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


function G = build_template_grid_from_points_step04_local(cfg, x, v, point_weight, lap_id, x_range)
G = struct('ok', false, 'x_grid', [], 'v_grid_raw', [], 'v_smooth', [], ...
    'count_grid', [], 'lap_count_grid', [], 'v_std_grid', [], ...
    'weight_grid_raw', [], 'valid_bin_mask', []);
if numel(x_range) ~= 2 || any(~isfinite(x_range)) || x_range(2) <= x_range(1)
    return;
end
x = x(:);
v = v(:);
point_weight = point_weight(:);
lap_id = lap_id(:);
finite = isfinite(x) & isfinite(v);
x = x(finite);
v = v(finite);
point_weight = point_weight(finite);
lap_id = lap_id(finite);
if numel(x) < cfg.step04_min_template_points
    return;
end
edges = x_range(1):cfg.step04_grid_dx_mm:x_range(2);
if numel(edges) < 5
    edges = linspace(x_range(1), x_range(2), 20);
end
x_grid = 0.5 * (edges(1:end-1) + edges(2:end));
x_grid = x_grid(:);
bin = discretize(x, edges);
n_bins = numel(x_grid);
v_grid_raw = nan(n_bins, 1);
weight_grid_raw = nan(n_bins, 1);
v_std_grid = nan(n_bins, 1);
lap_count_grid = zeros(n_bins, 1);
count_grid = zeros(n_bins, 1);
for ib = 1:n_bins
    mask = bin == ib;
    count_grid(ib) = nnz(mask);
    if count_grid(ib) > 0
        v_grid_raw(ib) = median(v(mask), 'omitnan');
        v_std_grid(ib) = std(v(mask), 'omitnan');
        weight_grid_raw(ib) = mean(point_weight(mask), 'omitnan');
        laps_here = lap_id(mask);
        lap_count_grid(ib) = numel(unique(laps_here(isfinite(laps_here))));
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
G.ok = true;
G.x_grid = x_grid;
G.v_grid_raw = v_grid_raw;
G.v_smooth = smooth_template_local(v_interp, cfg.step04_smooth_span_bins);
G.count_grid = count_grid;
G.lap_count_grid = lap_count_grid;
G.v_std_grid = v_std_grid;
G.weight_grid_raw = weight_grid_raw;
G.valid_bin_mask = valid_bin_mask;
end


function Domain = determine_trusted_domain_local(cfg, x_grid, v_grid, dv_dx, count_grid, lap_count_grid, v_std_grid, valid_bin_mask, baseline, amplitude, max_abs_grad)
%DETERMINE_TRUSTED_DOMAIN_LOCAL Select the calibration-grade x-domain.
%
% The domain is not a raw coverage interval.  It must be a waveform-response
% core that is covered by enough low-speed laps and repeatable enough for
% later direct-template fitting.  Follow the 20250527 foundation / legacy
% GradientXRange030 idea: gradient identifies the informative span; amplitude
% is an auxiliary guard, not a substitute that can expand the trusted domain
% into low-slope pulse edges.
x_grid = x_grid(:);
v_grid = v_grid(:);
dv_dx = dv_dx(:);
count_grid = count_grid(:);
lap_count_grid = lap_count_grid(:);
v_std_grid = v_std_grid(:);
valid_bin_mask = valid_bin_mask(:);

Domain = make_empty_domain_diagnostic_step04_local(size(x_grid));
if isempty(x_grid) || nnz(valid_bin_mask) < 5 || numel(v_grid) ~= numel(x_grid)
    Domain.status = "invalid_insufficient_grid";
    return;
end

edge_count = max(5, round(0.08 * nnz(valid_bin_mask)));
valid_idx = find(valid_bin_mask & isfinite(v_grid));
if numel(valid_idx) >= 2 * edge_count
    edge_values = [v_grid(valid_idx(1:edge_count)); v_grid(valid_idx(end-edge_count+1:end))];
else
    edge_values = v_grid(valid_idx);
end
if ~isfinite(baseline)
    baseline = median(edge_values, 'omitnan');
end
y = v_grid - baseline;
if abs(min(y(valid_bin_mask), [], 'omitnan')) > abs(max(y(valid_bin_mask), [], 'omitnan'))
    response = -y;
    polarity = -1;
else
    response = y;
    polarity = 1;
end
response(response < 0) = 0;
peak_amp = max(response(valid_bin_mask), [], 'omitnan');
if ~isfinite(peak_amp) || peak_amp <= 0
    if isfinite(amplitude) && amplitude > 0
        peak_amp = amplitude;
    else
        Domain.status = "invalid_no_response";
        return;
    end
end

noise_sigma = 1.4826 * mad(edge_values - median(edge_values, 'omitnan'), 1);
if ~isfinite(noise_sigma)
    noise_sigma = 0;
end
amp_threshold = max(cfg.step04_domain_noise_sigma_factor * noise_sigma, ...
    cfg.step04_amplitude_min_ratio * peak_amp);
support_mask = valid_bin_mask & ...
    count_grid >= cfg.step04_domain_min_points_per_bin & ...
    lap_count_grid >= cfg.step04_domain_min_laps_per_bin;
amplitude_mask = valid_bin_mask & response >= amp_threshold;
gradient_mask = false(size(x_grid));
if isfinite(max_abs_grad) && max_abs_grad > 0
    gradient_mask = valid_bin_mask & abs(dv_dx) >= cfg.step04_gradient_min_ratio * max_abs_grad;
end
if all(~isfinite(v_std_grid))
    repeat_mask = true(size(x_grid));
else
    repeat_limit = max(cfg.step04_domain_std_abs_limit_v, cfg.step04_domain_repeat_ratio * peak_amp);
    repeat_mask = valid_bin_mask & (v_std_grid <= repeat_limit | ~isfinite(v_std_grid));
end

if strcmpi(cfg.step04_domain_center_mode, 'peak')
    high_mask = valid_bin_mask & response >= 0.85 * peak_amp;
    if nnz(high_mask) >= 3 && sum(response(high_mask), 'omitnan') > eps
        x_center = sum(x_grid(high_mask) .* response(high_mask), 'omitnan') ./ ...
            sum(response(high_mask), 'omitnan');
    else
        response_for_peak = response;
        response_for_peak(~valid_bin_mask | ~isfinite(response_for_peak)) = -Inf;
        [~, i_peak] = max(response_for_peak);
        x_center = x_grid(i_peak);
    end
else
    x_center = 0;
end
if ~isfinite(x_center)
    Domain.status = "invalid_no_center";
    return;
end

% Gradient defines the informative core; amplitude only rejects noise-floor
% gradient artifacts.  Coverage/repeatability remain mandatory.
response_region = gradient_mask;
if isfinite(cfg.step04_amplitude_min_ratio) && cfg.step04_amplitude_min_ratio > 0
    response_region = response_region & amplitude_mask;
end
trusted_mask = support_mask & repeat_mask & response_region;

left_candidates = find(x_grid < x_center & trusted_mask);
right_candidates = find(x_grid > x_center & trusted_mask);

% Fallback only when one side has no informative gradient bins.  This keeps
% the route runnable on weak pulses, but still requires coverage/repeatability
% and a stricter amplitude floor.
if isempty(left_candidates) || isempty(right_candidates)
    amp_support_threshold = max(cfg.step04_amplitude_min_ratio * peak_amp, ...
        0.05 * peak_amp);
    amp_support = support_mask & repeat_mask & response >= amp_support_threshold;
    left_candidates = find(x_grid < x_center & amp_support);
    right_candidates = find(x_grid > x_center & amp_support);
end

if isempty(left_candidates) || isempty(right_candidates)
    Domain.status = "invalid_no_trusted_region";
    Domain.support_mask = support_mask;
    Domain.amplitude_mask = amplitude_mask;
    Domain.gradient_mask = gradient_mask;
    Domain.repeat_mask = repeat_mask;
    Domain.trusted_mask = trusted_mask;
    Domain.response_region_mask = response_region;
    return;
end

left_L = max(x_center - x_grid(left_candidates(1)), 0);
right_L = max(x_grid(right_candidates(end)) - x_center, 0);
left_L = max(left_L, cfg.step04_min_half_width_mm);
right_L = max(right_L, cfg.step04_min_half_width_mm);
if isfinite(cfg.step04_max_half_width_mm) && cfg.step04_max_half_width_mm > 0
    left_L = min(left_L, cfg.step04_max_half_width_mm);
    right_L = min(right_L, cfg.step04_max_half_width_mm);
end

coverage_left = min(x_grid(support_mask), [], 'omitnan');
coverage_right = max(x_grid(support_mask), [], 'omitnan');
left_x = max(x_center - left_L + cfg.step04_domain_boundary_margin_mm, coverage_left);
right_x = min(x_center + right_L - cfg.step04_domain_boundary_margin_mm, coverage_right);
if right_x <= left_x || (right_x - left_x) < cfg.step04_min_domain_width_mm
    Domain.status = "invalid_domain_too_narrow";
    Domain.support_mask = support_mask;
    Domain.amplitude_mask = amplitude_mask;
    Domain.gradient_mask = gradient_mask;
    Domain.repeat_mask = repeat_mask;
    Domain.trusted_mask = trusted_mask;
    return;
end

x_support = domain_from_mask_step04_local(x_grid, support_mask, x_center);
x_response = domain_from_mask_step04_local(x_grid, response_region & valid_bin_mask, x_center);
x_query_safe = [left_x + cfg.step04_query_guard_mm, right_x - cfg.step04_query_guard_mm];
if x_query_safe(2) <= x_query_safe(1)
    x_query_safe = [NaN NaN];
end

Domain.x_domain = [left_x, right_x];
Domain.x_query_safe_domain = x_query_safe;
Domain.x_support_domain_mm = x_support;
Domain.x_response_domain_mm = x_response;
Domain.domain_mask = x_grid >= left_x & x_grid <= right_x & valid_bin_mask;
Domain.support_mask = support_mask;
Domain.amplitude_mask = amplitude_mask;
Domain.gradient_mask = gradient_mask;
Domain.repeat_mask = repeat_mask;
Domain.trusted_mask = trusted_mask;
Domain.response_region_mask = response_region;
Domain.x_center = x_center;
Domain.baseline = baseline;
Domain.polarity = polarity;
Domain.peak_amp = peak_amp;
Domain.noise_sigma = noise_sigma;
Domain.amp_threshold = amp_threshold;
Domain.min_points_per_bin = cfg.step04_domain_min_points_per_bin;
Domain.min_laps_per_bin = cfg.step04_domain_min_laps_per_bin;
Domain.repeat_limit = max(cfg.step04_domain_std_abs_limit_v, cfg.step04_domain_repeat_ratio * peak_amp);
Domain.status = "ok";
end


function Domain = make_empty_domain_diagnostic_step04_local(sz)
Domain = struct();
Domain.x_domain = [NaN NaN];
Domain.x_query_safe_domain = [NaN NaN];
Domain.x_support_domain_mm = [NaN NaN];
Domain.x_response_domain_mm = [NaN NaN];
Domain.domain_mask = false(sz);
Domain.support_mask = false(sz);
Domain.amplitude_mask = false(sz);
Domain.gradient_mask = false(sz);
Domain.repeat_mask = false(sz);
Domain.trusted_mask = false(sz);
Domain.response_region_mask = false(sz);
Domain.x_center = NaN;
Domain.baseline = NaN;
Domain.polarity = NaN;
Domain.peak_amp = NaN;
Domain.noise_sigma = NaN;
Domain.amp_threshold = NaN;
Domain.min_points_per_bin = NaN;
Domain.min_laps_per_bin = NaN;
Domain.repeat_limit = NaN;
Domain.status = "empty";
end


function x_domain = domain_from_mask_step04_local(x, mask, x_center)
x_domain = [NaN NaN];
idx = find(mask);
if isempty(idx)
    return;
end
[~, k] = min(abs(x(idx) - x_center));
i_center = idx(k);
i_left = i_center;
while i_left > 1 && mask(i_left - 1)
    i_left = i_left - 1;
end
i_right = i_center;
while i_right < numel(x) && mask(i_right + 1)
    i_right = i_right + 1;
end
x_domain = [x(i_left), x(i_right)];
end

function weight_grid = build_template_weight_grid_local(count_grid, valid_bin_mask, domain_mask, dv_dx, max_abs_grad, selected_weight_grid)
count_grid = count_grid(:);
valid_bin_mask = valid_bin_mask(:);
domain_mask = domain_mask(:);
dv_dx = dv_dx(:);
selected_weight_grid = selected_weight_grid(:);
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
if numel(selected_weight_grid) ~= numel(count_grid)
    selected_weight_grid = ones(size(count_grid));
end
selected_weight_grid(~isfinite(selected_weight_grid)) = 1;
selected_weight_grid = max(selected_weight_grid, 0.05);
weight_grid(valid_bin_mask) = coverage_w(valid_bin_mask) .* ...
    gradient_w(valid_bin_mask) .* selected_weight_grid(valid_bin_mask);
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
metadata.schema_version = Template.SchemaVersion;
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
metadata.opr_timing_method = Template.OPR_Timing_Method;
metadata.opr_events_per_revolution = Template.OPR_Events_Per_Revolution;
metadata.opr_events_per_revolution_source = Template.OPR_Events_Per_Revolution_Source;
metadata.probe_arrival_method = Template.Probe_Arrival_Method;
metadata.blade_id_definition = Template.Blade_ID_Definition;
metadata.standard_angle_source = Template.Standard_Angle_Source;
metadata.standard_angle_field = Template.Standard_Angle_Field;
metadata.final_template_point_policy = Template.Final_Template_Point_Policy;
metadata.baseline_mode = cfg.step04_baseline_mode;
metadata.stable_window_selection_policy = cfg.step04_stable_window_selection_policy;
metadata.center_seed_source = cfg.step04_center_seed_source;
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
                'Color', [0.82 0.82 0.82], 'MarkerSize', 2, ...
                'DisplayName', 'wide waveform points');
        end

        yl_now = ylim(ax);
        if all(isfinite(entry.x_domain))
            shade_domain_local(ax, entry.x_domain, yl_now, colors(blade_id, :));
        end
        if ~isempty(entry.x_used_points_preview)
            plot(ax, entry.x_used_points_preview, entry.v_used_points_preview, '.', ...
                'Color', [0.28 0.28 0.28], 'MarkerSize', 2, ...
                'DisplayName', 'used calibration points');
        end
        plot(ax, entry.x_grid, entry.v_grid_raw, '-', 'LineWidth', 0.8, ...
            'Color', [0.35 0.35 0.35], 'DisplayName', 'final binned median');
        plot(ax, entry.x_grid, entry.v_grid, '-', 'LineWidth', 1.8, ...
            'Color', colors(blade_id, :), 'DisplayName', 'final smoothed template');
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


function value = parse_bool_env_local(name, default_value)
txt = lower(strtrim(getenv(name)));
if isempty(txt)
    value = default_value;
elseif ismember(txt, {'1', 'true', 'yes', 'on'})
    value = true;
elseif ismember(txt, {'0', 'false', 'no', 'off'})
    value = false;
else
    warning('Step04:InvalidBooleanEnv', '%s=%s is not logical; using default.', name, txt);
    value = default_value;
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


function value = get_numeric_field_local(S, name, default_value)
value = default_value;
if isstruct(S) && isfield(S, name)
    candidate = S.(name);
    if isnumeric(candidate) && isscalar(candidate) && isfinite(candidate)
        value = candidate;
    end
end
end


function value = get_field_or_default_local(S, field_name, default_value)
value = default_value;
if isstruct(S) && isfield(S, field_name)
    candidate = S.(field_name);
    if ~isempty(candidate)
        value = candidate;
    end
end
end


function validate_low_speed_template_contract_local(Template, cfg)
required_top = {'SchemaVersion', 'Dataset', 'Case_Name', 'Sensor_IDs', ...
    'Blades_Num', 'R_Tip_mm', 'OPR_Events_Per_Revolution', ...
    'Standard_Angle_Source', 'Standard_Angle_Field', 'Standard_Relative_Angles', ...
    'Final_Template_Point_Policy', 'SensorBlade', 'Summary_Table', 'Metadata'};
for i = 1:numel(required_top)
    if ~isfield(Template, required_top{i})
        error('Step04:TemplateContract', 'Template missing required top-level field: %s', required_top{i});
    end
end
if ~strcmp(char(string(Template.SchemaVersion)), 'LowSpeedOPRCenterStdTemplate/v1')
    error('Step04:TemplateContract', 'Unexpected Template.SchemaVersion: %s', char(string(Template.SchemaVersion)));
end
if ~strcmp(char(string(Template.Final_Template_Point_Policy)), 'trusted_domain_recalibration')
    error('Step04:TemplateContract', 'Template.Final_Template_Point_Policy must be trusted_domain_recalibration.');
end
if ~isnumeric(Template.OPR_Events_Per_Revolution) || ...
        ~isscalar(Template.OPR_Events_Per_Revolution) || ...
        Template.OPR_Events_Per_Revolution <= 0
    error('Step04:TemplateContract', 'Template.OPR_Events_Per_Revolution must be a positive scalar.');
end

required_entry = {'sensor_id', 'blade_id', 'x_grid', 'v_grid', 'dv_dx', ...
    'weight_grid', 'count_grid', 'lap_count_grid', 'v_std_grid', ...
    'valid_grid_mask', 'x_collect_domain', 'x_support_domain', ...
    'x_response_domain', 'x_domain', 'x_query_safe_domain', ...
    'xc', 'xc_mm', 'xc_detected', 'xc_detected_mm', ...
    'source_point_policy', 'final_template_point_policy', ...
    'quality_status', 'final_point_count', 'trusted_point_count'};

expected_count = numel(cfg.sensor_ids) * cfg.blades_num;
if numel(Template.SensorBlade) ~= expected_count
    error('Step04:TemplateContract', ...
        'Template.SensorBlade count is %d, expected %d.', numel(Template.SensorBlade), expected_count);
end

for i = 1:numel(Template.SensorBlade)
    entry = Template.SensorBlade(i);
    for j = 1:numel(required_entry)
        if ~isfield(entry, required_entry{j})
            error('Step04:TemplateContract', ...
                'SensorBlade(%d) missing required field: %s', i, required_entry{j});
        end
    end
    sid = entry.sensor_id;
    blade_id = entry.blade_id;
    if ~ismember(sid, cfg.sensor_ids) || blade_id < 1 || blade_id > cfg.blades_num
        error('Step04:TemplateContract', ...
            'SensorBlade(%d) has invalid sensor/blade id: S%d B%d.', i, sid, blade_id);
    end
    if ~strcmp(char(string(entry.final_template_point_policy)), char(string(entry.source_point_policy)))
        error('Step04:TemplateContract', ...
            'SensorBlade S%d B%d final_template_point_policy must match source_point_policy.', sid, blade_id);
    end
    quality_status = lower(strtrim(string(entry.quality_status)));
    if ~ismember(quality_status, ["good", "usable"])
        continue;
    end
    x = entry.x_grid(:);
    v = entry.v_grid(:);
    dv = entry.dv_dx(:);
    if numel(x) < 5 || numel(v) ~= numel(x) || numel(dv) ~= numel(x)
        error('Step04:TemplateContract', ...
            'SensorBlade S%d B%d has inconsistent x/v/dv grid lengths.', sid, blade_id);
    end
    finite_curve = isfinite(x) & isfinite(v);
    if nnz(finite_curve) < 5
        error('Step04:TemplateContract', ...
            'SensorBlade S%d B%d has too few finite template points.', sid, blade_id);
    end
    if any(diff(x(finite_curve)) <= 0)
        error('Step04:TemplateContract', ...
            'SensorBlade S%d B%d x_grid must be strictly increasing over finite points.', sid, blade_id);
    end
    x_domain = entry.x_domain(:).';
    x_query = entry.x_query_safe_domain(:).';
    if numel(x_domain) < 2 || any(~isfinite(x_domain(1:2))) || x_domain(2) <= x_domain(1)
        error('Step04:TemplateContract', ...
            'SensorBlade S%d B%d has invalid x_domain.', sid, blade_id);
    end
    if numel(x_query) >= 2 && all(isfinite(x_query(1:2)))
        if x_query(1) < x_domain(1) || x_query(2) > x_domain(2) || x_query(2) <= x_query(1)
            error('Step04:TemplateContract', ...
                'SensorBlade S%d B%d has invalid x_query_safe_domain.', sid, blade_id);
        end
    end
end

fprintf('>>> [Step04] Template contract OK: %s, %d sensor-blade entries.\n', ...
    char(string(Template.SchemaVersion)), numel(Template.SensorBlade));
end


function ang = wrap_to_180_local(ang)
ang = mod(ang + 180, 360) - 180;
end
