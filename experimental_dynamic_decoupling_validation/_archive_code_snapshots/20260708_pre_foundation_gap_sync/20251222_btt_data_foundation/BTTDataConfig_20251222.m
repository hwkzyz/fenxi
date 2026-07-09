function cfg = BTTDataConfig_20251222()
%BTTDATACONFIG_20251222 Data facts for the 20251222 BTT foundation route.
%
% This config is intentionally limited to stable data/project facts and
% output locations. Step-specific algorithm choices should stay at the top
% of each Step script.

route_dir = fileparts(mfilename('fullpath'));
validation_root = fileparts(route_dir);
rotating_route_dir = fullfile(validation_root, '20251222_low_speed_rotating_calibration');
gap_prior_dir = fullfile(validation_root, '20251222_low_speed_gap_prior_decoupling');

cfg = struct();
cfg.dataset = '20251222';
cfg.route_dir = route_dir;
cfg.validation_root = validation_root;

cfg.dataset_root = 'E:\试验数据\20251222\传感器数据';
cfg.low_speed_case = '1000rpm无振动';
cfg.dynamic_cases = {'1000_2500_3500'};
cfg.strain_root = 'E:\试验数据\20251222\应变片数据';

cfg.output_root = fullfile(route_dir, 'output');
cfg.step01_output_dir = fullfile(cfg.output_root, 'step01_low_speed_reference');
cfg.step02_output_dir = fullfile(cfg.output_root, 'step02_dynamic_btt');
cfg.step03_output_dir = fullfile(cfg.output_root, 'step03_observation_bundle');
cfg.step04_output_dir = fullfile(cfg.output_root, 'step04_low_speed_template');
cfg.step05_direct_template_output_dir = fullfile(cfg.output_root, 'step05_single_sync_direct_template');
cfg.optional_diagnostics_output_dir = fullfile(cfg.output_root, 'optional_diagnostics');
cfg.optional_template_guided_feature_baseline_output_dir = fullfile( ...
    cfg.optional_diagnostics_output_dir, 'template_guided_feature_baseline');
cfg.figure_root = fullfile(cfg.output_root, 'figures');
cfg.step01_figure_dir = fullfile(cfg.figure_root, 'step01_low_speed_reference');
cfg.step02_figure_dir = fullfile(cfg.figure_root, 'step02_dynamic_btt');
cfg.step03_figure_dir = fullfile(cfg.figure_root, 'step03_observation_bundle');
cfg.step04_figure_dir = fullfile(cfg.figure_root, 'step04_low_speed_template');
cfg.step05_direct_template_figure_dir = fullfile(cfg.figure_root, 'step05_single_sync_direct_template');
cfg.optional_diagnostics_figure_dir = fullfile(cfg.figure_root, 'optional_diagnostics');
cfg.optional_template_guided_feature_baseline_figure_dir = fullfile( ...
    cfg.optional_diagnostics_figure_dir, 'template_guided_feature_baseline');
cfg.save_figures = true;

% Physical raw-channel map from the 20251222 acquisition notes:
%   CH1 = optical fiber probe 3
%   CH2 = eddy-current probe 3
%   CH3 = optical fiber probe 4
%   CH4 = OPR
%   CH5 = capacitance probe 2
%   CH6 = strain gauge
%   CH7 = eddy-current probe 1
% The current foundation Step01-Step03 products were built for CH1-CH3.
cfg.candidate_sensor_ids = [1 2 3];
cfg.sensor_ids = [1 2 3];
cfg.optical_fiber_ids = [1 3];
cfg.fiber_sensor_ids = cfg.optical_fiber_ids;
cfg.btt_reference_sensor_type = 'optical_fiber';
cfg.capacitance_ids = [5];
cfg.eddy_current_ids = [2 7];
cfg.strain_channel_ids = [6];
cfg.unused_ids = [8];
cfg.opr_id = 4;

cfg.blades_num = 6;
cfg.opr_timing_reference = 'rising_edge';
cfg.opr_reference_time_column = 'opr_start_time_s';
cfg.opr_timing_method = 'threshold_rising_edge';
cfg.standard_angle_reference_field = 'Standard_Relative_Angles_StartEdge';
cfg.standard_angle_reference_name = 'opr_threshold_start_edge';
cfg.step04_standard_angle_source = 'startedge';
cfg.opr_events_per_revolution = 6;
cfg.opr_events_per_revolution_policy = 'fixed_known_from_experiment';
cfg.step04_opr_events_per_revolution = 6;
cfg.step04_opr_events_per_revolution_policy = 'fixed_known_from_experiment';
cfg.sample_rate_hz = 5e6;
cfg.pinlv = cfg.sample_rate_hz;
cfg.r_tip_mm = 62.0;

cfg.initial_trim_points = 20000;
cfg.gap_points = 1000;
cfg.opr_threshold = 1.5;
cfg.sensor_threshold_default = 0.5;
cfg.reference_sensor_id = 1;
cfg.reference_search_revs = 3;
cfg.max_alignment_search_starts = 24;
cfg.max_laps_process = 5000;

% Preferred first-round bootstrap sources. Step01 uses these only to adapt
% existing low-speed reference products into the foundation interface.
cfg.source = struct();
cfg.source.rotating_route_dir = rotating_route_dir;
cfg.source.gap_prior_dir = gap_prior_dir;
cfg.source.newflow_sensor_config_file = fullfile(rotating_route_dir, ...
    'output', 'new_flow', '00_low_speed_sensor_config', 'Sensor_Config_20251222.mat');
cfg.source.newflow_sensor_config_summary_file = fullfile(rotating_route_dir, ...
    'output', 'new_flow', '00_low_speed_sensor_config', 'Sensor_Config_20251222_Summary.csv');
cfg.source.newflow_opr_center_angle_file = fullfile(rotating_route_dir, ...
    'output', 'new_flow', '00_low_speed_sensor_config', 'Sensor_Config_20251222_OPRCenterAngles.csv');
cfg.source.legacy_sensor_config_file = fullfile(gap_prior_dir, 'legacy', 'output', ...
    [cfg.low_speed_case, '_reference'], 'Sensor_Config_20251222.mat');
cfg.source.legacy_sensor_config_summary_file = fullfile(gap_prior_dir, 'legacy', 'output', ...
    [cfg.low_speed_case, '_reference'], 'Sensor_Config_Summary_20251222.csv');
cfg.source.newflow_high_speed_peak_cache_dir = fullfile(rotating_route_dir, ...
    'output', 'new_flow', '04_high_speed_peak_cache');
cfg.source.newflow_high_speed_numbering_dir = fullfile(rotating_route_dir, ...
    'output', 'new_flow', '04_high_speed_numbering');
cfg.source.legacy_dynamic_output_root = fullfile(gap_prior_dir, 'legacy', 'output');
cfg.source.legacy_dynamic_output_dir = fullfile(cfg.source.legacy_dynamic_output_root, cfg.dynamic_cases{1});

cfg.sensor_thresholds = containers.Map('KeyType', 'double', 'ValueType', 'double');
for sid = cfg.candidate_sensor_ids
    cfg.sensor_thresholds(sid) = cfg.sensor_threshold_default;
end
end
