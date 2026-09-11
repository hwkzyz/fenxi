function cfg = BTTDataConfig_20241106()
%BTTDATACONFIG_20241106 Method-neutral BTT data-foundation config.
%
% This config is limited to data facts and extraction/reference settings.
% Direct-template, gap-library, eta, dx, and EO-search parameters belong in
% later method scripts, not in this data-foundation layer.

route_dir = fileparts(mfilename('fullpath'));
packageCfg = Config_20241106();

cfg = struct();
cfg.dataset = '20241106';
cfg.route_dir = route_dir;

cfg.dataset_root = fullfile('D:\', '博士-国科', '试验台数据', '新试验', '20241106_2');
cfg.low_speed_case = '900';
cfg.dynamic_cases = {'3000_3150'};
cfg.dataset_root = packageCfg.paths.rawDatasetRoot;

cfg.output_root = fullfile(route_dir, 'output');
cfg.step01_output_dir = fullfile(cfg.output_root, 'step01_low_speed_reference');
cfg.step02_output_dir = fullfile(cfg.output_root, 'step02_dynamic_btt');
cfg.step03_output_dir = fullfile(cfg.output_root, 'step03_observation_bundle');
cfg.step04_output_dir = fullfile(cfg.output_root, 'step04_low_speed_template');
cfg.step05_output_dir = fullfile(cfg.output_root, 'step05_single_sync_direct_template');
cfg.step06a_output_dir = fullfile(cfg.output_root, 'step06a_strain_rpm_resonance_evidence');
cfg.step06b_output_dir = fullfile(cfg.output_root, 'step06b_step05_strain_tarc_validation');
cfg.figure_root = fullfile(cfg.output_root, 'figures');
cfg.step01_figure_dir = fullfile(cfg.figure_root, 'step01_low_speed_reference');
cfg.step02_figure_dir = fullfile(cfg.figure_root, 'step02_dynamic_btt');
cfg.step03_figure_dir = fullfile(cfg.figure_root, 'step03_observation_bundle');
cfg.step04_figure_dir = fullfile(cfg.figure_root, 'step04_low_speed_template');
cfg.step05_figure_dir = fullfile(cfg.figure_root, 'step05_single_sync_direct_template');
cfg.step06a_figure_dir = fullfile(cfg.figure_root, 'step06a_strain_rpm_resonance_evidence');
cfg.step06b_figure_dir = fullfile(cfg.figure_root, 'step06b_step05_strain_tarc_validation');
cfg.output_root = packageCfg.paths.preparedFoundation;
cfg.step01_output_dir = fullfile(cfg.output_root, 'step01_low_speed_reference');
cfg.step02_output_dir = fullfile(cfg.output_root, 'step02_dynamic_btt');
cfg.step03_output_dir = fullfile(cfg.output_root, 'step03_observation_bundle');
cfg.step04_output_dir = fullfile(cfg.output_root, 'step04_low_speed_template');
cfg.step05_output_dir = packageCfg.paths.foundationResults;
cfg.step06a_output_dir = fullfile(cfg.output_root, 'step06a_strain_rpm_resonance_evidence');
cfg.step06b_output_dir = fullfile(cfg.output_root, 'step06b_step05_strain_tarc_validation');
cfg.figure_root = fullfile(packageCfg.paths.figures, 'foundation');
cfg.step01_figure_dir = fullfile(cfg.figure_root, 'step01_low_speed_reference');
cfg.step02_figure_dir = fullfile(cfg.figure_root, 'step02_dynamic_btt');
cfg.step03_figure_dir = fullfile(cfg.figure_root, 'step03_observation_bundle');
cfg.step04_figure_dir = fullfile(cfg.figure_root, 'step04_low_speed_template');
cfg.step05_figure_dir = fullfile(cfg.figure_root, 'step05_single_sync_direct_template');
cfg.step06a_figure_dir = fullfile(cfg.figure_root, 'step06a_strain_rpm_resonance_evidence');
cfg.step06b_figure_dir = fullfile(cfg.figure_root, 'step06b_step05_strain_tarc_validation');
cfg.save_figures = packageCfg.run.saveFigures;

cfg.candidate_sensor_ids = [2 3 4 5 6 7];
cfg.sensor_ids = [2 3 4 5 6 7];
cfg.low_speed_extract_sensor_ids = cfg.sensor_ids;
cfg.optional_sensor_ids = [4 6];
cfg.optical_fiber_ids = [4 6];
cfg.fiber_sensor_ids = cfg.optical_fiber_ids;
cfg.capacitance_ids = [5 7];
cfg.eddy_current_ids = [2 3];
cfg.unused_ids = [];
cfg.opr_id = 1;

cfg.blades_num = 6;
cfg.opr_pulses_per_rev = 1;
cfg.step04_opr_events_per_revolution = 1;
cfg.step04_opr_events_per_revolution_policy = 'fixed_known_from_experiment';
cfg.step04_low_speed_nominal_rpm = 900;
cfg.sample_rate_hz = 5e6;
cfg.pinlv = cfg.sample_rate_hz;
cfg.r_tip_mm = 65.0;

cfg.initial_trim_points = 20000;
cfg.gap_points = 10000;
cfg.opr_threshold = 2.0;
cfg.reference_sensor_id = 5;
cfg.reference_search_revs = 3;
cfg.max_alignment_search_starts = 24;
cfg.max_laps_process = 5000;

cfg.sensor_threshold_ids = [1 2 3 4 5 6 7];
cfg.sensor_threshold_values = [2.0 -0.2 -0.2 2.0 1.0 2.0 1.0];
cfg.sensor_threshold_default = 0.5;
cfg.sensor_thresholds = containers.Map('KeyType', 'double', 'ValueType', 'double');
for i = 1:numel(cfg.sensor_threshold_ids)
    cfg.sensor_thresholds(cfg.sensor_threshold_ids(i)) = cfg.sensor_threshold_values(i);
end

cfg.sensor_peak_fit_degree = 7;
cfg.sensor_peak_fit_refine_points = 200;

cfg.opr_timing_method = 'multi_threshold_center';
cfg.opr_center_correction_method = 'multi_threshold_center';
cfg.opr_center_level_ratios = [0.30 0.40 0.50 0.60 0.70];
cfg.opr_center_smooth_span_points = 16;
cfg.opr_center_pad_fraction = 0.25;
cfg.opr_center_min_pad_points = 8;
cfg.probe_arrival_method = 'polynomial_peak_fit';
cfg.probe_arrival_baseline_mode = 'edge_median';
cfg.probe_arrival_integral_signal = 'baseline_corrected';
cfg.probe_smooth_span_points = 21;
cfg.time_index_mode = 'local';
cfg.step02_gap_points = cfg.gap_points;
cfg.step02_segment_pad_fraction = 0.30;
cfg.step02_segment_pad_min_points = 4;
cfg.fingerprint_min_corr = 0.85;
cfg.fingerprint_min_corr_gap = 0.05;
cfg.fingerprint_min_shift_consistency = 0.60;
cfg.step02_fill_missing_blade_ids = false;
cfg.step02_revolution_vote_min_corr = -inf;
cfg.step02_progress_rev_interval = 200;
cfg.step03_target_blades = 1:cfg.blades_num;
cfg.step03_default_time_window_s = [];
cfg.step03_outlier_abs_mm = 20;
cfg.step03_histogram_bin_count = 120;
cfg.step03_plot_max_points_per_series = 6000;
cfg.step03_consistency_policy = 'warn';
cfg.step03_require_jiluopr_three_columns = true;
cfg.step03_use_current_opr_for_speed = true;
cfg.step03_use_robust_outlier_filter = false;
cfg.step03_robust_mad_k = 6.0;
cfg.step03_min_points_for_robust_filter = 12;
cfg.step03_save_legacy_vib_to_step02 = true;
cfg.step03_save_vib_to_step03 = true;
cfg.step03_standard_angle_preference = 'opr_center';
cfg.step03_wrap_angle_residual = true;
cfg.step04_standard_angle_source = 'opr_center';
cfg.step04_pulse_pad_fraction = 0.30;
cfg.step04_pulse_min_pad_points = 40;
cfg.step04_max_points_per_pulse = 450;
cfg.step04_max_pulses_per_sensor = inf;
cfg.step04_max_laps_per_sensor = cfg.max_laps_process;
cfg.step04_x_collect_abs_limit_mm = 12.0;
cfg.step04_grid_dx_mm = 0.02;
cfg.step04_min_bin_count = 5;
cfg.step04_smooth_span_bins = 9;
cfg.step04_gradient_min_ratio = 0.30;
cfg.step04_amplitude_min_ratio = 0.02;
cfg.step04_min_half_width_mm = 2.5;
cfg.step04_max_half_width_mm = 4.2;
cfg.step04_domain_center_mode = 'zero';
cfg.step04_domain_min_points_per_bin = 10;
cfg.step04_domain_min_laps_per_bin = 10;
cfg.step04_query_safe_margin_mm = 0.20;
cfg.step04_query_guard_mm = cfg.step04_query_safe_margin_mm;
cfg.step04_min_template_points = 1500;
cfg.step04_min_template_pulses = 20;
cfg.step04_min_domain_width_mm = 2.0 * cfg.step04_min_half_width_mm;
cfg.step04_plot_max_templates_per_figure = 6;
cfg.step04_initial_gap_points = cfg.gap_points;
cfg.step04_opr_center_level_ratios = cfg.opr_center_level_ratios;
% 20241106 legacy NewFlow main identification uses B4/S257 by default:
%   Step06_RunIdentificationByBlade_20241106: S06.blades = 4,
%   S06.analysisSensors = [2 5 7], S06.outputLabel = 'B4_only'.
% The identification windowing also follows that script: 20 target laps
% are evaluated as 3-lap windows with a 1-lap sliding step.
% Keep the all-blade list only as an explicit scan/audit option.
cfg.step05_target_blades = 4;
cfg.step04_target_blades = cfg.step05_target_blades;
cfg.step05_scan_blades = 1:cfg.blades_num;
cfg.step05_analysis_sensors = [2 5 7];
cfg.step06k_sensor_ids = cfg.optical_fiber_ids;
cfg.step06k_output_dir = fullfile(cfg.output_root, ...
    'step06k_fe_core_strain_btt_reference_fiber46');
cfg.step06k_figure_dir = fullfile(cfg.figure_root, ...
    'step06k_fe_core_strain_btt_reference_fiber46');
cfg.step05_analysis_start_time_s = 75.0;
cfg.step05_target_laps = 20;
cfg.step05_window_laps = 3;
cfg.step05_sliding_step_laps = 1;
cfg.step05_freq_search_hz = [300 1000];
cfg.step05_sensor_eta_limit_mm = 0.20;
cfg.step05_sensor_eta_reg_weight_v_per_mm = 0.02;
cfg.step05_dynamic_window_mode = 'legacy_row_bounds';
cfg.step05_pulse_pad_sec = 0;
cfg.step05_core_mask_mode = 'legacy_base_query_safe';
cfg.step05_phase_safe_margin_mm = 0.02;
cfg.step05_template_source = 'step04';
cfg.step05_legacy_template_library_file = '';
cfg.step05_timing_source = 'step02';
cfg.step05_legacy_pulse_dir = '';
cfg.strain_root = cfg.dataset_root;
cfg.step06_strain_files = {'AI1-01_20241106170148.mat', 'AI1-02_20241106170148.mat', ...
    'AI1-03_20241106170148.mat', 'AI1-04_20241106170148.mat'};
cfg.step06_strain_labels = {'AI1-01', 'AI1-02', 'AI1-03', 'AI1-04'};
cfg.step06_default_strain_channel = 'AI1-04';
cfg.step06_default_strain_file = fullfile(cfg.strain_root, 'AI1-04_20241106170148.mat');
cfg.step06_strain_sample_rate_hz = 10000;
cfg.step06_strain_to_btt_offset_s = -102.6;
cfg.step06_align_order_candidates = 6:18;
cfg.step06_fft_freq_band_hz = [0 1000];
cfg.step06_raw_fft_peak_band_hz = [20 800];
cfg.step06_target_blade = cfg.step05_target_blades(1);
cfg.standard_angle_reference = 'opr_pulse_center';
cfg.low_speed_blade1_rule = 'reference_sensor_opr_revolution_cyclic_match';
cfg.standard_angle_value = 'mean';
cfg.fingerprint_quality_policy = 'warn';
end
