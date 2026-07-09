function cfg = build_single_sync_experiment_config_20250527()
%BUILD_SINGLE_SYNC_EXPERIMENT_CONFIG_20250527
% Central config for the rewritten 20250527 single-sync SG pipeline.

script_dir = fileparts(mfilename('fullpath'));
base_cfg = Get_20250527_BTT_Config();

cfg = struct();
cfg.project_dir = script_dir;
cfg.base_cfg = base_cfg;

cfg.reference_case_name = base_cfg.low_speed_case;
cfg.dynamic_case_name = base_cfg.dynamic_cases{1};
cfg.static_data_dir = fullfile(base_cfg.dataset_root, cfg.reference_case_name);
cfg.dynamic_data_dir = fullfile(base_cfg.dataset_root, cfg.dynamic_case_name);

cfg.target_blade = 1;
cfg.target_blades = 1:base_cfg.blades_num;
cfg.num_blades = base_cfg.blades_num;
cfg.opr_pulses_per_rev = cfg.num_blades;
cfg.analysis_sensors = base_cfg.sensor_ids(:).';
cfg.analysis_start_time = 1.5;
cfg.target_laps = 20;
cfg.analysis_win_size = 3;
cfg.sliding_step = 1;
cfg.resonance_lap_range = [];
cfg.resonance_window_ids = [];

cfg.pinlv = base_cfg.pinlv;
cfg.r_tip_mm = base_cfg.r_tip_mm;
cfg.opr_channel = base_cfg.opr_id;
cfg.opr_threshold = base_cfg.opr_threshold;
cfg.gap_threshold = base_cfg.gap_points;
cfg.poly_order = 5;
cfg.pulse_window_sec = 6e-4;
cfg.pulse_pad_sec = 2 / cfg.pinlv;
cfg.dynamic_window_mode = 'peak_centered_fixed';
cfg.baseline_pad_points = 4000;
cfg.weight_floor = 0.05;
cfg.static_fit_max_points = 250000;
cfg.static_wide_max_points = 250000;
cfg.static_calibration_start_time = 1.0;
cfg.static_target_laps = 30;
cfg.static_calibration_time_range = [];

cfg.freq_search_hz = [500, 700];
cfg.reference_freq_hz = 580.0;
cfg.reference_eo = [];
cfg.freq_tolerance_hz = 5.0;
cfg.eo_pad = 2;
cfg.amplitude_limit_mm = 0.30;
cfg.offset_limit_mm = 0.20;

cfg.sensor_thresholds = zeros(1, max([base_cfg.sensor_ids, base_cfg.opr_id]));
for sid = base_cfg.sensor_ids
    cfg.sensor_thresholds(sid) = base_cfg.sensor_thresholds(sid);
end
cfg.sensor_thresholds(cfg.opr_channel) = cfg.opr_threshold;

cfg.sensor_config_file = fullfile(base_cfg.reference_output_dir, 'Sensor_Config_20250527.mat');
cfg.rebuild_case = false;
cfg.rebuild_calib = false;
cfg.save_experiment_case = false;
cfg.step5_case_mode = 'memory';
cfg.step5_force_rebuild_small_case = false;
cfg.show_figures = true;
cfg.step5_plot_case = true;
cfg.step5_plot_speed = true;
cfg.step5_plot_window_trends = true;
cfg.step5_plot_best_window = true;
cfg.step5_plot_seed_scan = true;
cfg.step5_save_figures = false;
cfg.step5_figure_dir = fullfile(cfg.project_dir, 'output', 'step5_diagnostics');
cfg.step4_save_figures = false;
cfg.step4_figure_dir = fullfile(cfg.project_dir, 'output', 'step4_diagnostics');
cfg.custom_calib_file = '';

cfg = resolve_single_sync_experiment_paths_20250527(cfg);
end
