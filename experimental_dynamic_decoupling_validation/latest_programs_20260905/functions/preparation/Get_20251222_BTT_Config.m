function cfg = Get_20251222_BTT_Config()
%GET_20251222_BTT_CONFIG Centralized path/channel config for 20251222 BTT work.

packageCfg = Config_20251222();

cfg.dataset_root = packageCfg.paths.rawDatasetRoot;
cfg.low_speed_case = packageCfg.case.lowSpeedCase;
cfg.dynamic_cases = { ...
    '1000_2500_3500'
    };

cfg.output_root = packageCfg.paths.preparedLegacy;
cfg.reference_output_dir = fullfile(cfg.output_root, [cfg.low_speed_case, '_reference']);
cfg.strain_root = packageCfg.paths.strainRoot;

cfg.sensor_ids = [1, 2, 3];
cfg.capacitance_ids = [1, 2, 3];
cfg.eddy_current_ids = [];
cfg.unused_ids = [5, 6, 7, 8];
cfg.opr_id = 4;

cfg.blades_num = 6;
cfg.pinlv = 5e6;
cfg.r_tip_mm = 62.0;
cfg.initial_trim_points = 20000;
cfg.gap_points = 1000;
cfg.opr_threshold = 1.5;
cfg.sensor_threshold_default = 0.5;

cfg.reference_sensor_id = 1;
cfg.reference_search_revs = 3;
cfg.max_alignment_search_starts = 24;
cfg.max_laps_process = 5000;
cfg.step2_match_laps = 4;
cfg.step2_min_match_laps = 2;
cfg.step2_max_candidate_starts = 12;
cfg.step2_corr_std_weight = 0.15;
cfg.step2_feature_corr_weight = 1.0;
cfg.step2_angle_err_weight = 0.015;
cfg.step2_shift_change_penalty = 0.08;
cfg.step2_progress_rev_interval = 200;
cfg.step3_default_btt_sensor_ids = [1, 2, 3];
cfg.step3_default_strain_channel = 'AI1-03';
cfg.step3_default_alignment_offset_sec = 0;
cfg.step3_default_sync_window_sec = [-1, 148];
cfg.step3_default_stft_xlim_sec = [0, 148];
cfg.step3_default_fft_window_sec = [50, 56];
cfg.step3_manual_fft_window_sec = [50, 56];
cfg.step3_activity_grid_dt_sec = 0.02;
cfg.step3_activity_smooth_sec = 0.80;
cfg.step3_activity_threshold_ratio = 0.45;
cfg.step3_activity_pad_sec = 2.0;
cfg.step3_activity_min_window_sec = 8.0;
cfg.step3_activity_fft_window_sec = 1.5;
cfg.step3_stability_derivative_weight = 0.8;
cfg.step3_despike_enable = true;
cfg.step3_despike_window_sec = 0.01;
cfg.step3_despike_sigma = 6.0;
cfg.step3_align_freq_band_hz = [520, 620];
cfg.step3_align_freq_plot_hz = [450, 700];
cfg.step3_align_order_candidates = 8:20;
cfg.step3_align_stft_window_sec = 0.20;
cfg.step3_align_stft_overlap_ratio = 0.50;
cfg.step3_align_tau_step_sec = 0.05;
cfg.step3_align_min_valid_points = 40;
cfg.step3_align_frot_margin_hz = 12;
cfg.step3_align_frot_grid_step_hz = 0.2;
cfg.step3_align_frot_transition_weight = 0.18;
cfg.step3_align_frot_min_harmonics = 4;
cfg.step3_display_window_mode = 'full_btt';
cfg.step3_align_use_activity_window = true;
cfg.step3_align_activity_pad_sec = 5.0;
cfg.step3_align_min_activity_window_sec = 18.0;
cfg.step3_align_frot_median_span = 7;
cfg.step3_align_frot_sgolay_span = 11;
cfg.step3_strain_baseline_remove = true;
cfg.step3_strain_baseline_method = 'movmedian';
cfg.step3_strain_baseline_window_points = 10001;

cfg.sensor_thresholds = containers.Map('KeyType', 'double', 'ValueType', 'double');
for sid = cfg.sensor_ids
    cfg.sensor_thresholds(sid) = cfg.sensor_threshold_default;
end
end

