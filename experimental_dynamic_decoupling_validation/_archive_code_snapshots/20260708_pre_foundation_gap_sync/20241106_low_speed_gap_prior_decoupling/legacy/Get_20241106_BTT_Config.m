function cfg = Get_20241106_BTT_Config()
%GET_20241106_BTT_CONFIG Centralized path/channel config for 20241106 BTT work.

script_dir = fileparts(mfilename('fullpath'));

cfg.dataset_root = 'D:\博士-国科\试验台数据\新试验\20241106_2';
cfg.low_speed_case = '900';
cfg.dynamic_cases = { ...
    '3000_3150'
    };

cfg.output_root = fullfile(script_dir, 'output');
cfg.reference_output_dir = fullfile(cfg.output_root, [cfg.low_speed_case, '_reference']);
cfg.strain_root = '';

cfg.sensor_ids = [2, 3, 4, 5, 6, 7];
cfg.step2_sensor_ids = [5, 7];
cfg.capacitance_ids = [5, 7];
cfg.eddy_current_ids = [2, 3];
cfg.unused_ids = [];
cfg.opr_id = 1;

cfg.blades_num = 6;
cfg.opr_pulses_per_rev = 1;
cfg.pinlv = 5e6;
cfg.r_tip_mm = 65.0;
cfg.initial_trim_points = 20000;
cfg.gap_points = 10000;
cfg.opr_threshold = 2.0;
cfg.sensor_threshold_default = 1.0;

cfg.reference_sensor_id = 5;
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
cfg.step3_default_btt_sensor_ids = [5, 7];
cfg.step3_default_strain_channel = '';
cfg.step3_default_alignment_offset_sec = 0;
cfg.step3_default_sync_window_sec = [0, 200];
cfg.step3_default_stft_xlim_sec = [0, 200];
cfg.step3_default_fft_window_sec = [55, 100];
cfg.step3_manual_fft_window_sec = [55, 100];
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
cfg.step3_align_freq_band_hz = [520, 700];
cfg.step3_align_freq_plot_hz = [450, 760];
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
cfg.step3_strain_baseline_remove = false;
cfg.step3_strain_baseline_method = 'movmedian';
cfg.step3_strain_baseline_window_points = 10001;

cfg.sensor_thresholds = containers.Map('KeyType', 'double', 'ValueType', 'double');
cfg.sensor_thresholds(1) = 2.0;
cfg.sensor_thresholds(2) = -0.2;
cfg.sensor_thresholds(3) = -0.2;
cfg.sensor_thresholds(4) = 2.0;
cfg.sensor_thresholds(5) = 1.0;
cfg.sensor_thresholds(6) = 2.0;
cfg.sensor_thresholds(7) = 1.0;
end
