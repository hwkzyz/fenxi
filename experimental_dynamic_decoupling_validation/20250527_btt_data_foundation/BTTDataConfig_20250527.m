function cfg = BTTDataConfig_20250527()
%BTTDATACONFIG_20250527 Method-neutral BTT data-foundation config.
%
% This config is intentionally limited to data facts and extraction settings.
% It should not contain static-gap, super-Gaussian, or direct-template
% identification parameters.
%
% 2026-07 revision notes:
%   1) OPR timing is configurable. The default is multi-threshold pulse
%      center so that later OPRCenterStd workflows use a consistent
%      reference frame.
%   2) Probe arrival timing is configurable and defaults to half-area timing.
%   3) Raw time-index interpretation is configurable: local, absolute, or auto.
%   4) Blade-1 definition and fingerprint/angle quality gates are explicit.

route_dir = fileparts(mfilename('fullpath'));

cfg = struct();
cfg.dataset = '20250527';
cfg.route_dir = route_dir;

cfg.dataset_root = 'E:\试验数据\20250527\试验20250527';
cfg.low_speed_case = '20250526_910';
cfg.dynamic_cases = {'20250526_2500-3500_t400'};

cfg.output_root = fullfile(route_dir, 'output');
cfg.step01_output_dir = fullfile(cfg.output_root, 'step01_low_speed_reference');
cfg.step02_output_dir = fullfile(cfg.output_root, 'step02_dynamic_btt');
cfg.step03_output_dir = fullfile(cfg.output_root, 'step03_observation_bundle');
cfg.step04_output_dir = fullfile(cfg.output_root, 'step04_low_speed_template');
cfg.figure_root = fullfile(cfg.output_root, 'figures');
cfg.step01_figure_dir = fullfile(cfg.figure_root, 'step01_low_speed_reference');
cfg.step02_figure_dir = fullfile(cfg.figure_root, 'step02_dynamic_btt');
cfg.step03_figure_dir = fullfile(cfg.figure_root, 'step03_observation_bundle');
cfg.step04_figure_dir = fullfile(cfg.figure_root, 'step04_low_speed_template');
cfg.save_figures = true;

% Candidate channels in this experiment. The default analysis set follows
% the current rotating-calibration route; extend to [1 3 6 7 8] if needed.
cfg.candidate_sensor_ids = [1 3 6 7 8];
cfg.sensor_ids = [1 3 6];
cfg.capacitance_ids = [1 3 6];
cfg.eddy_current_ids = [7 8];
cfg.unused_ids = [2 5];
cfg.opr_id = 4;

cfg.blades_num = 6;
cfg.sample_rate_hz = 5e6;
cfg.r_tip_mm = 62.0;

cfg.initial_trim_points = 20000;
cfg.gap_points = 1000;
cfg.opr_threshold = 1.5;
cfg.sensor_threshold_default = 0.5;
cfg.reference_sensor_id = 1;
cfg.reference_search_revs = 3;
cfg.max_alignment_search_starts = 24;
cfg.max_laps_process = 5000;

% Timing and coordinate-reference settings.
% - 'multi_threshold_center': OPR pulse center from several height levels.
% - 'rising_edge': threshold rising-edge crossing; kept for comparison.
cfg.opr_timing_method = 'multi_threshold_center';
cfg.opr_center_level_ratios = [0.30 0.40 0.50 0.60 0.70];
cfg.opr_center_smooth_span_points = 16;
cfg.opr_center_pad_fraction = 0.25;
cfg.opr_center_min_pad_points = 8;

% Probe-arrival method. Use the same value in low-speed reference and
% dynamic extraction if Step02 is rewritten around this data layer.
% Supported by the current low-speed extractor:
%   'half_area', 'threshold_centroid', 'peak', 'polynomial_centroid'.
cfg.probe_arrival_method = 'half_area';
cfg.probe_smooth_span_points = 21;

% Raw MAT first-column interpretation.
% - 'absolute': first column is an absolute sample index/time base.
% - 'local': each file restarts from a local sample index/time base.
% - 'auto': detect whether consecutive OPR files reset or continue.
cfg.time_index_mode = 'auto';

% Blade-1 convention. This is an internal relative blade label unless a
% physical blade offset is later supplied.
cfg.low_speed_blade1_rule = 'reference_sensor_max_peak_in_initial_revs';
cfg.blade_id_offset = 0;

% Fingerprint-alignment quality gates. Keep policy as 'warn' during method
% development; switch to 'error' for strict batch processing.
cfg.fingerprint_min_corr = 0.85;
cfg.fingerprint_min_corr_gap = 0.05;
cfg.fingerprint_quality_policy = 'warn';  % 'warn', 'error', or 'none'

% Standard-angle statistics. The compatibility field
% Sensor_Config.Standard_Relative_Angles uses this selected statistic.
cfg.standard_angle_value = 'mean';  % 'mean' or 'median'

cfg.sensor_thresholds = containers.Map('KeyType', 'double', 'ValueType', 'double');
for sid = cfg.candidate_sensor_ids
    cfg.sensor_thresholds(sid) = cfg.sensor_threshold_default;
end

% Compatibility alias used by some legacy helper-style code.
cfg.pinlv = cfg.sample_rate_hz;
end
