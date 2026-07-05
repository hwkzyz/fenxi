clc; clear; close all;

%STEP03_BUILD_DYNAMIC_OBSERVATION_BUNDLE_20250527 Build method-neutral BTT observations.
%
% Improved version for the 20250527 BTT data-foundation route.
%
% Step03 converts Step02 blade-arrival products into a reusable observation
% layer. It intentionally stops at the data/observation layer and does not
% perform low-speed template fitting, direct waveform inversion, decoupling,
% or modal identification.
%
% Main improvements over the previous draft:
%   1. Validate jiluOPR format and Step02 metadata before using jiluOPR(:,1)
%      as the OPR reference time.
%   2. Check consistency between Step01 standard-angle reference and Step02
%      dynamic OPR timing method.
%   3. Support the improved Step01 standard-angle statistics fields
%      (mean/median/std/IQR/count) and carry angle quality into observations.
%   4. Read Step02 dynamic-label quality diagnostics when available.
%   5. Use detailed quality flags rather than one generic invalid label.
%   6. Recompute the integration speed from the current OPR vector by default
%      and only use omega.mat for diagnostics/plotting.
%   7. Apply optional robust per sensor-blade outlier filtering.
%   8. Save full bundle, summary tables, compatibility vib files, and metadata
%      in the Step03 output folder. Legacy writing back to Step02 is optional.
%
% Usage:
%   Run this script directly after Step01 and Step02.

cfg = BTTDataConfig_20250527();
cfg = apply_step03_defaults_local(cfg);

%% Step03 run settings
show_plots = true;
target_cases = cfg.dynamic_cases;

%% Step03 observation-bundle settings
step03_target_blades = cfg.step03_target_blades;
step03_default_time_window_s = cfg.step03_default_time_window_s;
step03_outlier_abs_mm = cfg.step03_outlier_abs_mm;
step03_histogram_bin_count = cfg.step03_histogram_bin_count;
step03_plot_max_points_per_series = cfg.step03_plot_max_points_per_series;
step03_consistency_policy = cfg.step03_consistency_policy;
step03_require_jiluopr_three_columns = cfg.step03_require_jiluopr_three_columns;
step03_use_current_opr_for_speed = cfg.step03_use_current_opr_for_speed;
step03_use_robust_outlier_filter = cfg.step03_use_robust_outlier_filter;
step03_robust_mad_k = cfg.step03_robust_mad_k;
step03_min_points_for_robust_filter = cfg.step03_min_points_for_robust_filter;
step03_save_legacy_vib_to_step02 = cfg.step03_save_legacy_vib_to_step02;
step03_save_vib_to_step03 = cfg.step03_save_vib_to_step03;
step03_standard_angle_preference = cfg.step03_standard_angle_preference;
step03_wrap_angle_residual = cfg.step03_wrap_angle_residual;

cfg.step03_target_blades = step03_target_blades;
cfg.step03_default_time_window_s = step03_default_time_window_s;
cfg.step03_outlier_abs_mm = step03_outlier_abs_mm;
cfg.step03_histogram_bin_count = step03_histogram_bin_count;
cfg.step03_plot_max_points_per_series = step03_plot_max_points_per_series;
cfg.step03_consistency_policy = step03_consistency_policy;
cfg.step03_require_jiluopr_three_columns = step03_require_jiluopr_three_columns;
cfg.step03_use_current_opr_for_speed = step03_use_current_opr_for_speed;
cfg.step03_use_robust_outlier_filter = step03_use_robust_outlier_filter;
cfg.step03_robust_mad_k = step03_robust_mad_k;
cfg.step03_min_points_for_robust_filter = step03_min_points_for_robust_filter;
cfg.step03_save_legacy_vib_to_step02 = step03_save_legacy_vib_to_step02;
cfg.step03_save_vib_to_step03 = step03_save_vib_to_step03;
cfg.step03_standard_angle_preference = step03_standard_angle_preference;
cfg.step03_wrap_angle_residual = step03_wrap_angle_residual;

if ischar(target_cases) || isstring(target_cases)
    target_cases = cellstr(target_cases);
end

sensor_config_path = fullfile(cfg.step01_output_dir, 'Sensor_Config_20250527.mat');
if ~isfile(sensor_config_path)
    error('Missing Step01 sensor config. Please run Step01_Build_LowSpeed_Reference_20250527.m first:\n  %s', sensor_config_path);
end

loaded_cfg = load(sensor_config_path, 'Sensor_Config');
Sensor_Config = loaded_cfg.Sensor_Config;
check_step01_config_for_step03_local(cfg, Sensor_Config);

all_summary = table();
for iCase = 1:numel(target_cases)
    case_name = char(target_cases{iCase});
    step02_case_dir = fullfile(cfg.step02_output_dir, case_name);

    if ~isfolder(step02_case_dir) || ~isfile(fullfile(step02_case_dir, 'jiluOPR.mat'))
        error('Missing Step02 products for %s. Please run Step02_Extract_Dynamic_BTT_20250527.m first.', case_name);
    end

    fprintf('\n>>> [Step03][%s] Building observation bundle from:\n  %s\n', case_name, step02_case_dir);
    bundle = build_case_observation_bundle_local(case_name, cfg, Sensor_Config, sensor_config_path);

    case_output_dir = fullfile(cfg.step03_output_dir, case_name);
    if exist(case_output_dir, 'dir') ~= 7
        mkdir(case_output_dir);
    end
    bundle.Step03_Output_Dir = case_output_dir;
    bundle.Figure_Dir = fullfile(cfg.step03_figure_dir, case_name);

    observation_table = bundle.Observation_Table; %#ok<NASGU>
    summary_table = bundle.Summary_Table; %#ok<NASGU>
    metadata = bundle.Metadata; %#ok<NASGU>

    save(fullfile(case_output_dir, 'BTT_Observation_Bundle_20250527.mat'), 'bundle', 'metadata', '-v7.3');
    writetable(observation_table, fullfile(case_output_dir, 'BTT_Observation_LongTable_20250527.csv'));
    writetable(summary_table, fullfile(case_output_dir, 'BTT_Observation_Summary_20250527.csv'));

    save_compatibility_vib_files_local(bundle, case_output_dir, step02_case_dir, cfg);

    fprintf('>>> [Step03][%s] Observation rows = %d, valid rows = %d\n', ...
        case_name, height(bundle.Observation_Table), nnz(bundle.Observation_Table.is_valid));
    fprintf('>>> [Step03][%s] Output dir:\n  %s\n', case_name, case_output_dir);

    if show_plots || cfg.save_figures
        plot_step03_case_local(bundle, cfg, show_plots);
    end

    all_summary = [all_summary; summary_table]; %#ok<AGROW>
end

if ~isempty(all_summary)
    disp(all_summary);
end


function cfg = apply_step03_defaults_local(cfg)
if ~isfield(cfg, 'step03_target_blades') || isempty(cfg.step03_target_blades)
    cfg.step03_target_blades = 1:cfg.blades_num;
end
if ~isfield(cfg, 'step03_default_time_window_s')
    cfg.step03_default_time_window_s = [];
end
if ~isfield(cfg, 'step03_outlier_abs_mm') || isempty(cfg.step03_outlier_abs_mm)
    cfg.step03_outlier_abs_mm = 20;
end
if ~isfield(cfg, 'step03_histogram_bin_count') || isempty(cfg.step03_histogram_bin_count)
    cfg.step03_histogram_bin_count = 120;
end
if ~isfield(cfg, 'step03_plot_max_points_per_series') || isempty(cfg.step03_plot_max_points_per_series)
    cfg.step03_plot_max_points_per_series = 6000;
end

% New Step03 behavior controls.
if ~isfield(cfg, 'step03_consistency_policy')
    cfg.step03_consistency_policy = 'warn';  % 'warn', 'error', or 'none'
end
if ~isfield(cfg, 'step03_require_jiluopr_three_columns')
    cfg.step03_require_jiluopr_three_columns = false;
end
if ~isfield(cfg, 'step03_use_current_opr_for_speed')
    cfg.step03_use_current_opr_for_speed = true;
end
if ~isfield(cfg, 'step03_use_robust_outlier_filter')
    cfg.step03_use_robust_outlier_filter = false;
end
if ~isfield(cfg, 'step03_robust_mad_k')
    cfg.step03_robust_mad_k = 6.0;
end
if ~isfield(cfg, 'step03_min_points_for_robust_filter')
    cfg.step03_min_points_for_robust_filter = 12;
end
if ~isfield(cfg, 'step03_save_legacy_vib_to_step02')
    cfg.step03_save_legacy_vib_to_step02 = true;
end
if ~isfield(cfg, 'step03_save_vib_to_step03')
    cfg.step03_save_vib_to_step03 = true;
end
if ~isfield(cfg, 'step03_standard_angle_preference')
    cfg.step03_standard_angle_preference = 'selected'; % 'selected', 'mean', or 'median'
end
if ~isfield(cfg, 'step03_wrap_angle_residual')
    cfg.step03_wrap_angle_residual = true;
end
end


function check_step01_config_for_step03_local(cfg, Sensor_Config)
if ~isfield(Sensor_Config, 'Standard_Relative_Angles')
    error('Sensor_Config has no Standard_Relative_Angles field. Rebuild Step01 first.');
end
if isfield(Sensor_Config, 'Sensor_IDs')
    missing = setdiff(cfg.sensor_ids(:).', Sensor_Config.Sensor_IDs(:).');
    if ~isempty(missing)
        handle_consistency_issue_local(cfg, sprintf( ...
            'Step03 sensors %s are missing in Sensor_Config.Sensor_IDs.', mat2str(missing)));
    end
end
if isfield(Sensor_Config, 'OPR_Timing_Method') && isfield(cfg, 'opr_timing_method')
    if ~strcmpi(char(Sensor_Config.OPR_Timing_Method), char(cfg.opr_timing_method))
        handle_consistency_issue_local(cfg, sprintf( ...
            'Step01 OPR timing method (%s) differs from cfg.opr_timing_method (%s).', ...
            char(Sensor_Config.OPR_Timing_Method), char(cfg.opr_timing_method)));
    end
end
end


function bundle = build_case_observation_bundle_local(case_name, cfg, Sensor_Config, sensor_config_path)
step02_case_dir = fullfile(cfg.step02_output_dir, case_name);
step02_diag = load_step02_diagnostics_local(step02_case_dir);

[opr_times, opr_product, opr_metadata] = load_and_validate_opr_product_local(step02_case_dir, cfg, Sensor_Config);
[omega_diag_time_s, omega_diag_rpm, omega_diag_rad_s, omega_source] = ...
    load_or_compute_rpm_local(fullfile(step02_case_dir, 'omega.mat'), opr_times, cfg.blades_num);

% Always use speed derived from the exact OPR vector used for angular
% integration unless explicitly disabled.
if cfg.step03_use_current_opr_for_speed
    [omega_time_s, omega_rpm, omega_rad_s] = compute_rpm_from_opr_local(opr_times, cfg.blades_num);
    speed_source = 'recomputed_from_current_jiluOPR_column1';
else
    omega_time_s = omega_diag_time_s;
    omega_rpm = omega_diag_rpm;
    omega_rad_s = omega_diag_rad_s;
    speed_source = omega_source;
end

validate_omega_consistency_local(cfg, opr_times, omega_diag_time_s, omega_diag_rpm, omega_source);

[F_omega_deg_s, F_rpm] = build_speed_interpolants_local(opr_times, omega_time_s, omega_rpm, cfg.blades_num);
angle_ref = get_standard_angle_bundle_local(Sensor_Config, cfg);

all_tables = {};
per_sensor = repmat(struct( ...
    'sensor_id', NaN, ...
    'raw_jilublade_path', '', ...
    'vib_final_path_step03', '', ...
    'vib_final_path_legacy_step02', '', ...
    'observation_table', table(), ...
    'jilublade_vib', [], ...
    'step02_label_quality', struct()), 1, numel(cfg.sensor_ids));

for iSensor = 1:numel(cfg.sensor_ids)
    sid = cfg.sensor_ids(iSensor);
    probe_path = fullfile(step02_case_dir, sprintf('jilublade_probe%d.mat', sid));
    if ~isfile(probe_path)
        warning('Step03:MissingProbe', 'Missing probe file for case %s CH%d: %s', case_name, sid, probe_path);
        continue;
    end

    loaded_probe = load(probe_path);
    if isfield(loaded_probe, 'jilublade')
        jilublade = loaded_probe.jilublade;
    else
        fn = fieldnames(loaded_probe);
        jilublade = loaded_probe.(fn{1});
    end
    if size(jilublade, 2) < 4
        warning('Step03:InvalidProbe', 'Probe file has fewer than 4 columns for case %s: %s', case_name, probe_path);
        continue;
    end

    label_quality = get_step02_sensor_quality_local(step02_diag, sid);
    obs = convert_jilublade_to_observations_local( ...
        case_name, sid, jilublade, opr_times, F_omega_deg_s, F_rpm, angle_ref, label_quality, cfg);

    all_tables{end + 1} = obs.table; %#ok<AGROW>
    per_sensor(iSensor).sensor_id = sid;
    per_sensor(iSensor).raw_jilublade_path = probe_path;
    per_sensor(iSensor).vib_final_path_step03 = fullfile(cfg.step03_output_dir, case_name, sprintf('jilublade_probe%d_vib_final.mat', sid));
    per_sensor(iSensor).vib_final_path_legacy_step02 = fullfile(step02_case_dir, sprintf('jilublade_probe%d_vib_final.mat', sid));
    per_sensor(iSensor).observation_table = obs.table;
    per_sensor(iSensor).jilublade_vib = obs.jilublade_vib;
    per_sensor(iSensor).step02_label_quality = label_quality;
end

if isempty(all_tables)
    error('No valid Step02 probe files found for case %s.', case_name);
end

observation_table = vertcat(all_tables{:});
if ~isempty(cfg.step03_default_time_window_s)
    tw = cfg.step03_default_time_window_s(:).';
    keep = observation_table.arrival_time_s >= tw(1) & observation_table.arrival_time_s <= tw(2);
    observation_table = observation_table(keep, :);
end

observation_table = apply_robust_outlier_filter_local(observation_table, cfg);
summary_table = summarize_observations_local(case_name, observation_table, cfg, step02_diag);

metadata = build_step03_metadata_local(cfg, Sensor_Config, opr_metadata, opr_product, speed_source, omega_source, angle_ref, step02_diag);

bundle = struct();
bundle.Dataset = cfg.dataset;
bundle.Case_Name = char(case_name);
bundle.CreatedBy = mfilename;
bundle.CreatedOn = datestr(now, 31);
bundle.Route_Dir = cfg.route_dir;
bundle.Step01_Sensor_Config_Path = sensor_config_path;
bundle.Step02_Case_Dir = step02_case_dir;
bundle.Step03_Output_Dir = fullfile(cfg.step03_output_dir, case_name);
bundle.Figure_Dir = fullfile(cfg.step03_figure_dir, case_name);
bundle.Sensor_IDs = cfg.sensor_ids;
bundle.Blades_Num = cfg.blades_num;
bundle.R_Tip_mm = cfg.r_tip_mm;
bundle.OPR_Times_s = opr_times(:);
bundle.OPR_Product = opr_product;
bundle.Omega_Time_s = omega_time_s(:);
bundle.Omega_RPM = omega_rpm(:);
bundle.Omega_Rad_s = omega_rad_s(:);
bundle.Omega_Diagnostic_Time_s = omega_diag_time_s(:);
bundle.Omega_Diagnostic_RPM = omega_diag_rpm(:);
bundle.Omega_Diagnostic_Rad_s = omega_diag_rad_s(:);
bundle.Speed_Source = speed_source;
bundle.Standard_Relative_Angles = angle_ref.value;
bundle.Standard_Angle_Reference = angle_ref;
bundle.Step02_Diagnostics = step02_diag;
bundle.Metadata = metadata;
bundle.Observation_Table = observation_table;
bundle.Summary_Table = summary_table;
bundle.Per_Sensor = per_sensor;
end


function [opr_times, opr_product, metadata] = load_and_validate_opr_product_local(step02_case_dir, cfg, Sensor_Config)
opr_path = fullfile(step02_case_dir, 'jiluOPR.mat');
if ~isfile(opr_path)
    error('Missing Step02 OPR file: %s', opr_path);
end
loaded_opr = load(opr_path);
if ~isfield(loaded_opr, 'jiluOPR') || size(loaded_opr.jiluOPR, 2) < 1
    error('Invalid OPR product: %s', opr_path);
end
jiluOPR = loaded_opr.jiluOPR;
metadata = struct();
if isfield(loaded_opr, 'metadata')
    metadata = loaded_opr.metadata;
end

if size(jiluOPR, 2) < 3
    msg = sprintf(['jiluOPR in %s has only %d column(s). Improved Step02 should save ' ...
        '[opr_arrival_or_center_s, opr_start_s, opr_end_s]. Existing file may be an old product.'], ...
        opr_path, size(jiluOPR, 2));
    if cfg.step03_require_jiluopr_three_columns
        error(msg);
    else
        handle_consistency_issue_local(cfg, msg);
    end
end

if isfield(metadata, 'jiluOPR_columns')
    cols = metadata.jiluOPR_columns;
    if iscell(cols) && ~isempty(cols)
        first_col = char(cols{1});
        if isempty(regexpi(first_col, 'arrival|center', 'once'))
            handle_consistency_issue_local(cfg, sprintf( ...
                'jiluOPR first column metadata is "%s", not clearly an arrival/center reference.', first_col));
        end
    end
else
    handle_consistency_issue_local(cfg, 'jiluOPR metadata is missing. Cannot verify first-column meaning.');
end

if isfield(metadata, 'opr_timing_method') && isfield(Sensor_Config, 'OPR_Timing_Method')
    if ~strcmpi(char(metadata.opr_timing_method), char(Sensor_Config.OPR_Timing_Method))
        handle_consistency_issue_local(cfg, sprintf( ...
            'Step02 OPR timing method (%s) differs from Step01 standard-angle reference (%s).', ...
            char(metadata.opr_timing_method), char(Sensor_Config.OPR_Timing_Method)));
    end
elseif isfield(Sensor_Config, 'OPR_Timing_Method')
    handle_consistency_issue_local(cfg, sprintf( ...
        'Step02 OPR metadata missing; Step01 standard-angle reference is %s.', char(Sensor_Config.OPR_Timing_Method)));
end

opr_times = jiluOPR(:, 1);
opr_times = opr_times(isfinite(opr_times));
if numel(opr_times) <= cfg.blades_num
    error('Too few OPR pulses in %s.', opr_path);
end
if any(diff(opr_times) <= 0)
    error('OPR times in %s are not strictly increasing.', opr_path);
end

opr_product = struct();
opr_product.path = opr_path;
opr_product.jiluOPR = jiluOPR;
opr_product.column_count = size(jiluOPR, 2);
opr_product.metadata = metadata;
end


function step02_diag = load_step02_diagnostics_local(step02_case_dir)
step02_diag = struct();
step02_diag.available = false;
step02_diag.case_data = [];
step02_diag.metadata = struct();
step02_diag.summary_table = table();
step02_diag.source_file = '';

full_file = fullfile(step02_case_dir, 'Step02_Dynamic_BTT_Extraction_20250527.mat');
if isfile(full_file)
    S = load(full_file);
    step02_diag.available = true;
    step02_diag.source_file = full_file;
    if isfield(S, 'case_data')
        step02_diag.case_data = S.case_data;
    end
    if isfield(S, 'metadata')
        step02_diag.metadata = S.metadata;
    end
    if isfield(S, 'case_summary')
        step02_diag.summary_table = S.case_summary;
    end
end

summary_mat = fullfile(step02_case_dir, 'Step02_Summary_20250527.mat');
if isempty(step02_diag.summary_table) && isfile(summary_mat)
    S = load(summary_mat);
    if isfield(S, 'case_summary')
        step02_diag.summary_table = S.case_summary;
    end
    if isfield(S, 'metadata') && isempty(fieldnames(step02_diag.metadata))
        step02_diag.metadata = S.metadata;
    end
end

summary_csv = fullfile(step02_case_dir, 'Step02_Summary_20250527.csv');
if isempty(step02_diag.summary_table) && isfile(summary_csv)
    try
        step02_diag.summary_table = readtable(summary_csv);
    catch
        step02_diag.summary_table = table();
    end
end
end


function label_quality = get_step02_sensor_quality_local(step02_diag, sid)
label_quality = struct( ...
    'available', false, ...
    'best_corr', NaN, ...
    'second_best_corr', NaN, ...
    'corr_gap', NaN, ...
    'shift_consistency', NaN, ...
    'quality_status', 'unknown', ...
    'valid_revolution_count', NaN, ...
    'unlabeled_pulse_count', NaN);

T = step02_diag.summary_table;
if ~isempty(T) && any(strcmp(T.Properties.VariableNames, 'sensor_id'))
    idx = find(T.sensor_id == sid, 1, 'first');
    if ~isempty(idx)
        label_quality.available = true;
        label_quality.best_corr = get_table_value_local(T, idx, 'best_corr', NaN);
        label_quality.second_best_corr = get_table_value_local(T, idx, 'second_best_corr', NaN);
        label_quality.corr_gap = get_table_value_local(T, idx, 'corr_gap', NaN);
        label_quality.shift_consistency = get_table_value_local(T, idx, 'shift_consistency', NaN);
        label_quality.valid_revolution_count = get_table_value_local(T, idx, 'valid_revolution_count', NaN);
        label_quality.unlabeled_pulse_count = get_table_value_local(T, idx, 'unlabeled_pulse_count', NaN);
        if any(strcmp(T.Properties.VariableNames, 'quality_status'))
            label_quality.quality_status = char(T.quality_status(idx));
        end
        return;
    end
end

if step02_diag.available && isfield(step02_diag.case_data, 'channels') && sid <= numel(step02_diag.case_data.channels)
    ch = step02_diag.case_data.channels(sid);
    label_quality.available = true;
    label_quality.best_corr = get_optional_field_local(ch, 'best_corr', NaN);
    label_quality.second_best_corr = get_optional_field_local(ch, 'second_best_corr', NaN);
    label_quality.corr_gap = get_optional_field_local(ch, 'corr_gap', NaN);
    label_quality.shift_consistency = get_optional_field_local(ch, 'shift_consistency', NaN);
    label_quality.valid_revolution_count = get_optional_field_local(ch, 'valid_revolution_count', NaN);
    label_quality.unlabeled_pulse_count = get_optional_field_local(ch, 'unlabeled_pulse_count', NaN);
    label_quality.quality_status = get_optional_field_local(ch, 'quality_status', 'unknown');
end
end


function value = get_table_value_local(T, idx, name, default_value)
value = default_value;
if any(strcmp(T.Properties.VariableNames, name))
    tmp = T.(name)(idx);
    if isnumeric(tmp) || islogical(tmp)
        value = tmp(1);
    elseif iscell(tmp)
        value = tmp{1};
    elseif isstring(tmp)
        value = char(tmp(1));
    else
        value = tmp(1);
    end
end
end


function angle_ref = get_standard_angle_bundle_local(Sensor_Config, cfg)
preference = lower(strtrim(cfg.step03_standard_angle_preference));
value = [];
value_name = '';
if strcmp(preference, 'mean') && isfield(Sensor_Config, 'Standard_Relative_Angles_Mean')
    value = Sensor_Config.Standard_Relative_Angles_Mean;
    value_name = 'mean';
elseif strcmp(preference, 'median') && isfield(Sensor_Config, 'Standard_Relative_Angles_Median')
    value = Sensor_Config.Standard_Relative_Angles_Median;
    value_name = 'median';
elseif strcmp(preference, 'selected') && isfield(Sensor_Config, 'Standard_Relative_Angles')
    value = Sensor_Config.Standard_Relative_Angles;
    value_name = get_optional_field_local(Sensor_Config, 'Standard_Relative_Angles_SelectedStatistic', 'selected');
elseif isfield(Sensor_Config, 'Standard_Relative_Angles')
    value = Sensor_Config.Standard_Relative_Angles;
    value_name = 'legacy_or_selected';
else
    error('No usable standard relative angle matrix found in Sensor_Config.');
end

angle_ref = struct();
angle_ref.value = value;
angle_ref.value_name = char(value_name);
angle_ref.reference = get_optional_field_local(Sensor_Config, 'Standard_Relative_Angles_Reference', 'unknown');
angle_ref.mean = get_optional_field_local(Sensor_Config, 'Standard_Relative_Angles_Mean', nan(size(value)));
angle_ref.median = get_optional_field_local(Sensor_Config, 'Standard_Relative_Angles_Median', nan(size(value)));
angle_ref.std = get_optional_field_local(Sensor_Config, 'Standard_Relative_Angles_Std', nan(size(value)));
angle_ref.iqr = get_optional_field_local(Sensor_Config, 'Standard_Relative_Angles_IQR', nan(size(value)));
angle_ref.count = get_optional_field_local(Sensor_Config, 'Standard_Relative_Angles_Count', nan(size(value)));
angle_ref.min = get_optional_field_local(Sensor_Config, 'Standard_Relative_Angles_Min', nan(size(value)));
angle_ref.max = get_optional_field_local(Sensor_Config, 'Standard_Relative_Angles_Max', nan(size(value)));
end


function obs = convert_jilublade_to_observations_local( ...
    case_name, sid, jilublade, opr_times, F_omega_deg_s, F_rpm, angle_ref, label_quality, cfg)

n = size(jilublade, 1);
pulse_index = (1:n).';
start_time_s = jilublade(:, 1);
end_time_s = jilublade(:, 2);
arrival_time_s = jilublade(:, 3);
blade_id = round(jilublade(:, 4));

prev_opr_index = find_previous_opr_indices_local(opr_times, arrival_time_s);
revolution_index = nan(n, 1);
valid_prev = isfinite(prev_opr_index) & prev_opr_index >= 1;
revolution_index(valid_prev) = floor((prev_opr_index(valid_prev) - 1) / cfg.blades_num) + 1;

theta_actual_deg = nan(n, 1);
theta_standard_deg = nan(n, 1);
theta_error_deg = nan(n, 1);
displacement_mm = nan(n, 1);
rpm = nan(n, 1);
std_angle_std_deg = nan(n, 1);
std_angle_iqr_deg = nan(n, 1);
std_angle_count = nan(n, 1);
quality_flag = repmat({'valid'}, n, 1);

for i = 1:n
    t_meas = arrival_time_s(i);
    bid = blade_id(i);
    idx_prev = prev_opr_index(i);

    if ~isfinite(t_meas)
        quality_flag{i} = 'missing_arrival_time';
        continue;
    end
    if ~isfinite(idx_prev) || idx_prev < 1 || idx_prev > numel(opr_times)
        quality_flag{i} = 'missing_previous_opr';
        continue;
    end
    if ~isfinite(bid) || bid < 1 || bid > cfg.blades_num
        quality_flag{i} = 'invalid_blade_id';
        continue;
    end
    if sid > size(angle_ref.value, 1) || bid > size(angle_ref.value, 2)
        quality_flag{i} = 'standard_angle_index_out_of_range';
        continue;
    end

    theta_standard = angle_ref.value(sid, bid);
    if ~isfinite(theta_standard)
        quality_flag{i} = 'missing_standard_angle';
        continue;
    end

    t_ref = opr_times(idx_prev);
    if ~isfinite(t_ref) || t_meas < t_ref
        quality_flag{i} = 'negative_time_interval';
        continue;
    end

    t_grid = linspace(t_ref, t_meas, 10);
    theta_actual = trapz(t_grid, F_omega_deg_s(t_grid));
    theta_diff = theta_actual - theta_standard;
    if cfg.step03_wrap_angle_residual
        theta_diff = wrap_to_180_local(theta_diff);
    end

    theta_actual_deg(i) = theta_actual;
    theta_standard_deg(i) = theta_standard;
    theta_error_deg(i) = theta_diff;
    displacement_mm(i) = theta_diff * (pi / 180) * cfg.r_tip_mm;
    rpm(i) = F_rpm(t_meas);
    std_angle_std_deg(i) = get_angle_matrix_value_local(angle_ref.std, sid, bid);
    std_angle_iqr_deg(i) = get_angle_matrix_value_local(angle_ref.iqr, sid, bid);
    std_angle_count(i) = get_angle_matrix_value_local(angle_ref.count, sid, bid);

    if ~isfinite(displacement_mm(i))
        quality_flag{i} = 'nonfinite_displacement';
    elseif abs(displacement_mm(i)) > cfg.step03_outlier_abs_mm
        quality_flag{i} = 'absolute_outlier_displacement';
    end
end

is_valid = strcmp(quality_flag, 'valid');

case_col = repmat({char(case_name)}, n, 1);
sensor_col = repmat(sid, n, 1);
step02_best_corr = repmat(label_quality.best_corr, n, 1);
step02_corr_gap = repmat(label_quality.corr_gap, n, 1);
step02_shift_consistency = repmat(label_quality.shift_consistency, n, 1);
step02_label_quality = repmat({char(label_quality.quality_status)}, n, 1);

T = table(case_col, sensor_col, pulse_index, revolution_index, blade_id, ...
    start_time_s, end_time_s, arrival_time_s, prev_opr_index, ...
    theta_actual_deg, theta_standard_deg, theta_error_deg, displacement_mm, rpm, ...
    std_angle_std_deg, std_angle_iqr_deg, std_angle_count, ...
    step02_best_corr, step02_corr_gap, step02_shift_consistency, step02_label_quality, ...
    is_valid, quality_flag, ...
    'VariableNames', {'case_name', 'sensor_id', 'pulse_index', 'revolution_index', 'blade_id', ...
    'start_time_s', 'end_time_s', 'arrival_time_s', 'prev_opr_index', ...
    'theta_actual_deg', 'theta_standard_deg', 'theta_error_deg', 'displacement_mm', 'rpm', ...
    'std_angle_std_deg', 'std_angle_iqr_deg', 'std_angle_count', ...
    'step02_best_corr', 'step02_corr_gap', 'step02_shift_consistency', 'step02_label_quality', ...
    'is_valid', 'quality_flag'});

jilublade_vib = nan(n, 12);
jilublade_vib(:, 1:4) = jilublade(:, 1:4);
jilublade_vib(:, 5) = theta_actual_deg;
jilublade_vib(:, 6) = displacement_mm;
jilublade_vib(:, 7) = theta_standard_deg;
jilublade_vib(:, 8) = theta_error_deg;
jilublade_vib(:, 9) = prev_opr_index;
jilublade_vib(:, 10) = revolution_index;
jilublade_vib(:, 11) = rpm;
jilublade_vib(:, 12) = is_valid;

obs = struct();
obs.table = T;
obs.jilublade_vib = jilublade_vib;
end


function value = get_angle_matrix_value_local(M, sid, bid)
value = NaN;
if ~isempty(M) && sid <= size(M, 1) && bid <= size(M, 2)
    value = M(sid, bid);
end
end


function observation_table = apply_robust_outlier_filter_local(observation_table, cfg)
if ~cfg.step03_use_robust_outlier_filter || isempty(observation_table)
    return;
end
sensor_ids = unique(observation_table.sensor_id(:).');
blade_ids = unique(observation_table.blade_id(isfinite(observation_table.blade_id)).');
for sid = sensor_ids
    for bid = blade_ids
        mask = observation_table.sensor_id == sid & observation_table.blade_id == bid & observation_table.is_valid;
        idx = find(mask);
        if numel(idx) < cfg.step03_min_points_for_robust_filter
            continue;
        end
        y = observation_table.displacement_mm(idx);
        med_y = median(y, 'omitnan');
        mad_y = median(abs(y - med_y), 'omitnan');
        sigma = 1.4826 * mad_y;
        if ~isfinite(sigma) || sigma < eps
            continue;
        end
        bad_local = abs(y - med_y) > cfg.step03_robust_mad_k * sigma;
        if any(bad_local)
            bad_idx = idx(bad_local);
            observation_table.is_valid(bad_idx) = false;
            observation_table.quality_flag(bad_idx) = {'robust_outlier_displacement'};
        end
    end
end
end


function summary_table = summarize_observations_local(case_name, observation_table, cfg, step02_diag)
rows = repmat(struct( ...
    'case_name', '', ...
    'sensor_id', NaN, ...
    'blade_id', NaN, ...
    'point_count', NaN, ...
    'valid_count', NaN, ...
    'valid_ratio', NaN, ...
    'time_start_s', NaN, ...
    'time_end_s', NaN, ...
    'rpm_mean', NaN, ...
    'disp_mean_mm', NaN, ...
    'disp_median_mm', NaN, ...
    'disp_std_mm', NaN, ...
    'disp_rms_mm', NaN, ...
    'disp_min_mm', NaN, ...
    'disp_max_mm', NaN, ...
    'std_angle_std_deg', NaN, ...
    'std_angle_iqr_deg', NaN, ...
    'std_angle_count', NaN, ...
    'step02_best_corr', NaN, ...
    'step02_corr_gap', NaN, ...
    'step02_shift_consistency', NaN, ...
    'step02_label_quality', ''), numel(cfg.sensor_ids) * cfg.blades_num, 1);

row = 0;
for iSensor = 1:numel(cfg.sensor_ids)
    sid = cfg.sensor_ids(iSensor);
    label_quality = get_step02_sensor_quality_local(step02_diag, sid);
    for blade_id = 1:cfg.blades_num
        row = row + 1;
        mask_all = observation_table.sensor_id == sid & observation_table.blade_id == blade_id;
        mask_valid = mask_all & observation_table.is_valid;
        disp_v = observation_table.displacement_mm(mask_valid);
        rpm_v = observation_table.rpm(mask_valid);
        t_v = observation_table.arrival_time_s(mask_valid);

        rows(row).case_name = char(case_name);
        rows(row).sensor_id = sid;
        rows(row).blade_id = blade_id;
        rows(row).point_count = nnz(mask_all);
        rows(row).valid_count = nnz(mask_valid);
        if nnz(mask_all) > 0
            rows(row).valid_ratio = nnz(mask_valid) / nnz(mask_all);
        end
        if ~isempty(t_v)
            rows(row).time_start_s = min(t_v);
            rows(row).time_end_s = max(t_v);
        end
        if ~isempty(rpm_v)
            rows(row).rpm_mean = mean(rpm_v, 'omitnan');
        end
        if ~isempty(disp_v)
            rows(row).disp_mean_mm = mean(disp_v, 'omitnan');
            rows(row).disp_median_mm = median(disp_v, 'omitnan');
            rows(row).disp_std_mm = std(disp_v, 'omitnan');
            rows(row).disp_rms_mm = sqrt(mean(disp_v.^2, 'omitnan'));
            rows(row).disp_min_mm = min(disp_v, [], 'omitnan');
            rows(row).disp_max_mm = max(disp_v, [], 'omitnan');
        end
        if any(mask_all)
            rows(row).std_angle_std_deg = first_finite_local(observation_table.std_angle_std_deg(mask_all));
            rows(row).std_angle_iqr_deg = first_finite_local(observation_table.std_angle_iqr_deg(mask_all));
            rows(row).std_angle_count = first_finite_local(observation_table.std_angle_count(mask_all));
        end
        rows(row).step02_best_corr = label_quality.best_corr;
        rows(row).step02_corr_gap = label_quality.corr_gap;
        rows(row).step02_shift_consistency = label_quality.shift_consistency;
        rows(row).step02_label_quality = char(label_quality.quality_status);
    end
end

summary_table = struct2table(rows);
end


function x = first_finite_local(v)
x = NaN;
v = v(:);
idx = find(isfinite(v), 1, 'first');
if ~isempty(idx)
    x = v(idx);
end
end


function save_compatibility_vib_files_local(bundle, case_output_dir, step02_case_dir, cfg)
metadata = bundle.Metadata; %#ok<NASGU>
for i = 1:numel(bundle.Per_Sensor)
    sid = bundle.Per_Sensor(i).sensor_id;
    if ~isfinite(sid) || isempty(bundle.Per_Sensor(i).jilublade_vib)
        continue;
    end
    jilublade = bundle.Per_Sensor(i).jilublade_vib; %#ok<NASGU>
    if cfg.step03_save_vib_to_step03
        save(fullfile(case_output_dir, sprintf('jilublade_probe%d_vib_final.mat', sid)), ...
            'jilublade', 'metadata');
    end
    if cfg.step03_save_legacy_vib_to_step02
        save(fullfile(step02_case_dir, sprintf('jilublade_probe%d_vib_final.mat', sid)), ...
            'jilublade', 'metadata');
    end
end
end


function metadata = build_step03_metadata_local(cfg, Sensor_Config, opr_metadata, opr_product, speed_source, omega_source, angle_ref, step02_diag)
metadata = struct();
metadata.dataset = cfg.dataset;
metadata.created_by = mfilename;
metadata.created_on = datestr(now, 31);
metadata.description = 'Method-neutral dynamic observation bundle: arrival -> actual angle -> standard-angle residual -> equivalent tip displacement.';
metadata.jiluOPR_source = opr_product.path;
metadata.jiluOPR_columns = get_optional_field_local(opr_metadata, 'jiluOPR_columns', {});
metadata.step02_opr_timing_method = get_optional_field_local(opr_metadata, 'opr_timing_method', 'unknown');
metadata.step01_opr_timing_method = get_optional_field_local(Sensor_Config, 'OPR_Timing_Method', 'unknown');
metadata.standard_angle_reference = angle_ref.reference;
metadata.standard_angle_value_name = angle_ref.value_name;
metadata.speed_source_for_integration = speed_source;
metadata.omega_file_source = omega_source;
metadata.robust_outlier_filter_enabled = cfg.step03_use_robust_outlier_filter;
metadata.robust_mad_k = cfg.step03_robust_mad_k;
metadata.absolute_outlier_limit_mm = cfg.step03_outlier_abs_mm;
metadata.step02_diagnostics_available = step02_diag.available;
metadata.step02_diagnostics_source = step02_diag.source_file;
metadata.vib_final_columns = {'pulse_start_s','pulse_end_s','probe_arrival_s','blade_id', ...
    'theta_actual_deg','displacement_mm','theta_standard_deg','theta_error_deg', ...
    'prev_opr_index','revolution_index','rpm','is_valid'};
end


function [omega_time_s, omega_rpm, omega_rad_s, source] = load_or_compute_rpm_local(omega_path, opr_times, blades_num)
omega_time_s = [];
omega_rpm = [];
omega_rad_s = [];
source = 'recomputed_from_opr';

if isfile(omega_path)
    loaded = load(omega_path);
    if isfield(loaded, 'omega_time_s') && isfield(loaded, 'omega_rpm')
        omega_time_s = loaded.omega_time_s(:);
        omega_rpm = loaded.omega_rpm(:);
        if isfield(loaded, 'omega_rad_s')
            omega_rad_s = loaded.omega_rad_s(:);
        else
            omega_rad_s = omega_rpm * 2 * pi / 60;
        end
        source = 'omega_mat';
        return;
    end
end

[omega_time_s, omega_rpm, omega_rad_s] = compute_rpm_from_opr_local(opr_times, blades_num);
end


function [omega_time_s, omega_rpm, omega_rad_s] = compute_rpm_from_opr_local(opr_times, blades_num)
if numel(opr_times) <= blades_num
    omega_time_s = [];
    omega_rpm = [];
    omega_rad_s = [];
    return;
end
rev_period = opr_times((blades_num + 1):end) - opr_times(1:(end - blades_num));
omega_time_s = opr_times(1:(end - blades_num));
omega_rpm = 60 ./ rev_period;
omega_rad_s = 2 * pi ./ rev_period;
end


function validate_omega_consistency_local(cfg, opr_times, omega_time_s, omega_rpm, omega_source)
if isempty(omega_time_s) || isempty(omega_rpm) || ~strcmpi(omega_source, 'omega_mat')
    return;
end
expected_n = numel(opr_times) - cfg.blades_num;
if numel(omega_time_s) ~= expected_n
    handle_consistency_issue_local(cfg, sprintf( ...
        'omega.mat length (%d) differs from current OPR-derived length (%d). Step03 will use current OPR for integration if enabled.', ...
        numel(omega_time_s), expected_n));
end
end


function [F_omega_deg_s, F_rpm] = build_speed_interpolants_local(opr_times, omega_time_s, omega_rpm, blades_num)
speed_time = opr_times(1:(end - blades_num));
speed_deg_s = 360 ./ max(opr_times((blades_num + 1):end) - opr_times(1:(end - blades_num)), eps);
F_omega_deg_s = griddedInterpolant(speed_time(:), speed_deg_s(:), 'linear', 'nearest');

if isempty(omega_time_s) || isempty(omega_rpm)
    F_rpm = griddedInterpolant(speed_time(:), (speed_deg_s(:) / 360) * 60, 'linear', 'nearest');
else
    F_rpm = griddedInterpolant(omega_time_s(:), omega_rpm(:), 'linear', 'nearest');
end
end


function prev_idx = find_previous_opr_indices_local(opr_times, t)
prev_idx = nan(size(t));
if isempty(opr_times) || isempty(t)
    return;
end
edges = [-inf; opr_times(:); inf];
bin = discretize(t(:), edges);
idx = bin - 1;
idx(idx < 1 | idx > numel(opr_times)) = NaN;
prev_idx(:) = idx;
end


function ang = wrap_to_180_local(ang)
ang = mod(ang + 180, 360) - 180;
end


function handle_consistency_issue_local(cfg, msg)
policy = 'warn';
if isfield(cfg, 'step03_consistency_policy')
    policy = lower(strtrim(cfg.step03_consistency_policy));
end
switch policy
    case 'error'
        error('%s', msg);
    case 'warn'
        warning('Step03:Consistency', '%s', msg);
    case 'none'
        return;
    otherwise
        warning('Step03:Consistency', '%s', msg);
end
end


function value = get_optional_field_local(S, name, default_value)
value = default_value;
if isstruct(S) && isfield(S, name)
    value = S.(name);
end
end


function plot_step03_case_local(bundle, cfg, show_plots)
figure_dir = bundle.Figure_Dir;
if exist(figure_dir, 'dir') ~= 7
    mkdir(figure_dir);
end

T = bundle.Observation_Table;
S = bundle.Summary_Table;
sensor_ids = bundle.Sensor_IDs;
blades_num = bundle.Blades_Num;
colors = lines(max(numel(sensor_ids), 1));

fig1 = figure('Name', sprintf('20250527 Step03 Displacement Time Series - %s', bundle.Case_Name), ...
    'Color', 'w', 'Position', [70, 40, 1500, 930], 'NumberTitle', 'off', ...
    'Visible', visibility_state_local(show_plots));
tiledlayout(fig1, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for blade_id = 1:blades_num
    nexttile;
    hold on;
    grid on;
    box on;
    for iSensor = 1:numel(sensor_ids)
        sid = sensor_ids(iSensor);
        mask = T.is_valid & T.sensor_id == sid & T.blade_id == blade_id;
        if ~any(mask)
            continue;
        end
        [tt, yy] = decimate_series_for_plot_local(T.arrival_time_s(mask), T.displacement_mm(mask), cfg.step03_plot_max_points_per_series);
        plot(tt, yy, '.', 'Color', colors(iSensor, :), 'MarkerSize', 5, ...
            'DisplayName', sprintf('CH%d', sid));
    end
    yline(0, 'k:', 'HandleVisibility', 'off');
    xlabel('Time (s)');
    ylabel('Displacement (mm)');
    title(sprintf('Blade %d equivalent tangential displacement', blade_id));
    legend('Location', 'best');
end
save_step03_figure_local(fig1, figure_dir, cfg, 'Step03_Displacement_TimeSeries');

fig2 = figure('Name', sprintf('20250527 Step03 Sensor-Blade Summary - %s', bundle.Case_Name), ...
    'Color', 'w', 'Position', [100, 80, 1280, 650], 'NumberTitle', 'off', ...
    'Visible', visibility_state_local(show_plots));
tiledlayout(fig2, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
[std_mat, rms_mat, count_mat] = summary_to_matrices_local(S, sensor_ids, blades_num);

nexttile;
imagesc(1:blades_num, 1:numel(sensor_ids), std_mat);
set(gca, 'YTick', 1:numel(sensor_ids), 'YTickLabel', compose('CH%d', sensor_ids), 'XTick', 1:blades_num);
xlabel('Blade ID');
ylabel('Sensor');
title('Valid displacement STD (mm)');
colorbar;

nexttile;
imagesc(1:blades_num, 1:numel(sensor_ids), rms_mat);
set(gca, 'YTick', 1:numel(sensor_ids), 'YTickLabel', compose('CH%d', sensor_ids), 'XTick', 1:blades_num);
xlabel('Blade ID');
ylabel('Sensor');
title('Valid displacement RMS (mm)');
colorbar;

nexttile;
imagesc(1:blades_num, 1:numel(sensor_ids), count_mat);
set(gca, 'YTick', 1:numel(sensor_ids), 'YTickLabel', compose('CH%d', sensor_ids), 'XTick', 1:blades_num);
xlabel('Blade ID');
ylabel('Sensor');
title('Valid point count');
colorbar;
save_step03_figure_local(fig2, figure_dir, cfg, 'Step03_SensorBlade_SummaryHeatmap');

fig3 = figure('Name', sprintf('20250527 Step03 RPM and Coverage - %s', bundle.Case_Name), ...
    'Color', 'w', 'Position', [130, 100, 1300, 720], 'NumberTitle', 'off', ...
    'Visible', visibility_state_local(show_plots));
tiledlayout(fig3, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
hold on;
grid on;
box on;
plot(bundle.Omega_Time_s, bundle.Omega_RPM, 'k-', 'LineWidth', 1.2);
xlabel('Time (s)');
ylabel('RPM');
title(sprintf('OPR-derived RPM trend | %s | speed: %s', bundle.Case_Name, bundle.Speed_Source), 'Interpreter', 'none');

nexttile;
hold on;
grid on;
box on;
t_valid = T.arrival_time_s(T.is_valid);
if ~isempty(t_valid)
    edges = linspace(min(t_valid), max(t_valid), 80);
    for iSensor = 1:numel(sensor_ids)
        sid = sensor_ids(iSensor);
        counts = histcounts(T.arrival_time_s(T.is_valid & T.sensor_id == sid), edges);
        centers = 0.5 * (edges(1:end-1) + edges(2:end));
        plot(centers, counts, '-', 'LineWidth', 1.2, 'Color', colors(iSensor, :), ...
            'DisplayName', sprintf('CH%d', sid));
    end
end
xlabel('Time (s)');
ylabel('Valid arrivals / bin');
title('Observation coverage by sensor');
legend('Location', 'best');
save_step03_figure_local(fig3, figure_dir, cfg, 'Step03_RPM_And_Coverage');

fig4 = figure('Name', sprintf('20250527 Step03 Displacement Distribution - %s', bundle.Case_Name), ...
    'Color', 'w', 'Position', [160, 130, 1300, 760], 'NumberTitle', 'off', ...
    'Visible', visibility_state_local(show_plots));
tiledlayout(fig4, numel(sensor_ids), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for iSensor = 1:numel(sensor_ids)
    sid = sensor_ids(iSensor);
    nexttile;
    hold on;
    grid on;
    box on;
    for blade_id = 1:blades_num
        mask = T.is_valid & T.sensor_id == sid & T.blade_id == blade_id;
        if any(mask)
            histogram(T.displacement_mm(mask), cfg.step03_histogram_bin_count, ...
                'DisplayStyle', 'stairs', 'LineWidth', 1.0, 'DisplayName', sprintf('B%d', blade_id));
        end
    end
    xlabel('Displacement (mm)');
    ylabel('Count');
    title(sprintf('CH%d displacement distribution by blade', sid));
    legend('Location', 'eastoutside');
end
save_step03_figure_local(fig4, figure_dir, cfg, 'Step03_Displacement_Distribution');

fig5 = figure('Name', sprintf('20250527 Step03 Summary Table - %s', bundle.Case_Name), ...
    'Color', 'w', 'Position', [190, 150, 1450, 540], 'NumberTitle', 'off', ...
    'Visible', visibility_state_local(show_plots));
display_names = {'case_name', 'sensor_id', 'blade_id', 'point_count', 'valid_count', ...
    'valid_ratio', 'rpm_mean', 'disp_mean_mm', 'disp_std_mm', 'disp_rms_mm', ...
    'std_angle_std_deg', 'step02_best_corr', 'step02_shift_consistency'};
display_names = display_names(ismember(display_names, S.Properties.VariableNames));
table_for_display = S(:, display_names);
uitable('Data', sanitize_table_cells_local(table2cell(table_for_display)), ...
    'ColumnName', table_for_display.Properties.VariableNames, ...
    'Units', 'normalized', ...
    'Position', [0 0 1 1]);
save_step03_figure_local(fig5, figure_dir, cfg, 'Step03_Summary_Table');

fig6 = figure('Name', sprintf('20250527 Step03 Quality Flags - %s', bundle.Case_Name), ...
    'Color', 'w', 'Position', [220, 180, 1300, 520], 'NumberTitle', 'off', ...
    'Visible', visibility_state_local(show_plots));
plot_quality_flag_counts_local(T);
save_step03_figure_local(fig6, figure_dir, cfg, 'Step03_Quality_Flags');

if ~show_plots
    close([fig1 fig2 fig3 fig4 fig5 fig6]);
end
end


function plot_quality_flag_counts_local(T)
flags = T.quality_flag;
if iscell(flags)
    flags = string(flags);
end
[u, ~, ic] = unique(flags);
counts = accumarray(ic, 1);
bar(categorical(u), counts);
grid on;
box on;
xlabel('Quality flag');
ylabel('Count');
title('Observation quality flags');
xtickangle(30);
end


function [std_mat, rms_mat, count_mat] = summary_to_matrices_local(S, sensor_ids, blades_num)
std_mat = nan(numel(sensor_ids), blades_num);
rms_mat = nan(numel(sensor_ids), blades_num);
count_mat = nan(numel(sensor_ids), blades_num);
for iSensor = 1:numel(sensor_ids)
    sid = sensor_ids(iSensor);
    for blade_id = 1:blades_num
        idx = find(S.sensor_id == sid & S.blade_id == blade_id, 1, 'first');
        if isempty(idx)
            continue;
        end
        std_mat(iSensor, blade_id) = S.disp_std_mm(idx);
        rms_mat(iSensor, blade_id) = S.disp_rms_mm(idx);
        count_mat(iSensor, blade_id) = S.valid_count(idx);
    end
end
end


function [tt, yy] = decimate_series_for_plot_local(t, y, max_points)
t = t(:);
y = y(:);
valid = isfinite(t) & isfinite(y);
t = t(valid);
y = y(valid);
if numel(t) <= max_points
    tt = t;
    yy = y;
    return;
end
idx = unique(round(linspace(1, numel(t), max_points)));
tt = t(idx);
yy = y(idx);
end


function save_step03_figure_local(fig, figure_dir, cfg, tag)
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
    elseif iscell(value)
        if isempty(value)
            cells_out{i} = '';
        else
            cells_out{i} = value{1};
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
