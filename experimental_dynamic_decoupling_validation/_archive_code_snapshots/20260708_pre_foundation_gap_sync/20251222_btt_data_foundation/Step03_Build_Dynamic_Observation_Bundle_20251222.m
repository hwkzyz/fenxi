clc; clear; close all;

%STEP03_BUILD_DYNAMIC_OBSERVATION_BUNDLE_20251222 Build dynamic BTT observations.
%
% Step03 consumes only the new foundation products:
%   Step01: Sensor_Config_20251222.mat
%   Step02: DynamicBTTFeature_20251222.mat
%
% It converts dynamic arrivals into a method-neutral observation bundle:
%   observed angle - low-speed standard angle = angle deviation
%   angle deviation * tip radius = equivalent tip displacement
%
% This step does not run template matching, gap-prior fitting, eta/dx
% calibration, VP screening, or modal identification.

%% 1. Run settings
cfg = BTTDataConfig_20251222();
target_cases = cfg_cellstr_local(cfg.dynamic_cases);
show_tables = true;
show_figures = true;
force_rebuild = true;

%% 2. Step03 settings
Step03_Setting = struct();
Step03_Setting.source_mode = 'step01_step02_foundation_products';
Step03_Setting.standard_angle_source = sprintf('Step01.%s', ...
    get_field_or_default_local(cfg, 'standard_angle_reference_field', 'Standard_Relative_Angles'));
Step03_Setting.observed_angle_source = 'step02_relative_angle_with_oprtable_fallback';
Step03_Setting.rotation_direction = +1;
Step03_Setting.observed_angle_offset_deg = 0;
Step03_Setting.wrap_angle_deviation = true;
Step03_Setting.absolute_outlier_limit_mm = 20;
Step03_Setting.warning_outlier_limit_mm = 2.0;
Step03_Setting.use_robust_outlier_filter = false;
Step03_Setting.robust_mad_k = 6.0;
Step03_Setting.min_points_for_robust_filter = 12;
Step03_Setting.histogram_bin_count = 120;
Step03_Setting.plot_max_points_per_sensor = 8000;
Step03_Setting.save_figures = cfg.save_figures;

%% 3. Load Step01 standard-angle reference
sensor_config_file = fullfile(cfg.step01_output_dir, 'Sensor_Config_20251222.mat');
if ~isfile(sensor_config_file)
    error('Missing Step01 Sensor_Config. Run Step01_Build_LowSpeed_Reference_20251222 first:\n  %s', sensor_config_file);
end
S1 = load(sensor_config_file, 'Sensor_Config', 'AngleQualityTable');
Sensor_Config = S1.Sensor_Config;
AngleQualityTable = table();
if isfield(S1, 'AngleQualityTable')
    AngleQualityTable = S1.AngleQualityTable;
end
StandardAngle = resolve_standard_angles_local(Sensor_Config, cfg);
validate_standard_angles_local(StandardAngle, cfg);
Step03_Setting.standard_angle_source = ['Step01.', StandardAngle.SourceField];

%% 4. Main loop
all_case_summary = table();
for iCase = 1:numel(target_cases)
    case_name = char(target_cases{iCase});
    step02_case_dir = fullfile(cfg.step02_output_dir, case_name);
    step03_case_dir = fullfile(cfg.step03_output_dir, case_name);
    fig_dir = fullfile(cfg.step03_figure_dir, case_name);
    ensure_dir_local(step03_case_dir);
    if Step03_Setting.save_figures
        ensure_dir_local(fig_dir);
    end

    bundle_out = fullfile(step03_case_dir, 'BTT_Observation_Bundle_20251222.mat');
    if ~force_rebuild && isfile(bundle_out)
        fprintf('Step03 output already exists:\n  %s\n', bundle_out);
        continue;
    end

    fprintf('\n=== Step03 observation bundle [%s] ===\n', case_name);
    Step02 = load_step02_feature_local(step02_case_dir);
    DynamicPulseTable = Step02.DynamicPulseTable;
    OPRTable = Step02.OPRTable;
    OmegaTable = Step02.OmegaTable;
    SensorNumberingQuality = Step02.SensorNumberingQuality;
    [RPMReferenceTable, OmegaConsistencyTable] = build_rpm_reference_from_opr_local( ...
        case_name, OPRTable, OmegaTable);

    InputConsistencyReport = check_step03_input_consistency_local( ...
        cfg, Sensor_Config, StandardAngle, Step02);
    [ObservationTable, AngleSourceCheckTable, RevIDConsistencyTable] = build_observation_table_local( ...
        case_name, DynamicPulseTable, OPRTable, RPMReferenceTable, StandardAngle, ...
        AngleQualityTable, SensorNumberingQuality, cfg, Step03_Setting);
    ObservationTable = apply_optional_outlier_filter_local(ObservationTable, cfg, Step03_Setting);

    SummaryTable = summarize_observations_local(case_name, ObservationTable, cfg);
    QualityFlagTable = summarize_quality_flags_local(case_name, ObservationTable);
    Step03Meta = build_step03_metadata_local(cfg, Step03_Setting, case_name, ...
        sensor_config_file, step02_case_dir, StandardAngle, Step02);

    BTT_Observation_Bundle = struct();
    BTT_Observation_Bundle.Dataset = cfg.dataset;
    BTT_Observation_Bundle.CaseName = case_name;
    BTT_Observation_Bundle.SourceMode = Step03_Setting.source_mode;
    BTT_Observation_Bundle.Step01SensorConfigFile = sensor_config_file;
    BTT_Observation_Bundle.Step02CaseDir = step02_case_dir;
    BTT_Observation_Bundle.Step03OutputDir = step03_case_dir;
    BTT_Observation_Bundle.FigureDir = fig_dir;
    BTT_Observation_Bundle.SensorIDs = cfg.sensor_ids;
    BTT_Observation_Bundle.BladesNum = cfg.blades_num;
    BTT_Observation_Bundle.RTipMM = cfg.r_tip_mm;
    BTT_Observation_Bundle.StandardAngle = StandardAngle;
    BTT_Observation_Bundle.OPRTable = OPRTable;
    BTT_Observation_Bundle.OmegaTable = OmegaTable;
    BTT_Observation_Bundle.RPMReferenceTable = RPMReferenceTable;
    BTT_Observation_Bundle.Observation_Table = ObservationTable;
    BTT_Observation_Bundle.ObservationTable = ObservationTable;
    BTT_Observation_Bundle.Summary_Table = SummaryTable;
    BTT_Observation_Bundle.SummaryTable = SummaryTable;
    BTT_Observation_Bundle.QualityFlagTable = QualityFlagTable;
    BTT_Observation_Bundle.AngleSourceCheckTable = AngleSourceCheckTable;
    BTT_Observation_Bundle.RevIDConsistencyTable = RevIDConsistencyTable;
    BTT_Observation_Bundle.OmegaConsistencyTable = OmegaConsistencyTable;
    BTT_Observation_Bundle.InputConsistencyReport = InputConsistencyReport;
    BTT_Observation_Bundle.Step02SensorNumberingQuality = SensorNumberingQuality;
    BTT_Observation_Bundle.Metadata = Step03Meta;

    observation_table = ObservationTable; %#ok<NASGU>
    summary_table = SummaryTable; %#ok<NASGU>
    quality_flag_table = QualityFlagTable; %#ok<NASGU>
    angle_source_check_table = AngleSourceCheckTable; %#ok<NASGU>
    rev_id_consistency_table = RevIDConsistencyTable; %#ok<NASGU>
    omega_consistency_table = OmegaConsistencyTable; %#ok<NASGU>
    input_consistency_report = InputConsistencyReport; %#ok<NASGU>
    metadata = Step03Meta; %#ok<NASGU>

    save(bundle_out, 'BTT_Observation_Bundle', 'observation_table', ...
        'summary_table', 'quality_flag_table', 'angle_source_check_table', ...
        'rev_id_consistency_table', 'omega_consistency_table', ...
        'input_consistency_report', 'metadata', '-v7.3');
    save(fullfile(step03_case_dir, 'Step03_Metadata_20251222.mat'), 'Step03Meta');
    writetable(ObservationTable, fullfile(step03_case_dir, 'BTT_Observation_LongTable_20251222.csv'));
    writetable(SummaryTable, fullfile(step03_case_dir, 'BTT_Observation_Summary_20251222.csv'));
    writetable(QualityFlagTable, fullfile(step03_case_dir, 'BTT_Observation_QualityFlags_20251222.csv'));
    writetable(AngleSourceCheckTable, fullfile(step03_case_dir, 'Step03_AngleSourceCheck_20251222.csv'));
    writetable(RevIDConsistencyTable, fullfile(step03_case_dir, 'Step03_RevIDConsistency_20251222.csv'));
    writetable(OmegaConsistencyTable, fullfile(step03_case_dir, 'Step03_OmegaConsistency_20251222.csv'));
    writetable(InputConsistencyReport, fullfile(step03_case_dir, 'Step03_InputConsistency_20251222.csv'));
    save_compatibility_vib_files_local(step03_case_dir, ObservationTable, cfg, Step03Meta);

    fprintf('Observation rows = %d, valid rows = %d\n', height(ObservationTable), nnz(ObservationTable.is_valid));
    fprintf('Saved bundle:\n  %s\n', bundle_out);

    if show_figures || Step03_Setting.save_figures
        plot_step03_diagnostics_local(BTT_Observation_Bundle, cfg, Step03_Setting, ...
            fig_dir, show_figures, Step03_Setting.save_figures);
    end

    if show_tables
        disp(SummaryTable);
        disp(QualityFlagTable);
        disp(AngleSourceCheckTable);
        disp(RevIDConsistencyTable);
        disp(OmegaConsistencyTable);
    end
    all_case_summary = [all_case_summary; SummaryTable]; %#ok<AGROW>
end

if ~isempty(all_case_summary)
    fprintf('\n=== Step03 all-case summary ===\n');
    disp(all_case_summary);
end

function cases = cfg_cellstr_local(cases)
if ischar(cases) || isstring(cases)
    cases = cellstr(cases);
end
end

function ensure_dir_local(folder)
if exist(folder, 'dir') ~= 7
    mkdir(folder);
end
end

function Step02 = load_step02_feature_local(step02_case_dir)
feature_file = fullfile(step02_case_dir, 'DynamicBTTFeature_20251222.mat');
if ~isfile(feature_file)
    error('Missing Step02 DynamicBTTFeature. Run Step02_Extract_Dynamic_BTT_20251222 first:\n  %s', feature_file);
end
S = load(feature_file);
if isfield(S, 'DynamicBTTFeature')
    F = S.DynamicBTTFeature;
else
    F = struct();
end
Step02 = struct();
Step02.FeatureFile = feature_file;
Step02.DynamicBTTFeature = F;
Step02.DynamicPulseTable = get_feature_table_local(S, F, 'DynamicPulseTable');
Step02.OPRTable = get_feature_table_local(S, F, 'OPRTable');
Step02.OmegaTable = get_feature_table_local(S, F, 'OmegaTable');
Step02.SensorNumberingQuality = get_feature_table_local(S, F, 'SensorNumberingQuality');
Step02.Step02Summary = get_feature_table_local(S, F, 'Step02Summary');
Step02.Metadata = struct();
if isfield(F, 'Metadata')
    Step02.Metadata = F.Metadata;
elseif isfield(S, 'Step02Meta')
    Step02.Metadata = S.Step02Meta;
end
required_tables = {'DynamicPulseTable','OPRTable','OmegaTable'};
for i = 1:numel(required_tables)
    if isempty(Step02.(required_tables{i})) || ~istable(Step02.(required_tables{i}))
        error('Step02 feature lacks required table %s:\n  %s', required_tables{i}, feature_file);
    end
end
end

function T = get_feature_table_local(S, F, name)
T = table();
if isfield(S, name) && istable(S.(name))
    T = S.(name);
elseif isstruct(F) && isfield(F, name) && istable(F.(name))
    T = F.(name);
end
end

function StandardAngle = resolve_standard_angles_local(Sensor_Config, cfg)
preferred_field = get_field_or_default_local(cfg, ...
    'standard_angle_reference_field', 'Standard_Relative_Angles');
if isfield(Sensor_Config, preferred_field)
    A = Sensor_Config.(preferred_field);
    source_name = preferred_field;
elseif isfield(Sensor_Config, 'Standard_Relative_Angles_OPRCenter')
    A = Sensor_Config.Standard_Relative_Angles_OPRCenter;
    source_name = 'Standard_Relative_Angles_OPRCenter';
elseif isfield(Sensor_Config, 'Standard_Relative_Angles')
    A = Sensor_Config.Standard_Relative_Angles;
    source_name = 'Standard_Relative_Angles';
else
    error('Sensor_Config has no usable standard angle matrix.');
end

StandardAngle = struct();
StandardAngle.Value = A;
StandardAngle.SourceField = source_name;
StandardAngle.Reference = get_field_or_default_local(Sensor_Config, ...
    'Standard_Relative_Angles_Reference', 'unknown');
StandardAngle.OPRTimingMethod = get_field_or_default_local(Sensor_Config, ...
    'OPR_Timing_Method', 'unknown');
StandardAngle.SelectedStatistic = get_field_or_default_local(Sensor_Config, ...
    'Standard_Relative_Angles_SelectedStatistic', 'unknown');
StandardAngle.SensorIDs = get_field_or_default_local(Sensor_Config, 'Sensor_IDs', cfg.sensor_ids);
StandardAngle.BladesNum = get_field_or_default_local(Sensor_Config, 'Blades_Num', cfg.blades_num);
end

function validate_standard_angles_local(StandardAngle, cfg)
if size(StandardAngle.Value, 2) < cfg.blades_num
    error('Standard angle matrix has %d blade columns, expected at least %d.', ...
        size(StandardAngle.Value, 2), cfg.blades_num);
end
for sid = cfg.sensor_ids(:).'
    row_idx = standard_angle_row_index_local(StandardAngle, sid);
    if ~isfinite(row_idx) || row_idx < 1 || row_idx > size(StandardAngle.Value, 1)
        error('Standard angle matrix has no row for sensor %d.', sid);
    end
end
end

function T = check_step03_input_consistency_local(cfg, Sensor_Config, StandardAngle, Step02)
rows = [];
sensor_ids_in_config = get_field_or_default_local(Sensor_Config, 'Sensor_IDs', []);
rows = [rows; consistency_row_local('Sensor_Config.Sensor_IDs', ...
    mat2str(sensor_ids_in_config), ...
    mat2str(cfg.sensor_ids), isequal(sort(sensor_ids_in_config(:).'), sort(cfg.sensor_ids(:).')))]; %#ok<AGROW>
rows = [rows; consistency_row_local('Sensor_Config.Blades_Num', ...
    string(get_field_or_default_local(Sensor_Config, 'Blades_Num', NaN)), ...
    string(cfg.blades_num), get_field_or_default_local(Sensor_Config, 'Blades_Num', NaN) == cfg.blades_num)]; %#ok<AGROW>
rows = [rows; consistency_row_local('Sensor_Config.OPR_ID', ...
    string(get_field_or_default_local(Sensor_Config, 'OPR_ID', NaN)), ...
    string(cfg.opr_id), get_field_or_default_local(Sensor_Config, 'OPR_ID', NaN) == cfg.opr_id)]; %#ok<AGROW>
rows = [rows; consistency_row_local('Sensor_Config.Sample_Rate_Hz', ...
    string(get_field_or_default_local(Sensor_Config, 'Sample_Rate_Hz', NaN)), ...
    string(cfg.sample_rate_hz), abs(get_field_or_default_local(Sensor_Config, 'Sample_Rate_Hz', NaN) - cfg.sample_rate_hz) < eps(cfg.sample_rate_hz))]; %#ok<AGROW>

step01_ref = string(StandardAngle.Reference);
step02_timing_reference = string(get_field_or_default_local(Step02.Metadata, ...
    'OPRTimingReference', get_field_or_default_local(Step02.Metadata, ...
    'OPRCenterMethod', 'unknown')));
rows = [rows; consistency_row_local('Step01 standard-angle reference', ...
    step01_ref, 'known OPR timing reference', contains(lower(step01_ref), 'opr'))]; %#ok<AGROW>
rows = [rows; consistency_row_local('Step02 OPR timing reference', ...
    step02_timing_reference, 'known timing reference', step02_timing_reference ~= "unknown")]; %#ok<AGROW>

events_per_rev = unique(Step02.OPRTable.opr_events_per_revolution(isfinite(Step02.OPRTable.opr_events_per_revolution)));
rows = [rows; consistency_row_local('Step02 OPR events per revolution', ...
    mat2str(events_per_rev(:).'), 'single finite value', numel(events_per_rev) == 1)]; %#ok<AGROW>

T = struct2table(rows);
if any(~T.is_consistent)
    warning('Step03:InputConsistency', ...
        'Step03 input consistency report contains warnings. Inspect Step03_InputConsistency_20251222.csv.');
end
end

function row = consistency_row_local(check_name, actual_value, expected_value, is_consistent)
row = struct();
row.check_name = string(check_name);
row.actual_value = string(actual_value);
row.expected_value = string(expected_value);
row.is_consistent = logical(is_consistent);
if row.is_consistent
    row.status = "ok";
else
    row.status = "check";
end
end

function [RPMReferenceTable, OmegaConsistencyTable] = build_rpm_reference_from_opr_local(case_name, OPRTable, OmegaTable)
events_per_rev = first_finite_local(OPRTable.opr_events_per_revolution);
opr_t = resolve_opr_reference_times_from_table_local(OPRTable);
if ~isfinite(events_per_rev) || events_per_rev < 1 || numel(opr_t) <= events_per_rev
    error('Cannot recompute RPM from OPRTable: invalid opr_events_per_revolution.');
end
idx = (1:(numel(opr_t) - events_per_rev)).';
rev_period_s = opr_t(idx + events_per_rev) - opr_t(idx);
rpm = 60 ./ rev_period_s;
omega_rad_s = 2 * pi ./ rev_period_s;
RPMReferenceTable = table( ...
    repmat(string(case_name), numel(idx), 1), idx, opr_t(idx), rev_period_s, ...
    rpm, omega_rad_s, repmat(events_per_rev, numel(idx), 1), ...
    repmat(string('recomputed_from_step02_OPRTable'), numel(idx), 1), ...
    'VariableNames', {'case_name','rpm_index','omega_time_s','rev_period_s', ...
    'omega_rpm','omega_rad_s','opr_events_per_revolution','rpm_source'});

OmegaConsistencyTable = compare_opr_rpm_with_omega_local(case_name, RPMReferenceTable, OmegaTable);
end

function T = compare_opr_rpm_with_omega_local(case_name, RPMReferenceTable, OmegaTable)
if isempty(OmegaTable) || ~istable(OmegaTable) || height(OmegaTable) < 2
    T = table(string(case_name), 0, NaN, NaN, NaN, NaN, string("missing_omega_table"), ...
        'VariableNames', {'case_name','comparison_count','median_abs_rpm_diff', ...
        'max_abs_rpm_diff','median_signed_rpm_diff','rpm_diff_ratio_median','status'});
    return;
end
valid = isfinite(OmegaTable.omega_time_s) & isfinite(OmegaTable.omega_rpm);
t = OmegaTable.omega_time_s(valid);
rpm = OmegaTable.omega_rpm(valid);
[t, order] = sort(t(:));
rpm = rpm(order);
[t, ia] = unique(t, 'stable');
rpm = rpm(ia);
if numel(t) < 2
    T = table(string(case_name), 0, NaN, NaN, NaN, NaN, string("too_few_unique_omega_samples"), ...
        'VariableNames', {'case_name','comparison_count','median_abs_rpm_diff', ...
        'max_abs_rpm_diff','median_signed_rpm_diff','rpm_diff_ratio_median','status'});
    return;
end
F = griddedInterpolant(t, rpm(:), 'linear', 'nearest');
omega_on_ref = F(RPMReferenceTable.omega_time_s);
diff_rpm = RPMReferenceTable.omega_rpm - omega_on_ref;
valid_diff = isfinite(diff_rpm) & isfinite(omega_on_ref) & abs(omega_on_ref) > eps;
median_abs_diff = median_or_nan_local(abs(diff_rpm(valid_diff)));
ratio_med = median_or_nan_local(abs(diff_rpm(valid_diff)) ./ abs(omega_on_ref(valid_diff)));
status = "ok";
if ratio_med > 1e-3
    status = "check";
end
T = table(string(case_name), nnz(valid_diff), median_abs_diff, ...
    first_or_nan_local(abs(diff_rpm(valid_diff)), 'max'), median_or_nan_local(diff_rpm(valid_diff)), ...
    ratio_med, status, ...
    'VariableNames', {'case_name','comparison_count','median_abs_rpm_diff', ...
    'max_abs_rpm_diff','median_signed_rpm_diff','rpm_diff_ratio_median','status'});
end

function [ObservationTable, AngleSourceCheckTable, RevIDConsistencyTable] = build_observation_table_local(case_name, PulseTable, OPRTable, RPMReferenceTable, StandardAngle, AngleQualityTable, SensorNumberingQuality, cfg, Step03_Setting)
T = PulseTable;
n = height(T);
standard_angle_deg = nan(n, 1);
observed_angle_deg = nan(n, 1);
fallback_observed_angle_deg = nan(n, 1);
angle_deviation_deg = nan(n, 1);
displacement_mm = nan(n, 1);
rpm = nan(n, 1);
rev_id_from_opr = nan(n, 1);
event_in_rev_from_opr = nan(n, 1);
prev_opr_event_id_from_opr = nan(n, 1);
rev_id_mismatch = false(n, 1);
event_in_rev_mismatch = false(n, 1);
std_angle_quality = strings(n, 1);
quality_flag = strings(n, 1);
quality_flag(:) = "valid";

F_rpm = build_rpm_interpolant_local(RPMReferenceTable);
[prev_opr_event_id_from_opr, rev_id_from_opr, event_in_rev_from_opr, fallback_observed_angle_deg] = ...
    derive_opr_reference_for_arrivals_local(T.arrival_time_s, OPRTable, Step03_Setting);
if any(strcmp('rev_id', T.Properties.VariableNames))
    rev_id_mismatch = isfinite(T.rev_id) & isfinite(rev_id_from_opr) & T.rev_id ~= rev_id_from_opr;
end
if any(strcmp('event_in_rev', T.Properties.VariableNames))
    event_in_rev_mismatch = isfinite(T.event_in_rev) & isfinite(event_in_rev_from_opr) & ...
        T.event_in_rev ~= event_in_rev_from_opr;
end

for i = 1:n
    sid = T.sensor_id(i);
    bid = T.assigned_blade_id(i);
    if ~isfinite(T.arrival_time_s(i))
        quality_flag(i) = "missing_arrival_time";
        continue;
    end
    if ~isfinite(sid) || ~ismember(sid, cfg.sensor_ids)
        quality_flag(i) = "invalid_sensor_id";
        continue;
    end
    if ~isfinite(bid) || bid < 1 || bid > cfg.blades_num
        quality_flag(i) = "invalid_blade_id";
        continue;
    end

    row_idx = standard_angle_row_index_local(StandardAngle, sid);
    theta_std = StandardAngle.Value(row_idx, bid);
    if ~isfinite(theta_std)
        quality_flag(i) = "missing_standard_angle";
        continue;
    end

    theta_obs = NaN;
    if any(strcmp('relative_angle_deg', T.Properties.VariableNames))
        theta_obs = T.relative_angle_deg(i);
    end
    if ~isfinite(theta_obs)
        theta_obs = fallback_observed_angle_deg(i);
    end
    if ~isfinite(theta_obs)
        quality_flag(i) = "missing_observed_angle";
        continue;
    end
    theta_obs = normalize_observed_angle_local(theta_obs, Step03_Setting);

    theta_err = theta_obs - theta_std;
    if Step03_Setting.wrap_angle_deviation
        theta_err = wrap_to_180_local(theta_err);
    end
    disp_mm = theta_err * (pi / 180) * cfg.r_tip_mm;

    standard_angle_deg(i) = theta_std;
    observed_angle_deg(i) = theta_obs;
    angle_deviation_deg(i) = theta_err;
    displacement_mm(i) = disp_mm;
    rpm(i) = F_rpm(T.arrival_time_s(i));
    std_angle_quality(i) = lookup_angle_quality_local(AngleQualityTable, sid, bid);

    if any(strcmp('is_valid', T.Properties.VariableNames)) && ~T.is_valid(i)
        quality_flag(i) = "invalid_step02_assignment";
    elseif ~isfinite(disp_mm)
        quality_flag(i) = "nonfinite_displacement";
    elseif abs(disp_mm) > Step03_Setting.absolute_outlier_limit_mm
        quality_flag(i) = "absolute_outlier_displacement";
    end
end

is_valid = quality_flag == "valid";
if any(strcmp('is_valid', T.Properties.VariableNames))
    is_valid = is_valid & T.is_valid;
end

step02_quality = attach_sensor_quality_local(T, SensorNumberingQuality);
ObservationTable = table( ...
    repmat(string(case_name), n, 1), T.sensor_id, T.assigned_blade_id, ...
    T.rev_id, T.event_in_rev, T.prev_opr_event_id, T.source_row, ...
    rev_id_from_opr, event_in_rev_from_opr, prev_opr_event_id_from_opr, ...
    rev_id_mismatch, event_in_rev_mismatch, ...
    T.pulse_start_time_s, T.pulse_end_time_s, T.arrival_time_s, ...
    rpm, standard_angle_deg, observed_angle_deg, fallback_observed_angle_deg, ...
    angle_deviation_deg, displacement_mm, abs(displacement_mm) > Step03_Setting.warning_outlier_limit_mm, ...
    T.pulse_peak_value, T.pulse_feature_value, T.assignment_quality, ...
    step02_quality.sensor_best_corr_median, step02_quality.sensor_corr_gap_median, ...
    string(T.quality_flags), std_angle_quality, quality_flag, is_valid, ...
    repmat(string(Step03_Setting.observed_angle_source), n, 1), ...
    repmat(string(Step03_Setting.standard_angle_source), n, 1), ...
    'VariableNames', {'case_name','sensor_id','blade_id','rev_id','event_in_rev', ...
    'prev_opr_event_id','pulse_index','rev_id_from_opr','event_in_rev_from_opr', ...
    'prev_opr_event_id_from_opr','rev_id_mismatch','event_in_rev_mismatch', ...
    'pulse_start_time_s','pulse_end_time_s','arrival_time_s','rpm', ...
    'standard_angle_deg','observed_angle_deg','fallback_observed_angle_deg', ...
    'angle_deviation_deg','displacement_mm','large_displacement_warning', ...
    'pulse_peak_value','pulse_feature_value', ...
    'step02_assignment_quality','step02_sensor_best_corr_median', ...
    'step02_sensor_corr_gap_median','step02_quality_flag','standard_angle_quality', ...
    'quality_flag','is_valid','observed_angle_source','standard_angle_source'});

AngleSourceCheckTable = build_angle_source_check_table_local(case_name, T, ...
    fallback_observed_angle_deg, OPRTable, cfg, Step03_Setting);
RevIDConsistencyTable = build_rev_id_consistency_table_local(case_name, ObservationTable, cfg);
warn_if_diagnostics_bad_local(AngleSourceCheckTable, RevIDConsistencyTable);
end

function row_idx = standard_angle_row_index_local(StandardAngle, sid)
row_idx = NaN;
sensor_ids = StandardAngle.SensorIDs(:).';
idx = find(sensor_ids == sid, 1, 'first');
if ~isempty(idx) && idx <= size(StandardAngle.Value, 1)
    row_idx = idx;
elseif sid <= size(StandardAngle.Value, 1)
    row_idx = sid;
end
end

function F_rpm = build_rpm_interpolant_local(OmegaTable)
valid = isfinite(OmegaTable.omega_time_s) & isfinite(OmegaTable.omega_rpm);
t = OmegaTable.omega_time_s(valid);
rpm = OmegaTable.omega_rpm(valid);
[t, order] = sort(t(:));
rpm = rpm(order);
[t, ia] = unique(t, 'stable');
rpm = rpm(ia);
if numel(t) < 2
    error('OmegaTable has too few valid unique RPM samples.');
end
F_rpm = griddedInterpolant(t, rpm(:), 'linear', 'nearest');
end

function [prev_event_id, rev_id, event_in_rev, theta_obs] = derive_opr_reference_for_arrivals_local(arrival_time_s, OPRTable, Step03_Setting)
n = numel(arrival_time_s);
prev_event_id = nan(n, 1);
rev_id = nan(n, 1);
event_in_rev = nan(n, 1);
theta_obs = nan(n, 1);
opr_times = resolve_opr_reference_times_from_table_local(OPRTable);
for i = 1:n
    t = arrival_time_s(i);
    if ~isfinite(t)
        continue;
    end
    idx = find(opr_times <= t, 1, 'last');
    if isempty(idx) || idx < 1
        continue;
    end
    events_per_rev = OPRTable.opr_events_per_revolution(idx);
    idx_next = idx + events_per_rev;
    prev_event_id(i) = OPRTable.event_id(idx);
    rev_id(i) = OPRTable.rev_id(idx);
    event_in_rev(i) = OPRTable.event_in_rev(idx);
    if idx_next > height(OPRTable)
        continue;
    end
    rev_period = opr_times(idx_next) - opr_times(idx);
    if rev_period <= 0
        continue;
    end
    theta_local = 360 * (t - opr_times(idx)) / rev_period;
    theta_obs(i) = normalize_observed_angle_local(theta_local, Step03_Setting);
end
end

function opr_times = resolve_opr_reference_times_from_table_local(OPRTable)
if ismember('opr_reference_time_s', OPRTable.Properties.VariableNames)
    opr_times = OPRTable.opr_reference_time_s;
elseif ismember('opr_start_time_s', OPRTable.Properties.VariableNames)
    opr_times = OPRTable.opr_start_time_s;
else
    opr_times = OPRTable.opr_center_time_s;
end
end

function theta = normalize_observed_angle_local(theta, Step03_Setting)
theta = Step03_Setting.rotation_direction * theta + Step03_Setting.observed_angle_offset_deg;
theta = mod(theta, 360);
end

function T = build_angle_source_check_table_local(case_name, PulseTable, fallback_angle_deg, OPRTable, cfg, Step03_Setting)
rows = [];
has_relative = any(strcmp('relative_angle_deg', PulseTable.Properties.VariableNames));
events_per_rev = first_finite_local(OPRTable.opr_events_per_revolution);
for sid = cfg.sensor_ids(:).'
    mask = PulseTable.sensor_id == sid;
    rel = nan(nnz(mask), 1);
    if has_relative
        rel = PulseTable.relative_angle_deg(mask);
        rel(isfinite(rel)) = normalize_observed_angle_local(rel(isfinite(rel)), Step03_Setting);
    end
    fb = fallback_angle_deg(mask);
    both = isfinite(rel) & isfinite(fb);
    diff_deg = wrap_to_180_local(rel(both) - fb(both));
    rows = [rows; struct( ... %#ok<AGROW>
        'case_name', string(case_name), ...
        'sensor_id', sid, ...
        'opr_events_per_revolution', events_per_rev, ...
        'has_relative_angle_deg', has_relative, ...
        'valid_relative_angle_ratio', safe_ratio_local(nnz(isfinite(rel)), max(nnz(mask), 1)), ...
        'relative_angle_min', first_or_nan_local(rel, 'min'), ...
        'relative_angle_max', first_or_nan_local(rel, 'max'), ...
        'fallback_angle_valid_ratio', safe_ratio_local(nnz(isfinite(fb)), max(nnz(mask), 1)), ...
        'median_abs_relative_vs_fallback_deg', median_or_nan_local(abs(diff_deg)), ...
        'median_relative_minus_fallback_deg', median_or_nan_local(diff_deg), ...
        'angle_source_status', angle_source_status_local(has_relative, rel, fb, diff_deg))];
end
T = struct2table(rows);
end

function status = angle_source_status_local(has_relative, rel, fb, diff_deg)
if ~has_relative
    status = "fallback_only";
elseif safe_ratio_local(nnz(isfinite(rel)), numel(rel)) < 0.99
    status = "low_relative_angle_coverage";
elseif safe_ratio_local(nnz(isfinite(fb)), numel(fb)) < 0.99
    status = "low_fallback_angle_coverage";
elseif median_or_nan_local(abs(diff_deg)) > 0.5
    status = "relative_fallback_mismatch";
else
    status = "ok";
end
end

function T = build_rev_id_consistency_table_local(case_name, ObservationTable, cfg)
rows = [];
for sid = cfg.sensor_ids(:).'
    mask = ObservationTable.sensor_id == sid;
    rows = [rows; struct( ... %#ok<AGROW>
        'case_name', string(case_name), ...
        'sensor_id', sid, ...
        'row_count', nnz(mask), ...
        'rev_id_mismatch_count', nnz(ObservationTable.rev_id_mismatch(mask)), ...
        'rev_id_mismatch_ratio', safe_ratio_local(nnz(ObservationTable.rev_id_mismatch(mask)), nnz(mask)), ...
        'event_in_rev_mismatch_count', nnz(ObservationTable.event_in_rev_mismatch(mask)), ...
        'event_in_rev_mismatch_ratio', safe_ratio_local(nnz(ObservationTable.event_in_rev_mismatch(mask)), nnz(mask)))];
end
T = struct2table(rows);
end

function warn_if_diagnostics_bad_local(AngleSourceCheckTable, RevIDConsistencyTable)
bad_angle = AngleSourceCheckTable.angle_source_status ~= "ok";
if any(bad_angle)
    warning('Step03:AngleSourceCheck', ...
        'One or more sensors have weak relative-angle diagnostics. Inspect Step03_AngleSourceCheck_20251222.csv.');
end
if any(RevIDConsistencyTable.rev_id_mismatch_ratio > 0.01)
    warning('Step03:RevIDMismatch', ...
        'Step02 rev_id and OPR-derived rev_id mismatch ratio exceeds 1%%. Inspect Step03_RevIDConsistency_20251222.csv.');
end
end

function label = lookup_angle_quality_local(AngleQualityTable, sid, bid)
label = "unknown";
if isempty(AngleQualityTable) || ~istable(AngleQualityTable)
    return;
end
names = AngleQualityTable.Properties.VariableNames;
if all(ismember({'SensorID','BladeID'}, names))
    idx = find(AngleQualityTable.SensorID == sid & AngleQualityTable.BladeID == bid, 1, 'first');
elseif all(ismember({'sensor_id','blade_id'}, names))
    idx = find(AngleQualityTable.sensor_id == sid & AngleQualityTable.blade_id == bid, 1, 'first');
else
    idx = [];
end
if isempty(idx)
    return;
end
if any(strcmp('QualityStatus', names))
    label = string(AngleQualityTable.QualityStatus(idx));
elseif any(strcmp('quality_status', names))
    label = string(AngleQualityTable.quality_status(idx));
end
end

function Q = attach_sensor_quality_local(PulseTable, SensorNumberingQuality)
n = height(PulseTable);
Q = struct();
Q.sensor_best_corr_median = nan(n, 1);
Q.sensor_corr_gap_median = nan(n, 1);
if isempty(SensorNumberingQuality) || ~istable(SensorNumberingQuality)
    return;
end
for i = 1:height(SensorNumberingQuality)
    sid = SensorNumberingQuality.sensor_id(i);
    mask = PulseTable.sensor_id == sid;
    if any(strcmp('best_corr_median', SensorNumberingQuality.Properties.VariableNames))
        Q.sensor_best_corr_median(mask) = SensorNumberingQuality.best_corr_median(i);
    end
    if any(strcmp('corr_gap_median', SensorNumberingQuality.Properties.VariableNames))
        Q.sensor_corr_gap_median(mask) = SensorNumberingQuality.corr_gap_median(i);
    end
end
end

function ObservationTable = apply_optional_outlier_filter_local(ObservationTable, cfg, Step03_Setting)
if ~Step03_Setting.use_robust_outlier_filter
    return;
end
for sid = cfg.sensor_ids(:).'
    for bid = 1:cfg.blades_num
        mask = ObservationTable.sensor_id == sid & ObservationTable.blade_id == bid & ObservationTable.is_valid;
        idx = find(mask);
        if numel(idx) < Step03_Setting.min_points_for_robust_filter
            continue;
        end
        y = ObservationTable.displacement_mm(idx);
        med_y = median(y, 'omitnan');
        mad_y = median(abs(y - med_y), 'omitnan');
        sigma = 1.4826 * mad_y;
        if ~isfinite(sigma) || sigma < eps
            continue;
        end
        bad = abs(y - med_y) > Step03_Setting.robust_mad_k * sigma;
        if any(bad)
            bad_idx = idx(bad);
            ObservationTable.is_valid(bad_idx) = false;
            ObservationTable.quality_flag(bad_idx) = "robust_outlier_displacement";
        end
    end
end
end

function SummaryTable = summarize_observations_local(case_name, ObservationTable, cfg)
rows = [];
for sid = cfg.sensor_ids(:).'
    for bid = 1:cfg.blades_num
        mask_all = ObservationTable.sensor_id == sid & ObservationTable.blade_id == bid;
        mask_valid = mask_all & ObservationTable.is_valid;
        y = ObservationTable.displacement_mm(mask_valid);
        rpm = ObservationTable.rpm(mask_valid);
        t = ObservationTable.arrival_time_s(mask_valid);
        rows = [rows; struct( ... %#ok<AGROW>
            'case_name', string(case_name), ...
            'sensor_id', sid, ...
            'blade_id', bid, ...
            'point_count', nnz(mask_all), ...
            'valid_count', nnz(mask_valid), ...
            'valid_ratio', safe_ratio_local(nnz(mask_valid), nnz(mask_all)), ...
            'large_disp_warning_count', nnz(mask_all & ObservationTable.large_displacement_warning), ...
            'large_disp_warning_ratio', safe_ratio_local(nnz(mask_all & ObservationTable.large_displacement_warning), nnz(mask_all)), ...
            'time_start_s', first_or_nan_local(t, 'min'), ...
            'time_end_s', first_or_nan_local(t, 'max'), ...
            'rpm_mean', mean_or_nan_local(rpm), ...
            'rpm_min', first_or_nan_local(rpm, 'min'), ...
            'rpm_max', first_or_nan_local(rpm, 'max'), ...
            'disp_mean_mm', mean_or_nan_local(y), ...
            'disp_median_mm', median_or_nan_local(y), ...
            'disp_std_mm', std_or_nan_local(y), ...
            'disp_rms_mm', rms_or_nan_local(y), ...
            'disp_min_mm', first_or_nan_local(y, 'min'), ...
            'disp_max_mm', first_or_nan_local(y, 'max'))];
    end
end
SummaryTable = struct2table(rows);
end

function T = summarize_quality_flags_local(case_name, ObservationTable)
flags = unique(ObservationTable.quality_flag);
rows = [];
for i = 1:numel(flags)
    flag = flags(i);
    count = nnz(ObservationTable.quality_flag == flag);
    rows = [rows; struct( ... %#ok<AGROW>
        'case_name', string(case_name), ...
        'quality_flag', string(flag), ...
        'count', count, ...
        'ratio', count / max(height(ObservationTable), 1))];
end
T = struct2table(rows);
end

function Meta = build_step03_metadata_local(cfg, Step03_Setting, case_name, sensor_config_file, step02_case_dir, StandardAngle, Step02)
Meta = struct();
Meta.Dataset = cfg.dataset;
Meta.CaseName = case_name;
Meta.CreatedBy = mfilename;
Meta.CreatedAt = string(datetime('now'));
Meta.SourceMode = Step03_Setting.source_mode;
Meta.Step01SensorConfigFile = sensor_config_file;
Meta.Step02FeatureFile = Step02.FeatureFile;
Meta.Step02CaseDir = step02_case_dir;
Meta.StandardAngleSourceField = StandardAngle.SourceField;
Meta.StandardAngleReference = StandardAngle.Reference;
Meta.Step01OPRTimingMethod = StandardAngle.OPRTimingMethod;
Meta.Step02SourceMode = get_field_or_default_local(Step02.Metadata, 'SourceMode', 'unknown');
Meta.Step02OPRCenterMethod = get_field_or_default_local(Step02.Metadata, 'OPRCenterMethod', 'unknown');
Meta.ObservedAngleSource = Step03_Setting.observed_angle_source;
Meta.RPMSource = 'recomputed_from_step02_OPRTable';
Meta.FallbackAngleSemantics = 'local angle after the previous OPR event, using a one-revolution period from opr_events_per_revolution';
Meta.RotationDirection = Step03_Setting.rotation_direction;
Meta.ObservedAngleOffsetDeg = Step03_Setting.observed_angle_offset_deg;
Meta.WrapAngleDeviation = Step03_Setting.wrap_angle_deviation;
Meta.AbsoluteOutlierLimitMM = Step03_Setting.absolute_outlier_limit_mm;
Meta.WarningOutlierLimitMM = Step03_Setting.warning_outlier_limit_mm;
Meta.RTipMM = cfg.r_tip_mm;
Meta.DisplacementDefinition = 'displacement_mm = wrapTo180(observed_angle_deg - standard_angle_deg) * pi/180 * r_tip_mm';
Meta.DisplacementPhysicalMeaning = 'Tangential equivalent tip displacement from angular timing deviation, not radial clearance displacement.';
Meta.jilublade_vib_final_columns = { ...
    'pulse_start_time_s', 'pulse_end_time_s', 'arrival_time_s', 'blade_id', ...
    'observed_angle_deg', 'displacement_mm', 'standard_angle_deg', ...
    'angle_deviation_deg', 'prev_opr_event_id', 'rev_id', 'rpm', 'is_valid'};
Meta.Note = 'Step03 builds method-neutral observation rows from Step01 standard angles and Step02 raw dynamic arrivals.';
end

function save_compatibility_vib_files_local(step03_case_dir, ObservationTable, cfg, Step03Meta)
metadata = Step03Meta; %#ok<NASGU>
for sid = cfg.sensor_ids(:).'
    rows = ObservationTable(ObservationTable.sensor_id == sid, :);
    jilublade = [rows.pulse_start_time_s, rows.pulse_end_time_s, rows.arrival_time_s, rows.blade_id, ...
        rows.observed_angle_deg, rows.displacement_mm, rows.standard_angle_deg, rows.angle_deviation_deg, ...
        rows.prev_opr_event_id, rows.rev_id, rows.rpm, rows.is_valid]; %#ok<NASGU>
    save(fullfile(step03_case_dir, sprintf('jilublade_probe%d_vib_final.mat', sid)), ...
        'jilublade', 'metadata');
end
end

function plot_step03_diagnostics_local(bundle, cfg, Step03_Setting, fig_dir, show_figures, save_figures)
visibility = 'off';
if show_figures
    visibility = 'on';
end
T = bundle.Observation_Table;

fig = figure('Name', 'Step03 Displacement Trend', 'Color', 'w', 'Visible', visibility, 'Position', [80 80 1350 800]);
tiledlayout(fig, numel(cfg.sensor_ids), 1, 'Padding', 'compact', 'TileSpacing', 'compact');
colors = lines(cfg.blades_num);
for i = 1:numel(cfg.sensor_ids)
    sid = cfg.sensor_ids(i);
    ax = nexttile;
    hold(ax, 'on');
    rows = T(T.sensor_id == sid & T.is_valid, :);
    rows = select_evenly_for_plot_local(rows, Step03_Setting.plot_max_points_per_sensor);
    for bid = 1:cfg.blades_num
        mask = rows.blade_id == bid;
        scatter(ax, rows.arrival_time_s(mask), rows.displacement_mm(mask), 5, colors(bid, :), 'filled');
    end
    grid(ax, 'on');
    xlabel(ax, 'Arrival time (s)');
    ylabel(ax, sprintf('CH%d disp (mm)', sid));
    title(ax, sprintf('CH%d observed displacement by blade', sid));
end
save_figure_local(fig, fig_dir, 'Step03_Displacement_Trend_20251222', save_figures);

fig = figure('Name', 'Step03 Angle Deviation Trend', 'Color', 'w', 'Visible', visibility, 'Position', [100 100 1350 800]);
tiledlayout(fig, numel(cfg.sensor_ids), 1, 'Padding', 'compact', 'TileSpacing', 'compact');
for i = 1:numel(cfg.sensor_ids)
    sid = cfg.sensor_ids(i);
    rows = T(T.sensor_id == sid & T.is_valid, :);
    rows = select_evenly_for_plot_local(rows, Step03_Setting.plot_max_points_per_sensor);
    ax = nexttile;
    hold(ax, 'on');
    plot_handles = gobjects(cfg.blades_num, 1);
    for bid = 1:cfg.blades_num
        mask = rows.blade_id == bid;
        plot_handles(bid) = scatter(ax, rows.arrival_time_s(mask), rows.angle_deviation_deg(mask), ...
            5, colors(bid, :), 'filled');
    end
    grid(ax, 'on');
    xlabel(ax, 'Arrival time (s)');
    ylabel(ax, sprintf('CH%d angle error (deg)', sid));
    title(ax, sprintf('CH%d angle deviation: observed - standard', sid));
    if i == 1
        legend(ax, plot_handles, compose('B%d', 1:cfg.blades_num), 'Location', 'eastoutside');
    end
end
save_figure_local(fig, fig_dir, 'Step03_AngleDeviation_Trend_20251222', save_figures);

fig = figure('Name', 'Step03 RPM and Coverage', 'Color', 'w', 'Visible', visibility, 'Position', [120 120 1250 650]);
tiledlayout(fig, 2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
ax = nexttile;
plot(ax, bundle.RPMReferenceTable.omega_time_s, bundle.RPMReferenceTable.omega_rpm, 'k-', 'LineWidth', 1.1);
grid(ax, 'on');
xlabel(ax, 'Time (s)');
ylabel(ax, 'RPM');
title(ax, 'Step03 RPM reference recomputed from Step02 OPRTable');
ax = nexttile;
valid_t = T.arrival_time_s(T.is_valid);
if ~isempty(valid_t)
    edges = linspace(min(valid_t), max(valid_t), 80);
    centers = 0.5 * (edges(1:end-1) + edges(2:end));
    hold(ax, 'on');
    line_styles = {'-', '--', ':'};
    for k = 1:numel(cfg.sensor_ids)
        sid = cfg.sensor_ids(k);
        y = T.arrival_time_s(T.is_valid & T.sensor_id == sid);
        counts = histcounts(y, edges);
        stairs(ax, centers, counts, 'LineWidth', 1.2, ...
            'LineStyle', line_styles{1 + mod(k - 1, numel(line_styles))});
    end
    legend(ax, compose('CH%d', cfg.sensor_ids), 'Location', 'best');
end
grid(ax, 'on');
xlabel(ax, 'Arrival time (s)');
ylabel(ax, 'Valid observations/bin');
title(ax, 'Observation coverage by sensor');
save_figure_local(fig, fig_dir, 'Step03_RPM_And_Coverage_20251222', save_figures);

fig = figure('Name', 'Step03 Displacement Distribution', 'Color', 'w', 'Visible', visibility, 'Position', [140 140 1200 720]);
tiledlayout(fig, numel(cfg.sensor_ids), 1, 'Padding', 'compact', 'TileSpacing', 'compact');
for i = 1:numel(cfg.sensor_ids)
    sid = cfg.sensor_ids(i);
    ax = nexttile;
    y = T.displacement_mm(T.sensor_id == sid & T.is_valid);
    histogram(ax, y, Step03_Setting.histogram_bin_count, 'FaceColor', [0.2 0.45 0.75], 'EdgeColor', 'none');
    grid(ax, 'on');
    xlabel(ax, 'Displacement (mm)');
    ylabel(ax, 'Count');
    title(ax, sprintf('CH%d displacement distribution', sid));
end
save_figure_local(fig, fig_dir, 'Step03_Displacement_Distribution_20251222', save_figures);

fig = figure('Name', 'Step03 Quality Flags', 'Color', 'w', 'Visible', visibility, 'Position', [160 160 1050 520]);
ax = axes(fig);
Q = bundle.QualityFlagTable;
bar(ax, categorical(Q.quality_flag), Q.count);
grid(ax, 'on');
ylabel(ax, 'Rows');
title(ax, 'Step03 observation quality flags');
save_figure_local(fig, fig_dir, 'Step03_Quality_Flags_20251222', save_figures);

fig = figure('Name', 'Step03 Observed vs Standard Angle', 'Color', 'w', 'Visible', visibility, 'Position', [180 180 1200 760]);
tiledlayout(fig, numel(cfg.sensor_ids), 1, 'Padding', 'compact', 'TileSpacing', 'compact');
for i = 1:numel(cfg.sensor_ids)
    sid = cfg.sensor_ids(i);
    rows = T(T.sensor_id == sid & T.is_valid, :);
    rows = select_evenly_for_plot_local(rows, Step03_Setting.plot_max_points_per_sensor);
    ax = nexttile;
    scatter(ax, rows.standard_angle_deg, rows.observed_angle_deg, 5, rows.blade_id, 'filled');
    grid(ax, 'on');
    xlabel(ax, 'Standard angle (deg)');
    ylabel(ax, 'Observed angle (deg)');
    title(ax, sprintf('CH%d observed vs standard angle', sid));
end
save_figure_local(fig, fig_dir, 'Step03_Observed_vs_Standard_Angle_20251222', save_figures);

fig = figure('Name', 'Step03 Angle Source Check', 'Color', 'w', 'Visible', visibility, 'Position', [200 200 1200 760]);
tiledlayout(fig, numel(cfg.sensor_ids), 1, 'Padding', 'compact', 'TileSpacing', 'compact');
for i = 1:numel(cfg.sensor_ids)
    sid = cfg.sensor_ids(i);
    rows = T(T.sensor_id == sid & isfinite(T.fallback_observed_angle_deg), :);
    rows = select_evenly_for_plot_local(rows, Step03_Setting.plot_max_points_per_sensor);
    ax = nexttile;
    delta = wrap_to_180_local(rows.observed_angle_deg - rows.fallback_observed_angle_deg);
    plot(ax, rows.arrival_time_s, delta, '.', 'MarkerSize', 4);
    grid(ax, 'on');
    xlabel(ax, 'Arrival time (s)');
    ylabel(ax, sprintf('CH%d rel - fallback (deg)', sid));
    title(ax, sprintf('CH%d relative angle source consistency', sid));
end
save_figure_local(fig, fig_dir, 'Step03_AngleSource_Check_20251222', save_figures);

fig = figure('Name', 'Step03 Rev ID Mismatch', 'Color', 'w', 'Visible', visibility, 'Position', [220 220 1200 760]);
tiledlayout(fig, numel(cfg.sensor_ids), 1, 'Padding', 'compact', 'TileSpacing', 'compact');
for i = 1:numel(cfg.sensor_ids)
    sid = cfg.sensor_ids(i);
    rows = T(T.sensor_id == sid, :);
    rows = select_evenly_for_plot_local(rows, Step03_Setting.plot_max_points_per_sensor);
    ax = nexttile;
    plot(ax, rows.arrival_time_s, double(rows.rev_id_mismatch), '.', 'MarkerSize', 4);
    grid(ax, 'on');
    ylim(ax, [-0.1 1.1]);
    xlabel(ax, 'Arrival time (s)');
    ylabel(ax, sprintf('CH%d mismatch', sid));
    title(ax, sprintf('CH%d Step02 rev_id vs OPR-derived rev_id', sid));
end
save_figure_local(fig, fig_dir, 'Step03_RevID_Mismatch_20251222', save_figures);
end

function rows = select_evenly_for_plot_local(rows, max_points)
if height(rows) <= max_points
    return;
end
[~, order] = sort(rows.arrival_time_s);
rows = rows(order, :);
idx = unique(round(linspace(1, height(rows), max_points)));
rows = rows(idx, :);
end

function save_figure_local(fig, fig_dir, base_name, save_figures)
if ~save_figures
    return;
end
ensure_dir_local(fig_dir);
png_file = fullfile(fig_dir, [base_name, '.png']);
pdf_file = fullfile(fig_dir, [base_name, '.pdf']);
try
    exportgraphics(fig, png_file, 'Resolution', 300);
    exportgraphics(fig, pdf_file, 'ContentType', 'vector');
catch
    saveas(fig, png_file);
    saveas(fig, pdf_file);
end
end

function value = get_field_or_default_local(S, field_name, default_value)
if isstruct(S) && isfield(S, field_name) && ~isempty(S.(field_name))
    value = S.(field_name);
else
    value = default_value;
end
end

function x = wrap_to_180_local(x)
x = mod(x + 180, 360) - 180;
end

function r = safe_ratio_local(a, b)
if b == 0
    r = NaN;
else
    r = a / b;
end
end

function x = first_or_nan_local(v, mode_name)
v = v(isfinite(v));
if isempty(v)
    x = NaN;
elseif strcmp(mode_name, 'min')
    x = min(v);
elseif strcmp(mode_name, 'max')
    x = max(v);
else
    x = v(1);
end
end

function x = first_finite_local(v)
v = v(isfinite(v));
if isempty(v)
    x = NaN;
else
    x = v(1);
end
end

function x = mean_or_nan_local(v)
v = v(isfinite(v));
if isempty(v)
    x = NaN;
else
    x = mean(v);
end
end

function x = median_or_nan_local(v)
v = v(isfinite(v));
if isempty(v)
    x = NaN;
else
    x = median(v);
end
end

function x = std_or_nan_local(v)
v = v(isfinite(v));
if numel(v) < 2
    x = NaN;
else
    x = std(v);
end
end

function x = rms_or_nan_local(v)
v = v(isfinite(v));
if isempty(v)
    x = NaN;
else
    x = sqrt(mean(v.^2));
end
end
