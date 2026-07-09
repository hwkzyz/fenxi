function Analyze_TemplateCalibration_Impact_20250527(mode)
%ANALYZE_TEMPLATECALIBRATION_IMPACT_20250527
% Compare how two low-speed waveform calibrations affect Step05 results.
%
% Usage:
%   Analyze_TemplateCalibration_Impact_20250527
%   Analyze_TemplateCalibration_Impact_20250527('prepare')
%   Analyze_TemplateCalibration_Impact_20250527('summarize')
%   Analyze_TemplateCalibration_Impact_20250527('rerun')
%
% The script keeps the comparison fair by:
%   1. adapting the rotating-calibration low-speed templates to the
%      foundation Template.SensorBlade format;
%   2. rerunning the same foundation Step05 logic with each template;
%   3. comparing window-paired identification and strain-validation outputs.

if nargin < 1 || isempty(mode)
    mode = 'all';
end
mode = lower(strtrim(char(mode)));

clc;
fprintf('\n=== Template calibration impact analysis: %s ===\n', mode);

P = make_analysis_settings_local();
ensure_dir_local(P.out_dir);
ensure_dir_local(P.fig_dir);

switch mode
    case 'prepare'
        prepare_template_inputs_local(P);
    case 'summarize'
        summarize_available_results_local(P);
    case 'rerun'
        prepare_template_inputs_local(P);
        run_reidentification_children_local(P);
        summarize_available_results_local(P);
    case 'all'
        prepare_template_inputs_local(P);
        summarize_available_results_local(P);
    otherwise
        error('Unknown mode: %s. Use prepare, summarize, rerun, or all.', mode);
end
end


function P = make_analysis_settings_local()
script_dir = fileparts(mfilename('fullpath'));
cfg = BTTProjectConfig_20250527();

P = struct();
P.cfg = cfg;
P.script_dir = script_dir;
P.case_name = '20250526_2500-3500_t400';
P.dataset = cfg.dataset;
P.blade_id = 1;
P.sensors = [1 3 6];
P.sensor_tag = ['S', sprintf('%d', P.sensors)];

P.foundation_tag = 'impF';
P.calibration_tag = 'impC';

P.out_dir = fullfile(cfg.output_root, 'template_calibration_impact', P.case_name);
P.fig_dir = fullfile(cfg.figure_root, 'template_calibration_impact', P.case_name);

P.foundation_template_file = fullfile(cfg.step04_output_dir, ...
    sprintf('Template_OPRCenterStd_LowSpeed_AllBlades_%s_%s.mat', ...
    P.sensor_tag, P.dataset));

P.calibration_root = fullfile(script_dir, '..', ...
    '20250527_low_speed_rotating_calibration');
P.calibration_template_library_dir = fullfile(P.calibration_root, ...
    'output', 'template_library');

P.adapted_template_file = fullfile(P.out_dir, sprintf( ...
    'Template_CalibrationAdapted_ForFoundationStep05_B%d_%s_%s.mat', ...
    P.blade_id, P.sensor_tag, P.dataset));

P.template_comparison_csv = fullfile(P.out_dir, sprintf( ...
    'TemplateCalibration_TemplateCurveComparison_B%d_%s_%s.csv', ...
    P.blade_id, P.sensor_tag, P.dataset));
P.coordinate_check_csv = fullfile(P.out_dir, sprintf( ...
    'TemplateCalibration_CoordinateGaugeCheck_B%d_%s_%s.csv', ...
    P.blade_id, P.sensor_tag, P.dataset));
P.window_delta_csv = fullfile(P.out_dir, sprintf( ...
    'TemplateCalibration_WindowDelta_B%d_%s_%s.csv', ...
    P.blade_id, P.sensor_tag, P.dataset));
P.step06b_delta_csv = fullfile(P.out_dir, sprintf( ...
    'TemplateCalibration_StrainValidationDelta_B%d_%s_%s.csv', ...
    P.blade_id, P.sensor_tag, P.dataset));
P.summary_csv = fullfile(P.out_dir, sprintf( ...
    'TemplateCalibration_ImpactSummary_B%d_%s_%s.csv', ...
    P.blade_id, P.sensor_tag, P.dataset));
P.commands_txt = fullfile(P.out_dir, sprintf( ...
    'TemplateCalibration_RerunCommands_B%d_%s_%s.txt', ...
    P.blade_id, P.sensor_tag, P.dataset));

P.step05_foundation_root = fullfile(P.cfg.output_root, ...
    ['step05_single_sync_direct_template_', P.foundation_tag]);
P.step05_calibration_root = fullfile(P.cfg.output_root, ...
    ['step05_single_sync_direct_template_', P.calibration_tag]);
P.step06b_foundation_root = fullfile(P.cfg.output_root, ...
    ['step06b_step05_strain_tarc_validation_', P.foundation_tag]);
P.step06b_calibration_root = fullfile(P.cfg.output_root, ...
    ['step06b_step05_strain_tarc_validation_', P.calibration_tag]);
end


function prepare_template_inputs_local(P)
fprintf('\n[1] Preparing adapted template and template-level comparison.\n');

if ~isfile(P.foundation_template_file)
    error('Foundation template not found: %s', P.foundation_template_file);
end

loaded = load(P.foundation_template_file, 'Template');
FoundationTemplate = loaded.Template;

[CalibrationBySensor, calibration_files] = load_calibration_sensor_templates_local(P);
AdaptedTemplate = adapt_calibration_template_local( ...
    FoundationTemplate, CalibrationBySensor, calibration_files, P);

Template = AdaptedTemplate; %#ok<NASGU>
metadata = AdaptedTemplate.Metadata; %#ok<NASGU>
save(P.adapted_template_file, 'Template', 'metadata', '-v7.3');
fprintf('Adapted calibration template saved:\n  %s\n', P.adapted_template_file);

TemplateComparison = compare_template_curves_local( ...
    FoundationTemplate, AdaptedTemplate, P);
writetable(TemplateComparison, P.template_comparison_csv);
fprintf('Template comparison saved:\n  %s\n', P.template_comparison_csv);

CoordinateGaugeCheck = build_coordinate_gauge_check_local( ...
    FoundationTemplate, AdaptedTemplate, P);
writetable(CoordinateGaugeCheck, P.coordinate_check_csv);
fprintf('Coordinate/gauge check saved:\n  %s\n', P.coordinate_check_csv);

write_rerun_commands_local(P);
plot_template_comparison_local(FoundationTemplate, AdaptedTemplate, ...
    TemplateComparison, P);
end


function [CalibrationBySensor, files] = load_calibration_sensor_templates_local(P)
CalibrationBySensor = struct();
files = strings(numel(P.sensors), 1);

for i = 1:numel(P.sensors)
    sid = P.sensors(i);
    file = fullfile(P.calibration_template_library_dir, ...
        sprintf('Template_%s_B%d_S%d.mat', P.dataset, P.blade_id, sid));

    if ~isfile(file)
        error('Calibration sensor template not found: %s', file);
    end

    D = load(file, 'Template');
    if ~isfield(D, 'Template') || ~isfield(D.Template, 'Sensor')
        error('Calibration template must contain Template.Sensor: %s', file);
    end

    S = D.Template.Sensor;
    if numel(S) ~= 1
        idx = find([S.sensor_id] == sid, 1, 'first');
        if isempty(idx)
            error('Cannot find CH%d in calibration template: %s', sid, file);
        end
        S = S(idx);
    end

    key = sensor_key_local(sid);
    CalibrationBySensor.(key) = S;
    files(i) = string(file);
end
end


function Adapted = adapt_calibration_template_local(Foundation, CalBySensor, files, P)
Adapted = Foundation;
Adapted.CreatedBy = mfilename;
Adapted.CreatedOn = datestr(now, 31);
Adapted.CalibrationImpact_Source = '20250527_low_speed_rotating_calibration';
Adapted.CalibrationImpact_SourceFiles = files;
Adapted.CalibrationImpact_Adapter = ...
    'Template.Sensor single-blade records converted to foundation SensorBlade records';

for i = 1:numel(P.sensors)
    sid = P.sensors(i);
    key = sensor_key_local(sid);
    cal = CalBySensor.(key);
    idx = find_sensorblade_index_local(Adapted, sid, P.blade_id);
    if isempty(idx)
        error('Foundation template lacks CH%d B%d entry.', sid, P.blade_id);
    end

    entry = Adapted.SensorBlade(idx);
    x_grid = cal.x_grid(:);
    v_grid = cal.v_grid(:);
    dv_dx = cal.dv_dx(:);

    if isfield(cal, 'bin_weight') && ~isempty(cal.bin_weight)
        weight_grid = cal.bin_weight(:);
    else
        weight_grid = ones(size(x_grid));
    end

    if isfield(cal, 'x_domain') && numel(cal.x_domain) >= 2
        x_domain = double(cal.x_domain(:)).';
    else
        x_domain = [min(x_grid, [], 'omitnan'), max(x_grid, [], 'omitnan')];
    end

    domain_mask = x_grid >= x_domain(1) & x_grid <= x_domain(2) & ...
        isfinite(x_grid) & isfinite(v_grid);
    valid_grid_mask = isfinite(x_grid) & isfinite(v_grid) & isfinite(dv_dx);

    entry.x_grid = x_grid;
    entry.v_grid = v_grid;
    entry.v_grid_raw = v_grid;
    entry.v_grid_baseline_removed = v_grid - get_numeric_field_local(cal, 'baseline', 0);
    entry.dv_dx = dv_dx;
    entry.weight_grid = weight_grid;
    entry.count_grid = max(1, round(100 * weight_grid));
    entry.valid_grid_mask = valid_grid_mask;
    entry.domain_effective_mask = domain_mask;
    entry.domain_mask = domain_mask;
    entry.x_domain = x_domain;
    entry.baseline = get_numeric_field_local(cal, 'baseline', NaN);
    entry.amplitude = max(v_grid, [], 'omitnan') - min(v_grid, [], 'omitnan');
    entry.pulse_count = get_numeric_field_local(cal, 'lap_count', NaN);
    entry.point_count = get_numeric_field_local(cal, 'point_count', numel(x_grid));
    entry.quality_status = 'calibration_adapted';
    entry.x_points_preview = [];
    entry.v_points_preview = [];

    Adapted.SensorBlade(idx) = entry;
end

Adapted.Summary_Table = build_summary_table_from_template_local(Adapted, P);
Adapted.Metadata = struct();
Adapted.Metadata.dataset = P.dataset;
Adapted.Metadata.created_by = mfilename;
Adapted.Metadata.created_on = Adapted.CreatedOn;
Adapted.Metadata.source = Adapted.CalibrationImpact_Source;
Adapted.Metadata.source_files = files;
Adapted.Metadata.blade_id = P.blade_id;
Adapted.Metadata.sensor_ids = P.sensors;
end


function T = compare_template_curves_local(F, C, P)
rows = {};

for i = 1:numel(P.sensors)
    sid = P.sensors(i);
    f = get_template_entry_local(F, sid, P.blade_id);
    c = get_template_entry_local(C, sid, P.blade_id);

    overlap = [max(f.x_domain(1), c.x_domain(1)), ...
        min(f.x_domain(2), c.x_domain(2))];
    has_overlap = all(isfinite(overlap)) && overlap(2) > overlap(1);

    n_grid = 1200;
    if has_overlap
        xq = linspace(overlap(1), overlap(2), n_grid).';
        vf = interp1(f.x_grid(:), f.v_grid(:), xq, 'linear');
        vc = interp1(c.x_grid(:), c.v_grid(:), xq, 'linear');
        gf = interp1(f.x_grid(:), f.dv_dx(:), xq, 'linear');
        gc = interp1(c.x_grid(:), c.dv_dx(:), xq, 'linear');
        keep = isfinite(vf) & isfinite(vc);
        rmse_v = sqrt(mean((vf(keep) - vc(keep)).^2, 'omitnan'));
        mae_v = mean(abs(vf(keep) - vc(keep)), 'omitnan');
        max_abs_diff_v = max(abs(vf(keep) - vc(keep)), [], 'omitnan');
        corr_v = corr_pair_local(vf(keep), vc(keep));
        grad_scale_median = median(abs(gc) ./ max(abs(gf), eps), 'omitnan');
    else
        rmse_v = NaN;
        mae_v = NaN;
        max_abs_diff_v = NaN;
        corr_v = NaN;
        grad_scale_median = NaN;
    end

    rows(end+1, :) = {sid, P.blade_id, ...
        f.x_domain(1), f.x_domain(2), diff(f.x_domain), ...
        c.x_domain(1), c.x_domain(2), diff(c.x_domain), ...
        safe_divide_local(diff(c.x_domain), diff(f.x_domain)), ...
        max(f.v_grid, [], 'omitnan') - min(f.v_grid, [], 'omitnan'), ...
        max(c.v_grid, [], 'omitnan') - min(c.v_grid, [], 'omitnan'), ...
        max(abs(f.dv_dx), [], 'omitnan'), ...
        max(abs(c.dv_dx), [], 'omitnan'), ...
        overlap(1), overlap(2), diff(overlap), ...
        rmse_v, mae_v, max_abs_diff_v, corr_v, grad_scale_median}; %#ok<AGROW>
end

T = cell2table(rows, 'VariableNames', { ...
    'sensor_id', 'blade_id', ...
    'foundation_domain_left_mm', 'foundation_domain_right_mm', ...
    'foundation_domain_width_mm', ...
    'calibration_domain_left_mm', 'calibration_domain_right_mm', ...
    'calibration_domain_width_mm', 'domain_width_ratio_cal_over_foundation', ...
    'foundation_p2p_v', 'calibration_p2p_v', ...
    'foundation_max_abs_dv_dx', 'calibration_max_abs_dv_dx', ...
    'common_domain_left_mm', 'common_domain_right_mm', 'common_domain_width_mm', ...
    'curve_rmse_v_common', 'curve_mae_v_common', ...
    'curve_max_abs_diff_v_common', 'curve_corr_common', ...
    'median_abs_gradient_ratio_cal_over_foundation'});
end


function T = build_coordinate_gauge_check_local(F, C, P)
rows = {};

for i = 1:numel(P.sensors)
    sid = P.sensors(i);
    f = get_template_entry_local(F, sid, P.blade_id);
    c = get_template_entry_local(C, sid, P.blade_id);

    theta_f = NaN;
    theta_c = NaN;
    if isfield(F, 'Standard_Relative_Angles')
        theta_f = F.Standard_Relative_Angles(sid, P.blade_id);
    end
    if isfield(C, 'Standard_Relative_Angles')
        theta_c = C.Standard_Relative_Angles(sid, P.blade_id);
    end

    rows(end+1, :) = {sid, P.blade_id, theta_f, theta_c, ...
        wrap_to_180_local(theta_c - theta_f), ...
        f.baseline, c.baseline, c.baseline - f.baseline, ...
        get_numeric_field_local(c, 'pulse_count', NaN), ...
        get_numeric_field_local(c, 'point_count', NaN), ...
        string(getfield_default_local(c, 'quality_status', 'unknown'))}; %#ok<AGROW>
end

T = cell2table(rows, 'VariableNames', { ...
    'sensor_id', 'blade_id', ...
    'foundation_theta_std_deg', 'adapted_theta_std_deg', ...
    'theta_std_delta_deg', ...
    'foundation_baseline_v', 'calibration_baseline_v', ...
    'baseline_delta_v_cal_minus_foundation', ...
    'calibration_lap_count_or_pulse_count', ...
    'calibration_point_count', 'adapted_quality_status'});
end


function write_rerun_commands_local(P)
lines = {};
lines{end+1} = 'Run from PowerShell in the workspace root:';
lines{end+1} = '';
lines{end+1} = sprintf("cd '%s'", P.script_dir);
lines{end+1} = '';
lines{end+1} = 'Foundation-template rerun:';
lines{end+1} = sprintf("$env:STEP05_OUTPUT_TAG='%s'; $env:STEP05_TEMPLATE_OVERRIDE_FILE=''; $env:STEP05_SHOW_PLOTS='0'; $env:STEP05_SAVE_FIGURES='0'; matlab -batch ""cd(''%s''); Step05_SingleSync_DirectTemplate_Identification_20250527""", P.foundation_tag, P.script_dir);
lines{end+1} = sprintf("$env:STEP06B_OUTPUT_TAG='%s'; $env:STEP06B_STEP05_ROOT='%s'; $env:STEP06B_FORCE_REBUILD='1'; $env:STEP06B_SHOW_PLOTS='0'; matlab -batch ""cd(''%s''); Step06B_Validate_Step05_With_Strain_TARC_20250527""", P.foundation_tag, P.step05_foundation_root, P.script_dir);
lines{end+1} = '';
lines{end+1} = 'Calibration-template rerun:';
lines{end+1} = sprintf("$env:STEP05_OUTPUT_TAG='%s'; $env:STEP05_TEMPLATE_OVERRIDE_FILE='%s'; $env:STEP05_SHOW_PLOTS='0'; $env:STEP05_SAVE_FIGURES='0'; matlab -batch ""cd(''%s''); Step05_SingleSync_DirectTemplate_Identification_20250527""", P.calibration_tag, P.adapted_template_file, P.script_dir);
lines{end+1} = sprintf("$env:STEP06B_OUTPUT_TAG='%s'; $env:STEP06B_STEP05_ROOT='%s'; $env:STEP06B_FORCE_REBUILD='1'; $env:STEP06B_SHOW_PLOTS='0'; matlab -batch ""cd(''%s''); Step06B_Validate_Step05_With_Strain_TARC_20250527""", P.calibration_tag, P.step05_calibration_root, P.script_dir);
lines{end+1} = '';
lines{end+1} = 'Summarize after reruns:';
lines{end+1} = sprintf("matlab -batch ""cd(''%s''); Analyze_TemplateCalibration_Impact_20250527(''summarize'')""", P.script_dir);

    fid = fopen(P.commands_txt, 'w');
    cleanup = onCleanup(@() fclose(fid));
    for i = 1:numel(lines)
        fprintf(fid, '%s\n', lines{i});
    end
    clear cleanup;
fprintf('Rerun command note saved:\n  %s\n', P.commands_txt);
end


function run_reidentification_children_local(P)
fprintf('\n[2] Running full Step05/Step06B child processes.\n');

run_child_step05_local(P, P.foundation_tag, '', P.step05_foundation_root);
run_child_step06b_local(P, P.foundation_tag, P.step05_foundation_root);
run_child_step05_local(P, P.calibration_tag, P.adapted_template_file, ...
    P.step05_calibration_root);
run_child_step06b_local(P, P.calibration_tag, P.step05_calibration_root);
end


function run_child_step05_local(P, tag, template_file, expected_root)
old_tag = getenv('STEP05_OUTPUT_TAG');
old_tpl = getenv('STEP05_TEMPLATE_OVERRIDE_FILE');
old_show = getenv('STEP05_SHOW_PLOTS');
old_save = getenv('STEP05_SAVE_FIGURES');
cleanup = onCleanup(@() restore_env4_local( ...
    'STEP05_OUTPUT_TAG', old_tag, ...
    'STEP05_TEMPLATE_OVERRIDE_FILE', old_tpl, ...
    'STEP05_SHOW_PLOTS', old_show, ...
    'STEP05_SAVE_FIGURES', old_save));

setenv('STEP05_OUTPUT_TAG', tag);
setenv('STEP05_TEMPLATE_OVERRIDE_FILE', template_file);
setenv('STEP05_SHOW_PLOTS', '0');
setenv('STEP05_SAVE_FIGURES', '0');

cmd = sprintf('matlab -batch "cd(''%s''); Step05_SingleSync_DirectTemplate_Identification_20250527"', ...
    P.script_dir);
[status, txt] = system(cmd);
log_file = fullfile(P.out_dir, sprintf('log_step05_%s.txt', tag));
write_text_local(log_file, txt);
if status ~= 0
    error('Step05 child process failed for %s. See %s', tag, log_file);
end
if ~isfolder(expected_root)
    error('Expected Step05 output root not created: %s', expected_root);
end
clear cleanup;
end


function run_child_step06b_local(P, tag, step05_root)
old_tag = getenv('STEP06B_OUTPUT_TAG');
old_root = getenv('STEP06B_STEP05_ROOT');
old_force = getenv('STEP06B_FORCE_REBUILD');
old_show = getenv('STEP06B_SHOW_PLOTS');
cleanup = onCleanup(@() restore_env4_local( ...
    'STEP06B_OUTPUT_TAG', old_tag, ...
    'STEP06B_STEP05_ROOT', old_root, ...
    'STEP06B_FORCE_REBUILD', old_force, ...
    'STEP06B_SHOW_PLOTS', old_show));

setenv('STEP06B_OUTPUT_TAG', tag);
setenv('STEP06B_STEP05_ROOT', step05_root);
setenv('STEP06B_FORCE_REBUILD', '1');
setenv('STEP06B_SHOW_PLOTS', '0');

cmd = sprintf('matlab -batch "cd(''%s''); Step06B_Validate_Step05_With_Strain_TARC_20250527"', ...
    P.script_dir);
[status, txt] = system(cmd);
log_file = fullfile(P.out_dir, sprintf('log_step06b_%s.txt', tag));
write_text_local(log_file, txt);
if status ~= 0
    error('Step06B child process failed for %s. See %s', tag, log_file);
end
clear cleanup;
end


function summarize_available_results_local(P)
fprintf('\n[3] Summarizing available Step05/Step06B results.\n');

summary_rows = {};

if has_step05_pair_local(P)
    WDelta = compare_step05_window_results_local(P);
    writetable(WDelta, P.window_delta_csv);
    fprintf('Window-level Step05 delta saved:\n  %s\n', P.window_delta_csv);

    Stats05 = build_step05_summary_stats_local(WDelta);
    summary_rows = [summary_rows; table_to_rows_local(Stats05)]; %#ok<AGROW>
    plot_identification_delta_local(WDelta, P);
else
    fprintf('Step05 paired rerun results are not complete yet; skipping Step05 delta.\n');
end

if has_step06b_pair_local(P)
    SDelta = compare_step06b_validation_results_local(P);
    writetable(SDelta, P.step06b_delta_csv);
    fprintf('Step06B strain-validation delta saved:\n  %s\n', P.step06b_delta_csv);

    Stats06 = build_step06b_summary_stats_local(SDelta);
    summary_rows = [summary_rows; table_to_rows_local(Stats06)]; %#ok<AGROW>
    plot_step06b_delta_local(SDelta, P);
else
    fprintf('Step06B paired validation results are not complete yet; skipping Step06B delta.\n');
end

if ~isempty(summary_rows)
    ImpactSummary = cell2table(summary_rows, 'VariableNames', { ...
        'block', 'metric', 'n', 'foundation_mean', 'calibration_mean', ...
        'delta_mean_cal_minus_foundation', 'delta_median', ...
        'delta_std', 'delta_abs_mean', 'paired_t', 'paired_p_approx', ...
        'corr_foundation_calibration'});
    writetable(ImpactSummary, P.summary_csv);
    fprintf('Impact summary saved:\n  %s\n', P.summary_csv);
else
    fprintf('No paired rerun outputs found. Use mode rerun or commands in:\n  %s\n', ...
        P.commands_txt);
end
end


function tf = has_step05_pair_local(P)
tf = isfile(step05_trend_file_local(P, P.step05_foundation_root)) && ...
    isfile(step05_trend_file_local(P, P.step05_calibration_root));
end


function tf = has_step06b_pair_local(P)
tf = isfile(step06b_response_file_local(P, P.step06b_foundation_root)) && ...
    isfile(step06b_response_file_local(P, P.step06b_calibration_root));
end


function file = step05_trend_file_local(P, root)
file = fullfile(root, P.case_name, sprintf( ...
    'Trend_Step05_SingleSyncDirectTemplate_B%d_%s_%s.csv', ...
    P.blade_id, P.sensor_tag, P.dataset));
end


function file = step06b_response_file_local(P, root)
file = fullfile(root, P.case_name, sprintf( ...
    'Step06B_ResponseAmp_Comparison_B%d_%s_%s.csv', ...
    P.blade_id, P.sensor_tag, P.dataset));
end


function D = compare_step05_window_results_local(P)
F = readtable(step05_trend_file_local(P, P.step05_foundation_root));
C = readtable(step05_trend_file_local(P, P.step05_calibration_root));

common = intersect(F.window_id, C.window_id, 'stable');
n = numel(common);
D = table();
D.window_id = common(:);

for i = 1:n
    wf = F(F.window_id == common(i), :);
    wc = C(C.window_id == common(i), :);
    if i == 1
        D.time_center_s = nan(n,1);
        D.rpm_mean = nan(n,1);
        D.A_foundation_mm = nan(n,1);
        D.A_calibration_mm = nan(n,1);
        D.fn_foundation_hz = nan(n,1);
        D.fn_calibration_hz = nan(n,1);
        D.EO_foundation = nan(n,1);
        D.EO_calibration = nan(n,1);
        D.rmse_foundation_mv = nan(n,1);
        D.rmse_calibration_mv = nan(n,1);
    end
    D.time_center_s(i) = table_value_any_local(wf, ...
        {'time_center_s', 'window_center_time_s'}, NaN);
    D.rpm_mean(i) = table_value_any_local(wf, ...
        {'rpm_mean', 'rot_rpm_mean'}, NaN);
    D.A_foundation_mm(i) = abs(table_value_any_local(wf, ...
        {'A_id_mm', 'A_id'}, NaN));
    D.A_calibration_mm(i) = abs(table_value_any_local(wc, ...
        {'A_id_mm', 'A_id'}, NaN));
    D.fn_foundation_hz(i) = table_value_any_local(wf, ...
        {'fn_id_hz', 'fn_id'}, NaN);
    D.fn_calibration_hz(i) = table_value_any_local(wc, ...
        {'fn_id_hz', 'fn_id'}, NaN);
    D.EO_foundation(i) = table_value_local(wf, 'EO_id', NaN);
    D.EO_calibration(i) = table_value_local(wc, 'EO_id', NaN);
    D.rmse_foundation_mv(i) = rmse_mv_from_trend_local(wf);
    D.rmse_calibration_mv(i) = rmse_mv_from_trend_local(wc);
end

D.A_delta_mm = D.A_calibration_mm - D.A_foundation_mm;
D.A_ratio_cal_over_foundation = safe_divide_local(D.A_calibration_mm, D.A_foundation_mm);
D.A_delta_percent = 100 * safe_divide_local(D.A_delta_mm, D.A_foundation_mm);
D.fn_delta_hz = D.fn_calibration_hz - D.fn_foundation_hz;
D.EO_changed = D.EO_calibration ~= D.EO_foundation;
D.rmse_delta_mv = D.rmse_calibration_mv - D.rmse_foundation_mv;
D.rmse_ratio_cal_over_foundation = safe_divide_local(D.rmse_calibration_mv, D.rmse_foundation_mv);
D.A_bland_altman_mean_mm = 0.5 * (D.A_calibration_mm + D.A_foundation_mm);
D.A_bland_altman_diff_mm = D.A_delta_mm;
end


function D = compare_step06b_validation_results_local(P)
F = readtable(step06b_response_file_local(P, P.step06b_foundation_root));
C = readtable(step06b_response_file_local(P, P.step06b_calibration_root));

common = intersect(F.window_id, C.window_id, 'stable');
n = numel(common);
D = table();
D.window_id = common(:);
D.time_center_s = nan(n,1);
D.Step05_A_foundation_mm = nan(n,1);
D.Step05_A_calibration_mm = nan(n,1);
D.Strain_EqDisp_foundation_mm = nan(n,1);
D.Strain_EqDisp_calibration_mm = nan(n,1);
D.EqDisp_AbsError_foundation_mm = nan(n,1);
D.EqDisp_AbsError_calibration_mm = nan(n,1);
D.phase_opt_TARC_foundation = nan(n,1);
D.phase_opt_TARC_calibration = nan(n,1);

for i = 1:n
    wf = F(F.window_id == common(i), :);
    wc = C(C.window_id == common(i), :);
    D.time_center_s(i) = table_value_local(wf, 'time_center_s', NaN);
    D.Step05_A_foundation_mm(i) = table_value_local(wf, 'Step05_A_mm', NaN);
    D.Step05_A_calibration_mm(i) = table_value_local(wc, 'Step05_A_mm', NaN);
    D.Strain_EqDisp_foundation_mm(i) = table_value_local(wf, 'Strain_EqDispAmp_mm', NaN);
    D.Strain_EqDisp_calibration_mm(i) = table_value_local(wc, 'Strain_EqDispAmp_mm', NaN);
    D.EqDisp_AbsError_foundation_mm(i) = table_value_local(wf, 'EqDisp_AbsError_mm', NaN);
    D.EqDisp_AbsError_calibration_mm(i) = table_value_local(wc, 'EqDisp_AbsError_mm', NaN);
    D.phase_opt_TARC_foundation(i) = table_value_local(wf, 'phase_opt_TARC', NaN);
    D.phase_opt_TARC_calibration(i) = table_value_local(wc, 'phase_opt_TARC', NaN);
end

D.Step05_A_delta_mm = D.Step05_A_calibration_mm - D.Step05_A_foundation_mm;
D.Step05_A_ratio_cal_over_foundation = safe_divide_local( ...
    D.Step05_A_calibration_mm, D.Step05_A_foundation_mm);
D.EqDisp_AbsError_delta_mm = ...
    D.EqDisp_AbsError_calibration_mm - D.EqDisp_AbsError_foundation_mm;
D.phase_opt_TARC_delta = ...
    D.phase_opt_TARC_calibration - D.phase_opt_TARC_foundation;
end


function Stats = build_step05_summary_stats_local(D)
Stats = [
    paired_metric_stats_local('Step05', 'A_mm', ...
    D.A_foundation_mm, D.A_calibration_mm)
    paired_metric_stats_local('Step05', 'fn_hz', ...
    D.fn_foundation_hz, D.fn_calibration_hz)
    paired_metric_stats_local('Step05', 'weighted_rmse_mv', ...
    D.rmse_foundation_mv, D.rmse_calibration_mv)];
end


function Stats = build_step06b_summary_stats_local(D)
Stats = [
    paired_metric_stats_local('Step06B', 'Step05_A_mm', ...
    D.Step05_A_foundation_mm, D.Step05_A_calibration_mm)
    paired_metric_stats_local('Step06B', 'Strain_EqDispAmp_mm', ...
    D.Strain_EqDisp_foundation_mm, D.Strain_EqDisp_calibration_mm)
    paired_metric_stats_local('Step06B', 'EqDisp_AbsError_mm', ...
    D.EqDisp_AbsError_foundation_mm, D.EqDisp_AbsError_calibration_mm)
    paired_metric_stats_local('Step06B', 'phase_opt_TARC', ...
    D.phase_opt_TARC_foundation, D.phase_opt_TARC_calibration)];
end


function S = paired_metric_stats_local(block, metric, foundation, calibration)
foundation = foundation(:);
calibration = calibration(:);
keep = isfinite(foundation) & isfinite(calibration);
f = foundation(keep);
c = calibration(keep);
d = c - f;
n = numel(d);

S = struct();
S.block = string(block);
S.metric = string(metric);
S.n = n;
S.foundation_mean = mean(f, 'omitnan');
S.calibration_mean = mean(c, 'omitnan');
S.delta_mean_cal_minus_foundation = mean(d, 'omitnan');
S.delta_median = median(d, 'omitnan');
S.delta_std = std(d, 'omitnan');
S.delta_abs_mean = mean(abs(d), 'omitnan');
S.paired_t = NaN;
S.paired_p_approx = NaN;
S.corr_foundation_calibration = corr_pair_local(f, c);

if n >= 2 && S.delta_std > 0
    S.paired_t = S.delta_mean_cal_minus_foundation / (S.delta_std / sqrt(n));
    S.paired_p_approx = 2 * normal_tail_approx_local(abs(S.paired_t));
end
end


function rows = table_to_rows_local(S)
rows = cell(numel(S), 12);
for i = 1:numel(S)
    rows(i, :) = {char(S(i).block), char(S(i).metric), S(i).n, ...
        S(i).foundation_mean, S(i).calibration_mean, ...
        S(i).delta_mean_cal_minus_foundation, S(i).delta_median, ...
        S(i).delta_std, S(i).delta_abs_mean, S(i).paired_t, ...
        S(i).paired_p_approx, S(i).corr_foundation_calibration};
end
end


function plot_template_comparison_local(F, C, T, P)
fig = figure('Color', 'w', 'Name', 'Template calibration comparison', ...
    'Visible', 'off');
tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

colors = lines(numel(P.sensors));

ax = nexttile(tl, 1);
hold(ax, 'on'); grid(ax, 'on');
for i = 1:numel(P.sensors)
    sid = P.sensors(i);
    f = get_template_entry_local(F, sid, P.blade_id);
    c = get_template_entry_local(C, sid, P.blade_id);
    plot(ax, f.x_grid, f.v_grid, '-', 'Color', colors(i,:), 'LineWidth', 1.4);
    plot(ax, c.x_grid, c.v_grid, '--', 'Color', colors(i,:), 'LineWidth', 1.4);
end
xlabel(ax, 'x (mm)'); ylabel(ax, 'Voltage (V)');
title(ax, 'Template curves');

ax = nexttile(tl, 2);
bar(ax, categorical("CH" + string(T.sensor_id)), ...
    [T.foundation_domain_width_mm, T.calibration_domain_width_mm]);
grid(ax, 'on'); ylabel(ax, 'Domain width (mm)');
legend(ax, {'Foundation', 'Calibration'}, 'Location', 'best');
title(ax, 'Trusted-domain width');

ax = nexttile(tl, 3);
bar(ax, categorical("CH" + string(T.sensor_id)), ...
    [T.foundation_p2p_v, T.calibration_p2p_v]);
grid(ax, 'on'); ylabel(ax, 'Peak-to-peak voltage (V)');
legend(ax, {'Foundation', 'Calibration'}, 'Location', 'best');
title(ax, 'Template voltage span');

ax = nexttile(tl, 4);
bar(ax, categorical("CH" + string(T.sensor_id)), T.curve_rmse_v_common);
grid(ax, 'on'); ylabel(ax, 'RMSE on common domain (V)');
title(ax, 'Curve difference');

saveas(fig, fullfile(P.fig_dir, sprintf( ...
    'TemplateCalibration_Fig01_TemplateComparison_B%d_%s.png', ...
    P.blade_id, P.sensor_tag)));
close(fig);
end


function plot_identification_delta_local(D, P)
fig = figure('Color', 'w', 'Name', 'Step05 template impact', ...
    'Visible', 'off');
tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile(tl, 1);
plot(ax, D.window_id, D.A_foundation_mm, 'o-', 'LineWidth', 1.2); hold(ax, 'on');
plot(ax, D.window_id, D.A_calibration_mm, 's-', 'LineWidth', 1.2);
grid(ax, 'on'); xlabel(ax, 'Window'); ylabel(ax, 'A (mm)');
legend(ax, {'Foundation', 'Calibration'}, 'Location', 'best');
title(ax, 'Identified amplitude');

ax = nexttile(tl, 2);
bar(ax, D.window_id, D.A_delta_mm);
grid(ax, 'on'); xlabel(ax, 'Window'); ylabel(ax, '\DeltaA (mm)');
title(ax, 'Calibration minus foundation');

ax = nexttile(tl, 3);
scatter(ax, D.A_bland_altman_mean_mm, D.A_bland_altman_diff_mm, 35, 'filled');
hold(ax, 'on'); grid(ax, 'on');
m = mean(D.A_bland_altman_diff_mm, 'omitnan');
s = std(D.A_bland_altman_diff_mm, 'omitnan');
yline(ax, m, '-');
yline(ax, m + 1.96*s, '--');
yline(ax, m - 1.96*s, '--');
xlabel(ax, 'Mean A (mm)'); ylabel(ax, 'A diff (mm)');
title(ax, 'Bland-Altman: A');

ax = nexttile(tl, 4);
plot(ax, D.window_id, D.rmse_foundation_mv, 'o-', 'LineWidth', 1.2); hold(ax, 'on');
plot(ax, D.window_id, D.rmse_calibration_mv, 's-', 'LineWidth', 1.2);
grid(ax, 'on'); xlabel(ax, 'Window'); ylabel(ax, 'Weighted RMSE (mV)');
legend(ax, {'Foundation', 'Calibration'}, 'Location', 'best');
title(ax, 'Waveform fit residual');

saveas(fig, fullfile(P.fig_dir, sprintf( ...
    'TemplateCalibration_Fig02_Step05Impact_B%d_%s.png', ...
    P.blade_id, P.sensor_tag)));
close(fig);
end


function plot_step06b_delta_local(D, P)
fig = figure('Color', 'w', 'Name', 'Step06B validation impact', ...
    'Visible', 'off');
tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile(tl, 1);
plot(ax, D.window_id, D.Step05_A_foundation_mm, 'o-', 'LineWidth', 1.2); hold(ax, 'on');
plot(ax, D.window_id, D.Step05_A_calibration_mm, 's-', 'LineWidth', 1.2);
plot(ax, D.window_id, D.Strain_EqDisp_foundation_mm, 'k--', 'LineWidth', 1.2);
grid(ax, 'on'); xlabel(ax, 'Window'); ylabel(ax, 'Amplitude (mm)');
legend(ax, {'Foundation Step05', 'Calibration Step05', 'Strain eq.'}, ...
    'Location', 'best');
title(ax, 'Amplitude against strain equivalent displacement');

ax = nexttile(tl, 2);
bar(ax, D.window_id, D.EqDisp_AbsError_delta_mm);
grid(ax, 'on'); xlabel(ax, 'Window'); ylabel(ax, '\Delta abs error (mm)');
title(ax, 'Eq-disp validation error change');

ax = nexttile(tl, 3);
plot(ax, D.window_id, D.phase_opt_TARC_foundation, 'o-', 'LineWidth', 1.2); hold(ax, 'on');
plot(ax, D.window_id, D.phase_opt_TARC_calibration, 's-', 'LineWidth', 1.2);
grid(ax, 'on'); xlabel(ax, 'Window'); ylabel(ax, 'Phase-opt TARC');
legend(ax, {'Foundation', 'Calibration'}, 'Location', 'best');
title(ax, 'TARC consistency');

ax = nexttile(tl, 4);
scatter(ax, D.Step05_A_foundation_mm, D.Step05_A_calibration_mm, 35, 'filled');
hold(ax, 'on'); grid(ax, 'on');
lims = [min([xlim(ax), ylim(ax)]), max([xlim(ax), ylim(ax)])];
plot(ax, lims, lims, 'k--'); xlim(ax, lims); ylim(ax, lims);
xlabel(ax, 'Foundation A (mm)'); ylabel(ax, 'Calibration A (mm)');
title(ax, 'Paired amplitude');

saveas(fig, fullfile(P.fig_dir, sprintf( ...
    'TemplateCalibration_Fig03_StrainValidationImpact_B%d_%s.png', ...
    P.blade_id, P.sensor_tag)));
close(fig);
end


function T = build_summary_table_from_template_local(Template, P)
rows = {};
for i = 1:numel(P.sensors)
    sid = P.sensors(i);
    e = get_template_entry_local(Template, sid, P.blade_id);
    rows(end+1, :) = {sid, P.blade_id, e.pulse_count, e.point_count, ...
        e.x_domain(1), e.x_domain(2), diff(e.x_domain), ...
        e.baseline, e.amplitude, max(abs(e.dv_dx), [], 'omitnan'), ...
        string(e.quality_status)}; %#ok<AGROW>
end
T = cell2table(rows, 'VariableNames', { ...
    'sensor_id', 'blade_id', 'pulse_count', 'point_count', ...
    'domain_left_mm', 'domain_right_mm', 'domain_width_mm', ...
    'baseline_v', 'amplitude_v', 'max_abs_gradient_v_per_mm', ...
    'quality_status'});
end


function idx = find_sensorblade_index_local(Template, sid, blade_id)
idx = [];
for i = 1:numel(Template.SensorBlade)
    if Template.SensorBlade(i).sensor_id == sid && ...
            Template.SensorBlade(i).blade_id == blade_id
        idx = i;
        return;
    end
end
end


function entry = get_template_entry_local(Template, sid, blade_id)
idx = find_sensorblade_index_local(Template, sid, blade_id);
if isempty(idx)
    error('Template entry not found: CH%d B%d.', sid, blade_id);
end
entry = Template.SensorBlade(idx);
end


function v = table_value_local(T, name, default_value)
if ismember(name, T.Properties.VariableNames)
    x = T.(name);
    if isempty(x)
        v = default_value;
    else
        v = x(1);
    end
else
    v = default_value;
end
end


function v = table_value_any_local(T, names, default_value)
v = default_value;
for i = 1:numel(names)
    if ismember(names{i}, T.Properties.VariableNames)
        v = table_value_local(T, names{i}, default_value);
        return;
    end
end
end


function v = rmse_mv_from_trend_local(T)
if ismember('weighted_voltage_rmse_mv', T.Properties.VariableNames)
    v = table_value_local(T, 'weighted_voltage_rmse_mv', NaN);
elseif ismember('weighted_voltage_rmse', T.Properties.VariableNames)
    v = 1000 * table_value_local(T, 'weighted_voltage_rmse', NaN);
else
    v = NaN;
end
end


function c = corr_pair_local(a, b)
a = a(:);
b = b(:);
keep = isfinite(a) & isfinite(b);
if nnz(keep) < 2
    c = NaN;
    return;
end
R = corrcoef(a(keep), b(keep));
c = R(1,2);
end


function y = safe_divide_local(a, b)
y = a ./ b;
y(~isfinite(y)) = NaN;
end


function key = sensor_key_local(sid)
key = sprintf('CH%d', sid);
end


function ensure_dir_local(d)
if exist(d, 'dir') ~= 7
    mkdir(d);
end
end


function value = get_numeric_field_local(S, name, default_value)
if isfield(S, name) && ~isempty(S.(name))
    value = double(S.(name));
    if numel(value) > 1
        value = value(1);
    end
else
    value = default_value;
end
end


function value = getfield_default_local(S, name, default_value)
if isfield(S, name) && ~isempty(S.(name))
    value = S.(name);
else
    value = default_value;
end
end


function p = normal_tail_approx_local(z)
p = 0.5 * erfc(z ./ sqrt(2));
end


function ang = wrap_to_180_local(ang)
ang = mod(ang + 180, 360) - 180;
end


function write_text_local(file, txt)
fid = fopen(file, 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s', txt);
clear cleanup;
end


function restore_env4_local(k1, v1, k2, v2, k3, v3, k4, v4)
setenv(k1, v1);
setenv(k2, v2);
setenv(k3, v3);
setenv(k4, v4);
end
