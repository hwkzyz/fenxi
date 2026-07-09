function Summary = Diag_Step05_RouteAblation_20251222()
%DIAG_STEP05_ROUTEABLATION_20251222
% One-factor-at-a-time diagnosis for the Step05 eta-zero route change.
%
% The question is why eta=0 does not recover the official NoEta result.
% This harness creates temporary worker scripts under output/diagnostics,
% changes only a small set of route constants, runs the first two windows,
% and summarizes EO/RMSE/boundary behavior.

clc; close all;

this_dir = fileparts(mfilename('fullpath'));
cfg = BTTProjectConfig_20251222();
case_name = cfg.dynamic_cases{1};
out_dir = fullfile(cfg.output_root, 'diagnostics', ...
    'step05_route_ablation_20251222');
worker_dir = fullfile(out_dir, 'workers');
if exist(worker_dir, 'dir') ~= 7
    mkdir(worker_dir);
end
zero_eta_file = fullfile(out_dir, 'RouteAblation_ZeroEta_B1_S123.mat');
write_zero_eta_file_local(zero_eta_file, cfg.sensor_ids(:).');

official_file = fullfile(this_dir, ...
    'Step05_SingleSync_DirectTemplate_Identification_20251222.m');
legacy_file = fullfile(this_dir, 'legacy_opr_center_flow_20260707', ...
    'Step05_SingleSync_DirectTemplate_Identification_20251222.m');

experiments = build_experiment_set_local();
rows = repmat(empty_result_row_local(), numel(experiments), 1);

old_env = capture_env_local();
cleanup_obj = onCleanup(@() restore_env_local(old_env)); %#ok<NASGU>

old_dir = pwd;
dir_cleanup = onCleanup(@() cd(old_dir)); %#ok<NASGU>
cd(this_dir);
addpath(this_dir);
addpath(fileparts(this_dir));

for i = 1:numel(experiments)
    ex = experiments(i);
    fprintf('\n=== Route ablation %d/%d: %s ===\n', ...
        i, numel(experiments), ex.tag);
    worker_file = fullfile(worker_dir, ...
        sprintf('Worker_%02d_%s.m', i, ex.tag));

    if strcmp(ex.source, "official")
        source_file = official_file;
    else
        source_file = legacy_file;
    end

    make_worker_local(source_file, worker_file, ex);

    setenv('STEP05_OUTPUT_TAG', ['abl_', char(ex.tag)]);
    setenv('STEP05_DEBUG_MAX_WINDOWS', '2');
    setenv('STEP05_TARGET_BLADES', '1');
    setenv('STEP05_SHOW_PLOTS', '0');
    setenv('STEP05_SAVE_FIGURES', '0');
    setenv('STEP05_ANALYSIS_START_TIME', char(ex.start_env));
    setenv('STEP05_TARGET_LAPS', char(ex.target_laps_env));
    if ex.static_eta_enable
        setenv('STEP05_STATIC_ETA_FILE', zero_eta_file);
    else
        setenv('STEP05_STATIC_ETA_FILE', '');
    end

    try
        run(worker_file);
        result_file = resolve_result_file_local(cfg, case_name, ex);
        rows(i) = summarize_result_local(ex, worker_file, result_file);
    catch ME
        rows(i) = summarize_failure_local(ex, worker_file, ME);
    end
end

ResultTable = struct2table(rows);
summary_csv = fullfile(out_dir, ...
    'Step05_RouteAblation_B1_S123_20251222.csv');
writetable(ResultTable, summary_csv);

Summary = struct();
Summary.OutputDir = out_dir;
Summary.ResultTable = ResultTable;
Summary.SummaryCSV = summary_csv;

summary_mat = fullfile(out_dir, ...
    'Step05_RouteAblation_B1_S123_20251222.mat');
save(summary_mat, 'Summary', '-v7.3');

fprintf('\n=== Step05 route ablation summary ===\n');
disp(ResultTable(:, {'tag', 'source', 'start_time_s', 'core_mask_mode', ...
    'main_pulse_mode', 'phase_safe_expansion', 'active_eta_model', 'eo_sequence', ...
    'median_rmse_v', 'median_A_mm', 'median_dx_c_mm', ...
    'A_upper_hit_fraction', 'dx_lower_hit_fraction', 'status'}));
fprintf('Saved summary CSV:\n  %s\n', summary_csv);
end


function experiments = build_experiment_set_local()
base = struct( ...
    'tag', "", ...
    'source', "", ...
    'start_env', "50.0", ...
    'target_laps_env', "20", ...
    'start_time_s', 50.0, ...
    'target_laps', 20, ...
    'core_mask_mode', "", ...
    'main_pulse_mode', "", ...
    'phase_safe_expansion', false, ...
    'active_eta_model', "", ...
    'static_eta_enable', false, ...
    'description', "");
experiments = repmat(base, 0, 1);

experiments(end+1) = make_experiment_local( ...
    "official_noeta", "official", "50.0", "20", ...
    "sg_threshold_mainpulse_trust", "legacy_threshold", false, "none", false, ...
    "Official route reference.");

experiments(end+1) = make_experiment_local( ...
    "official_eta0", "official", "50.0", "20", ...
    "sg_threshold_mainpulse_trust", "legacy_threshold", false, "fixed_joint_static_eta", true, ...
    "Official settings but force fixed eta model with zero eta prior.");

experiments(end+1) = make_experiment_local( ...
    "official_start502", "official", "50.2", "10", ...
    "sg_threshold_mainpulse_trust", "legacy_threshold", false, "none", false, ...
    "Only move official NoEta to the legacy start/lap subset.");

experiments(end+1) = make_experiment_local( ...
    "legacy_eta0_raw", "legacy", "50.2", "10", ...
    "per_pass_main_pulse_query_safe", "adaptive20", true, "fixed_joint_static_eta", true, ...
    "Raw archived fixed-eta-zero route.");

experiments(end+1) = make_experiment_local( ...
    "legacy_start500", "legacy", "50.0", "20", ...
    "per_pass_main_pulse_query_safe", "adaptive20", true, "fixed_joint_static_eta", true, ...
    "Legacy fixed-eta route with official start/lap settings.");

experiments(end+1) = make_experiment_local( ...
    "legacy_nophase", "legacy", "50.2", "10", ...
    "per_pass_main_pulse_query_safe", "adaptive20", false, "fixed_joint_static_eta", true, ...
    "Legacy fixed-eta route with phase-safe expansion disabled.");

experiments(end+1) = make_experiment_local( ...
    "legacy_sgcore", "legacy", "50.2", "10", ...
    "sg_threshold_mainpulse_trust", "adaptive20", true, "fixed_joint_static_eta", true, ...
    "Legacy fixed-eta route with official core-mask mode.");

experiments(end+1) = make_experiment_local( ...
    "legacy_legpulse", "legacy", "50.2", "10", ...
    "per_pass_main_pulse_query_safe", "legacy_threshold", true, ...
    "fixed_joint_static_eta", true, ...
    "Legacy fixed-eta route with official main-pulse mode.");

experiments(end+1) = make_experiment_local( ...
    "legacy_allofficial", "legacy", "50.0", "20", ...
    "sg_threshold_mainpulse_trust", "legacy_threshold", false, ...
    "fixed_joint_static_eta", true, ...
    "Legacy code with official start/core/main-pulse/phase settings.");
end


function ex = make_experiment_local(tag, source, start_env, target_laps_env, ...
    core_mask_mode, main_pulse_mode, phase_safe_expansion, active_eta_model, ...
    static_eta_enable, description)
ex = struct();
ex.tag = string(tag);
ex.source = string(source);
ex.start_env = string(start_env);
ex.target_laps_env = string(target_laps_env);
ex.start_time_s = str2double(start_env);
ex.target_laps = str2double(target_laps_env);
ex.core_mask_mode = string(core_mask_mode);
ex.main_pulse_mode = string(main_pulse_mode);
ex.phase_safe_expansion = logical(phase_safe_expansion);
ex.active_eta_model = string(active_eta_model);
ex.static_eta_enable = logical(static_eta_enable);
ex.description = string(description);
end


function make_worker_local(source_file, worker_file, ex)
txt = fileread(source_file);

txt = regexprep(txt, ...
    "S\.method_name = '[^']*';", ...
    sprintf("S.method_name = 'route_ablation_%s';", ex.tag), 'once');
txt = regexprep(txt, ...
    "S\.method_file_tag = '[^']*';", ...
    sprintf("S.method_file_tag = 'RouteAblation_%s';", ex.tag), 'once');
txt = regexprep(txt, ...
    "S\.core_mask_mode = '[^']*';", ...
    sprintf("S.core_mask_mode = '%s';", ex.core_mask_mode), 'once');
txt = regexprep(txt, ...
    "S\.main_pulse_mode = '[^']*';", ...
    sprintf("S.main_pulse_mode = '%s';", ex.main_pulse_mode), 'once');
txt = regexprep(txt, ...
    "S\.phase_safe_expansion = (true|false);", ...
    sprintf("S.phase_safe_expansion = %s;", bool_text_local(ex.phase_safe_expansion)), 'once');
txt = regexprep(txt, ...
    "S\.active_eta_model = '[^']*';", ...
    sprintf("S.active_eta_model = '%s';", ex.active_eta_model), 'once');
txt = regexprep(txt, ...
    "S\.eta_models = \{S\.active_eta_model\};", ...
    "S.eta_models = {S.active_eta_model};", 'once');
txt = regexprep(txt, ...
    "S\.static_eta_source = '[^']*';", ...
    "S.static_eta_source = 'route_ablation_zero_eta';", 'once');
txt = regexprep(txt, ...
    "S\.static_eta_enable = (true|false);", ...
    sprintf("S.static_eta_enable = %s;", bool_text_local(ex.static_eta_enable)), 'once');
txt = regexprep(txt, ...
    "S\.sensor_eta_prior_mm = zeros\(1, numel\(S\.analysis_sensors\)\);", ...
    "S.sensor_eta_prior_mm = zeros(1, numel(S.analysis_sensors));", 'once');

% Workers are diagnostic-only; allow environment variables to drive the
% ablation matrix without tripping official route guards.
guard_patterns = {
    "if ~isempty\(core_mask_mode_env\)\s+error\([^\n]*(\n\s*'[^']*'\s*\.\.\.)?[\s\S]*?official foundation Step05 uses[^\n]*\);\s+end"
    "if ~isempty\(static_eta_env\)\s+error\([^\n]*\);\s+end"
    "if ~isempty\(static_eta_file_env\)\s+error\([^\n]*\);\s+end"
    "if ~isempty\(phase_safe_env\)\s+error\([^\n]*(\n\s*'[^']*'\s*\.\.\.)?[\s\S]*?official foundation Step05[^\n]*\);\s+end"
    };
for i = 1:numel(guard_patterns)
    txt = regexprep(txt, guard_patterns{i}, '');
end

fid = fopen(worker_file, 'w');
if fid < 0
    error('Cannot write worker file: %s', worker_file);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s', txt);
delete(cleanup);
end


function write_zero_eta_file_local(zero_eta_file, sensors)
Summary = struct();
Summary.TargetBlade = 1;
Summary.SensorIDs = sensors;
Summary.ReferenceSensorID = sensors(1);
Summary.Consensus = struct();
Summary.Consensus.status = "ok";
Summary.Consensus.diagnostic_status = "ok";
Summary.Consensus.eta_median_mm = zeros(1, numel(sensors));
Summary.Consensus.eta_plan_iqr_mm = zeros(1, numel(sensors));
Summary.Consensus.max_abs_plan_iqr_mm = 0;
Summary.Consensus.reference_sensor_id = sensors(1);
Summary.Consensus.eta_source_policy = "route_ablation_zero_eta";
Summary.Consensus.primary_plan_name = "zero";
save(zero_eta_file, 'Summary');
end


function txt = bool_text_local(tf)
if tf
    txt = 'true';
else
    txt = 'false';
end
end


function result_file = resolve_result_file_local(cfg, case_name, ex)
sensor_tag = ['S', sprintf('%d', cfg.sensor_ids)];
out_tag = ['abl_', char(ex.tag)];
if strcmp(ex.source, "official")
    result_name = sprintf( ...
        'Result_Step05_RouteAblation_%s_B1_%s_20251222.mat', ...
        char(ex.tag), sensor_tag);
else
    result_name = sprintf( ...
        'Result_Step05_RouteAblation_%s_B1_%s_20251222.mat', ...
        char(ex.tag), sensor_tag);
end

candidate = fullfile(cfg.output_root, ...
    ['step05_single_sync_direct_template_', out_tag], ...
    case_name, result_name);
if isfile(candidate)
    result_file = candidate;
    return;
end

files = dir(fullfile(cfg.output_root, ...
    ['step05_single_sync_direct_template_', out_tag], ...
    case_name, 'Result_Step05_*.mat'));
if isempty(files)
    result_file = candidate;
else
    result_file = fullfile(files(1).folder, files(1).name);
end
end


function row = empty_result_row_local()
row = struct( ...
    'tag', "", ...
    'source', "", ...
    'description', "", ...
    'worker_file', "", ...
    'result_file', "", ...
    'status', "", ...
    'failure_reason', "", ...
    'start_time_s', NaN, ...
    'target_laps', NaN, ...
    'core_mask_mode', "", ...
    'main_pulse_mode', "", ...
    'phase_safe_expansion', false, ...
    'active_eta_model', "", ...
    'static_eta_enable', false, ...
    'window_count', NaN, ...
    'eo_sequence', "", ...
    'median_rmse_v', NaN, ...
    'median_A_mm', NaN, ...
    'median_dx_c_mm', NaN, ...
    'A_upper_hit_fraction', NaN, ...
    'dx_lower_hit_fraction', NaN, ...
    'time_start_sequence', "");
end


function row = summarize_result_local(ex, worker_file, result_file)
row = base_row_local(ex, worker_file);
row.result_file = string(result_file);
if ~isfile(result_file)
    row.status = "missing_result_file";
    row.failure_reason = "Worker completed but expected result file is missing";
    return;
end

loaded = load(result_file, 'Result');
if ~isfield(loaded, 'Result') || ~isfield(loaded.Result, 'Trend')
    row.status = "invalid_result_file";
    row.failure_reason = "Result.Trend missing";
    return;
end

T = loaded.Result.Trend;
ok = strcmpi(string(T.status), 'ok');
if ~any(ok)
    row.status = "no_ok_windows";
    row.window_count = 0;
    return;
end

eo = T.EO_id(ok);
A = T.A_id(ok);
dx = T.dx_c_id(ok);
rmse = T.weighted_voltage_rmse(ok);
row.status = "ok";
row.window_count = nnz(ok);
row.eo_sequence = string(strjoin(compose('%d', eo(:).'), '|'));
row.median_rmse_v = median(rmse, 'omitnan');
row.median_A_mm = median(A, 'omitnan');
row.median_dx_c_mm = median(dx, 'omitnan');
row.A_upper_hit_fraction = mean(abs(A - 0.50) <= 1e-6, 'omitnan');
row.dx_lower_hit_fraction = mean(abs(dx + 0.35) <= 1e-6, 'omitnan');

starts = nan(numel(loaded.Result.WindowResult), 1);
for i = 1:numel(loaded.Result.WindowResult)
    tw = loaded.Result.WindowResult(i).time_window_s;
    if numel(tw) >= 1
        starts(i) = tw(1);
    end
end
starts = starts(isfinite(starts));
row.time_start_sequence = string(strjoin(compose('%.4f', starts(:).'), '|'));
end


function row = summarize_failure_local(ex, worker_file, ME)
row = base_row_local(ex, worker_file);
row.status = "failed";
row.failure_reason = string(ME.message);
end


function row = base_row_local(ex, worker_file)
row = empty_result_row_local();
row.tag = ex.tag;
row.source = ex.source;
row.description = ex.description;
row.worker_file = string(worker_file);
row.start_time_s = ex.start_time_s;
row.target_laps = ex.target_laps;
row.core_mask_mode = ex.core_mask_mode;
row.main_pulse_mode = ex.main_pulse_mode;
row.phase_safe_expansion = ex.phase_safe_expansion;
row.active_eta_model = ex.active_eta_model;
row.static_eta_enable = ex.static_eta_enable;
end


function env = capture_env_local()
names = {'STEP05_OUTPUT_TAG', 'STEP05_DEBUG_MAX_WINDOWS', ...
    'STEP05_TARGET_BLADES', 'STEP05_SHOW_PLOTS', 'STEP05_SAVE_FIGURES', ...
    'STEP05_ANALYSIS_START_TIME', 'STEP05_TARGET_LAPS', ...
    'STEP05_STATIC_ETA_FILE'};
env = struct();
for i = 1:numel(names)
    env.(names{i}) = getenv(names{i});
end
end


function restore_env_local(env)
names = fieldnames(env);
for i = 1:numel(names)
    setenv(names{i}, env.(names{i}));
end
end
