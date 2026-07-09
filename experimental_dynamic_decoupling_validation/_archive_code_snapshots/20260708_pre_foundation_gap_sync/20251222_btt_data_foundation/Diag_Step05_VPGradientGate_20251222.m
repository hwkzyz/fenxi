function Summary = Diag_Step05_VPGradientGate_20251222()
%DIAG_STEP05_VPGRADIENTGATE_20251222
% Test whether low-gradient points are responsible for the legacy fixed-eta
% VP seed pathology. Workers are generated under output/diagnostics only.

clc; close all;

this_dir = fileparts(mfilename('fullpath'));
cfg = BTTProjectConfig_20251222();
case_name = cfg.dynamic_cases{1};
out_dir = fullfile(cfg.output_root, 'diagnostics', ...
    'step05_vp_gradient_gate_20251222');
worker_dir = fullfile(out_dir, 'workers');
if exist(worker_dir, 'dir') ~= 7
    mkdir(worker_dir);
end

legacy_file = fullfile(this_dir, 'legacy_opr_center_flow_20260707', ...
    'Step05_SingleSync_DirectTemplate_Identification_20251222.m');
zero_eta_file = fullfile(out_dir, 'VPGradientGate_ZeroEta_B1_S123.mat');
write_zero_eta_file_local(zero_eta_file, cfg.sensor_ids(:).');

experiments = build_experiment_set_local();
rows = repmat(empty_result_row_local(), numel(experiments), 1);
seed_rows = repmat(empty_seed_row_local(), 0, 1);

old_env = capture_env_local();
cleanup_obj = onCleanup(@() restore_env_local(old_env)); %#ok<NASGU>

old_dir = pwd;
dir_cleanup = onCleanup(@() cd(old_dir)); %#ok<NASGU>
cd(this_dir);
addpath(this_dir);
addpath(fileparts(this_dir));

for i = 1:numel(experiments)
    ex = experiments(i);
    fprintf('\n=== VP gradient gate %d/%d: %s ===\n', ...
        i, numel(experiments), ex.tag);

    worker_file = fullfile(worker_dir, sprintf('Worker_%02d_%s.m', ...
        i, ex.tag));
    make_worker_local(legacy_file, worker_file, ex);

    setenv('STEP05_OUTPUT_TAG', char(ex.output_tag));
    setenv('STEP05_DEBUG_MAX_WINDOWS', '2');
    setenv('STEP05_TARGET_BLADES', '1');
    setenv('STEP05_SHOW_PLOTS', '0');
    setenv('STEP05_SAVE_FIGURES', '0');
    setenv('STEP05_ANALYSIS_START_TIME', '50.0');
    setenv('STEP05_TARGET_LAPS', '20');
    setenv('STEP05_STATIC_ETA_FILE', zero_eta_file);

    try
        run(worker_file);
        result_file = resolve_result_file_local(cfg, case_name, ex);
        rows(i) = summarize_result_local(ex, worker_file, result_file);
        seed_rows = [seed_rows; summarize_seed_rows_local(ex, result_file)]; %#ok<AGROW>
    catch ME
        rows(i) = summarize_failure_local(ex, worker_file, ME);
    end
end

ResultTable = struct2table(rows);
SeedTable = struct2table(seed_rows);

summary_csv = fullfile(out_dir, 'Step05_VPGradientGate_Summary_B1_S123_20251222.csv');
seed_csv = fullfile(out_dir, 'Step05_VPGradientGate_Seeds_B1_S123_20251222.csv');
writetable(ResultTable, summary_csv);
writetable(SeedTable, seed_csv);

Summary = struct();
Summary.OutputDir = out_dir;
Summary.ResultTable = ResultTable;
Summary.SeedTable = SeedTable;
Summary.SummaryCSV = summary_csv;
Summary.SeedCSV = seed_csv;
save(fullfile(out_dir, 'Step05_VPGradientGate_B1_S123_20251222.mat'), ...
    'Summary', '-v7.3');

fprintf('\n=== Step05 VP gradient-gate summary ===\n');
disp(ResultTable(:, {'tag', 'gradient_gate_enable', 'gradient_min_ratio', ...
    'gradient_weight_power', 'seed_clamp_enable', 'vp_rank_score_mode', ...
    'eo_sequence', 'median_rmse_v', 'median_A_mm', 'median_dx_c_mm', ...
    'A_upper_hit_fraction', 'dx_lower_hit_fraction', 'status'}));
fprintf('Saved:\n  %s\n  %s\n', summary_csv, seed_csv);
end


function experiments = build_experiment_set_local()
base = struct( ...
    'tag', "", ...
    'output_tag', "", ...
    'gradient_gate_enable', false, ...
    'gradient_min_ratio', 0.10, ...
    'gradient_reference_quantile', 95, ...
    'gradient_weight_power', 0, ...
    'seed_clamp_enable', false, ...
    'vp_rank_score_mode', "linear_residual", ...
    'description', "");
experiments = repmat(base, 0, 1);

experiments(end+1) = make_experiment_local( ...
    "legacy_allofficial_baseline", "vpgg_base", false, 0.10, 95, 0, false, ...
    "linear_residual", "Baseline legacy all-official fixed eta=0.");
experiments(end+1) = make_experiment_local( ...
    "fullseed_only", "vpgg_full", false, 0.10, 95, 0, false, ...
    "full_seed", "Only switch ranking from linear residual to full-seed voltage RMSE.");
experiments(end+1) = make_experiment_local( ...
    "clamp_fullseed", "vpgg_cfull", false, 0.10, 95, 0, true, ...
    "full_seed", "Seed clamp plus full-seed voltage RMSE ranking.");
experiments(end+1) = make_experiment_local( ...
    "gate_only", "vpgg_gate", true, 0.10, 95, 0, false, ...
    "linear_residual", "VP-only low-gradient removal.");
experiments(end+1) = make_experiment_local( ...
    "gate_weight", "vpgg_wgt", true, 0.10, 95, 2, false, ...
    "linear_residual", "VP-only low-gradient removal plus gradient weighting.");
experiments(end+1) = make_experiment_local( ...
    "gate_weight_clamp", "vpgg_clamp", true, 0.10, 95, 2, true, ...
    "linear_residual", "Gradient gate/weight plus seed parameter clamp.");
experiments(end+1) = make_experiment_local( ...
    "gate_weight_clamp_fullseed", "vpgg_gwcf", true, 0.10, 95, 2, true, ...
    "full_seed", "Gate/weight/clamp plus full-seed voltage RMSE ranking.");
experiments(end+1) = make_experiment_local( ...
    "gate25_weight_clamp_fullseed", "vpgg_g25f", true, 0.25, 95, 2, true, ...
    "full_seed", "Stronger VP gradient gate plus weight/clamp/full-seed ranking.");
end


function ex = make_experiment_local(tag, output_tag, gate, ratio, quantile, power, clamp, rank_mode, description)
ex = struct();
ex.tag = string(tag);
ex.output_tag = string(output_tag);
ex.gradient_gate_enable = logical(gate);
ex.gradient_min_ratio = ratio;
ex.gradient_reference_quantile = quantile;
ex.gradient_weight_power = power;
ex.seed_clamp_enable = logical(clamp);
ex.vp_rank_score_mode = string(rank_mode);
ex.description = string(description);
end


function make_worker_local(source_file, worker_file, ex)
txt = fileread(source_file);

txt = regexprep(txt, "S\.method_name = '[^']*';", ...
    sprintf("S.method_name = 'vp_gradient_gate_%s';", ex.tag), 'once');
txt = regexprep(txt, "S\.method_file_tag = '[^']*';", ...
    sprintf("S.method_file_tag = 'VPGradientGate_%s';", ex.tag), 'once');
txt = regexprep(txt, "S\.core_mask_mode = '[^']*';", ...
    "S.core_mask_mode = 'sg_threshold_mainpulse_trust';", 'once');
txt = regexprep(txt, "S\.main_pulse_mode = '[^']*';", ...
    "S.main_pulse_mode = 'legacy_threshold';", 'once');
txt = regexprep(txt, "S\.phase_safe_expansion = (true|false);", ...
    "S.phase_safe_expansion = false;", 'once');
txt = regexprep(txt, "S\.active_eta_model = '[^']*';", ...
    "S.active_eta_model = 'fixed_joint_static_eta';", 'once');
txt = regexprep(txt, "S\.static_eta_source = '[^']*';", ...
    "S.static_eta_source = 'vp_gradient_gate_zero_eta';", 'once');
txt = regexprep(txt, "S\.static_eta_enable = (true|false);", ...
    "S.static_eta_enable = true;", 'once');
txt = regexprep(txt, "S\.vp_rank_score_mode = '[^']*';", ...
    sprintf('S.vp_rank_score_mode = ''%s'';', char(ex.vp_rank_score_mode)), 'once');

settings = sprintf([ ...
    '\nS.vp_gradient_gate_enable = %s;' ...
    '\nS.vp_gradient_min_ratio = %.12g;' ...
    '\nS.vp_gradient_reference_quantile = %.12g;' ...
    '\nS.vp_gradient_weight_power = %.12g;' ...
    '\nS.vp_seed_clamp_enable = %s;\n'], ...
    bool_text_local(ex.gradient_gate_enable), ex.gradient_min_ratio, ...
    ex.gradient_reference_quantile, ex.gradient_weight_power, ...
    bool_text_local(ex.seed_clamp_enable));
txt = regexprep(txt, "(S\.vp_rank_score_mode = '[^']*';)", ...
    "$1" + string(settings), 'once');

gate_snip = sprintf([ ...
    'grad_abs_vp = abs(Tp);\n' ...
    'grad_ref_vp = prctile(grad_abs_vp(isfinite(grad_abs_vp)), ' ...
    'getfield_default_local(S, ''vp_gradient_reference_quantile'', 95));\n' ...
    'if isfield(S, ''vp_gradient_gate_enable'') && S.vp_gradient_gate_enable && ' ...
    'isfinite(grad_ref_vp) && grad_ref_vp > 0\n' ...
    '    grad_ok_vp = grad_abs_vp >= max(0, S.vp_gradient_min_ratio) * grad_ref_vp;\n' ...
    '    v = v(grad_ok_vp);\n' ...
    '    theta = theta(grad_ok_vp);\n' ...
    '    T0 = T0(grad_ok_vp);\n' ...
    '    Tp = Tp(grad_ok_vp);\n' ...
    '    w = w(grad_ok_vp);\n' ...
    '    sensor_id = sensor_id(grad_ok_vp);\n' ...
    'end\n']);
txt = regexprep(txt, "sensor_id = sensor_id\(valid\);\s*", ...
    sprintf('sensor_id = sensor_id(valid);\n%s', gate_snip), 'once');

weight_snip = sprintf([ ...
    'grad_power_vp = getfield_default_local(S, ''vp_gradient_weight_power'', 0);\n' ...
    'if grad_power_vp > 0\n' ...
    '    grad_abs_w = abs(Tp);\n' ...
    '    grad_ref_w = prctile(grad_abs_w(isfinite(grad_abs_w)), ' ...
    'getfield_default_local(S, ''vp_gradient_reference_quantile'', 95));\n' ...
    '    if isfinite(grad_ref_w) && grad_ref_w > 0\n' ...
    '        grad_weight_vp = min(grad_abs_w ./ max(grad_ref_w, eps), 1) .^ grad_power_vp;\n' ...
    '        wq = max(S.weight_floor, wq .* grad_weight_vp);\n' ...
    '    end\n' ...
    'end\n']);
txt = regexprep(txt, "wq = max\(wq, S\.weight_floor\);\s*", ...
    sprintf('wq = max(wq, S.weight_floor);\n%s', weight_snip), 'once');

clamp_snip = sprintf([ ...
    'if isfield(S, ''vp_seed_clamp_enable'') && S.vp_seed_clamp_enable\n' ...
    '    if isfinite(S.amplitude_limit_mm)\n' ...
    '        A = min(abs(A), abs(S.amplitude_limit_mm));\n' ...
    '    end\n' ...
    '    dx_c = clamp_scalar_local(dx_c, S.dx_c_limit_mm);\n' ...
    '    eta = clamp_vector_local(eta, S.sensor_eta_limit_mm);\n' ...
    'end\n']);
txt = regexprep(txt, ...
    "eta\(~isfinite\(eta\)\) = 0;\s*\n\s*\[wrmse, prmse\] = evaluate_full_template_model_local", ...
    sprintf('eta(~isfinite(eta)) = 0;\n%s\n    [wrmse, prmse] = evaluate_full_template_model_local', clamp_snip), ...
    'once');

fid = fopen(worker_file, 'w');
if fid < 0
    error('Cannot write worker file: %s', worker_file);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s', txt);
delete(cleanup);
end


function row = summarize_result_local(ex, worker_file, result_file)
row = base_row_local(ex, worker_file);
row.result_file = string(result_file);
if ~isfile(result_file)
    row.status = "missing_result_file";
    row.failure_reason = "Expected result file missing";
    return;
end
loaded = load(result_file, 'Result');
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
end


function rows = summarize_seed_rows_local(ex, result_file)
rows = repmat(empty_seed_row_local(), 0, 1);
if ~isfile(result_file)
    return;
end
loaded = load(result_file, 'Result');
target_eo = [3 10 12 14 23];
for iw = 1:min(2, numel(loaded.Result.WindowResult))
    if ~isfield(loaded.Result.WindowResult(iw), 'Result')
        continue;
    end
    R = loaded.Result.WindowResult(iw).Result;
    if ~isfield(R, 'VPSeedTable') || isempty(R.VPSeedTable)
        continue;
    end
    T = R.VPSeedTable;
    for eo = target_eo
        idx = find(T.EO == eo, 1, 'first');
        if isempty(idx)
            continue;
        end
        row = empty_seed_row_local();
        row.tag = ex.tag;
        row.window_id = iw;
        row.EO = eo;
        row.A_seed_mm = table_value_local(T, idx, 'A');
        row.dx_seed_mm = table_value_local(T, idx, 'dx_c');
        row.rank_score = table_value_local(T, idx, 'rank_score');
        row.vp_linear_rmse = table_value_local(T, idx, 'vp_linear_rmse');
        row.weighted_voltage_rmse = table_value_local(T, idx, 'weighted_voltage_rmse');
        rows(end+1, 1) = row; %#ok<AGROW>
    end
end
end


function value = table_value_local(T, idx, name)
value = NaN;
if istable(T) && ismember(name, T.Properties.VariableNames) && idx <= height(T)
    value = T.(name)(idx);
end
end


function row = empty_result_row_local()
row = struct( ...
    'tag', "", ...
    'description', "", ...
    'worker_file', "", ...
    'result_file', "", ...
    'status', "", ...
    'failure_reason', "", ...
    'gradient_gate_enable', false, ...
    'gradient_min_ratio', NaN, ...
    'gradient_weight_power', NaN, ...
    'seed_clamp_enable', false, ...
    'vp_rank_score_mode', "", ...
    'window_count', NaN, ...
    'eo_sequence', "", ...
    'median_rmse_v', NaN, ...
    'median_A_mm', NaN, ...
    'median_dx_c_mm', NaN, ...
    'A_upper_hit_fraction', NaN, ...
    'dx_lower_hit_fraction', NaN);
end


function row = empty_seed_row_local()
row = struct( ...
    'tag', "", ...
    'window_id', NaN, ...
    'EO', NaN, ...
    'A_seed_mm', NaN, ...
    'dx_seed_mm', NaN, ...
    'rank_score', NaN, ...
    'vp_linear_rmse', NaN, ...
    'weighted_voltage_rmse', NaN);
end


function row = summarize_failure_local(ex, worker_file, ME)
row = base_row_local(ex, worker_file);
row.status = "failed";
row.failure_reason = string(ME.message);
end


function row = base_row_local(ex, worker_file)
row = empty_result_row_local();
row.tag = ex.tag;
row.description = ex.description;
row.worker_file = string(worker_file);
row.gradient_gate_enable = ex.gradient_gate_enable;
row.gradient_min_ratio = ex.gradient_min_ratio;
row.gradient_weight_power = ex.gradient_weight_power;
row.seed_clamp_enable = ex.seed_clamp_enable;
row.vp_rank_score_mode = ex.vp_rank_score_mode;
end


function result_file = resolve_result_file_local(cfg, case_name, ex)
sensor_tag = ['S', sprintf('%d', cfg.sensor_ids)];
out_tag = char(ex.output_tag);
result_name = sprintf( ...
    'Result_Step05_VPGradientGate_%s_B1_%s_20251222.mat', ...
    char(ex.tag), sensor_tag);
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
Summary.Consensus.eta_source_policy = "vp_gradient_gate_zero_eta";
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
