function Summary = Compare_Step05_FreeEta_20250527()
%COMPARE_STEP05_FREEETA_20250527
% Separate FreeEta pilot for deriving fixed eta outside the formal Step05.

clc; close all;

this_dir = fileparts(mfilename('fullpath'));
addpath(fileparts(this_dir));
cfg = BTTProjectConfig_20250527();
case_name = cfg.dynamic_cases{1};
sensors = cfg.sensor_ids(:).';
target_blade = resolve_single_target_blade_local(cfg);
sensor_tag = ['S', sprintf('%d', sensors)];

out_dir = fullfile(cfg.output_root, 'diagnostics', ...
    'step05_free_eta_compare_20250527');
if exist(out_dir, 'dir') ~= 7
    mkdir(out_dir);
end

worker_file = fullfile(out_dir, 'DiagWorker_Step05_FreeEta_20250527.m');
build_free_eta_worker_local(this_dir, worker_file);

max_windows = strtrim(getenv('STEP05_COMPARE_FREE_MAX_WINDOWS'));
if isempty(max_windows)
    max_windows = 'inf';
end

old_env = capture_env_local();
cleanup_obj = onCleanup(@() restore_env_local(old_env)); %#ok<NASGU>
old_dir = pwd;
dir_cleanup = onCleanup(@() cd(old_dir)); %#ok<NASGU>

setenv('STEP05_OUTPUT_TAG', 'freeeta_pilot');
setenv('STEP05_DEBUG_MAX_WINDOWS', max_windows);
setenv('STEP05_TARGET_BLADES', num2str(target_blade));
setenv('STEP05_SHOW_PLOTS', '0');
setenv('STEP05_SAVE_FIGURES', '0');
setenv('STEP05_STATIC_ETA_FILE', '');
setenv('STEP05_ETA_MODEL', '');

cd(this_dir);
run(worker_file);

result_file = fullfile(cfg.output_root, ...
    'step05_single_sync_direct_template_freeeta_pilot', ...
    case_name, sprintf( ...
    'Result_Step05_SingleSyncDirectTemplate_B%d_%s_20250527.mat', ...
    target_blade, sensor_tag));

Summary = summarize_free_eta_result_local( ...
    result_file, sensors, out_dir, max_windows);

fprintf('\n=== Step05 FreeEta pilot (20250527) ===\n');
fprintf('FreeEta result:\n  %s\n', result_file);
fprintf('Free eta median: %s mm\n', mat2str(Summary.free_eta_median_mm, 8));
fprintf('Free eta IQR:    %s mm\n', mat2str(Summary.free_eta_iqr_mm, 8));
fprintf('EO sequence: %s\n', Summary.eo_sequence);
fprintf('Median RMSE: %.6g V\n', Summary.median_rmse_v);
end


function build_free_eta_worker_local(this_dir, worker_file)
main_file = fullfile(this_dir, ...
    'Step05_SingleSync_DirectTemplate_Identification_20250527.m');
txt = fileread(main_file);

replacements = {
    "S.method_name = 'foundation_mainpulse_adaptive_fixed_joint_static_eta';", ...
    "S.method_name = 'foundation_mainpulse_adaptive_freeeta_pilot';"
    "S.version = 'FOUNDATION_MAINPULSE_ADAPTIVE_FIXEDJOINTETA_20260707';", ...
    "S.version = 'FOUNDATION_MAINPULSE_ADAPTIVE_FREEETA_PILOT_20260708';"
    "S.fit_sensor_eta = false;", ...
    "S.fit_sensor_eta = true;"
    "S.active_eta_model = 'fixed_joint_static_eta';", ...
    "S.active_eta_model = 'window_free_eta';"
    "S.sensor_eta_mode = 'fixed_static';", ...
    "S.sensor_eta_mode = 'window_free';"
    "S.static_eta_source = 'joint_dynamic_residual_eta_preview';", ...
    "S.static_eta_source = 'none';"
    };

for i = 1:size(replacements, 1)
    txt = strrep(txt, replacements{i, 1}, replacements{i, 2});
end

fid = fopen(worker_file, 'w');
if fid < 0
    error('Cannot create FreeEta worker: %s', worker_file);
end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, '%s', txt);
end


function Summary = summarize_free_eta_result_local( ...
    result_file, sensors, out_dir, max_windows)
if ~isfile(result_file)
    error('FreeEta result file not found: %s', result_file);
end

loaded_result = load(result_file, 'Result');
Result = loaded_result.Result;
T = Result.Trend;
ok = strcmpi(string(T.status), 'ok');

eta = nan(height(T), numel(sensors));
for iw = 1:numel(Result.WindowResult)
    if isfield(Result.WindowResult(iw), 'Result') && ...
            isfield(Result.WindowResult(iw).Result, 'sensor_eta_id')
        e = Result.WindowResult(iw).Result.sensor_eta_id(:).';
        n = min(numel(e), numel(sensors));
        eta(iw, 1:n) = e(1:n);
    end
end

free_median = median(eta(ok, :), 1, 'omitnan');
free_iqr = iqr(eta(ok, :), 1);

WindowTable = table();
WindowTable.window_id = T.window_id;
WindowTable.status = string(T.status);
WindowTable.EO_id = T.EO_id;
WindowTable.A_mm = T.A_id;
WindowTable.dx_c_mm = T.dx_c_id;
WindowTable.weighted_voltage_rmse_v = T.weighted_voltage_rmse;
for is = 1:numel(sensors)
    WindowTable.(sprintf('eta_CH%d_mm', sensors(is))) = eta(:, is);
end

target_blade = double(Result.TargetBlade);
sensor_tag = ['S', sprintf('%d', sensors)];
summary_csv = fullfile(out_dir, sprintf( ...
    'Step05_FreeEta_Compare_Windows_B%d_%s_20250527.csv', ...
    target_blade, sensor_tag));
writetable(WindowTable, summary_csv);

Summary = struct();
Summary.result_file = result_file;
Summary.max_windows = max_windows;
Summary.target_blade = target_blade;
Summary.sensor_ids = sensors;
Summary.free_eta_median_mm = free_median;
Summary.free_eta_iqr_mm = free_iqr;
Summary.window_count = nnz(ok);
Summary.eo_sequence = strjoin(compose('%d', T.EO_id(ok).'), '|');
Summary.median_rmse_v = median(T.weighted_voltage_rmse(ok), 'omitnan');
Summary.median_A_mm = median(T.A_id(ok), 'omitnan');
Summary.median_dx_c_mm = median(T.dx_c_id(ok), 'omitnan');
Summary.WindowTable = WindowTable;
Summary.summary_csv = summary_csv;
Summary.summary_mat = fullfile(out_dir, sprintf( ...
    'Step05_FreeEta_Compare_B%d_%s_20250527.mat', ...
    target_blade, sensor_tag));
save(Summary.summary_mat, 'Summary', '-v7.3');
end


function target_blade = resolve_single_target_blade_local(cfg)
target_blade = 1;
env_text = strtrim(getenv('STEP05_TARGET_BLADES'));
if ~isempty(env_text)
    tokens = regexp(env_text, '[,;\s]+', 'split');
    vals = [];
    for i = 1:numel(tokens)
        if isempty(tokens{i})
            continue;
        end
        v = str2double(tokens{i});
        if isfinite(v)
            vals(end+1) = round(v); %#ok<AGROW>
        end
    end
    vals = unique(vals, 'stable');
    if numel(vals) ~= 1
        error('Compare_Step05_FreeEta_20250527 requires exactly one target blade.');
    end
    target_blade = vals(1);
end
if target_blade < 1 || target_blade > cfg.blades_num
    error('Target blade B%d is outside 1..%d.', target_blade, cfg.blades_num);
end
end


function env = capture_env_local()
names = {'STEP05_OUTPUT_TAG', 'STEP05_DEBUG_MAX_WINDOWS', ...
    'STEP05_TARGET_BLADES', 'STEP05_SHOW_PLOTS', 'STEP05_SAVE_FIGURES', ...
    'STEP05_STATIC_ETA_FILE', 'STEP05_ETA_MODEL'};
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
