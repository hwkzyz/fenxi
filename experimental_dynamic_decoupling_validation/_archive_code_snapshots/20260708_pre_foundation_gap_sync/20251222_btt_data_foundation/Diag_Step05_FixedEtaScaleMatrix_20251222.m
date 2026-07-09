function Summary = Diag_Step05_FixedEtaScaleMatrix_20251222()
%DIAG_STEP05_FIXEDETASCALEMATRIX_20251222
% Diagnostic-only eta falsification matrix for the 20251222 foundation case.
%
% This script does not modify the official NoEta Step05 route. It reuses the
% archived fixed-joint-static-eta Step05 entry and tests whether the residual
% eta behaves like a physical sensor static bias under sign, scale, and
% channel-ablation perturbations.

clc; close all;

this_dir = fileparts(mfilename('fullpath'));
repo_dir = fileparts(this_dir);
legacy_dir = fullfile(this_dir, 'legacy_opr_center_flow_20260707');
legacy_step05 = fullfile(legacy_dir, ...
    'Step05_SingleSync_DirectTemplate_Identification_20251222.m');

addpath(repo_dir);
addpath(this_dir);

cfg = BTTProjectConfig_20251222();
case_name = cfg.dynamic_cases{1};
sensors = cfg.sensor_ids(:).';
target_blade = 1;
sensor_tag = ['S', sprintf('%d', sensors)];

out_dir = fullfile(cfg.output_root, 'diagnostics', ...
    'step05_fixed_eta_scale_matrix_20251222');
eta_dir = fullfile(out_dir, 'eta_variants');
if exist(out_dir, 'dir') ~= 7
    mkdir(out_dir);
end
if exist(eta_dir, 'dir') ~= 7
    mkdir(eta_dir);
end

worker_step05 = fullfile(out_dir, ...
    'DiagWorker_Step05_FixedJointEta_20251222.m');
copyfile(legacy_step05, worker_step05, 'f');

base_eta_file = resolve_base_eta_file_local(cfg, target_blade, sensor_tag);
base_eta = load_eta_vector_local(base_eta_file, sensors);

max_windows = strtrim(getenv('STEP05_DIAG_MAX_WINDOWS'));
if isempty(max_windows)
    max_windows = '6';
end

variants = build_eta_variants_local(base_eta);
rows = repmat(empty_summary_row_local(), numel(variants), 1);

old_env = capture_env_local();
cleanup_obj = onCleanup(@() restore_env_local(old_env)); %#ok<NASGU>

old_dir = pwd;
dir_cleanup = onCleanup(@() cd(old_dir)); %#ok<NASGU>
cd(this_dir);

for iv = 1:numel(variants)
    v = variants(iv);
    eta_file = write_eta_summary_file_local( ...
        eta_dir, target_blade, sensors, sensor_tag, v.tag, v.eta);

    setenv('STEP05_STATIC_ETA_FILE', eta_file);
    setenv('STEP05_OUTPUT_TAG', ['etadiag_', v.tag]);
    setenv('STEP05_DEBUG_MAX_WINDOWS', max_windows);
    setenv('STEP05_TARGET_BLADES', num2str(target_blade));
    setenv('STEP05_SHOW_PLOTS', '0');
    setenv('STEP05_SAVE_FIGURES', '0');

    fprintf('\n=== Eta variant %d/%d: %s ===\n', iv, numel(variants), v.tag);
    fprintf('eta = %s mm\n', mat2str(v.eta, 8));

    run(worker_step05);

    result_file = fullfile(cfg.output_root, ...
        ['step05_single_sync_direct_template_etadiag_', v.tag], ...
        case_name, sprintf( ...
        'Result_Step05_FoundationMainPulseAdaptiveFixedJointEta_B%d_%s_20251222.mat', ...
        target_blade, sensor_tag));
    rows(iv) = summarize_variant_result_local(v, eta_file, result_file);
end

VariantTable = struct2table(rows);
summary_csv = fullfile(out_dir, ...
    sprintf('Step05_FixedEtaScaleMatrix_B%d_%s_20251222.csv', ...
    target_blade, sensor_tag));
writetable(VariantTable, summary_csv);

Summary = struct();
Summary.output_dir = out_dir;
Summary.base_eta_file = base_eta_file;
Summary.base_eta_mm = base_eta;
Summary.max_windows = max_windows;
Summary.VariantTable = VariantTable;

summary_mat = fullfile(out_dir, ...
    sprintf('Step05_FixedEtaScaleMatrix_B%d_%s_20251222.mat', ...
    target_blade, sensor_tag));
save(summary_mat, 'Summary', '-v7.3');

fprintf('\n=== Fixed eta scale/sign diagnostic complete ===\n');
fprintf('Base eta: %s mm\n', mat2str(base_eta, 8));
fprintf('Max windows: %s\n', max_windows);
fprintf('Summary CSV:\n  %s\n', summary_csv);
disp(VariantTable(:, {'variant_tag', 'eta_mm_text', 'window_count', ...
    'eo_sequence', 'median_rmse_v', 'median_A_mm', 'median_dx_c_mm', ...
    'A_upper_hit_fraction', 'dx_lower_hit_fraction'}));
end


function eta_file = resolve_base_eta_file_local(cfg, target_blade, sensor_tag)
preferred = {
    fullfile(cfg.output_root, 'step05_joint_static_eta_preview_current_noeta_voltage', ...
    sprintf('JointStaticEtaPreview_B%d_%s.mat', target_blade, sensor_tag))
    fullfile(cfg.output_root, 'step05_joint_static_eta_preview', ...
    sprintf('JointStaticEtaPreview_B%d_%s.mat', target_blade, sensor_tag))
    };

for i = 1:numel(preferred)
    if isfile(preferred{i})
        eta_file = preferred{i};
        return;
    end
end

error('No base eta preview file found. Run Calibrate_Step05_JointStaticEta_FromPreview_20251222 first.');
end


function eta = load_eta_vector_local(eta_file, sensors)
loaded = load(eta_file, 'Summary');
if ~isfield(loaded, 'Summary') || ~isfield(loaded.Summary, 'Consensus') || ...
        ~isfield(loaded.Summary.Consensus, 'eta_median_mm')
    error('Invalid eta preview file: %s', eta_file);
end
eta = double(loaded.Summary.Consensus.eta_median_mm(:).');
if numel(eta) ~= numel(sensors) || any(~isfinite(eta))
    error('Eta vector in %s does not match sensors %s.', ...
        eta_file, mat2str(sensors));
end
eta(1) = 0;
end


function variants = build_eta_variants_local(base_eta)
z = zeros(size(base_eta));
variants = repmat(struct('tag', '', 'eta', []), 0, 1);
variants(end+1) = struct('tag', 'zero', 'eta', z);
variants(end+1) = struct('tag', 'base', 'eta', base_eta);
variants(end+1) = struct('tag', 'flip', 'eta', -base_eta);
variants(end+1) = struct('tag', 's025', 'eta', 0.25 .* base_eta);
variants(end+1) = struct('tag', 's050', 'eta', 0.50 .* base_eta);

eta_ch2 = z;
eta_ch2(2) = base_eta(2);
variants(end+1) = struct('tag', 'ch2', 'eta', eta_ch2);

eta_ch3 = z;
if numel(base_eta) >= 3
    eta_ch3(3) = base_eta(3);
end
variants(end+1) = struct('tag', 'ch3', 'eta', eta_ch3);
end


function eta_file = write_eta_summary_file_local( ...
    eta_dir, target_blade, sensors, sensor_tag, tag, eta)
Summary = struct();
Summary.TargetBlade = target_blade;
Summary.SensorIDs = sensors;
Summary.ReferenceSensorID = sensors(1);
Summary.Consensus = struct();
Summary.Consensus.status = "ok";
Summary.Consensus.diagnostic_status = "ok";
Summary.Consensus.eta_median_mm = eta;
Summary.Consensus.eta_plan_iqr_mm = zeros(size(eta));
Summary.Consensus.max_abs_plan_iqr_mm = 0;
Summary.Consensus.reference_sensor_id = sensors(1);
Summary.Consensus.eta_source_policy = "diagnostic_eta_variant";
Summary.Consensus.primary_plan_name = string(tag);
Summary.Options = struct('diagnostic_variant_tag', tag);

eta_file = fullfile(eta_dir, sprintf( ...
    'JointStaticEtaPreview_B%d_%s_%s.mat', target_blade, sensor_tag, tag));
save(eta_file, 'Summary');
end


function row = empty_summary_row_local()
row = struct( ...
    'variant_tag', "", ...
    'eta_mm_text', "", ...
    'eta_file', "", ...
    'result_file', "", ...
    'status', "", ...
    'failure_reason', "", ...
    'window_count', NaN, ...
    'eo_sequence', "", ...
    'eo_mode', NaN, ...
    'eo14_fraction', NaN, ...
    'median_rmse_v', NaN, ...
    'median_A_mm', NaN, ...
    'median_dx_c_mm', NaN, ...
    'A_upper_hit_fraction', NaN, ...
    'dx_lower_hit_fraction', NaN);
end


function row = summarize_variant_result_local(v, eta_file, result_file)
row = empty_summary_row_local();
row.variant_tag = string(v.tag);
row.eta_mm_text = string(mat2str(v.eta, 8));
row.eta_file = string(eta_file);
row.result_file = string(result_file);

if ~isfile(result_file)
    row.status = "missing_result_file";
    row.failure_reason = "Step05 run did not produce expected result file";
    return;
end

loaded = load(result_file, 'Result');
if ~isfield(loaded, 'Result') || ~isfield(loaded.Result, 'Trend') || ...
        ~istable(loaded.Result.Trend)
    row.status = "invalid_result_file";
    row.failure_reason = "Result.Trend missing";
    return;
end

T = loaded.Result.Trend;
ok = strcmpi(string(T.status), 'ok');
if ~any(ok)
    row.status = "no_ok_windows";
    row.failure_reason = "No ok windows";
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
row.eo_mode = mode(round(eo(isfinite(eo))));
row.eo14_fraction = mean(round(eo) == 14, 'omitnan');
row.median_rmse_v = median(rmse, 'omitnan');
row.median_A_mm = median(A, 'omitnan');
row.median_dx_c_mm = median(dx, 'omitnan');
row.A_upper_hit_fraction = mean(abs(A - 0.50) <= 1e-6, 'omitnan');
row.dx_lower_hit_fraction = mean(abs(dx + 0.35) <= 1e-6, 'omitnan');
end


function env = capture_env_local()
names = {'STEP05_STATIC_ETA_FILE', 'STEP05_OUTPUT_TAG', ...
    'STEP05_DEBUG_MAX_WINDOWS', 'STEP05_TARGET_BLADES', ...
    'STEP05_SHOW_PLOTS', 'STEP05_SAVE_FIGURES'};
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
