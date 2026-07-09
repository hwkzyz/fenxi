function SummarySet = Build_Step05_FixedEta_FromFreeEtaPilot_20251222()
%BUILD_STEP05_FIXEDETA_FROMFREEETAPILOT_20251222
% Estimate per-blade fixed eta from the separate FreeEta pilot route.
%
% Formal Step05 remains fixed-eta. This script runs the diagnostic FreeEta
% worker for each requested blade, summarizes per-window sensor_eta_id, and
% writes the JointStaticEtaPreview_B*_S*.mat files consumed by formal Step05.
%
% Optional environment variables:
%   STEP05_TARGET_BLADES              e.g. "1" or "1 2 3 4 5 6"
%   STEP05_COMPARE_FREE_MAX_WINDOWS   default inherited by FreeEta compare

clc; close all;

this_dir = fileparts(mfilename('fullpath'));
cfg = BTTProjectConfig_20251222();
sensors = cfg.sensor_ids(:).';
sensor_tag = ['S', sprintf('%d', sensors)];
target_blades = resolve_target_blades_local(cfg);

preview_dir = fullfile(cfg.output_root, 'step05_joint_static_eta_preview');
if exist(preview_dir, 'dir') ~= 7
    mkdir(preview_dir);
end

old_env = capture_env_local();
cleanup_obj = onCleanup(@() restore_env_local(old_env)); %#ok<NASGU>
old_dir = pwd;
dir_cleanup = onCleanup(@() cd(old_dir)); %#ok<NASGU>
cd(this_dir);

SummarySet = repmat(empty_summary_row_local(sensors), numel(target_blades), 1);

fprintf('\n=== Build fixed eta from FreeEta pilot (20251222) ===\n');
fprintf('Target blades: %s; sensors: %s\n', ...
    mat2str(target_blades), mat2str(sensors));

for ib = 1:numel(target_blades)
    blade_id = target_blades(ib);
    setenv('STEP05_TARGET_BLADES', num2str(blade_id));

    Free = Compare_Step05_FreeEta_20251222();
    [Preview, PlanTable, Reference] = ...
        build_preview_from_free_summary_local(Free, cfg, sensors);

    preview_mat = fullfile(preview_dir, sprintf( ...
        'JointStaticEtaPreview_B%d_%s.mat', blade_id, sensor_tag));
    plans_csv = fullfile(preview_dir, sprintf( ...
        'JointStaticEtaPreview_Plans_B%d_%s.csv', blade_id, sensor_tag));
    reference_csv = fullfile(preview_dir, sprintf( ...
        'JointStaticEtaPreview_Reference_B%d_%s.csv', blade_id, sensor_tag));

    Summary = Preview; %#ok<NASGU>
    save(preview_mat, 'Summary', '-v7.3');
    writetable(PlanTable, plans_csv);
    writetable(Reference, reference_csv);

    SummarySet(ib) = make_summary_row_local( ...
        blade_id, sensors, Preview, Free, preview_mat, plans_csv, reference_csv);

    fprintf('B%d fixed eta: %s mm; IQR: %s mm; EO seq: %s\n', ...
        blade_id, mat2str(Preview.Consensus.eta_median_mm, 8), ...
        mat2str(Preview.Consensus.eta_window_iqr_mm, 8), ...
        Free.eo_sequence);
end

SummaryTable = struct2table(SummarySet);
summary_csv = fullfile(preview_dir, ...
    'Step05_FixedEta_FromFreeEtaPilot_20251222.csv');
writetable(SummaryTable, summary_csv);

fprintf('Saved fixed-eta summary:\n  %s\n', summary_csv);
end


function [Preview, PlanTable, Reference] = ...
    build_preview_from_free_summary_local(Free, cfg, sensors)
WT = Free.WindowTable;
eta_cols = strings(1, numel(sensors));
eta = nan(height(WT), numel(sensors));
for is = 1:numel(sensors)
    eta_cols(is) = sprintf('eta_CH%d_mm', sensors(is));
    eta(:, is) = WT.(eta_cols(is));
end

ok = strcmpi(string(WT.status), "ok");
all_plan = summarize_plan_local("all_ok", WT, eta, ok, sensors);

eo_ok = round(WT.EO_id(ok));
eo_ok = eo_ok(isfinite(eo_ok));
dominant_eo = NaN;
if ~isempty(eo_ok)
    dominant_eo = mode(eo_ok);
end
dominant_mask = ok & round(WT.EO_id) == dominant_eo;
dominant_plan = summarize_plan_local("dominant_eo_ok", ...
    WT, eta, dominant_mask, sensors);

PlanTable = struct2table([all_plan; dominant_plan]);

consensus_eta = all_plan.eta_mm;
consensus_iqr = all_plan.eta_iqr_mm;
if isempty(consensus_eta) || numel(consensus_eta) ~= numel(sensors)
    consensus_eta = nan(1, numel(sensors));
end
if isempty(consensus_iqr) || numel(consensus_iqr) ~= numel(sensors)
    consensus_iqr = nan(1, numel(sensors));
end
consensus_eta(1) = 0;

Consensus = struct();
Consensus.status = "ok";
Consensus.diagnostic_status = "ok";
Consensus.reference_sensor_id = sensors(1);
Consensus.eta_source_policy = "free_eta_pilot_window_median";
Consensus.primary_plan_name = "all_ok";
Consensus.eta_median_mm = consensus_eta;
Consensus.eta_plan_iqr_mm = consensus_iqr;
Consensus.eta_window_iqr_mm = consensus_iqr;
Consensus.max_abs_plan_iqr_mm = max(abs(consensus_iqr), [], 'omitnan');

Preview = struct();
Preview.DatasetDir = cfg.route_dir;
Preview.NoEtaResultFile = "";
Preview.FreeEtaResultFile = string(Free.result_file);
Preview.Step04TemplateFile = fullfile(cfg.step04_output_dir, ...
    'Template_OPRCenterStd_LowSpeed_AllBlades_S123_20251222.mat');
Preview.OutputDir = fullfile(cfg.output_root, 'step05_joint_static_eta_preview');
Preview.TargetBlade = double(Free.target_blade);
Preview.SensorIDs = sensors;
Preview.ReferenceSensorID = sensors(1);
Preview.ReportedEO = dominant_eo;
Preview.Model = ['fixed eta from per-window FreeEta pilot: ' ...
    'V = T(x - dx_c - eta_s - A*sin(EO*theta+phi))'];
Preview.SignConvention = ...
    'Step05 forward model uses x_query = x - dx_c - eta_s - A*sin(EO*theta+phi).';
Preview.Options = struct( ...
    'source', 'Compare_Step05_FreeEta_20251222', ...
    'max_windows', string(Free.max_windows), ...
    'consensus_policy', 'all_ok_window_median');
Preview.PlanTable = PlanTable;
Preview.Reference = table(sensors(:), nan(numel(sensors), 1), ...
    Free.free_eta_median_mm(:), Free.free_eta_iqr_mm(:), ...
    'VariableNames', {'sensor_id', 'step04_xc_eta_mm', ...
    'free_window_eta_median_mm', 'free_window_eta_iqr_mm'});
Preview.Consensus = Consensus;

Reference = Preview.Reference;
end


function plan = summarize_plan_local(plan_name, WT, eta, mask, sensors)
mask = mask(:) & all(isfinite(eta), 2);
eta_plan = nan(1, numel(sensors));
eta_iqr = nan(1, numel(sensors));
if any(mask)
    eta_plan = median(eta(mask, :), 1, 'omitnan');
    eta_iqr = iqr(eta(mask, :), 1);
    eta_plan(1) = 0;
end

plan = struct();
plan.plan_name = string(plan_name);
plan.status = string("ok");
if ~any(mask)
    plan.status = string("failed");
end
plan.failure_reason = string("");
if ~any(mask)
    plan.failure_reason = string("no usable FreeEta windows");
end
plan.window_count = nnz(mask);
plan.used_windows = string(mat2str(WT.window_id(mask).'));
plan.eo_sequence = string(strjoin(compose('%d', WT.EO_id(mask).'), '|'));
plan.weighted_voltage_rmse_median_v = ...
    median(WT.weighted_voltage_rmse_v(mask), 'omitnan');
plan.A_median_mm = median(WT.A_mm(mask), 'omitnan');
plan.dx_c_median_mm = median(WT.dx_c_mm(mask), 'omitnan');
plan.eta_mm = eta_plan;
plan.eta_iqr_mm = eta_iqr;
for is = 1:numel(sensors)
    plan.(sprintf('eta_CH%d_mm', sensors(is))) = eta_plan(is);
    plan.(sprintf('eta_CH%d_iqr_mm', sensors(is))) = eta_iqr(is);
end
end


function target_blades = resolve_target_blades_local(cfg)
env_text = strtrim(getenv('STEP05_TARGET_BLADES'));
if isempty(env_text)
    target_blades = 1:cfg.blades_num;
    return;
end

tokens = regexp(env_text, '[,;\s]+', 'split');
target_blades = [];
for i = 1:numel(tokens)
    if isempty(tokens{i})
        continue;
    end
    v = str2double(tokens{i});
    if isfinite(v)
        target_blades(end+1) = round(v); %#ok<AGROW>
    end
end
target_blades = unique(target_blades, 'stable');
if isempty(target_blades)
    error('STEP05_TARGET_BLADES does not contain any valid blade id.');
end
if any(target_blades < 1 | target_blades > cfg.blades_num)
    error('STEP05_TARGET_BLADES must stay inside 1..%d.', cfg.blades_num);
end
end


function row = empty_summary_row_local(sensors)
row = struct();
row.target_blade = NaN;
row.sensor_tag = string("");
row.eta_median_mm = nan(1, numel(sensors));
row.eta_iqr_mm = nan(1, numel(sensors));
row.free_eta_result_file = string("");
row.preview_mat = string("");
row.plans_csv = string("");
row.reference_csv = string("");
row.eo_sequence = string("");
row.median_rmse_v = NaN;
end


function row = make_summary_row_local( ...
    blade_id, sensors, Preview, Free, preview_mat, plans_csv, reference_csv)
row = empty_summary_row_local(sensors);
row.target_blade = blade_id;
row.sensor_tag = string(['S', sprintf('%d', sensors)]);
row.eta_median_mm = Preview.Consensus.eta_median_mm;
row.eta_iqr_mm = Preview.Consensus.eta_window_iqr_mm;
row.free_eta_result_file = string(Free.result_file);
row.preview_mat = string(preview_mat);
row.plans_csv = string(plans_csv);
row.reference_csv = string(reference_csv);
row.eo_sequence = string(Free.eo_sequence);
row.median_rmse_v = Free.median_rmse_v;
end


function env = capture_env_local()
names = {'STEP05_TARGET_BLADES'};
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
