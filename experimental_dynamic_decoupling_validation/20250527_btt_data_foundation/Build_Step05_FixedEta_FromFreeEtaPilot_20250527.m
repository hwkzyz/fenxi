function SummarySet = Build_Step05_FixedEta_FromFreeEtaPilot_20250527()
%BUILD_STEP05_FIXEDETA_FROMFREEETAPILOT_20250527
% Run per-blade FreeEta pilot and write fixed-eta preview files for Step05.

clc; close all;

this_dir = fileparts(mfilename('fullpath'));
addpath(fileparts(this_dir));
cfg = BTTProjectConfig_20250527();
target_blades = resolve_target_blades_local(cfg);
preview_dir = fullfile(cfg.output_root, 'step05_joint_static_eta_preview');
template_file = fullfile(cfg.output_root, 'step04_low_speed_template', ...
    'Template_OPRCenterStd_LowSpeed_AllBlades_S136_20250527.mat');

old_env = capture_env_local();
cleanup_obj = onCleanup(@() restore_env_local(old_env)); %#ok<NASGU>
old_dir = pwd;
dir_cleanup = onCleanup(@() cd(old_dir)); %#ok<NASGU>
cd(this_dir);

SummarySet = repmat(struct( ...
    'target_blade', NaN, ...
    'eta_median_mm', [], ...
    'eta_iqr_mm', [], ...
    'preview_mat', string("")), numel(target_blades), 1);

fprintf('\n=== Build fixed eta from FreeEta pilot (20250527) ===\n');
fprintf('Target blades: %s\n', mat2str(target_blades));

for ib = 1:numel(target_blades)
    blade_id = target_blades(ib);
    setenv('STEP05_TARGET_BLADES', num2str(blade_id));
    Free = Compare_Step05_FreeEta_20250527();
    Preview = Build_Step05_FixedEta_FromFreeEtaResult( ...
        Free.result_file, preview_dir, ...
        'DatasetDir', this_dir, ...
        'Step04TemplateFile', template_file, ...
        'SourceScript', mfilename);

    SummarySet(ib).target_blade = blade_id;
    SummarySet(ib).eta_median_mm = Preview.Consensus.eta_median_mm;
    SummarySet(ib).eta_iqr_mm = Preview.Consensus.eta_window_iqr_mm;
    SummarySet(ib).preview_mat = Preview.PreviewMatFile;
end

summary_csv = fullfile(preview_dir, ...
    'Step05_FixedEta_FromFreeEtaPilot_20250527.csv');
writetable(struct2table(SummarySet), summary_csv);
fprintf('Saved fixed-eta summary:\n  %s\n', summary_csv);
end


function target_blades = resolve_target_blades_local(cfg)
env_text = strtrim(getenv('STEP05_TARGET_BLADES'));
if isempty(env_text)
    target_blades = 1;
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
