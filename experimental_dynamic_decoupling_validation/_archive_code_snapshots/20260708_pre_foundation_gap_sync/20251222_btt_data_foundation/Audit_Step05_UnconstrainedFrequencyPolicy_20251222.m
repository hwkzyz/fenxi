%% Audit_Step05_UnconstrainedFrequencyPolicy_20251222
% Verify that the official Step05 route does not apply frequency priors or
% hard parameter constraints, and that it never reports a wrong final
% frequency for the currently audited legacy-comparable blades.

clear; clc;

this_file = mfilename('fullpath');
foundation_dir = fileparts(this_file);
repo_root = fileparts(foundation_dir);
cfg = BTTDataConfig_20251222();

diag_dir = fullfile(cfg.output_root, 'diagnostics');
if exist(diag_dir, 'dir') ~= 7
    mkdir(diag_dir);
end

step05_file = fullfile(foundation_dir, ...
    'Step05_SingleSync_DirectTemplate_Identification_20251222.m');

policy_rows = audit_step05_source_policy_local(step05_file);
result_rows = audit_frequency_results_local(repo_root, cfg);

PolicyAudit = struct2table(policy_rows);
ResultAudit = struct2table(result_rows);

overall_pass = all(strcmpi(string(PolicyAudit.status), "PASS")) && ...
    all(strcmpi(string(ResultAudit.status), "PASS"));

OverallAudit = table( ...
    "Step05 unconstrained frequency policy", ...
    overall_pass, ...
    pass_fail_text_local(overall_pass), ...
    'VariableNames', {'audit_name','pass','status'});

policy_file = fullfile(diag_dir, ...
    sprintf('Step05_UnconstrainedPolicy_Audit_%s.csv', cfg.dataset));
result_file = fullfile(diag_dir, ...
    sprintf('Step05_UnconstrainedFrequency_ResultAudit_%s.csv', cfg.dataset));
overall_file = fullfile(diag_dir, ...
    sprintf('Step05_UnconstrainedFrequency_OverallAudit_%s.csv', cfg.dataset));

writetable(PolicyAudit, policy_file);
writetable(ResultAudit, result_file);
writetable(OverallAudit, overall_file);

fprintf('\n=== Step05 unconstrained frequency audit ===\n');
disp(PolicyAudit);
disp(ResultAudit);
disp(OverallAudit);
fprintf('Saved audit files:\n  %s\n  %s\n  %s\n', ...
    policy_file, result_file, overall_file);

function rows = audit_step05_source_policy_local(step05_file)
txt = fileread(step05_file);
checks = {
    'objective_corr', "S.objective_mode = 'corr'"
    'all_eo_no_topk', 'S.top_k_eo = inf'
    'amplitude_unbounded', 'S.amplitude_limit_mm = inf'
    'dx_unbounded', 'S.dx_c_limit_mm = inf'
    'eta_unbounded', 'S.sensor_eta_limit_mm = inf'
    'eta_reg_zero', 'S.sensor_eta_reg_weight_v_per_mm = 0'
    'overshoot_penalty_zero', 'S.overshoot_penalty_weight = 0'
    'joint_gate_present', 'S.joint_eo_gap_ratio_threshold = 0.03'
    'frequency_policy_present', 'frequency_output_policy'
    };

rows = repmat(struct('check_name', "", 'status', "", 'detail', ""), size(checks, 1), 1);
for i = 1:size(checks, 1)
    check_name = string(checks{i, 1});
    pattern = checks{i, 2};
    ok = contains(txt, pattern);
    rows(i).check_name = check_name;
    rows(i).status = pass_fail_text_local(ok);
    rows(i).detail = string(pattern);
end
end

function rows = audit_frequency_results_local(repo_root, cfg)
legacy_files = containers.Map('KeyType', 'double', 'ValueType', 'char');
legacy_files(1) = fullfile(repo_root, '20251222_low_speed_rotating_calibration', ...
    'output', 'identification', ...
    'Result_Step03_Main_VPTop3SynchronousWaveform_B1_S123_Main_GradientXRange030_20251222.mat');
legacy_files(2) = fullfile(repo_root, '20251222_low_speed_rotating_calibration', ...
    'output', 'identification', ...
    'Result_Step03_Main_VPTop3SynchronousWaveform_B2_S123_Main_DirectTemplate_OPRCenterStd_Dx035_NoEta_PrevWinPhaseSafe_20251222.mat');
legacy_files(3) = fullfile(repo_root, '20251222_low_speed_rotating_calibration', ...
    'output', 'identification', ...
    'Result_Step03_Main_VPTop3SynchronousWaveform_B3_S123_Main_DirectTemplate_OPRCenterStd_Dx035_NoEta_PrevWinPhaseSafe_20251222.mat');

foundation_files = containers.Map('KeyType', 'double', 'ValueType', 'char');
foundation_files(1) = fullfile(cfg.output_root, ...
    'step05_single_sync_direct_template_official_default_unconstrained_w18', ...
    cfg.dynamic_cases{1}, ...
    'Result_Step05_SingleSyncDirectTemplate_B1_S123_20251222.mat');
foundation_files(2) = fullfile(cfg.output_root, ...
    'step05_single_sync_direct_template_b2w18', cfg.dynamic_cases{1}, ...
    'Result_Step05_SingleSyncDirectTemplate_B2_S123_20251222.mat');
foundation_files(3) = fullfile(cfg.output_root, ...
    'step05_single_sync_direct_template_b3w18', cfg.dynamic_cases{1}, ...
    'Result_Step05_SingleSyncDirectTemplate_B3_S123_20251222.mat');

freq_tolerance_hz = 0.05;
blade_ids = [1 2 3];
rows = repmat(struct('blade_id', NaN, 'legacy_eo', NaN, ...
    'legacy_mean_freq_hz', NaN, 'foundation_joint_eo', NaN, ...
    'foundation_mean_freq_hz', NaN, 'foundation_status', "", ...
    'gap_ratio', NaN, 'freq_delta_hz', NaN, 'status', "", ...
    'detail', ""), numel(blade_ids), 1);

for i = 1:numel(blade_ids)
    blade_id = blade_ids(i);
    old = load(legacy_files(blade_id), 'Result');
    new = load(foundation_files(blade_id), 'Result');

    old_trend = old.Result.Trend;
    old_freq = mean(old_trend.fn_id(isfinite(old_trend.fn_id)), 'omitnan');
    old_eo = mode(old_trend.EO_id(isfinite(old_trend.EO_id)));

    joint = new.Result.JointEOSummary;
    resonance = new.Result.ResonanceSummary;
    new_status = string(getfield_default_local(joint, 'status', ''));
    new_freq = getfield_default_local(resonance, 'mean_freq_hz', NaN);
    new_eo = getfield_default_local(joint, 'dominant_eo', NaN);
    gap_ratio = getfield_default_local(joint, 'gap_ratio', NaN);
    freq_delta = new_freq - old_freq;

    if strcmpi(new_status, "ok")
        ok = isfinite(new_freq) && abs(freq_delta) <= freq_tolerance_hz;
        detail = sprintf('ok status requires frequency within %.3f Hz of legacy.', ...
            freq_tolerance_hz);
    else
        ok = ~isfinite(new_freq);
        detail = 'ambiguous status must not report a final frequency.';
    end

    rows(i).blade_id = blade_id;
    rows(i).legacy_eo = old_eo;
    rows(i).legacy_mean_freq_hz = old_freq;
    rows(i).foundation_joint_eo = new_eo;
    rows(i).foundation_mean_freq_hz = new_freq;
    rows(i).foundation_status = new_status;
    rows(i).gap_ratio = gap_ratio;
    rows(i).freq_delta_hz = freq_delta;
    rows(i).status = string(pass_fail_text_local(ok));
    rows(i).detail = string(detail);
end
end

function txt = pass_fail_text_local(ok)
if ok
    txt = "PASS";
else
    txt = "FAIL";
end
end

function v = getfield_default_local(S, field_name, default_value)
if isstruct(S) && isfield(S, field_name) && ~isempty(S.(field_name))
    v = S.(field_name);
else
    v = default_value;
end
end
