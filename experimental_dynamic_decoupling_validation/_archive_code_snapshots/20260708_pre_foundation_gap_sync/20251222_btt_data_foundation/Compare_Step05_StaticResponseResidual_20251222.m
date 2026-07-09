function Summary = Compare_Step05_StaticResponseResidual_20251222()
%COMPARE_STEP05_STATICRESPONSERESIDUAL_20251222
% Controlled comparison for the 20251222 Step05 failure mode.
%
% This is not the formal Step05 entry.  It reuses the saved foundation
% dynamic bundles and compares three static-response interpretations:
%
%   1) formal direct non-parametric template voltage fitting;
%   2) Step04 non-parametric template inversion;
%   3) parametric static-response replay using the known-success SG Step4.
%
% The purpose is to separate OPR/dynamic-window errors from static
% calibration/model errors before changing the official Step05 route.

cfg = BTTProjectConfig_20251222();

result_file = resolve_compare_result_file_local(cfg);
old_diag_env = getenv('STEP05_DIAG_RESULT_FILE');
cleanup_obj = onCleanup(@() setenv('STEP05_DIAG_RESULT_FILE', old_diag_env)); %#ok<NASGU>
setenv('STEP05_DIAG_RESULT_FILE', result_file);

Diag = Diagnose_Step05_TemplateInversion_20251222();
W = Diag.WindowSummary;

out_dir = fullfile(cfg.output_root, 'diagnostics', ...
    'step05_static_response_residual_compare_20251222');
if exist(out_dir, 'dir') ~= 7
    mkdir(out_dir);
end

ok = strcmpi(W.status, "ok");
Compare = table();
Compare.window_id = W.window_id;
Compare.direct_template_eo = W.direct_fit_eo;
Compare.direct_template_rmse_v = W.direct_fit_rmse_v;
Compare.nonparam_inverse_common_eo = W.inv_common_eo;
Compare.nonparam_inverse_common_eo14_rank = W.inv_common_eo14_rank;
Compare.nonparam_inverse_sensor_eta_eo = W.inv_eta_eo;
Compare.nonparam_inverse_sensor_eta_eo14_rank = W.inv_eta_eo14_rank;
Compare.sg_success_static_eo = W.sg_replay_eo;
Compare.sg_success_static_rmse_v = W.sg_replay_rmse_v;
Compare.sg_success_static_eo14_rank = W.sg_replay_eo14_rank;
Compare.step04_parametric_sg_eo = W.step04_sg_eo;
Compare.step04_parametric_sg_eo14_rank = W.step04_sg_eo14_rank;
Compare.is_valid_window = ok;

summary_csv = fullfile(out_dir, ...
    'Step05_StaticResponseResidual_Compare_B1_S123_20251222.csv');
writetable(Compare, summary_csv);

Stats = struct();
Stats.source_result_file = result_file;
Stats.output_dir = out_dir;
Stats.valid_window_count = nnz(ok);
Stats.direct_template_eo_sequence = join_eo_sequence_local(W.direct_fit_eo(ok));
Stats.sg_success_static_eo_sequence = join_eo_sequence_local(W.sg_replay_eo(ok));
Stats.step04_parametric_sg_eo_sequence = join_eo_sequence_local(W.step04_sg_eo(ok));
Stats.direct_template_eo14_window_fraction = mean(W.direct_fit_eo(ok) == 14, 'omitnan');
Stats.sg_success_static_eo14_window_fraction = mean(W.sg_replay_eo(ok) == 14, 'omitnan');
Stats.step04_parametric_sg_eo14_window_fraction = mean(W.step04_sg_eo(ok) == 14, 'omitnan');
Stats.sg_success_static_eo14_median_rank = median(W.sg_replay_eo14_rank(ok), 'omitnan');
Stats.step04_parametric_sg_eo14_median_rank = median(W.step04_sg_eo14_rank(ok), 'omitnan');
Stats.nonparam_inverse_common_eo14_median_rank = median(W.inv_common_eo14_rank(ok), 'omitnan');

summary_mat = fullfile(out_dir, ...
    'Step05_StaticResponseResidual_Compare_B1_S123_20251222.mat');
Summary = struct();
Summary.Stats = Stats;
Summary.CompareTable = Compare;
Summary.Diagnostic = Diag;
save(summary_mat, 'Summary', '-v7.3');

fprintf('\n=== Step05 static-response residual comparison (20251222) ===\n');
fprintf('Source result:\n  %s\n', result_file);
fprintf('Comparison CSV:\n  %s\n', summary_csv);
fprintf('Valid windows: %d\n', Stats.valid_window_count);
fprintf('Direct-template EO sequence:      %s\n', Stats.direct_template_eo_sequence);
fprintf('SG-success static EO sequence:    %s\n', Stats.sg_success_static_eo_sequence);
fprintf('Step04-parametric SG EO sequence: %s\n', Stats.step04_parametric_sg_eo_sequence);
fprintf('EO14 window fraction: direct %.3f, SG-success %.3f, Step04-paramSG %.3f\n', ...
    Stats.direct_template_eo14_window_fraction, ...
    Stats.sg_success_static_eo14_window_fraction, ...
    Stats.step04_parametric_sg_eo14_window_fraction);
fprintf('EO14 median rank: nonparam inverse %.2f, SG-success %.2f, Step04-paramSG %.2f\n', ...
    Stats.nonparam_inverse_common_eo14_median_rank, ...
    Stats.sg_success_static_eo14_median_rank, ...
    Stats.step04_parametric_sg_eo14_median_rank);
end


function result_file = resolve_compare_result_file_local(cfg)
env_file = strtrim(getenv('STEP05_COMPARE_RESULT_FILE'));
if ~isempty(env_file)
    if ~isfile(env_file)
        error('STEP05_COMPARE_RESULT_FILE does not exist: %s', env_file);
    end
    result_file = env_file;
    return;
end

case_name = cfg.dynamic_cases{1};
preferred = {
    fullfile(cfg.output_root, ...
        'step05_single_sync_direct_template_main_noeta_sgwin50w20', ...
        case_name, ...
        'Result_Step05_FoundationMainPulseAdaptiveFixedJointEta_B1_S123_20251222.mat')
    fullfile(cfg.output_root, ...
        'step05_single_sync_direct_template', ...
        case_name, ...
        'Result_Step05_FoundationMainPulseAdaptiveFixedJointEta_B1_S123_20251222.mat')
    fullfile(cfg.output_root, ...
        'step05_single_sync_direct_template', ...
        case_name, ...
        'Result_Step05_FoundationMainPulseAdaptiveNoEta_B1_S123_20251222.mat')
    };

for i = 1:numel(preferred)
    if isfile(preferred{i})
        result_file = preferred{i};
        return;
    end
end

error(['No saved foundation Step05 result was found for comparison. ' ...
    'Run Step05 first, or set STEP05_COMPARE_RESULT_FILE.']);
end


function text = join_eo_sequence_local(eo)
if isempty(eo)
    text = "";
    return;
end
text = strjoin(compose('%d', eo(:).'), '|');
end
