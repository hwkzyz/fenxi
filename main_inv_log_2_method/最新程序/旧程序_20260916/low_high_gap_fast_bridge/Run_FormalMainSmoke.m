function Result = Run_FormalMainSmoke(useInformationSelection,snrDb)
%RUN_FORMALMAINSMOKE Regression test for the promoted Route-30 main entry.
if nargin<1,useInformationSelection=false;end
if nargin<2||isempty(snrDb),snrDb=20;end

root = fileparts(mfilename('fullpath'));
mainDir = fileparts(root);
addpath(mainDir, '-begin');
addpath(fullfile(mainDir, 'local_func'), '-begin');

ctx = load_inv_log_2_project_context(mainDir);
trust = load_fixed_trust_domain(mainDir);
lib = make_fixed_trust_template_library(ctx.gapList, ctx.xCell, ctx.yCell, ...
    NaN, ctx.cfgAna.xGridN, get_inv_log_2_model_def(), trust);
base = make_inv_log_2_demo_case(ctx, 0.2);
cfg = base.cfgCase;
cfg.RPM_low = min(cfg.RPM_high, 300);
cfg.NumRevs_low = 8;
cfg.snrDb = snrDb;
cfg.route30BaseSamples = 600;
cfg.route30EnableRhoFallback = false;
cfg.route30ForwardModel = "absolute";
cfg.route30UseInformationSelection = logical(useInformationSelection);

gLow = 0.8;
gHigh = 0.5;
lowData = simulate_low_speed_template(@(z) eval_gap_template(lib, gLow, z), ...
    cfg, lib.domain, cfg.snrDb, 20260806);

caseNames = ["single_531", "dual_700_1200"];
frequencyTruth = {[531], [700, 1200]};
amplitudeTruth = {[0.25], [0.25, 0.15]};
phaseTruth = {[pi/4], [pi/4, -pi/3]};
rows = repmat(empty_row(), numel(caseNames), 1);

for i = 1:numel(caseNames)
    data = simulate_rotating_waveform_from_template( ...
        @(z) eval_gap_template(lib, gHigh, z), cfg.RPM_high, ...
        cfg.NumRevs_high, cfg.fs, cfg.R_tip, cfg.alpha_k, 0, lib.domain, ...
        amplitudeTruth{i}, frequencyTruth{i}, phaseTruth{i}, 'noise_ratio');
    signalRms = rms(data.V_clean-min(data.V_clean));
    rng(20260806+i, 'twister');
    data.V_cap = data.V_clean + signalRms/10^(cfg.snrDb/20)*randn(size(data.V_clean));
    highMap = map_highspeed_to_space(data, cfg.alpha_k, cfg.R_tip, lib.domain, 0.02);

    fit = run_inv_log_2_low_high_main(lowData, highMap, lib, cfg);
    if numel(fit.VFit) ~= numel(highMap.V_a)
        error('The formal main result must predict the complete high-speed waveform.');
    end
    fEst = sort(fit.f_id(:).');
    fTrue = sort(frequencyTruth{i});
    rows(i).case_name = caseNames(i);
    rows(i).true_order = numel(fTrue);
    rows(i).estimated_order = fit.model_order;
    rows(i).g_low_hat = fit.low_speed_state.gHat;
    rows(i).g_high_hat = fit.g_used;
    rows(i).max_frequency_error_hz = max(abs(fEst-fTrue));
    rows(i).rmse_V = fit.rmse;
    rows(i).elapsed_s = fit.elapsed_s;
    rows(i).used_samples = fit.used_sample_count;
    rows(i).second_frequency_hz = fit.second_frequency_test.second_frequency_hz;
    rows(i).second_delta_bic = fit.second_frequency_test.delta_bic;
    rows(i).second_amplitude_z = fit.second_frequency_test.amplitude_z;
    rows(i).strong_second_frequency = fit.second_frequency_test.strong_second_frequency;
    if isfield(fit.information_selection,'spectral_error')
        rows(i).spectral_error=fit.information_selection.spectral_error;
    end
    rows(i).success = fit.model_order == numel(fTrue) && ...
        abs(fit.g_used-gHigh) <= 0.05 && rows(i).max_frequency_error_hz <= 2;
end

Audit = struct2table(rows);
outDir = fullfile(root, 'output', 'formal_main_smoke');
if ~exist(outDir, 'dir'), mkdir(outDir); end
writetable(Audit, fullfile(outDir, 'formal_main_smoke.csv'));
save(fullfile(outDir, 'formal_main_smoke.mat'), 'Audit');
disp(Audit);
if ~all(Audit.success)
    error('Route-30 formal main smoke test failed.');
end
Result = struct('Audit', Audit, 'outputDir', outDir);
end

function row = empty_row()
row = struct('case_name', "", 'true_order', NaN, 'estimated_order', NaN, ...
    'g_low_hat', NaN, 'g_high_hat', NaN, 'max_frequency_error_hz', NaN, ...
    'rmse_V', NaN, 'elapsed_s', NaN,'used_samples',NaN,...
    'second_frequency_hz',NaN,'second_delta_bic',NaN,...
    'second_amplitude_z',NaN,'strong_second_frequency',false,...
    'spectral_error',NaN,'success', false);
end
