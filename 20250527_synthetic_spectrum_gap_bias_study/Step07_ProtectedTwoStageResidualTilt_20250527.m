%% Step07: protected two-stage residual tilt diagnostic.
%
% Stage 1 is the paper-facing model:
%   q = [A, phi, d0, dg_s], dmu_s = 0.
%
% Stage 2 is diagnostic only:
%   fit residual dmu_s on the voltage residual after projecting J_mu away
%   from the protected q-subspace. The reported A remains the Stage-1 A.

cfg = study_config_20250527();
priorPath = fullfile(cfg.outputDir, 'Step01_synthetic_strain_prior.mat');
if ~isfile(priorPath)
    error('Run Step01 first. Missing file: %s', priorPath);
end
if ~isfile(cfg.resultFileS136)
    error('Missing Step5 S136 result: %s', cfg.resultFileS136);
end

P = load(priorPath, 'prior');
prior = P.prior;
S = load(cfg.resultFileS136, 'Result_Struct');
baseBundle = S.Result_Struct.BestWindow.bundle;

sensorSet = [1 3 6];
sensorTag = "S136";
bundleSet = synthetic_gap_bias_utils_20250527('filter_bundle_by_sensors', ...
    baseBundle, sensorSet);

caseNames = ["noise_only", "ch6_static_dg", "ch6_tilt_dmu", "ch6_dg_dmu"];
rows = table();

for iCase = 1:numel(caseNames)
    caseName = caseNames(iCase);
    caseCfg = cfg.ch6ShiftCases.(caseName);
    bundleSyn = synthetic_gap_bias_utils_20250527('make_synthetic_voltage_bundle', ...
        bundleSet, cfg, prior, caseCfg, cfg.syntheticVoltageSeed + 7000 + iCase);
    fitBundle = synthetic_gap_bias_utils_20250527('thin_bundle', ...
        bundleSyn, cfg.maxFitPoints);

    fitGapOnly = synthetic_gap_bias_utils_20250527('fit_nested_sg_model', ...
        fitBundle, fitBundle, "gap_only", cfg, prior, []);
    fitGapTilt = synthetic_gap_bias_utils_20250527('fit_nested_sg_model', ...
        fitBundle, fitBundle, "gap_tilt", cfg, prior, []);

    rows = [rows; make_nonlinear_row_local(caseName, sensorTag, ...
        "stage1_gap_fixed_tilt", "main", prior, fitGapOnly, fitGapOnly, NaN)]; %#ok<AGROW>
    rows = [rows; make_nonlinear_row_local(caseName, sensorTag, ...
        "current_gap_tilt_independent", "diagnostic_unprotected", prior, ...
        fitGapTilt, fitGapOnly, synthetic_gap_bias_utils_20250527( ...
        'get_sensor_value', fitGapTilt.sensor_ids, fitGapTilt.dmu_mm_per_mm, 6))]; %#ok<AGROW>

    baseState = struct();
    baseState.A_mm = fitGapOnly.A_mm;
    baseState.phi_rad = fitGapOnly.phi_rad;
    baseState.d0_mm = fitGapOnly.d0_mm;
    baseState.dg_mm = fitGapOnly.dg_mm;
    baseState.dmu_mm_per_mm = fitGapOnly.dmu_mm_per_mm;
    D = synthetic_gap_bias_utils_20250527('build_jacobian_columns', ...
        fitBundle, cfg, prior, baseState);

    residual = fitBundle.V(:) - D.VBase(:);
    weight = D.weight(:);
    sensorIds = D.sensorIds;
    freeIdx = 2:numel(sensorIds);
    protected = [D.J_A, D.J_phi, D.J_d0, D.J_dg(:, freeIdx)];
    rawDmu = D.J_dmu(:, freeIdx);
    gapOrthDmu = rawDmu;
    for j = 1:numel(freeIdx)
        Jg = D.J_dg(:, freeIdx(j));
        alpha = weighted_inner_local(Jg, rawDmu(:, j), weight) ./ ...
            max(weighted_inner_local(Jg, Jg, weight), eps);
        gapOrthDmu(:, j) = rawDmu(:, j) - alpha .* Jg;
    end
    safeDmu = synthetic_gap_bias_utils_20250527('residualize_columns', ...
        rawDmu, protected, weight);

    rows = [rows; make_post_dmu_row_local(caseName, sensorTag, ...
        "stage2_raw_dmu_A_locked", "raw", prior, fitGapOnly, residual, ...
        rawDmu, weight, sensorIds, freeIdx)]; %#ok<AGROW>
    rows = [rows; make_post_dmu_row_local(caseName, sensorTag, ...
        "stage2_gaporth_dmu_A_locked", "gap_orthogonal", prior, fitGapOnly, ...
        residual, gapOrthDmu, weight, sensorIds, freeIdx)]; %#ok<AGROW>
    rows = [rows; make_post_dmu_row_local(caseName, sensorTag, ...
        "stage2_safe_dmu_A_locked", "vib_gap_orthogonal", prior, fitGapOnly, ...
        residual, safeDmu, weight, sensorIds, freeIdx)]; %#ok<AGROW>

    rows = [rows; make_linear_update_row_local(caseName, sensorTag, ...
        "linear_raw_refit_can_move_A", "raw", prior, fitGapOnly, residual, ...
        protected, rawDmu, weight, sensorIds, freeIdx)]; %#ok<AGROW>
    rows = [rows; make_linear_update_row_local(caseName, sensorTag, ...
        "linear_safe_refit_can_move_A", "vib_gap_orthogonal", prior, ...
        fitGapOnly, residual, protected, safeDmu, weight, sensorIds, freeIdx)]; %#ok<AGROW>
end

csvPath = fullfile(cfg.outputDir, 'Step07_protected_two_stage_residual_tilt.csv');
writetable(rows, csvPath);

disp('Step07 protected two-stage residual tilt:');
disp(rows(:, {'case_name','model_name','A_policy','A_reported_mm', ...
    'amp_error_percent','fit_dg6_mm','fit_dmu6_mm_per_mm', ...
    'weighted_rmse_v','rmse_improvement_percent'}));
fprintf('Step07 completed.\n  %s\n', csvPath);

function row = make_nonlinear_row_local(caseName, sensorTag, modelName, policy, ...
    prior, fit, mainFit, dmu6)
row = base_row_local(caseName, sensorTag, modelName, policy, prior);
row.A_reported_mm = fit.A_mm;
row.amp_error_mm = fit.A_mm - prior.main_amp_tip_mm;
row.amp_error_percent = 100 * row.amp_error_mm ./ max(prior.main_amp_tip_mm, eps);
row.fit_dg6_mm = synthetic_gap_bias_utils_20250527('get_sensor_value', ...
    fit.sensor_ids, fit.dg_mm, 6);
row.fit_dmu6_mm_per_mm = dmu6;
row.weighted_rmse_v = fit.weighted_rmse_v;
row.baseline_gap_only_rmse_v = mainFit.weighted_rmse_v;
row.rmse_improvement_percent = 100 * (mainFit.weighted_rmse_v - fit.weighted_rmse_v) ./ ...
    max(mainFit.weighted_rmse_v, eps);
row.linear_delta_A_mm = NaN;
end

function row = make_post_dmu_row_local(caseName, sensorTag, modelName, dmuKind, ...
    prior, mainFit, residual, dmuBasis, weight, sensorIds, freeIdx)
[beta, detail] = synthetic_gap_bias_utils_20250527('weighted_linear_fit', ...
    dmuBasis, residual, weight);
row = base_row_local(caseName, sensorTag, modelName, "A_locked_" + dmuKind, prior);
row.A_reported_mm = mainFit.A_mm;
row.amp_error_mm = mainFit.A_mm - prior.main_amp_tip_mm;
row.amp_error_percent = 100 * row.amp_error_mm ./ max(prior.main_amp_tip_mm, eps);
row.fit_dg6_mm = synthetic_gap_bias_utils_20250527('get_sensor_value', ...
    mainFit.sensor_ids, mainFit.dg_mm, 6);
row.fit_dmu6_mm_per_mm = get_free_sensor_beta_local(sensorIds, freeIdx, beta, 6);
row.weighted_rmse_v = detail.weightedRmse;
row.baseline_gap_only_rmse_v = mainFit.weighted_rmse_v;
row.rmse_improvement_percent = 100 * (mainFit.weighted_rmse_v - detail.weightedRmse) ./ ...
    max(mainFit.weighted_rmse_v, eps);
row.linear_delta_A_mm = 0;
end

function row = make_linear_update_row_local(caseName, sensorTag, modelName, dmuKind, ...
    prior, mainFit, residual, protected, dmuBasis, weight, sensorIds, freeIdx)
X = [protected, dmuBasis];
[beta, detail] = synthetic_gap_bias_utils_20250527('weighted_linear_fit', X, residual, weight);
nProtected = size(protected, 2);
deltaA = beta(1);
dg6Delta = get_free_sensor_beta_local(sensorIds, freeIdx, beta(4:(3 + numel(freeIdx))), 6);
dmu6 = get_free_sensor_beta_local(sensorIds, freeIdx, beta((nProtected + 1):end), 6);

row = base_row_local(caseName, sensorTag, modelName, "linear_update_" + dmuKind, prior);
row.A_reported_mm = mainFit.A_mm + deltaA;
row.amp_error_mm = row.A_reported_mm - prior.main_amp_tip_mm;
row.amp_error_percent = 100 * row.amp_error_mm ./ max(prior.main_amp_tip_mm, eps);
row.fit_dg6_mm = synthetic_gap_bias_utils_20250527('get_sensor_value', ...
    mainFit.sensor_ids, mainFit.dg_mm, 6) + dg6Delta;
row.fit_dmu6_mm_per_mm = dmu6;
row.weighted_rmse_v = detail.weightedRmse;
row.baseline_gap_only_rmse_v = mainFit.weighted_rmse_v;
row.rmse_improvement_percent = 100 * (mainFit.weighted_rmse_v - detail.weightedRmse) ./ ...
    max(mainFit.weighted_rmse_v, eps);
row.linear_delta_A_mm = deltaA;
end

function row = base_row_local(caseName, sensorTag, modelName, policy, prior)
row = table();
row.case_name = string(caseName);
row.sensors = sensorTag;
row.model_name = modelName;
row.A_policy = policy;
row.A_true_mm = prior.main_amp_tip_mm;
end

function val = get_free_sensor_beta_local(sensorIds, freeIdx, beta, sensorId)
val = NaN;
for i = 1:numel(freeIdx)
    if sensorIds(freeIdx(i)) == sensorId && numel(beta) >= i
        val = beta(i);
        return;
    end
end
end

function val = weighted_inner_local(a, b, weight)
a = a(:);
b = b(:);
weight = max(weight(:), 0);
valid = isfinite(a) & isfinite(b) & isfinite(weight);
val = sum(weight(valid) .* a(valid) .* b(valid));
end
