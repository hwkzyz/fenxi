%% Step05: first-order orthogonalization of residual dmu.
%
% The protected subspace is [A, phi, d0, dg_s].  A residual dmu column is
% projected onto the orthogonal complement of this subspace before fitting.
% Therefore, to first order, residual dmu can only explain waveform content
% that cannot already be represented by vibration and static gap offset.

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
bundleSet = synthetic_gap_bias_utils_20250527('filter_bundle_by_sensors', baseBundle, sensorSet);

caseList = [
    struct('name', "dg_clean", 'dg6_mm', 0.05, 'dmu6_mm_per_mm', 0.000, 'noise_scale', 0.0)
    struct('name', "dg_noise", 'dg6_mm', 0.05, 'dmu6_mm_per_mm', 0.000, 'noise_scale', 1.0)
    struct('name', "dg_dmu_noise", 'dg6_mm', 0.03, 'dmu6_mm_per_mm', 0.030, 'noise_scale', 1.0)
    ];

rows = table();
for iCase = 1:numel(caseList)
    caseCfg = caseList(iCase);
    bundleSyn = synthetic_gap_bias_utils_20250527('make_synthetic_voltage_bundle', ...
        bundleSet, cfg, prior, caseCfg, cfg.syntheticVoltageSeed + 5000 + iCase);
    fitBundle = synthetic_gap_bias_utils_20250527('thin_bundle', ...
        bundleSyn, cfg.maxFitPoints);

    D = synthetic_gap_bias_utils_20250527('build_jacobian_columns', ...
        fitBundle, cfg, prior, []);
    y = fitBundle.V(:) - D.VBase(:);
    weight = D.weight(:);

    sensorIds = D.sensorIds;
    freeIdx = 2:numel(sensorIds);
    protected = [D.J_A, D.J_phi, D.J_d0, D.J_dg(:, freeIdx)];
    rawDmu = D.J_dmu(:, freeIdx);
    orthDmu = synthetic_gap_bias_utils_20250527('residualize_columns', ...
        rawDmu, protected, weight);

    [betaGap, detailGap] = synthetic_gap_bias_utils_20250527('weighted_linear_fit', ...
        protected, y, weight);
    rows = [rows; make_linear_row_local(caseCfg, sensorTag, "linear_gap_only", ...
        "none", prior, sensorIds, freeIdx, betaGap, detailGap, protected, ...
        rawDmu, orthDmu, weight)]; %#ok<AGROW>

    [betaRaw, detailRaw] = synthetic_gap_bias_utils_20250527('weighted_linear_fit', ...
        [protected, rawDmu], y, weight);
    rows = [rows; make_linear_row_local(caseCfg, sensorTag, "linear_gap_tilt_raw", ...
        "physical_dmu", prior, sensorIds, freeIdx, betaRaw, detailRaw, protected, ...
        rawDmu, orthDmu, weight)]; %#ok<AGROW>

    [betaOrth, detailOrth] = synthetic_gap_bias_utils_20250527('weighted_linear_fit', ...
        [protected, orthDmu], y, weight);
    rows = [rows; make_linear_row_local(caseCfg, sensorTag, "linear_gap_tilt_orth", ...
        "orthogonal_residual_dmu", prior, sensorIds, freeIdx, betaOrth, detailOrth, ...
        protected, rawDmu, orthDmu, weight)]; %#ok<AGROW>
end

csvPath = fullfile(cfg.outputDir, 'Step05_orthogonalized_residual_tilt.csv');
writetable(rows, csvPath);

disp('Step05 orthogonalized residual tilt:');
disp(rows(:, {'case_name','fit_name','dmu_column_meaning','A_est_mm', ...
    'amp_error_percent','fit_dg6_mm','fit_dmu6_mm_per_mm', ...
    'weighted_rmse_v','proj_dmu6_on_protected','dmu6_perp_norm_fraction'}));
fprintf('Step05 completed.\n  %s\n', csvPath);

function row = make_linear_row_local(caseCfg, sensorTag, fitName, dmuMeaning, ...
    prior, sensorIds, freeIdx, beta, detail, protected, rawDmu, orthDmu, weight)
row = table();
row.case_name = caseCfg.name;
row.sensors = sensorTag;
row.fit_name = fitName;
row.dmu_column_meaning = dmuMeaning;
row.A_true_mm = prior.main_amp_tip_mm;
row.A_est_mm = prior.main_amp_tip_mm + beta(1);
row.delta_A_mm = beta(1);
row.amp_error_percent = 100 * beta(1) ./ max(prior.main_amp_tip_mm, eps);
row.delta_phi_rad = beta(2);
row.delta_d0_mm = beta(3);
row.true_dg6_mm = caseCfg.dg6_mm;
row.true_dmu6_mm_per_mm = caseCfg.dmu6_mm_per_mm;
row.weighted_rmse_v = detail.weightedRmse;
row.linear_rank = detail.rank;
row.linear_cond = detail.cond;

nProtected = size(protected, 2);
betaDg = nan(1, numel(sensorIds));
betaDmu = nan(1, numel(sensorIds));
for i = 1:numel(freeIdx)
    sensorId = sensorIds(freeIdx(i));
    betaDg(freeIdx(i)) = beta(3 + i);
    idxDmu = nProtected + i;
    if numel(beta) >= idxDmu
        betaDmu(freeIdx(i)) = beta(idxDmu);
    end
    if sensorId == 6
        row.fit_dg6_mm = betaDg(freeIdx(i));
        row.fit_dmu6_mm_per_mm = betaDmu(freeIdx(i));
        row.proj_dmu6_on_protected = synthetic_gap_bias_utils_20250527( ...
            'projection_ratio', rawDmu(:, i), protected, weight);
        row.dmu6_perp_norm_fraction = synthetic_gap_bias_utils_20250527( ...
            'weighted_norm', orthDmu(:, i), weight) ./ ...
            max(synthetic_gap_bias_utils_20250527('weighted_norm', rawDmu(:, i), weight), eps);
    end
end

if ~ismember('fit_dg6_mm', string(row.Properties.VariableNames))
    row.fit_dg6_mm = NaN;
    row.fit_dmu6_mm_per_mm = NaN;
    row.proj_dmu6_on_protected = NaN;
    row.dmu6_perp_norm_fraction = NaN;
end
end
