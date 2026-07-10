%% Step04: pure static CH6 shift should not require residual tilt.
%
% This step separates three questions:
%   1) Does a clean pure dg6 case have a valid zero-dmu solution?
%   2) Does the independent gap_tilt optimizer find that solution?
%   3) How much do the dmu bound and noise make the optimizer drift?
%
% The warm-start and truth-start rows are diagnostics only. They are not a
% proposed paper-facing prior from gap_only to gap_tilt.

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
dmuLimitList = [0, 0.001, 0.003, 0.005, 0.010, 0.030, 0.080];

rows = table();
for iCase = 1:numel(caseList)
    caseCfg = caseList(iCase);
    bundleSyn = synthetic_gap_bias_utils_20250527('make_synthetic_voltage_bundle', ...
        bundleSet, cfg, prior, caseCfg, cfg.syntheticVoltageSeed + 4000 + iCase);
    fitBundle = synthetic_gap_bias_utils_20250527('thin_bundle', ...
        bundleSyn, cfg.maxFitPoints);

    cfgGap = cfg;
    fitGapOnly = synthetic_gap_bias_utils_20250527('fit_nested_sg_model', ...
        fitBundle, bundleSyn, "gap_only", cfgGap, prior, []);
    rows = [rows; make_result_row_local(caseCfg, sensorTag, "gap_only", ...
        "independent", NaN, prior, fitGapOnly, cfgGap, NaN)]; %#ok<AGROW>

    sensorIds = unique(fitBundle.S(:)).';
    freeSensorIds = sensorIds(2:end);
    nFree = numel(freeSensorIds);

    pGapWarm = [fitGapOnly.p(:).', zeros(1, nFree)];
    pTruth = build_truth_gap_tilt_seed_local(prior, cfg, caseCfg, freeSensorIds);

    for iLimit = 1:numel(dmuLimitList)
        cfgTilt = cfg;
        cfgTilt.fitDmuLimitMmPerMm = dmuLimitList(iLimit);

        fitInd = synthetic_gap_bias_utils_20250527('fit_nested_sg_model', ...
            fitBundle, bundleSyn, "gap_tilt", cfgTilt, prior, []);
        rows = [rows; make_result_row_local(caseCfg, sensorTag, "gap_tilt", ...
            "independent", dmuLimitList(iLimit), prior, fitInd, cfgTilt, NaN)]; %#ok<AGROW>

        fitWarm = synthetic_gap_bias_utils_20250527('fit_nested_sg_model', ...
            fitBundle, bundleSyn, "gap_tilt", cfgTilt, prior, pGapWarm);
        rows = [rows; make_result_row_local(caseCfg, sensorTag, "gap_tilt", ...
            "gap_only_zero_dmu_start", dmuLimitList(iLimit), prior, fitWarm, cfgTilt, ...
            fitGapOnly.A_mm)]; %#ok<AGROW>

        fitTruth = synthetic_gap_bias_utils_20250527('fit_nested_sg_model', ...
            fitBundle, bundleSyn, "gap_tilt", cfgTilt, prior, pTruth);
        rows = [rows; make_result_row_local(caseCfg, sensorTag, "gap_tilt", ...
            "truth_start", dmuLimitList(iLimit), prior, fitTruth, cfgTilt, NaN)]; %#ok<AGROW>
    end
end

csvPath = fullfile(cfg.outputDir, 'Step04_static_dg_identifiability_sweep.csv');
writetable(rows, csvPath);

disp('Step04 static-dg identifiability sweep:');
disp(rows(:, {'case_name','method','seed_policy','dmu_limit_mm_per_mm', ...
    'A_id_mm','amp_error_percent','fit_dg6_mm','fit_dmu6_mm_per_mm', ...
    'weighted_rmse_v'}));
fprintf('Step04 completed.\n  %s\n', csvPath);

function p0 = build_truth_gap_tilt_seed_local(prior, cfg, caseCfg, freeSensorIds)
dgFree = zeros(1, numel(freeSensorIds));
dmuFree = zeros(1, numel(freeSensorIds));
for i = 1:numel(freeSensorIds)
    if freeSensorIds(i) == 6
        dgFree(i) = caseCfg.dg6_mm;
        dmuFree(i) = caseCfg.dmu6_mm_per_mm;
    end
end
p0 = [prior.main_amp_tip_mm, cfg.truePhaseRad, cfg.trueD0Mm, dgFree, dmuFree];
end

function row = make_result_row_local(caseCfg, sensorTag, methodName, seedPolicy, ...
    dmuLimit, prior, fit, cfgFit, warmStartA)
row = table();
row.case_name = caseCfg.name;
row.sensors = sensorTag;
row.method = methodName;
row.seed_policy = seedPolicy;
row.dmu_limit_mm_per_mm = dmuLimit;
row.A_true_mm = prior.main_amp_tip_mm;
row.A_id_mm = fit.A_mm;
row.amp_error_mm = fit.A_mm - prior.main_amp_tip_mm;
row.amp_error_percent = 100 * (fit.A_mm - prior.main_amp_tip_mm) ./ ...
    max(prior.main_amp_tip_mm, eps);
row.A_minus_gap_only_warm_start_mm = fit.A_mm - warmStartA;
row.phi_id_rad = fit.phi_rad;
row.d0_id_mm = fit.d0_mm;
row.weighted_rmse_v = fit.weighted_rmse_v;
row.true_dg6_mm = caseCfg.dg6_mm;
row.true_dmu6_mm_per_mm = caseCfg.dmu6_mm_per_mm;
row.fit_dg6_mm = synthetic_gap_bias_utils_20250527('get_sensor_value', ...
    fit.sensor_ids, fit.dg_mm, 6);
row.fit_dmu6_mm_per_mm = synthetic_gap_bias_utils_20250527('get_sensor_value', ...
    fit.sensor_ids, fit.dmu_mm_per_mm, 6);
row.fit_dmu_limit_cfg = cfgFit.fitDmuLimitMmPerMm;
end
