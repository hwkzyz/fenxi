%% Fixed-trust inv_log_2: representative key cases.
close all;
clc;

mainDir = fileparts(fileparts(mfilename('fullpath')));
addpath(mainDir, '-begin');
ctx = load_inv_log_2_project_context(mainDir);
cfg = ctx.cfg;
cfgAna = ctx.cfgAna;
gapList = ctx.gapList;
xCell = ctx.xCell;
yCell = ctx.yCell;
outDir = ctx.outDir;
trustInfo = load_fixed_trust_domain(ctx.thisDir);

idxHoldout = find(abs(gapList - cfg.g_holdout) < 1e-12, 1);
FHighTrue = griddedInterpolant(xCell{idxHoldout}, yCell{idxHoldout}, ...
    'pchip', 'nearest');

templateFull = build_gap_template_library(gapList, xCell, yCell, cfg.g_holdout, cfgAna.xGridN);
modelDefs = [get_baseline_model_def(); get_inv_log_2_model_def()];
caseDefs = make_key_case_defs();
allRows = table();

for ic = 1:numel(caseDefs)
    caseCfg = caseDefs(ic);
    cfgCase = cfgAna;
    cfgCase.noiseMode = 'noise_ratio';
    cfgCase.snrDb = caseCfg.snrDb;
    cfgCase.A_true = caseCfg.A_true;
    cfgCase.f_true = caseCfg.f_true;
    cfgCase.phi_true = caseCfg.phi_true;
    cfgCase.noiseRatio = calibrate_full_waveform_snr_ratio(cfgCase, ...
        gapList, xCell, yCell, templateFull.domain, cfgCase.snrDb);

    rng(caseCfg.seed, 'twister');
    dataHigh = simulate_rotating_waveform_from_template(FHighTrue, cfgCase.RPM_high, ...
        cfgCase.NumRevs_high, cfgCase.fs, cfgCase.R_tip, cfgCase.alpha_k, ...
        cfgCase.noiseRatio, templateFull.domain, cfgCase.A_true, ...
        cfgCase.f_true, cfgCase.phi_true, cfgCase.noiseMode);

    truthCase = struct('g_high', cfg.g_holdout, 'A', cfgCase.A_true, ...
        'f', cfgCase.f_true, 'phi', cfgCase.phi_true);

    for im = 1:height(modelDefs)
        templateLib = make_fixed_trust_template_library(gapList, xCell, yCell, ...
            cfg.g_holdout, cfgAna.xGridN, modelDefs(im, :), trustInfo);
        highMap = map_highspeed_to_space(dataHigh, cfgCase.alpha_k, cfgCase.R_tip, ...
            templateLib.domain, cfgCase.fitActiveLevel);
        staticState = estimate_highspeed_static_gap_raw(highMap, templateLib, cfgCase);
        staticState.gHat = min(max(staticState.gHat + caseCfg.gapBiasMm, ...
            min(templateLib.gapTrain) - cfgCase.rawGapSearchMargin), ...
            max(templateLib.gapTrain) + cfgCase.rawGapSearchMargin);

        timerMethod = tic;
        result = run_two_vp_gap_vib_joint_method(highMap, templateLib, cfgCase, staticState);
        result.elapsed_s = toc(timerMethod);

        row = make_three_method_summary(result, truthCase);
        row.case_name = repmat(caseCfg.name, height(row), 1);
        row.model_label = repmat(modelDefs.model_label(im), height(row), 1);
        row.response_model = repmat(modelDefs.response_model(im), height(row), 1);
        row.basis_name = repmat(modelDefs.basis_name(im), height(row), 1);
        row.parameter = repmat(modelDefs.parameter(im), height(row), 1);
        row.basis_order = repmat(modelDefs.basis_order(im), height(row), 1);
        row.snr_dB = repmat(cfgCase.snrDb, height(row), 1);
        row.gap_bias_mm = repmat(caseCfg.gapBiasMm, height(row), 1);
        row.true_f1_Hz = repmat(cfgCase.f_true(1), height(row), 1);
        row.true_f2_Hz = repmat(cfgCase.f_true(2), height(row), 1);
        row.true_A1_mm = repmat(cfgCase.A_true(1), height(row), 1);
        row.true_A2_mm = repmat(cfgCase.A_true(2), height(row), 1);
        row.dx_used_mm = repmat(result.dx_used, height(row), 1);
        row.static_dx0_mm = repmat(staticState.dx0, height(row), 1);
        row.dx_update_mm = repmat(result.dx_used - staticState.dx0, height(row), 1);
        row.trust_x_min_mm = repmat(templateLib.domain(1), height(row), 1);
        row.trust_x_max_mm = repmat(templateLib.domain(2), height(row), 1);
        row.trust_mode = repmat("fixed_saved_trust_domain", height(row), 1);
        allRows = [allRows; row]; %#ok<AGROW>
    end
end

writetable(allRows, fullfile(outDir, 'inv_log_2_fixed_trust_key_cases.csv'));
save(fullfile(outDir, 'inv_log_2_fixed_trust_key_cases.mat'), ...
    'allRows', 'modelDefs', 'caseDefs', 'trustInfo', '-v7.3');
disp(allRows);

function caseDefs = make_key_case_defs()
caseDefs = [
    build_case("default_10dB", [0.25, 0.15], [500, 1300], [pi/4, -pi/3], 10, 0, 20260911)
    build_case("phase_sensitive_20dB", [0.25, 0.15], [500, 1200], [0, pi/2], 20, 0, 20260912)
    build_case("weak_mode2_5dB", [0.25, 0.03], [500, 1300], [pi/4, -pi/3], 5, 0, 20260913)
    build_case("gap_bias_minus005_10dB", [0.25, 0.15], [500, 1300], [pi/4, -pi/3], 10, -0.05, 20260914)
    ];
end

function caseDef = build_case(name, ATrue, fTrue, phiTrue, snrDb, gapBiasMm, seed)
caseDef = struct('name', string(name), 'A_true', ATrue, 'f_true', fTrue, ...
    'phi_true', phiTrue, 'snrDb', snrDb, 'gapBiasMm', gapBiasMm, 'seed', seed);
end

function noiseRatio = calibrate_full_waveform_snr_ratio(cfgCase, gapList, xCell, yCell, domain, snrDb)
idxHoldout = find(abs(gapList - cfgCase.g_holdout) < 1e-12, 1);
Fx = griddedInterpolant(xCell{idxHoldout}, yCell{idxHoldout}, 'pchip', 'nearest');
dataClean = simulate_rotating_waveform_from_template(Fx, cfgCase.RPM_high, ...
    cfgCase.NumRevs_high, cfgCase.fs, cfgCase.R_tip, cfgCase.alpha_k, ...
    0, domain, cfgCase.A_true, cfgCase.f_true, cfgCase.phi_true, 'noise_ratio');
signal = dataClean.V_clean(:) - min(dataClean.V_clean(:));
signalRms = sqrt(mean(signal.^2));
probe = interp1(xCell{idxHoldout}, yCell{idxHoldout}, ...
    linspace(domain(1), domain(2), 400)', 'pchip');
templateRange = range(probe);
noiseRatio = signalRms ./ (templateRange .* 10^(snrDb / 20));
end
