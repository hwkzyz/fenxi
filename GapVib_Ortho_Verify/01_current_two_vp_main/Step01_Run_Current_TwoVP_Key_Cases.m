%% Current two-VP main: key-case run.
close all;
clc;

thisDir = fileparts(mfilename('fullpath'));
projectDir = fileparts(thisDir);
rootDir = fileparts(projectDir);

addpath(fullfile(projectDir, 'func'));
addpath(thisDir, '-begin');

oldDir = resolve_legacy_config_dir(rootDir);
outDir = fullfile(projectDir, 'results', '01_current_two_vp_main');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

load(fullfile(oldDir, 'stage0_config.mat'), 'cfg');
cfg.dataFile = resolve_project_data_file(rootDir, cfg.dataFile);
cfgAna = make_multicase_cfg(cfg);
cfgAna.twoVpTopKFreqCandidates = 5;
cfgAna.twoVpJointCandidateKeep = 2;
cfgAna.twoVpJointMuList = [0, 0.15, 0.75];
cfgAna.twoVpDxRefineHalfWidth = 0.05;

[gapList, xCell, yCell] = load_stacked_gap_curves(cfg.dataFile);
templateLib = build_gap_template_library(gapList, xCell, yCell, ...
    cfg.g_holdout, cfgAna.xGridN);

idxHoldout = find(abs(gapList - cfg.g_holdout) < 1e-12, 1);
if isempty(idxHoldout)
    error('Holdout gap %.6g mm not found in data file.', cfg.g_holdout);
end
FHighTrue = griddedInterpolant(xCell{idxHoldout}, yCell{idxHoldout}, ...
    'pchip', 'nearest');

caseDefs = {
    'default_10dB', struct('A_true', [0.25, 0.15], 'f_true', [500, 1300], 'phi_true', [pi/4, -pi/3], 'snrDb', 10, 'gapBiasMm', 0, 'seed', 20260901)
    'phase_sensitive_500_1200_90deg_20dB', struct('A_true', [0.25, 0.15], 'f_true', [500, 1200], 'phi_true', [0, pi/2], 'snrDb', 20, 'gapBiasMm', 0, 'seed', 20260902)
    'weakA2_003_5dB', struct('A_true', [0.25, 0.03], 'f_true', [500, 1300], 'phi_true', [pi/4, -pi/3], 'snrDb', 5, 'gapBiasMm', 0, 'seed', 20260903)
    'gapbias_m005_10dB', struct('A_true', [0.25, 0.15], 'f_true', [500, 1300], 'phi_true', [pi/4, -pi/3], 'snrDb', 10, 'gapBiasMm', -0.05, 'seed', 20260904)
    };

allRows = table();
for ic = 1:size(caseDefs, 1)
    caseName = caseDefs{ic, 1};
    caseCfg = caseDefs{ic, 2};

    cfgCase = cfgAna;
    cfgCase.noiseMode = 'noise_ratio';
    cfgCase.snrDb = caseCfg.snrDb;
    cfgCase.A_true = caseCfg.A_true;
    cfgCase.f_true = caseCfg.f_true;
    cfgCase.phi_true = caseCfg.phi_true;
    cfgCase.noiseRatio = calibrate_full_waveform_snr_ratio(cfgCase, ...
        gapList, xCell, yCell, templateLib.domain, cfgCase.snrDb);

    rng(caseCfg.seed, 'twister');
    dataHigh = simulate_rotating_waveform_from_template(FHighTrue, cfgCase.RPM_high, ...
        cfgCase.NumRevs_high, cfgCase.fs, cfgCase.R_tip, cfgCase.alpha_k, ...
        cfgCase.noiseRatio, templateLib.domain, cfgCase.A_true, ...
        cfgCase.f_true, cfgCase.phi_true, cfgCase.noiseMode);
    highMap = map_highspeed_to_space(dataHigh, cfgCase.alpha_k, cfgCase.R_tip, ...
        templateLib.domain, cfgCase.fitActiveLevel);

    staticState = estimate_highspeed_static_gap_raw(highMap, templateLib, cfgCase);
    staticState.gHat = min(max(staticState.gHat + caseCfg.gapBiasMm, ...
        min(templateLib.gapTrain) - cfgCase.rawGapSearchMargin), ...
        max(templateLib.gapTrain) + cfgCase.rawGapSearchMargin);

    truthCase = struct('g_high', cfg.g_holdout, 'A', cfgCase.A_true, ...
        'f', cfgCase.f_true, 'phi', cfgCase.phi_true);

    timerMethod = tic;
    result = run_two_vp_gap_vib_joint_method(highMap, templateLib, cfgCase, staticState);
    result.elapsed_s = toc(timerMethod);

    row = make_three_method_summary(result, truthCase);
    row.case_name = repmat(string(caseName), height(row), 1);
    row.snr_dB = repmat(caseCfg.snrDb, height(row), 1);
    row.gap_bias_mm = repmat(caseCfg.gapBiasMm, height(row), 1);
    row.true_f1_Hz = repmat(cfgCase.f_true(1), height(row), 1);
    row.true_f2_Hz = repmat(cfgCase.f_true(2), height(row), 1);
    row.true_A1_mm = repmat(cfgCase.A_true(1), height(row), 1);
    row.true_A2_mm = repmat(cfgCase.A_true(2), height(row), 1);
    row.dx_used_mm = repmat(result.dx_used, height(row), 1);
    row.static_dx0_mm = repmat(staticState.dx0, height(row), 1);
    row.dx_update_mm = repmat(result.dx_used - staticState.dx0, height(row), 1);
    allRows = [allRows; row]; %#ok<AGROW>
end

writetable(allRows, fullfile(outDir, 'new_method_key_cases.csv'));
save(fullfile(outDir, 'new_method_key_cases.mat'), 'allRows', '-v7.3');
disp(allRows);

function noiseRatio = calibrate_full_waveform_snr_ratio(cfgCase, gapList, xCell, yCell, domain, snrDb)
idxHoldout = find(abs(gapList - cfgCase.g_holdout) < 1e-12, 1);
Fx = griddedInterpolant(xCell{idxHoldout}, yCell{idxHoldout}, ...
    'pchip', 'nearest');
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
