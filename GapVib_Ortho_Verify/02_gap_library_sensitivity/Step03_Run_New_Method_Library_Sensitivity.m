%% Step 03: run the new method under different gap-library definitions.

close all;
clc;

caseDir = fileparts(mfilename('fullpath'));
projectDir = fileparts(caseDir);
rootDir = fileparts(projectDir);
addpath(fullfile(projectDir, 'func'));
addpath(fullfile(projectDir, '01_current_two_vp_main'), '-begin');

oldDir = resolve_legacy_config_dir(rootDir);
outDir = fullfile(projectDir, 'results', '02_gap_library_sensitivity');
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
if isfield(cfg, 'g_holdout') && any(abs(gapList - cfg.g_holdout) < 1e-12)
    gTrue = cfg.g_holdout;
else
    gTrue = 1.0;
end
idxTrue = find(abs(gapList - gTrue) < 1e-12, 1);
if isempty(idxTrue)
    error('Truth gap %.6g mm not found in the gap library.', gTrue);
end

FTrue = griddedInterpolant(xCell{idxTrue}, yCell{idxTrue}, 'pchip', 'nearest');
templateTruth = build_gap_template_library(gapList, xCell, yCell, NaN, cfgAna.xGridN);

libraryDefs = {
    'ideal_in_library', build_gap_template_library(gapList, xCell, yCell, NaN, cfgAna.xGridN)
    'leave_one_out_inv_g_linear', build_gap_template_library(gapList, xCell, yCell, gTrue, cfgAna.xGridN)
    };

sparseMask = ismember(round(gapList, 12), round([0.2; 0.6; 1.0; 1.4], 12));
sparseMask(idxTrue) = false;
if nnz(sparseMask) >= 2
    libraryDefs(end+1, :) = {'sparse_leave_one_out_inv_g_linear', ...
        build_gap_template_library(gapList(sparseMask), xCell(sparseMask), ...
        yCell(sparseMask), NaN, cfgAna.xGridN)}; %#ok<SAGROW>
end

snrList = [20, 10, 5];
allRows = table();
for isnr = 1:numel(snrList)
    cfgCase = cfgAna;
    cfgCase.g_holdout = gTrue;
    cfgCase.A_true = [0.25, 0.15];
    cfgCase.f_true = [500, 1300];
    cfgCase.phi_true = [pi/4, -pi/3];
    cfgCase.noiseMode = 'noise_ratio';
    cfgCase.snrDb = snrList(isnr);
    cfgCase.noiseRatio = calibrate_full_waveform_snr_ratio(cfgCase, ...
        gapList, xCell, yCell, templateTruth.domain, cfgCase.snrDb);

    rng(20261020 + isnr, 'twister');
    dataHigh = simulate_rotating_waveform_from_template(FTrue, cfgCase.RPM_high, ...
        cfgCase.NumRevs_high, cfgCase.fs, cfgCase.R_tip, cfgCase.alpha_k, ...
        cfgCase.noiseRatio, templateTruth.domain, cfgCase.A_true, ...
        cfgCase.f_true, cfgCase.phi_true, cfgCase.noiseMode);
    highMap = map_highspeed_to_space(dataHigh, cfgCase.alpha_k, cfgCase.R_tip, ...
        templateTruth.domain, cfgCase.fitActiveLevel);

    for ilib = 1:size(libraryDefs, 1)
        libraryMode = string(libraryDefs{ilib, 1});
        templateLib = libraryDefs{ilib, 2};
        staticState = estimate_highspeed_static_gap_raw(highMap, templateLib, cfgCase);

        timerMethod = tic;
        result = run_two_vp_gap_vib_joint_method(highMap, templateLib, cfgCase, staticState);
        result.elapsed_s = toc(timerMethod);

        truthCase = struct('g_high', gTrue, 'A', cfgCase.A_true, ...
            'f', cfgCase.f_true, 'phi', cfgCase.phi_true);
        row = make_three_method_summary(result, truthCase);
        row.library_mode = repmat(libraryMode, height(row), 1);
        row.snr_dB = repmat(cfgCase.snrDb, height(row), 1);
        row.num_gap_templates = repmat(numel(templateLib.gapTrain), height(row), 1);
        row.static_g0_mm = repmat(staticState.gHat, height(row), 1);
        row.static_dx0_mm = repmat(staticState.dx0, height(row), 1);
        row.dx_used_mm = repmat(result.dx_used, height(row), 1);
        row.dx_update_mm = repmat(result.dx_used - staticState.dx0, height(row), 1);
        allRows = [allRows; row]; %#ok<AGROW>
    end
end

writetable(allRows, fullfile(outDir, 'new_method_library_sensitivity.csv'));
save(fullfile(outDir, 'new_method_library_sensitivity.mat'), 'allRows', '-v7.3');
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
