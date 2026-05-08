%% Step 04: sweep truth gap values for the new method library sensitivity.

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
templateTruth = build_gap_template_library(gapList, xCell, yCell, NaN, cfgAna.xGridN);

runMode = "smoke";  % "smoke", "quick", or "full".
switch runMode
    case "full"
        gapSweepList = gapList(:).';
        snrList = [20, 10, 5];
        libraryModeList = ["ideal_in_library", "leave_one_out_inv_g_linear", ...
            "leave_one_out_g_linear", "leave_one_out_g_pchip", ...
            "sparse_leave_one_out_inv_g_linear"];
    case "quick"
        gapSweepList = intersect([0.2, 0.4, 1.0, 1.4], gapList(:).', 'stable');
        snrList = 10;
        libraryModeList = ["ideal_in_library", "leave_one_out_inv_g_linear", ...
            "sparse_leave_one_out_inv_g_linear"];
    otherwise
        gapSweepList = intersect(0.2, gapList(:).', 'stable');
        snrList = 10;
        libraryModeList = ["ideal_in_library", "leave_one_out_inv_g_linear", ...
            "leave_one_out_g_linear"];
end

csvPath = fullfile(outDir, 'new_method_gap_sweep_library_sensitivity.csv');
matPath = fullfile(outDir, 'new_method_gap_sweep_library_sensitivity.mat');
if exist(csvPath, 'file')
    delete(csvPath);
end

allRows = table();
for isweep = 1:numel(gapSweepList)
    gTrue = gapSweepList(isweep);
    ig = find(abs(gapList - gTrue) < 1e-12, 1);
    FTrue = griddedInterpolant(xCell{ig}, yCell{ig}, 'pchip', 'nearest');
    libraryDefs = make_library_defs(gapList, xCell, yCell, gTrue, ...
        cfgAna.xGridN, libraryModeList);

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

        rng(20261040 + 1000 * ig + isnr, 'twister');
        dataHigh = simulate_rotating_waveform_from_template(FTrue, cfgCase.RPM_high, ...
            cfgCase.NumRevs_high, cfgCase.fs, cfgCase.R_tip, cfgCase.alpha_k, ...
            cfgCase.noiseRatio, templateTruth.domain, cfgCase.A_true, ...
            cfgCase.f_true, cfgCase.phi_true, cfgCase.noiseMode);
        highMap = map_highspeed_to_space(dataHigh, cfgCase.alpha_k, cfgCase.R_tip, ...
            templateTruth.domain, cfgCase.fitActiveLevel);

        for ilib = 1:size(libraryDefs, 1)
            libraryMode = string(libraryDefs{ilib, 1});
            templateLib = libraryDefs{ilib, 2};
            libraryErr = analyze_template_error(templateLib, gTrue, xCell{ig}, yCell{ig});
            staticState = estimate_highspeed_static_gap_raw(highMap, templateLib, cfgCase);

            timerMethod = tic;
            result = run_two_vp_gap_vib_joint_method(highMap, templateLib, cfgCase, staticState);
            result.elapsed_s = toc(timerMethod);

            truthCase = struct('g_high', gTrue, 'A', cfgCase.A_true, ...
                'f', cfgCase.f_true, 'phi', cfgCase.phi_true);
            row = make_three_method_summary(result, truthCase);
            row.library_mode = repmat(libraryMode, height(row), 1);
            row.snr_dB = repmat(cfgCase.snrDb, height(row), 1);
            row.g_true_mm = repmat(gTrue, height(row), 1);
            row.num_gap_templates = repmat(numel(templateLib.gapTrain), height(row), 1);
            row.template_rmse = repmat(libraryErr.template_rmse, height(row), 1);
            row.derivative_rmse = repmat(libraryErr.derivative_rmse, height(row), 1);
            row.static_g0_mm = repmat(staticState.gHat, height(row), 1);
            row.static_gap_error_mm = repmat(abs(staticState.gHat - gTrue), height(row), 1);
            row.static_dx0_mm = repmat(staticState.dx0, height(row), 1);
            row.dx_used_mm = repmat(result.dx_used, height(row), 1);
            row.dx_update_mm = repmat(result.dx_used - staticState.dx0, height(row), 1);
            allRows = [allRows; row]; %#ok<AGROW>
            writetable(allRows, csvPath);
        end
    end
end

save(matPath, 'allRows', '-v7.3');
disp(allRows);

function libraryDefs = make_library_defs(gapList, xCell, yCell, gTrue, xGridN, libraryModeList)
libraryDefs = {
    'ideal_in_library', build_gap_template_library(gapList, xCell, yCell, NaN, xGridN)
    'leave_one_out_inv_g_linear', build_gap_template_library(gapList, xCell, yCell, gTrue, xGridN)
    'leave_one_out_g_linear', build_custom_gap_template_library(gapList, xCell, yCell, gTrue, xGridN, "g_linear")
    'leave_one_out_g_pchip', build_custom_gap_template_library(gapList, xCell, yCell, gTrue, xGridN, "g_pchip")
    };

sparseMask = ismember(round(gapList, 12), round([0.2; 0.6; 1.0; 1.4], 12));
sparseMask(abs(gapList - gTrue) < 1e-12) = false;
if nnz(sparseMask) >= 2
    libraryDefs(end+1, :) = {'sparse_leave_one_out_inv_g_linear', ...
        build_gap_template_library(gapList(sparseMask), xCell(sparseMask), ...
        yCell(sparseMask), NaN, xGridN)}; %#ok<AGROW>
end

keepMask = ismember(string(libraryDefs(:, 1)), libraryModeList);
libraryDefs = libraryDefs(keepMask, :);
end

function err = analyze_template_error(templateLib, gTrue, xTrue, yTrue)
xEval = templateLib.xGrid(:);
yTruth = interp1(xTrue(:), yTrue(:), xEval, 'pchip');
yPred = eval_gap_template_generic(templateLib, gTrue, xEval);
dyTruth = gradient(yTruth, mean(diff(xEval)));
dyPred = eval_gap_derivative_generic(templateLib, gTrue, xEval);
err = struct();
err.template_rmse = sqrt(mean((yPred - yTruth).^2));
err.derivative_rmse = sqrt(mean((dyPred - dyTruth).^2));
end

function templateLib = build_custom_gap_template_library(gapList, xCell, yCell, gHoldout, xGridN, gapInterpMode)
if isnan(gHoldout)
    trainIdx = true(size(gapList));
else
    trainIdx = abs(gapList - gHoldout) > 1e-12;
end
trainIds = find(trainIdx);
gTrain = gapList(trainIds);
xMin = -inf;
xMax = inf;
for ii = 1:numel(trainIds)
    ig = trainIds(ii);
    xMin = max(xMin, min(xCell{ig}));
    xMax = min(xMax, max(xCell{ig}));
end
xPad = 0.05 * (xMax - xMin);
xGrid = linspace(xMin + xPad, xMax - xPad, xGridN)';
S = zeros(numel(gTrain), numel(xGrid));
for ii = 1:numel(trainIds)
    ig = trainIds(ii);
    S(ii, :) = interp1(xCell{ig}, yCell{ig}, xGrid, 'pchip');
end
dx = mean(diff(xGrid));
FxMat = zeros(size(S));
for ii = 1:size(S, 1)
    FxMat(ii, :) = gradient(S(ii, :), dx);
end
templateLib = struct('gapTrain', gTrain(:), 'xGrid', xGrid(:), ...
    'S', S, 'FxMat', FxMat, 'domain', [min(xGrid), max(xGrid)], ...
    'gapInterpMode', string(gapInterpMode));
end

function y = eval_gap_template_generic(templateLib, g, xq)
[curveGrid, ~] = eval_gap_grid_generic(templateLib, g);
y = interp1(templateLib.xGrid, curveGrid(:), xq, 'pchip', 'extrap');
end

function y = eval_gap_derivative_generic(templateLib, g, xq)
[~, dGrid] = eval_gap_grid_generic(templateLib, g);
y = interp1(templateLib.xGrid, dGrid(:), xq, 'pchip', 'extrap');
end

function [curveGrid, dGrid] = eval_gap_grid_generic(templateLib, g)
if isfield(templateLib, 'gapInterpMode') && templateLib.gapInterpMode == "g_linear"
    [gSort, order] = sort(templateLib.gapTrain(:), 'ascend');
    curveGrid = interp1(gSort, templateLib.S(order, :), g, 'linear', 'extrap');
    dGrid = interp1(gSort, templateLib.FxMat(order, :), g, 'linear', 'extrap');
elseif isfield(templateLib, 'gapInterpMode') && templateLib.gapInterpMode == "g_pchip"
    [gSort, order] = sort(templateLib.gapTrain(:), 'ascend');
    curveGrid = interp1(gSort, templateLib.S(order, :), g, 'pchip', 'extrap');
    dGrid = interp1(gSort, templateLib.FxMat(order, :), g, 'pchip', 'extrap');
else
    [curveGrid, dGrid] = eval_gap_grid(templateLib, g);
end
end

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
