%% Step 05: closed-loop test for equivalent static-template response models.
%
% This script checks whether the static-library models supported by
% 03_static_template_library_theory also improve high-speed gap-vibration
% identification. It keeps the current method unchanged and only changes the
% continuous static-template reconstruction model.

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

if ~exist('runMode', 'var') || isempty(runMode)
    runMode = "smoke";
end
runMode = string(runMode);

% "smoke" runs g=0.2 at 10 dB; "quick" adds middle gaps;
% "full" sweeps all measured gaps at [20, 10, 5] dB.
switch runMode
    case "full"
        gapSweepList = gapList(:).';
        snrList = [20, 10, 5];
    case "quick"
        gapSweepList = intersect([0.2, 0.4, 1.0, 1.4], gapList(:).', 'stable');
        snrList = 10;
    otherwise
        gapSweepList = intersect(0.2, gapList(:).', 'stable');
        snrList = 10;
end

modelDefs = make_closed_loop_model_defs();
allRows = table();
fileStem = make_output_stem(runMode);
csvPath = fullfile(outDir, fileStem + ".csv");
summaryCsvPath = fullfile(outDir, fileStem + "_summary.csv");
matPath = fullfile(outDir, fileStem + ".mat");
if exist(csvPath, 'file')
    delete(csvPath);
end

for igap = 1:numel(gapSweepList)
    gTrue = gapSweepList(igap);
    idxTrue = find(abs(gapList - gTrue) < 1e-12, 1);
    FTrue = griddedInterpolant(xCell{idxTrue}, yCell{idxTrue}, 'pchip', 'nearest');

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

        rng(20261105 + 1000 * idxTrue + isnr, 'twister');
        dataHigh = simulate_rotating_waveform_from_template(FTrue, cfgCase.RPM_high, ...
            cfgCase.NumRevs_high, cfgCase.fs, cfgCase.R_tip, cfgCase.alpha_k, ...
            cfgCase.noiseRatio, templateTruth.domain, cfgCase.A_true, ...
            cfgCase.f_true, cfgCase.phi_true, cfgCase.noiseMode);
        highMap = map_highspeed_to_space(dataHigh, cfgCase.alpha_k, cfgCase.R_tip, ...
            templateTruth.domain, cfgCase.fitActiveLevel);

        for imodel = 1:height(modelDefs)
            templateLib = build_response_template_library(gapList, xCell, yCell, ...
                gTrue, cfgAna.xGridN, modelDefs(imodel, :));
            libraryErr = analyze_template_error(templateLib, gTrue, xCell{idxTrue}, yCell{idxTrue});
            staticState = estimate_highspeed_static_gap_raw(highMap, templateLib, cfgCase);

            timerMethod = tic;
            result = run_two_vp_gap_vib_joint_method(highMap, templateLib, cfgCase, staticState);
            result.elapsed_s = toc(timerMethod);

            truthCase = struct('g_high', gTrue, 'A', cfgCase.A_true, ...
                'f', cfgCase.f_true, 'phi', cfgCase.phi_true);
            row = make_three_method_summary(result, truthCase);
            row.response_model = repmat(modelDefs.response_model(imodel), height(row), 1);
            row.model_label = repmat(modelDefs.model_label(imodel), height(row), 1);
            row.parameter = repmat(modelDefs.parameter(imodel), height(row), 1);
            row.basis_order = repmat(modelDefs.basis_order(imodel), height(row), 1);
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

summary = summarize_closed_loop(allRows);
writetable(summary, summaryCsvPath);
save(matPath, 'allRows', 'summary', 'modelDefs', 'runMode', ...
    'gapSweepList', 'snrList', '-v7.3');
disp(summary);

function fileStem = make_output_stem(runMode)
if runMode == "smoke"
    fileStem = "equivalent_response_closed_loop";
else
    fileStem = "equivalent_response_closed_loop_" + runMode;
end
end

function modelDefs = make_closed_loop_model_defs()
modelLabel = [
    "baseline_1_over_g"
    "power_law_p075_order1"
    "power_law_p025_order2"
    "field_basis_inv_log_2"
    "exponential_lambda030"
    ];
responseModel = [
    "inv_g_linear"
    "power_law"
    "power_law"
    "field_basis"
    "exponential"
    ];
parameter = [NaN; 0.75; 0.25; NaN; 0.30];
basisOrder = [1; 1; 2; 3; 1];
basisName = [
    ""
    ""
    ""
    "inv_log_2"
    ""
    ];
modelDefs = table(modelLabel, responseModel, parameter, basisOrder, basisName, ...
    'VariableNames', {'model_label','response_model','parameter','basis_order','basis_name'});
end

function templateLib = build_response_template_library(gapList, xCell, yCell, ...
    gHoldout, xGridN, modelDef)
templateLib = build_gap_template_library(gapList, xCell, yCell, gHoldout, xGridN);
switch string(modelDef.response_model)
    case "inv_g_linear"
        return;
    case "power_law"
        templateLib.gapInterpMode = "power_law";
        templateLib.gapInterpParameter = modelDef.parameter;
        templateLib.gapInterpBasisOrder = modelDef.basis_order;
    case "exponential"
        templateLib.gapInterpMode = "exponential";
        templateLib.gapInterpParameter = modelDef.parameter;
        templateLib.gapInterpBasisOrder = modelDef.basis_order;
    case "field_basis"
        templateLib.gapInterpMode = "field_basis";
        templateLib.gapInterpBasisName = string(modelDef.basis_name);
        templateLib.gapInterpBasisOrder = modelDef.basis_order;
end
end

function err = analyze_template_error(templateLib, gTrue, xTrue, yTrue)
xEval = templateLib.xGrid(:);
yTruth = interp1(xTrue(:), yTrue(:), xEval, 'pchip');
yPred = eval_gap_template(templateLib, gTrue, xEval);
dyTruth = gradient(yTruth, mean(diff(xEval)));
dyPred = eval_gap_derivative(templateLib, gTrue, xEval);
err = struct();
err.template_rmse = sqrt(mean((yPred - yTruth).^2));
err.derivative_rmse = sqrt(mean((dyPred - dyTruth).^2));
end

function summary = summarize_closed_loop(allRows)
parameterForGroup = allRows.parameter;
parameterForGroup(isnan(parameterForGroup)) = realmax;
[groupId, modelLabel, responseModel, parameter, basisOrder] = findgroups( ...
    allRows.model_label, allRows.response_model, parameterForGroup, allRows.basis_order);
parameter(parameter == realmax) = NaN;
summary = table(modelLabel, responseModel, parameter, basisOrder, ...
    splitapply(@mean, allRows.gap_error_mm, groupId), ...
    splitapply(@mean, allRows.mean_freq_error_Hz, groupId), ...
    splitapply(@mean, allRows.mean_amp_error_mm, groupId), ...
    splitapply(@mean, allRows.rmse, groupId), ...
    splitapply(@mean, allRows.template_rmse, groupId), ...
    splitapply(@mean, allRows.derivative_rmse, groupId), ...
    'VariableNames', {'model_label','response_model','parameter','basis_order', ...
    'mean_gap_error_mm','mean_freq_error_Hz','mean_amp_error_mm','mean_wave_rmse', ...
    'mean_template_rmse','mean_derivative_rmse'});
[~, order] = sort(summary.mean_freq_error_Hz + 1000 * summary.mean_gap_error_mm, 'ascend');
summary = summary(order, :);
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
