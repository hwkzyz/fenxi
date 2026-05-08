%% Fixed-trust inv_log_2: compact closed-loop comparison.
%
% This script reads the saved trusted domain and does not repeat the
% leave-one-gap window scan.
close all;
clc;

mainDir = fileparts(fileparts(mfilename('fullpath')));
addpath(mainDir, '-begin');
ctx = load_inv_log_2_project_context(mainDir);
cfgAna = ctx.cfgAna;
gapList = ctx.gapList;
xCell = ctx.xCell;
yCell = ctx.yCell;
outDir = ctx.outDir;
trustInfo = load_fixed_trust_domain(ctx.thisDir);

templateTruth = build_gap_template_library(gapList, xCell, yCell, NaN, cfgAna.xGridN);

if ~exist('runMode', 'var') || isempty(runMode)
    runMode = "quick";
end
runMode = string(runMode);

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

modelDefs = get_inv_log_2_model_suite();
allRows = table();
fileStem = make_output_stem(runMode);
csvPath = fullfile(outDir, fileStem + ".csv");
summaryCsvPath = fullfile(outDir, fileStem + "_summary.csv");
matPath = fullfile(outDir, fileStem + ".mat");

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

        rng(20261135 + 1000 * idxTrue + isnr, 'twister');
        dataHigh = simulate_rotating_waveform_from_template(FTrue, cfgCase.RPM_high, ...
            cfgCase.NumRevs_high, cfgCase.fs, cfgCase.R_tip, cfgCase.alpha_k, ...
            cfgCase.noiseRatio, templateTruth.domain, cfgCase.A_true, ...
            cfgCase.f_true, cfgCase.phi_true, cfgCase.noiseMode);

        for im = 1:height(modelDefs)
            templateLib = make_fixed_trust_template_library(gapList, xCell, yCell, ...
                gTrue, cfgAna.xGridN, modelDefs(im, :), trustInfo);
            highMap = map_highspeed_to_space(dataHigh, cfgCase.alpha_k, cfgCase.R_tip, ...
                templateLib.domain, cfgCase.fitActiveLevel);
            libraryErr = analyze_template_error(templateLib, gTrue, xCell{idxTrue}, yCell{idxTrue});
            staticState = estimate_highspeed_static_gap_raw(highMap, templateLib, cfgCase);

            timerMethod = tic;
            result = run_two_vp_gap_vib_joint_method(highMap, templateLib, cfgCase, staticState);
            result.elapsed_s = toc(timerMethod);

            truthCase = struct('g_high', gTrue, 'A', cfgCase.A_true, ...
                'f', cfgCase.f_true, 'phi', cfgCase.phi_true);
            row = make_three_method_summary(result, truthCase);
            row.model_label = repmat(modelDefs.model_label(im), height(row), 1);
            row.response_model = repmat(modelDefs.response_model(im), height(row), 1);
            row.basis_name = repmat(modelDefs.basis_name(im), height(row), 1);
            row.parameter = repmat(modelDefs.parameter(im), height(row), 1);
            row.basis_order = repmat(modelDefs.basis_order(im), height(row), 1);
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
            row.trust_x_min_mm = repmat(templateLib.domain(1), height(row), 1);
            row.trust_x_max_mm = repmat(templateLib.domain(2), height(row), 1);
            row.trust_half_width_mm = repmat(0.5 * diff(templateLib.domain), height(row), 1);
            row.trust_mode = repmat("fixed_saved_trust_domain", height(row), 1);
            allRows = [allRows; row]; %#ok<AGROW>
        end
    end
end

summary = summarize_closed_loop(allRows);
writetable(allRows, csvPath);
writetable(summary, summaryCsvPath);
save(matPath, 'allRows', 'summary', 'modelDefs', 'runMode', ...
    'gapSweepList', 'snrList', 'trustInfo', '-v7.3');
disp(summary);
fprintf('Fixed trust domain used: [%.4f, %.4f] mm\n', ...
    trustInfo.domain(1), trustInfo.domain(2));

function fileStem = make_output_stem(runMode)
if runMode == "smoke"
    fileStem = "inv_log_2_fixed_trust_closed_loop";
else
    fileStem = "inv_log_2_fixed_trust_closed_loop_" + runMode;
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
[groupId, modelLabel, responseModel, basisName, parameter, basisOrder] = findgroups( ...
    allRows.model_label, allRows.response_model, allRows.basis_name, ...
    parameterForGroup, allRows.basis_order);
parameter(parameter == realmax) = NaN;
summary = table(modelLabel, responseModel, basisName, parameter, basisOrder, ...
    splitapply(@mean, allRows.gap_error_mm, groupId), ...
    splitapply(@mean, allRows.mean_freq_error_Hz, groupId), ...
    splitapply(@mean, allRows.mean_amp_error_mm, groupId), ...
    splitapply(@mean, allRows.rmse, groupId), ...
    splitapply(@mean, allRows.template_rmse, groupId), ...
    splitapply(@mean, allRows.derivative_rmse, groupId), ...
    splitapply(@mean, allRows.trust_half_width_mm, groupId), ...
    'VariableNames', {'model_label','response_model','basis_name','parameter','basis_order', ...
    'mean_gap_error_mm','mean_freq_error_Hz','mean_amp_error_mm','mean_wave_rmse', ...
    'mean_template_rmse','mean_derivative_rmse','trust_half_width_mm'});
summary = sortrows(summary, {'mean_freq_error_Hz', 'mean_gap_error_mm'}, {'ascend', 'ascend'});
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
