%% Step 04: compare field-inspired static-template basis models.
%
% This script keeps the existing main method untouched. It tests whether
% field-mechanism-inspired basis functions such as [1, 1/g, log(g), 1/g^2]
% can match or improve the current best calibrated power-law model.

close all;
clc;

caseDir = fileparts(mfilename('fullpath'));
projectDir = fileparts(caseDir);
rootDir = fileparts(projectDir);
addpath(fullfile(projectDir, 'func'));

oldDir = resolve_legacy_config_dir(rootDir);
outDir = fullfile(projectDir, 'results', '03_static_template_library_theory');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

load(fullfile(oldDir, 'stage0_config.mat'), 'cfg');
cfg.dataFile = resolve_project_data_file(rootDir, cfg.dataFile);
cfgAna = make_multicase_cfg(cfg);
[gapList, xCell, yCell] = load_stacked_gap_curves(cfg.dataFile);

modelDefs = make_field_model_defs();
rows = table();
for im = 1:height(modelDefs)
    for ih = 1:numel(gapList)
        rows = [rows; evaluate_one_holdout(gapList, xCell, yCell, gapList(ih), ...
            cfgAna.xGridN, modelDefs(im, :))]; %#ok<AGROW>
    end
end

summary = summarize_field_model_errors(rows);
writetable(rows, fullfile(outDir, 'field_basis_model_leave_one_out.csv'));
writetable(summary, fullfile(outDir, 'field_basis_model_summary.csv'));
save(fullfile(outDir, 'field_basis_model_analysis.mat'), ...
    'rows', 'summary', 'modelDefs', '-v7.3');

disp(summary);

function modelDefs = make_field_model_defs()
modelLabel = [
    "baseline_1_over_g"
    "power_law_p025_order2"
    "cap_edge_1"
    "cap_edge_2"
    "cap_edge_3"
    "inv_poly_3"
    "inv_log_2"
    ];
modelType = [
    "inv_g_linear"
    "power_law"
    "field_basis"
    "field_basis"
    "field_basis"
    "field_basis"
    "field_basis"
    ];
basisName = [
    ""
    ""
    "cap_edge_1"
    "cap_edge_2"
    "cap_edge_3"
    "inv_poly_3"
    "inv_log_2"
    ];
parameter = [NaN; 0.25; NaN; NaN; NaN; NaN; NaN];
basisOrder = [1; 2; 4; 4; 4; 4; 3];
modelDefs = table(modelLabel, modelType, basisName, parameter, basisOrder, ...
    'VariableNames', {'model_label','model_type','basis_name','parameter','basis_order'});
end

function row = evaluate_one_holdout(gapList, xCell, yCell, gHoldout, xGridN, modelDef)
templateLib = build_gap_template_library(gapList, xCell, yCell, gHoldout, xGridN);
templateLib = configure_template_model(templateLib, modelDef);

idxHoldout = find(abs(gapList - gHoldout) < 1e-12, 1);
xGrid = templateLib.xGrid(:);
yTruth = interp1(xCell{idxHoldout}(:), yCell{idxHoldout}(:), xGrid, 'pchip');
dyTruth = gradient(yTruth, mean(diff(xGrid)));

yPred = eval_gap_template(templateLib, gHoldout, xGrid);
dyPred = eval_gap_derivative(templateLib, gHoldout, xGrid);

res = yPred(:) - yTruth(:);
dres = dyPred(:) - dyTruth(:);
if gHoldout < min(templateLib.gapTrain) || gHoldout > max(templateLib.gapTrain)
    reconstructionType = "boundary_extrapolation";
else
    reconstructionType = "interpolation";
end

row = table(string(modelDef.model_label), string(modelDef.model_type), ...
    string(modelDef.basis_name), modelDef.parameter, modelDef.basis_order, ...
    gHoldout, string(reconstructionType), numel(templateLib.gapTrain), ...
    sqrt(mean(res.^2)), mean(abs(res)), max(abs(res)), ...
    sqrt(mean(dres.^2)), mean(abs(dres)), max(abs(dres)), ...
    'VariableNames', {'model_label','model_type','basis_name','parameter', ...
    'basis_order','g_holdout_mm','reconstruction_type','num_gap_templates', ...
    'template_rmse','template_mean_abs','template_max_abs', ...
    'derivative_rmse','derivative_mean_abs','derivative_max_abs'});
end

function templateLib = configure_template_model(templateLib, modelDef)
switch string(modelDef.model_type)
    case "inv_g_linear"
        return;
    case "power_law"
        templateLib.gapInterpMode = "power_law";
        templateLib.gapInterpParameter = modelDef.parameter;
        templateLib.gapInterpBasisOrder = modelDef.basis_order;
    case "field_basis"
        templateLib.gapInterpMode = "field_basis";
        templateLib.gapInterpBasisName = string(modelDef.basis_name);
end
end

function summary = summarize_field_model_errors(rows)
parameterForGroup = rows.parameter;
parameterForGroup(isnan(parameterForGroup)) = realmax;
[groupId, modelLabel, modelType, basisName, parameter, basisOrder] = findgroups( ...
    rows.model_label, rows.model_type, rows.basis_name, parameterForGroup, rows.basis_order);
parameter(parameter == realmax) = NaN;
nGroup = max(groupId);

meanTemplate = splitapply(@mean, rows.template_rmse, groupId);
meanDerivative = splitapply(@mean, rows.derivative_rmse, groupId);
maxTemplate = splitapply(@max, rows.template_rmse, groupId);
maxDerivative = splitapply(@max, rows.derivative_rmse, groupId);

boundaryTemplate = nan(nGroup, 1);
boundaryDerivative = nan(nGroup, 1);
for ig = 1:nGroup
    mask = groupId == ig & rows.reconstruction_type == "boundary_extrapolation";
    boundaryTemplate(ig) = mean(rows.template_rmse(mask), 'omitnan');
    boundaryDerivative(ig) = mean(rows.derivative_rmse(mask), 'omitnan');
end

templateScale = median(meanTemplate, 'omitnan');
derivativeScale = median(meanDerivative, 'omitnan');
combinedScore = meanTemplate ./ max(templateScale, eps) + ...
    meanDerivative ./ max(derivativeScale, eps) + ...
    0.5 * boundaryTemplate ./ max(median(boundaryTemplate, 'omitnan'), eps);

summary = table(modelLabel, modelType, basisName, parameter, basisOrder, ...
    meanTemplate, meanDerivative, maxTemplate, maxDerivative, ...
    boundaryTemplate, boundaryDerivative, combinedScore, ...
    'VariableNames', {'model_label','model_type','basis_name','parameter', ...
    'basis_order','mean_template_rmse','mean_derivative_rmse', ...
    'max_template_rmse','max_derivative_rmse', ...
    'boundary_template_rmse','boundary_derivative_rmse','combined_score'});
summary = sortrows(summary, 'combined_score', 'ascend');
end
