%% Step 01: analyze physically motivated gap-coordinate models.
%
% This script evaluates how well different gap-direction coordinates
% reconstruct the finite static template library under leave-one-gap-out
% validation. The goal is to test whether the current 1/g linear coordinate is
% physically adequate, or whether a calibrated power-law/exponential response
% coordinate better reconstructs both F_g(x) and F'_g(x).

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

pGrid = 0.25:0.25:3.00;
lambdaGrid = 0.10:0.05:2.00;

modelDefs = make_model_defs(pGrid, lambdaGrid);
rows = table();
for im = 1:height(modelDefs)
    for ih = 1:numel(gapList)
        gHoldout = gapList(ih);
        result = evaluate_one_holdout(gapList, xCell, yCell, gHoldout, ...
            cfgAna.xGridN, modelDefs.model(im), modelDefs.parameter(im), ...
            modelDefs.basis_order(im));
        rows = [rows; result]; %#ok<AGROW>
    end
end

summary = summarize_model_errors(rows);
[~, order] = sort(summary.combined_score, 'ascend');
summary = summary(order, :);

writetable(rows, fullfile(outDir, 'gap_coordinate_model_leave_one_out.csv'));
writetable(summary, fullfile(outDir, 'gap_coordinate_model_summary.csv'));
save(fullfile(outDir, 'gap_coordinate_model_analysis.mat'), ...
    'rows', 'summary', 'pGrid', 'lambdaGrid', '-v7.3');

disp(summary);

function modelDefs = make_model_defs(pGrid, lambdaGrid)
model = strings(0, 1);
parameter = zeros(0, 1);
basisOrder = zeros(0, 1);

[model, parameter, basisOrder] = append_model(model, parameter, basisOrder, ...
    "g_linear", NaN, 1);
[model, parameter, basisOrder] = append_model(model, parameter, basisOrder, ...
    "inv_g_linear", NaN, 1);

for ip = 1:numel(pGrid)
    p = pGrid(ip);
    [model, parameter, basisOrder] = append_model(model, parameter, basisOrder, ...
        "power_law", p, 1);
    [model, parameter, basisOrder] = append_model(model, parameter, basisOrder, ...
        "power_law", p, 2);
end

for il = 1:numel(lambdaGrid)
    lambda = lambdaGrid(il);
    [model, parameter, basisOrder] = append_model(model, parameter, basisOrder, ...
        "exponential", lambda, 1);
end

modelDefs = table(model, parameter, basisOrder, ...
    'VariableNames', {'model','parameter','basis_order'});
end

function [model, parameter, basisOrder] = append_model(model, parameter, basisOrder, ...
    modelName, paramValue, orderValue)
model(end+1, 1) = modelName;
parameter(end+1, 1) = paramValue;
basisOrder(end+1, 1) = orderValue;
end

function row = evaluate_one_holdout(gapList, xCell, yCell, gHoldout, xGridN, ...
    modelName, parameter, basisOrder)
trainMask = abs(gapList - gHoldout) > 1e-12;
trainIds = find(trainMask);
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

idxHoldout = find(abs(gapList - gHoldout) < 1e-12, 1);
yTruth = interp1(xCell{idxHoldout}, yCell{idxHoldout}, xGrid, 'pchip');
dyTruth = gradient(yTruth, dx);

[yPred, dyPred] = reconstruct_template(gTrain, S, FxMat, gHoldout, ...
    modelName, parameter, basisOrder);

res = yPred(:) - yTruth(:);
dres = dyPred(:) - dyTruth(:);
if gHoldout < min(gTrain) || gHoldout > max(gTrain)
    reconstructionType = "boundary_extrapolation";
else
    reconstructionType = "interpolation";
end

row = table(string(modelName), parameter, basisOrder, gHoldout, ...
    string(reconstructionType), numel(gTrain), ...
    sqrt(mean(res.^2)), mean(abs(res)), max(abs(res)), ...
    sqrt(mean(dres.^2)), mean(abs(dres)), max(abs(dres)), ...
    'VariableNames', {'model','parameter','basis_order','g_holdout_mm', ...
    'reconstruction_type','num_gap_templates', ...
    'template_rmse','template_mean_abs','template_max_abs', ...
    'derivative_rmse','derivative_mean_abs','derivative_max_abs'});
end

function [curveGrid, dGrid] = reconstruct_template(gTrain, S, FxMat, gQuery, ...
    modelName, parameter, basisOrder)
Phi = build_gap_basis(gTrain, modelName, parameter, basisOrder);
phiQuery = build_gap_basis(gQuery, modelName, parameter, basisOrder);

coefCurve = Phi \ S;
coefDeriv = Phi \ FxMat;
curveGrid = phiQuery * coefCurve;
dGrid = phiQuery * coefDeriv;
end

function Phi = build_gap_basis(g, modelName, parameter, basisOrder)
g = g(:);
switch string(modelName)
    case "g_linear"
        eta = g;
        Phi = [ones(size(g)), eta];
    case "inv_g_linear"
        eta = 1 ./ g;
        Phi = [ones(size(g)), eta];
    case "power_law"
        p = parameter;
        if basisOrder == 1
            Phi = [ones(size(g)), g.^(-p)];
        else
            Phi = [ones(size(g)), g.^(-p), g.^(-(p + 1))];
        end
    case "exponential"
        lambda = parameter;
        eta = exp(-g ./ lambda);
        Phi = [ones(size(g)), eta];
    otherwise
        error('Unknown gap-coordinate model: %s', modelName);
end
end

function summary = summarize_model_errors(rows)
parameterForGroup = rows.parameter;
parameterForGroup(isnan(parameterForGroup)) = realmax;
[groupId, model, parameter, basisOrder] = findgroups(rows.model, ...
    parameterForGroup, rows.basis_order);
parameter(parameter == realmax) = NaN;
nGroup = max(groupId);

meanTemplate = splitapply(@mean, rows.template_rmse, groupId);
meanDerivative = splitapply(@mean, rows.derivative_rmse, groupId);
maxTemplate = splitapply(@max, rows.template_rmse, groupId);
maxDerivative = splitapply(@max, rows.derivative_rmse, groupId);

boundaryTemplate = zeros(nGroup, 1);
boundaryDerivative = zeros(nGroup, 1);
for ig = 1:nGroup
    mask = groupId == ig & rows.reconstruction_type == "boundary_extrapolation";
    boundaryTemplate(ig) = mean(rows.template_rmse(mask), 'omitnan');
    boundaryDerivative(ig) = mean(rows.derivative_rmse(mask), 'omitnan');
end

templateScale = median(meanTemplate, 'omitnan');
derivativeScale = median(meanDerivative, 'omitnan');
combinedScore = meanTemplate ./ max(templateScale, eps) + ...
    meanDerivative ./ max(derivativeScale, eps);

summary = table(model, parameter, basisOrder, meanTemplate, meanDerivative, ...
    maxTemplate, maxDerivative, boundaryTemplate, boundaryDerivative, ...
    combinedScore, ...
    'VariableNames', {'model','parameter','basis_order', ...
    'mean_template_rmse','mean_derivative_rmse', ...
    'max_template_rmse','max_derivative_rmse', ...
    'boundary_template_rmse','boundary_derivative_rmse', ...
    'combined_score'});
end
