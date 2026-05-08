%% Step 01: analyze static gap library reconstruction error.

close all;
clc;

caseDir = fileparts(mfilename('fullpath'));
projectDir = fileparts(caseDir);
rootDir = fileparts(projectDir);
addpath(fullfile(projectDir, 'func'));

oldDir = resolve_legacy_config_dir(rootDir);
outDir = fullfile(projectDir, 'results', '02_gap_library_sensitivity');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

load(fullfile(oldDir, 'stage0_config.mat'), 'cfg');
cfg.dataFile = resolve_project_data_file(rootDir, cfg.dataFile);
cfgAna = make_multicase_cfg(cfg);
[gapList, xCell, yCell] = load_stacked_gap_curves(cfg.dataFile);

rows = table();
for ih = 1:numel(gapList)
    gTrue = gapList(ih);
    xTrue = xCell{ih}(:);
    yTrue = yCell{ih}(:);

    templateCurrent = build_gap_template_library(gapList, xCell, yCell, ...
        gTrue, cfgAna.xGridN);
    rows = [rows; analyze_one_library("leave_one_out_inv_g_linear", ...
        templateCurrent, gTrue, xTrue, yTrue)]; %#ok<AGROW>

    templateIdeal = build_gap_template_library(gapList, xCell, yCell, ...
        NaN, cfgAna.xGridN);
    rows = [rows; analyze_one_library("ideal_in_library", ...
        templateIdeal, gTrue, xTrue, yTrue)]; %#ok<AGROW>

    sparseMask = ismember(round(gapList, 12), round([0.2; 0.6; 1.0; 1.4], 12));
    sparseMask(ih) = false;
    if nnz(sparseMask) >= 2
        templateSparse = build_gap_template_library(gapList(sparseMask), ...
            xCell(sparseMask), yCell(sparseMask), NaN, cfgAna.xGridN);
        rows = [rows; analyze_one_library("sparse_leave_one_out_inv_g_linear", ...
            templateSparse, gTrue, xTrue, yTrue)]; %#ok<AGROW>
    end

    templateGLin = build_custom_gap_template_library(gapList, xCell, yCell, ...
        gTrue, cfgAna.xGridN, "g_linear");
    rows = [rows; analyze_one_library("leave_one_out_g_linear", ...
        templateGLin, gTrue, xTrue, yTrue)]; %#ok<AGROW>

    templateGPchip = build_custom_gap_template_library(gapList, xCell, yCell, ...
        gTrue, cfgAna.xGridN, "g_pchip");
    rows = [rows; analyze_one_library("leave_one_out_g_pchip", ...
        templateGPchip, gTrue, xTrue, yTrue)]; %#ok<AGROW>
end

writetable(rows, fullfile(outDir, 'template_library_error.csv'));
save(fullfile(outDir, 'template_library_error.mat'), 'rows', '-v7.3');
disp(rows);

function row = analyze_one_library(libraryMode, templateLib, gTrue, xTrue, yTrue)
xEval = templateLib.xGrid(:);
yTruth = interp1(xTrue, yTrue, xEval, 'pchip');
yPred = eval_gap_template_generic(templateLib, gTrue, xEval);
dyTruth = gradient(yTruth, mean(diff(xEval)));
dyPred = eval_gap_derivative_generic(templateLib, gTrue, xEval);
res = yPred - yTruth;
dres = dyPred - dyTruth;
row = table(string(libraryMode), gTrue, numel(templateLib.gapTrain), ...
    sqrt(mean(res.^2)), max(abs(res)), mean(abs(res)), ...
    sqrt(mean(dres.^2)), max(abs(dres)), mean(abs(dres)), ...
    'VariableNames', {'library_mode','g_true_mm','num_gap_templates', ...
    'template_rmse','template_max_abs','template_mean_abs', ...
    'derivative_rmse','derivative_max_abs','derivative_mean_abs'});
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
