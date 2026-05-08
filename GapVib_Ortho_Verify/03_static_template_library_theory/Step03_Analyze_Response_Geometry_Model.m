%% Step 03: analyze the response-geometry model F_g(x) = H(g + s(x)).
%
% This script does not modify any existing main-method code. It introduces
% an independent theoretical analysis branch that fits a nonparametric
% monotone response H(d) together with a smooth geometry-offset profile s(x)
% from the static template library, then compares it against the current
% 1/g baseline and the best calibrated power-law model.

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

if ~exist('xGridNAnalyze', 'var') || isempty(xGridNAnalyze)
    xGridNAnalyze = min(cfgAna.xGridN, 401);
end
if ~exist('opts', 'var') || isempty(opts)
    opts = default_response_geometry_options();
end
if ~exist('analyzeGapList', 'var') || isempty(analyzeGapList)
    analyzeGapList = gapList(:).';
else
    analyzeGapList = intersect(analyzeGapList(:).', gapList(:).', 'stable');
end
allRows = table();

for ig = 1:numel(analyzeGapList)
    gHoldout = analyzeGapList(ig);
    idxHoldout = find(abs(gapList - gHoldout) < 1e-12, 1);
    [xGrid, yTruth] = build_holdout_truth_curve(gapList, xCell, yCell, gHoldout, xGridNAnalyze);

    baselineLib = build_gap_template_library(gapList, xCell, yCell, gHoldout, xGridNAnalyze);
    if numel(xGrid) ~= numel(baselineLib.xGrid) || any(abs(xGrid - baselineLib.xGrid) > 1e-12)
        xGrid = baselineLib.xGrid(:);
        yTruth = interp1(xCell{idxHoldout}(:), yCell{idxHoldout}(:), xGrid, 'pchip');
    end

    fprintf('Step03 holdout g = %.4g mm (%d/%d)\n', gHoldout, ig, numel(analyzeGapList));
    responseGeom = fit_response_geometry_holdout(gapList, xCell, yCell, gHoldout, xGridNAnalyze, opts);

    methods = make_method_defs();
    for im = 1:height(methods)
        methodName = string(methods.method_name(im));
        switch methodName
            case "baseline_1_over_g"
                yPred = eval_gap_template(baselineLib, gHoldout, xGrid);
                dyPred = eval_gap_derivative(baselineLib, gHoldout, xGrid);
            case "power_law_p025_order2"
                modelLib = baselineLib;
                modelLib.gapInterpMode = "power_law";
                modelLib.gapInterpParameter = 0.25;
                modelLib.gapInterpBasisOrder = 2;
                yPred = eval_gap_template(modelLib, gHoldout, xGrid);
                dyPred = eval_gap_derivative(modelLib, gHoldout, xGrid);
            case "response_geometry_joint"
                yPred = responseGeom.yPred(:);
                dyPred = responseGeom.dyPred(:);
        end

        dx = mean(diff(xGrid));
        dyTruth = gradient(yTruth, dx);
        row = table();
        row.method_name = methodName;
        row.g_holdout_mm = gHoldout;
        row.reconstruction_type = classify_holdout_type(gHoldout, gapList);
        row.template_rmse = sqrt(mean((yPred(:) - yTruth(:)).^2));
        row.derivative_rmse = sqrt(mean((dyPred(:) - dyTruth(:)).^2));
        if methodName == "response_geometry_joint"
            row.full_fit_rmse = responseGeom.trainFitRmse;
            row.refine_shift_rms_mm = sqrt(mean(responseGeom.sHat(:).^2));
            row.ref_col_index = responseGeom.refCol;
        else
            row.full_fit_rmse = NaN;
            row.refine_shift_rms_mm = NaN;
            row.ref_col_index = NaN;
        end
        allRows = [allRows; row]; %#ok<AGROW>
    end
end

summary = summarize_response_geometry_results(allRows);
writetable(allRows, fullfile(outDir, 'response_geometry_holdout_comparison.csv'));
writetable(summary, fullfile(outDir, 'response_geometry_summary.csv'));

fullFit = fit_response_geometry_holdout(gapList, xCell, yCell, NaN, cfgAna.xGridN, opts);
profileTable = table(fullFit.xGrid(:), fullFit.sHat(:), ...
    'VariableNames', {'x_mm', 's_hat_mm'});
writetable(profileTable, fullfile(outDir, 'response_geometry_fullfit_profile.csv'));

save(fullfile(outDir, 'response_geometry_analysis.mat'), ...
    'allRows', 'summary', 'profileTable', 'fullFit', 'opts', '-v7.3');

disp(summary);

function opts = default_response_geometry_options()
opts = struct();
opts.initShiftMax = 0.60;
opts.initShiftCount = 161;
opts.gDenseCount = 201;
opts.smoothWindow = 41;
opts.refineIters = 2;
opts.refineHalfWidthList = [0.06, 0.03];
opts.refineShiftCount = 61;
end

function methods = make_method_defs()
methods = table(["baseline_1_over_g"; "power_law_p025_order2"; "response_geometry_joint"], ...
    'VariableNames', {'method_name'});
end

function [xGrid, yTruth] = build_holdout_truth_curve(gapList, xCell, yCell, gHoldout, xGridN)
libAll = build_gap_template_library(gapList, xCell, yCell, NaN, xGridN);
xGrid = libAll.xGrid(:);
idxHoldout = find(abs(gapList - gHoldout) < 1e-12, 1);
yTruth = interp1(xCell{idxHoldout}(:), yCell{idxHoldout}(:), xGrid, 'pchip');
end

function label = classify_holdout_type(gHoldout, gapList)
if abs(gHoldout - min(gapList)) < 1e-12 || abs(gHoldout - max(gapList)) < 1e-12
    label = "boundary_extrapolation";
else
    label = "interpolation";
end
end

function model = fit_response_geometry_holdout(gapList, xCell, yCell, gHoldout, xGridN, opts)
templateLib = build_gap_template_library(gapList, xCell, yCell, gHoldout, xGridN);
gTrain = templateLib.gapTrain(:);
xGrid = templateLib.xGrid(:);
Y = templateLib.S;
dx = mean(diff(xGrid));

refCol = pick_reference_column(Y);
sHat = estimate_initial_shifts(gTrain, Y, refCol, opts);

for iter = 1:opts.refineIters
    HModel = fit_monotone_response(gTrain, sHat, Y, refCol);
    halfWidth = opts.refineHalfWidthList(min(iter, numel(opts.refineHalfWidthList)));
    sHat = refine_shifts_against_response(gTrain, Y, sHat, HModel, halfWidth, opts.refineShiftCount);
    sHat = smooth_and_center_shifts(sHat, opts.smoothWindow);
end

HModel = fit_monotone_response(gTrain, sHat, Y, refCol);
YFit = evaluate_response_geometry(HModel, sHat, gTrain);
trainFitRmse = sqrt(mean((YFit(:) - Y(:)).^2));

model = struct();
model.templateLib = templateLib;
model.gTrain = gTrain;
model.xGrid = xGrid;
model.sHat = sHat(:);
model.refCol = refCol;
model.HModel = HModel;
model.trainFitRmse = trainFitRmse;
if ~isnan(gHoldout)
    model.yPred = evaluate_response_geometry(HModel, sHat, gHoldout);
    model.dyPred = gradient(model.yPred(:), dx);
end
end

function refCol = pick_reference_column(Y)
[~, refCol] = max(std(Y, 0, 1));
end

function sHat = estimate_initial_shifts(gTrain, Y, refCol, opts)
gDense = linspace(min(gTrain), max(gTrain), opts.gDenseCount)';
shiftGrid = linspace(-opts.initShiftMax, opts.initShiftMax, opts.initShiftCount);
yRefDense = interp1(gTrain, Y(:, refCol), gDense, 'pchip');
nX = size(Y, 2);
sHat = zeros(nX, 1);

for j = 1:nX
    yColDense = interp1(gTrain, Y(:, j), gDense, 'pchip');
    errBest = inf;
    shiftBest = 0;
    for ishift = 1:numel(shiftGrid)
        delta = shiftGrid(ishift);
        gShift = gDense + delta;
        valid = gShift >= min(gTrain) & gShift <= max(gTrain);
        if nnz(valid) < max(10, round(0.25 * numel(gDense)))
            continue;
        end
        yRefShift = interp1(gDense, yRefDense, gShift(valid), 'pchip');
        err = mean((yColDense(valid) - yRefShift).^2);
        if err < errBest
            errBest = err;
            shiftBest = delta;
        end
    end
    sHat(j) = shiftBest;
end

sHat = smooth_and_center_shifts(sHat, opts.smoothWindow);
end

function sHat = refine_shifts_against_response(gTrain, Y, sHat0, HModel, halfWidth, nShift)
nX = size(Y, 2);
sHat = sHat0(:);
shiftGrid = linspace(-halfWidth, halfWidth, nShift);

for j = 1:nX
    errBest = inf;
    shiftBest = sHat0(j);
    for ishift = 1:numel(shiftGrid)
        sCand = sHat0(j) + shiftGrid(ishift);
        yPred = eval_monotone_response(HModel, gTrain + sCand);
        err = mean((Y(:, j) - yPred(:)).^2);
        if err < errBest
            errBest = err;
            shiftBest = sCand;
        end
    end
    sHat(j) = shiftBest;
end
end

function sHat = smooth_and_center_shifts(sHat, smoothWindow)
sHat = smoothdata(sHat(:), 'sgolay', min(smoothWindow, numel(sHat)));
sHat = sHat - mean(sHat);
end

function HModel = fit_monotone_response(gTrain, sHat, Y, refCol)
Z = gTrain(:) + sHat(:).';
zVec = Z(:);
yVec = Y(:);
[zSort, order] = sort(zVec, 'ascend');
ySort = yVec(order);

dirScore = corr(gTrain(:), Y(:, refCol), 'rows', 'complete');
if isnan(dirScore) || dirScore >= 0
    signDir = 1;
else
    signDir = -1;
end

yIso = signDir * isotonic_increasing(zSort, signDir * ySort);
[zUnique, ~, idxUnique] = unique(zSort, 'stable');
yUnique = accumarray(idxUnique, yIso, [], @mean);

HModel = struct();
HModel.zSupport = zUnique(:);
HModel.ySupport = yUnique(:);
HModel.signDir = signDir;
end

function yFit = evaluate_response_geometry(HModel, sHat, gValue)
zQuery = gValue + sHat(:).';
yFit = eval_monotone_response(HModel, zQuery);
end

function y = eval_monotone_response(HModel, zQuery)
y = interp1(HModel.zSupport, HModel.ySupport, zQuery, 'linear', 'extrap');
end

function yIso = isotonic_increasing(x, y)
% x is sorted. Standard pooled-adjacent-violators on y with unit weights.
n = numel(y);
level = zeros(n, 1);
weight = zeros(n, 1);
startIdx = zeros(n, 1);
nBlock = 0;

for i = 1:n
    nBlock = nBlock + 1;
    level(nBlock) = y(i);
    weight(nBlock) = 1;
    startIdx(nBlock) = i;
    while nBlock > 1 && level(nBlock - 1) > level(nBlock)
        newWeight = weight(nBlock - 1) + weight(nBlock);
        newLevel = (weight(nBlock - 1) * level(nBlock - 1) + ...
            weight(nBlock) * level(nBlock)) / newWeight;
        level(nBlock - 1) = newLevel;
        weight(nBlock - 1) = newWeight;
        nBlock = nBlock - 1;
    end
end

level = level(1:nBlock);
startIdx = startIdx(1:nBlock);
endIdx = [startIdx(2:end) - 1; n];
yIso = zeros(size(y));
for ib = 1:nBlock
    yIso(startIdx(ib):endIdx(ib)) = level(ib);
end
end

function summary = summarize_response_geometry_results(allRows)
[groupId, methodName] = findgroups(allRows.method_name);
summary = table(methodName, ...
    splitapply(@mean, allRows.template_rmse, groupId), ...
    splitapply(@max, allRows.template_rmse, groupId), ...
    splitapply(@mean, allRows.derivative_rmse, groupId), ...
    splitapply(@max, allRows.derivative_rmse, groupId), ...
    splitapply(@mean, allRows.full_fit_rmse, groupId), ...
    'VariableNames', {'method_name', 'mean_template_rmse', 'max_template_rmse', ...
    'mean_derivative_rmse', 'max_derivative_rmse', 'mean_full_fit_rmse'});

boundaryMask = allRows.reconstruction_type == "boundary_extrapolation";
boundaryRows = allRows(boundaryMask, :);
boundaryTemplate = nan(height(summary), 1);
boundaryDerivative = nan(height(summary), 1);
for i = 1:height(summary)
    mask = boundaryRows.method_name == summary.method_name(i);
    if any(mask)
        boundaryTemplate(i) = mean(boundaryRows.template_rmse(mask));
        boundaryDerivative(i) = mean(boundaryRows.derivative_rmse(mask));
    end
end
summary.boundary_mean_template_rmse = boundaryTemplate;
summary.boundary_mean_derivative_rmse = boundaryDerivative;
summary = sortrows(summary, {'mean_template_rmse', 'mean_derivative_rmse'}, {'ascend', 'ascend'});
end
