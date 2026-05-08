%% Step 02: verify equivalent response assumptions against the 1/g baseline.
%
% Step 01 ranks many candidate gap coordinates. This script turns that ranking
% into a stricter verification table:
%   1) select representative equivalent-response models;
%   2) compare each model with the current 1/g baseline on matched holdouts;
%   3) quantify full-library linearity with R^2 across x for F_g and F'_g.

close all;
clc;

caseDir = fileparts(mfilename('fullpath'));
projectDir = fileparts(caseDir);
rootDir = fileparts(projectDir);
addpath(fullfile(projectDir, 'func'));

outDir = fullfile(projectDir, 'results', '03_static_template_library_theory');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

summaryPath = fullfile(outDir, 'gap_coordinate_model_summary.csv');
rowsPath = fullfile(outDir, 'gap_coordinate_model_leave_one_out.csv');
if ~exist(summaryPath, 'file') || ~exist(rowsPath, 'file')
    run(fullfile(caseDir, 'Step01_Analyze_Gap_Coordinate_Models.m'));
end

summary = readtable(summaryPath, 'TextType', 'string');
rows = readtable(rowsPath, 'TextType', 'string');
selectedModels = select_representative_models(summary);

comparisonRows = compare_to_inv_g_baseline(rows, selectedModels);
linearityRows = evaluate_full_library_linearity(rootDir, selectedModels);
validationSummary = summarize_validation(comparisonRows, linearityRows);

writetable(selectedModels, fullfile(outDir, 'equivalent_response_selected_models.csv'));
writetable(comparisonRows, fullfile(outDir, 'equivalent_response_baseline_comparison.csv'));
writetable(linearityRows, fullfile(outDir, 'equivalent_response_full_library_linearity.csv'));
writetable(validationSummary, fullfile(outDir, 'equivalent_response_validation_summary.csv'));
save(fullfile(outDir, 'equivalent_response_validation.mat'), ...
    'selectedModels', 'comparisonRows', 'linearityRows', 'validationSummary', '-v7.3');

disp(validationSummary);

function selected = select_representative_models(summary)
selected = table();

selected = [selected; pick_exact(summary, "inv_g_linear", NaN, 1, "current_1_over_g_baseline")]; %#ok<AGROW>
selected = [selected; pick_best(summary, "power_law", 1, "best_power_law_first_order")]; %#ok<AGROW>
selected = [selected; pick_best(summary, "power_law", 2, "best_power_law_second_order")]; %#ok<AGROW>
selected = [selected; pick_best(summary, "exponential", 1, "best_exponential_first_order")]; %#ok<AGROW>
selected = [selected; add_label(summary(1, :), "best_overall")]; %#ok<AGROW>

[~, ia] = unique(compose_key(selected.model, selected.parameter, selected.basis_order), 'stable');
selected = selected(ia, :);
end

function row = pick_exact(summary, modelName, parameter, basisOrder, label)
if isnan(parameter)
    mask = summary.model == modelName & isnan(summary.parameter) & ...
        summary.basis_order == basisOrder;
else
    mask = summary.model == modelName & abs(summary.parameter - parameter) < 1e-12 & ...
        summary.basis_order == basisOrder;
end
row = add_label(summary(find(mask, 1), :), label);
end

function row = pick_best(summary, modelName, basisOrder, label)
mask = summary.model == modelName & summary.basis_order == basisOrder;
sub = summary(mask, :);
[~, idx] = min(sub.combined_score);
row = add_label(sub(idx, :), label);
end

function row = add_label(row, label)
if isempty(row)
    error('Could not select representative model: %s', label);
end
row.selection_label = string(label);
row = movevars(row, 'selection_label', 'Before', 1);
end

function comparisonRows = compare_to_inv_g_baseline(rows, selectedModels)
baselineRows = rows(rows.model == "inv_g_linear" & rows.basis_order == 1, :);
comparisonRows = table();

for im = 1:height(selectedModels)
    modelRow = selectedModels(im, :);
    modelRows = filter_rows(rows, modelRow.model, modelRow.parameter, modelRow.basis_order);
    for ir = 1:height(modelRows)
        g = modelRows.g_holdout_mm(ir);
        base = baselineRows(abs(baselineRows.g_holdout_mm - g) < 1e-12, :);
        if isempty(base)
            continue;
        end
        templateImprovement = 100 * (base.template_rmse - modelRows.template_rmse(ir)) ./ ...
            max(base.template_rmse, eps);
        derivativeImprovement = 100 * (base.derivative_rmse - modelRows.derivative_rmse(ir)) ./ ...
            max(base.derivative_rmse, eps);
        row = table(modelRow.selection_label, modelRow.model, modelRow.parameter, ...
            modelRow.basis_order, g, modelRows.reconstruction_type(ir), ...
            modelRows.template_rmse(ir), base.template_rmse, templateImprovement, ...
            modelRows.derivative_rmse(ir), base.derivative_rmse, derivativeImprovement, ...
            'VariableNames', {'selection_label','model','parameter','basis_order', ...
            'g_holdout_mm','reconstruction_type', ...
            'template_rmse','baseline_template_rmse','template_improvement_pct', ...
            'derivative_rmse','baseline_derivative_rmse','derivative_improvement_pct'});
        comparisonRows = [comparisonRows; row]; %#ok<AGROW>
    end
end
end

function modelRows = filter_rows(rows, modelName, parameter, basisOrder)
if isnan(parameter)
    mask = rows.model == modelName & isnan(rows.parameter) & rows.basis_order == basisOrder;
else
    mask = rows.model == modelName & abs(rows.parameter - parameter) < 1e-12 & ...
        rows.basis_order == basisOrder;
end
modelRows = rows(mask, :);
end

function linearityRows = evaluate_full_library_linearity(rootDir, selectedModels)
oldDir = resolve_legacy_config_dir(rootDir);
load(fullfile(oldDir, 'stage0_config.mat'), 'cfg');
cfg.dataFile = resolve_project_data_file(rootDir, cfg.dataFile);
cfgAna = make_multicase_cfg(cfg);
[gapList, xCell, yCell] = load_stacked_gap_curves(cfg.dataFile);

[xGrid, S, FxMat] = build_full_static_matrix(gapList, xCell, yCell, cfgAna.xGridN);
linearityRows = table();
for im = 1:height(selectedModels)
    modelRow = selectedModels(im, :);
    Phi = build_gap_basis(gapList, modelRow.model, modelRow.parameter, modelRow.basis_order);
    curveFit = Phi * (Phi \ S);
    derivFit = Phi * (Phi \ FxMat);
    curveR2 = columnwise_r2(S, curveFit);
    derivR2 = columnwise_r2(FxMat, derivFit);
    curveRmse = sqrt(mean((curveFit(:) - S(:)).^2));
    derivRmse = sqrt(mean((derivFit(:) - FxMat(:)).^2));
    row = table(modelRow.selection_label, modelRow.model, modelRow.parameter, ...
        modelRow.basis_order, numel(xGrid), ...
        mean(curveR2, 'omitnan'), local_percentile(curveR2, 5), curveRmse, ...
        mean(derivR2, 'omitnan'), local_percentile(derivR2, 5), derivRmse, ...
        'VariableNames', {'selection_label','model','parameter','basis_order', ...
        'x_grid_count','mean_template_r2','p05_template_r2','template_fit_rmse', ...
        'mean_derivative_r2','p05_derivative_r2','derivative_fit_rmse'});
    linearityRows = [linearityRows; row]; %#ok<AGROW>
end
end

function [xGrid, S, FxMat] = build_full_static_matrix(gapList, xCell, yCell, xGridN)
xMin = -inf;
xMax = inf;
for ig = 1:numel(gapList)
    xMin = max(xMin, min(xCell{ig}));
    xMax = min(xMax, max(xCell{ig}));
end
xPad = 0.05 * (xMax - xMin);
xGrid = linspace(xMin + xPad, xMax - xPad, xGridN)';
S = zeros(numel(gapList), numel(xGrid));
for ig = 1:numel(gapList)
    S(ig, :) = interp1(xCell{ig}, yCell{ig}, xGrid, 'pchip');
end
dx = mean(diff(xGrid));
FxMat = zeros(size(S));
for ig = 1:size(S, 1)
    FxMat(ig, :) = gradient(S(ig, :), dx);
end
end

function Phi = build_gap_basis(g, modelName, parameter, basisOrder)
g = g(:);
switch string(modelName)
    case "g_linear"
        Phi = [ones(size(g)), g];
    case "inv_g_linear"
        Phi = [ones(size(g)), 1 ./ g];
    case "power_law"
        if basisOrder == 1
            Phi = [ones(size(g)), g.^(-parameter)];
        else
            Phi = [ones(size(g)), g.^(-parameter), g.^(-(parameter + 1))];
        end
    case "exponential"
        Phi = [ones(size(g)), exp(-g ./ parameter)];
    otherwise
        error('Unknown model: %s', modelName);
end
end

function r2 = columnwise_r2(y, yFit)
yMean = mean(y, 1);
sst = sum((y - yMean).^2, 1);
sse = sum((y - yFit).^2, 1);
r2 = 1 - sse ./ max(sst, eps);
r2 = r2(:);
end

function validationSummary = summarize_validation(comparisonRows, linearityRows)
parameterForGroup = comparisonRows.parameter;
parameterForGroup(isnan(parameterForGroup)) = realmax;
[groupId, label, model, parameter, basisOrder] = findgroups(comparisonRows.selection_label, ...
    comparisonRows.model, parameterForGroup, comparisonRows.basis_order);
parameter(parameter == realmax) = NaN;

meanTemplateImprovement = splitapply(@mean, comparisonRows.template_improvement_pct, groupId);
meanDerivativeImprovement = splitapply(@mean, comparisonRows.derivative_improvement_pct, groupId);
minTemplateImprovement = splitapply(@min, comparisonRows.template_improvement_pct, groupId);
minDerivativeImprovement = splitapply(@min, comparisonRows.derivative_improvement_pct, groupId);

boundaryTemplateImprovement = zeros(max(groupId), 1);
boundaryDerivativeImprovement = zeros(max(groupId), 1);
for ig = 1:max(groupId)
    mask = groupId == ig & comparisonRows.reconstruction_type == "boundary_extrapolation";
    boundaryTemplateImprovement(ig) = mean(comparisonRows.template_improvement_pct(mask), 'omitnan');
    boundaryDerivativeImprovement(ig) = mean(comparisonRows.derivative_improvement_pct(mask), 'omitnan');
end

validationSummary = table(label, model, parameter, basisOrder, ...
    meanTemplateImprovement, meanDerivativeImprovement, ...
    minTemplateImprovement, minDerivativeImprovement, ...
    boundaryTemplateImprovement, boundaryDerivativeImprovement, ...
    'VariableNames', {'selection_label','model','parameter','basis_order', ...
    'mean_template_improvement_pct','mean_derivative_improvement_pct', ...
    'min_template_improvement_pct','min_derivative_improvement_pct', ...
    'boundary_template_improvement_pct','boundary_derivative_improvement_pct'});

linearityForJoin = removevars(linearityRows, {'model','parameter','basis_order'});
validationSummary = innerjoin(validationSummary, linearityForJoin, ...
    'Keys', 'selection_label');
[~, order] = sort(validationSummary.mean_template_improvement_pct + ...
    validationSummary.mean_derivative_improvement_pct, 'descend');
validationSummary = validationSummary(order, :);
end

function key = compose_key(model, parameter, basisOrder)
paramText = strings(size(parameter));
for ii = 1:numel(parameter)
    if isnan(parameter(ii))
        paramText(ii) = "NaN";
    else
        paramText(ii) = sprintf('%.12g', parameter(ii));
    end
end
key = model + "|" + paramText + "|" + string(basisOrder);
end
