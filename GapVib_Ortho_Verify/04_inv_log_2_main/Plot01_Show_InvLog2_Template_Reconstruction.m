%% Plot 01: static-template reconstruction comparison for inv_log_2.
close all;
clc;

ctx = load_inv_log_2_project_context(fileparts(mfilename('fullpath')));
gapList = ctx.gapList;
xCell = ctx.xCell;
yCell = ctx.yCell;
cfgAna = ctx.cfgAna;

if ~exist('gHoldout', 'var') || isempty(gHoldout)
    gHoldout = 0.2;
end

idxHoldout = find(abs(gapList - gHoldout) < 1e-12, 1);
if isempty(idxHoldout)
    error('Requested holdout gap %.6g mm is not available.', gHoldout);
end

modelDefs = table( ...
    ["baseline_1_over_g"; "field_basis_inv_log_2"; "power_law_p025_order2"], ...
    ["inv_g_linear"; "field_basis"; "power_law"], ...
    [""; "inv_log_2"; ""], ...
    [NaN; NaN; 0.25], ...
    [1; 3; 2], ...
    'VariableNames', {'model_label','response_model','basis_name','parameter','basis_order'});

xTrue = xCell{idxHoldout}(:);
yTrue = yCell{idxHoldout}(:);

fig = figure('Name', 'inv_log_2 template reconstruction', 'Color', 'w');
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile;
plot(xTrue, yTrue, 'k-', 'LineWidth', 1.6, 'DisplayName', sprintf('truth, g = %.1f mm', gHoldout));
hold(ax1, 'on');

ax2 = nexttile;
dyTruthRaw = gradient(yTrue, mean(diff(xTrue)));
plot(xTrue, dyTruthRaw, 'k-', 'LineWidth', 1.6, 'DisplayName', 'truth derivative');
hold(ax2, 'on');

colors = [0.15 0.15 0.15; 0.00 0.45 0.74; 0.85 0.33 0.10];
summaryRows = table();

for im = 1:height(modelDefs)
    trustOptions = make_gap_trust_options(cfgAna);
    templateLib = make_response_template_library(gapList, xCell, yCell, ...
        gHoldout, cfgAna.xGridN, modelDefs(im, :), trustOptions);
    xEval = templateLib.xGrid(:);
    yPred = eval_gap_template(templateLib, gHoldout, xEval);
    dyPred = eval_gap_derivative(templateLib, gHoldout, xEval);
    yTruthEval = interp1(xTrue, yTrue, xEval, 'pchip');
    dyTruthEval = gradient(yTruthEval, mean(diff(xEval)));

    plot(ax1, xEval, yPred, '-', 'Color', colors(im, :), 'LineWidth', 1.3, ...
        'DisplayName', char(modelDefs.model_label(im)));
    plot(ax2, xEval, dyPred, '-', 'Color', colors(im, :), 'LineWidth', 1.3, ...
        'DisplayName', char(modelDefs.model_label(im)));

    row = table(modelDefs.model_label(im), ...
        sqrt(mean((yPred - yTruthEval).^2)), ...
        sqrt(mean((dyPred - dyTruthEval).^2)), ...
        'VariableNames', {'model_label','template_rmse','derivative_rmse'});
    summaryRows = [summaryRows; row]; %#ok<AGROW>
end

title(ax1, '(a) Static template reconstruction', 'FontWeight', 'normal', 'FontSize', 9);
xlabel(ax1, 'Spatial coordinate x (mm)', 'FontSize', 9);
ylabel(ax1, 'Voltage (V)', 'FontSize', 9);
legend(ax1, 'Location', 'best');

title(ax2, '(b) Template derivative comparison', 'FontWeight', 'normal', 'FontSize', 9);
xlabel(ax2, 'Spatial coordinate x (mm)', 'FontSize', 9);
ylabel(ax2, 'dV/dx', 'FontSize', 9);
legend(ax2, 'Location', 'best');

annotationText = sprintf([ ...
    'RMSE summary\n', ...
    'baseline: %.4g / %.4g\n', ...
    'inv\\_log\\_2: %.4g / %.4g\n', ...
    'power\\_law: %.4g / %.4g'], ...
    summaryRows.template_rmse(1), summaryRows.derivative_rmse(1), ...
    summaryRows.template_rmse(2), summaryRows.derivative_rmse(2), ...
    summaryRows.template_rmse(3), summaryRows.derivative_rmse(3));
annotation(fig, 'textbox', [0.68, 0.13, 0.25, 0.18], ...
    'String', annotationText, ...
    'FitBoxToText', 'on', ...
    'BackgroundColor', 'w', ...
    'EdgeColor', [0.7 0.7 0.7], ...
    'Interpreter', 'tex', ...
    'FontName', 'Times New Roman', ...
    'FontSize', 8);

apply_inv_log_2_figure_style(ax1, 170, 120);
