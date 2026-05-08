%% Plot 07: trusted fitting-domain selection for inv_log_2.
close all;
clc;

ctx = load_inv_log_2_project_context(fileparts(mfilename('fullpath')));
gapList = ctx.gapList;
xCell = ctx.xCell;
yCell = ctx.yCell;
cfgAna = ctx.cfgAna;

if ~exist('gTrue', 'var') || isempty(gTrue)
    if any(abs(gapList - 0.2) < 1e-12)
        gTrue = 0.2;
    else
        gTrue = ctx.cfg.g_holdout;
    end
end

idxTrue = find(abs(gapList - gTrue) < 1e-12, 1);
modelInv = table("field_basis_inv_log_2", "field_basis", "inv_log_2", NaN, 3, ...
    'VariableNames', {'model_label','response_model','basis_name','parameter','basis_order'});
trustOptions = make_gap_trust_options(cfgAna);
templateFull = make_response_template_library(gapList, xCell, yCell, ...
    gTrue, cfgAna.xGridN, modelInv, struct('enable', false));
templateTrust = make_response_template_library(gapList, xCell, yCell, ...
    gTrue, cfgAna.xGridN, modelInv, trustOptions);
metrics = templateTrust.trustInfo.metrics;

fig = figure('Name', 'inv_log_2 trusted domain selection', 'Color', 'w');
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile;
hold(ax1, 'on');
trainGaps = templateTrust.gapTrain(:);
colors = turbo(numel(trainGaps));
for ig = 1:numel(trainGaps)
    idxTrain = find(abs(gapList(:) - trainGaps(ig)) < 1e-12, 1);
    plot(ax1, xCell{idxTrain}(:), yCell{idxTrain}(:), '-', ...
        'Color', colors(ig, :), 'LineWidth', 0.8);
end
plot(ax1, xCell{idxTrue}(:), yCell{idxTrue}(:), 'k--', 'LineWidth', 1.1, ...
    'DisplayName', 'held-out gap');
xline(ax1, templateTrust.domain(1), ':', 'Color', [0.10 0.10 0.10], ...
    'LineWidth', 1.0, 'DisplayName', 'trusted boundary');
xline(ax1, templateTrust.domain(2), ':', 'Color', [0.10 0.10 0.10], ...
    'LineWidth', 1.0, 'HandleVisibility', 'off');
title(ax1, '(a) Static library and selected window', 'FontWeight', 'normal', 'FontSize', 9);
xlabel(ax1, 'Spatial coordinate x (mm)', 'FontSize', 9);
ylabel(ax1, 'Voltage (V)', 'FontSize', 9);
legend(ax1, 'Location', 'best');

ax2 = nexttile;
yyaxis(ax2, 'left');
plot(ax2, metrics.half_width_mm, metrics.template_rmse, 'o-', ...
    'LineWidth', 1.2, 'MarkerSize', 4, 'DisplayName', 'template error');
ylabel(ax2, 'Normalized waveform RMSE', 'FontSize', 9);
yyaxis(ax2, 'right');
plot(ax2, metrics.half_width_mm, metrics.derivative_rmse, 's-', ...
    'LineWidth', 1.2, 'MarkerSize', 4, 'DisplayName', 'derivative error');
ylabel(ax2, 'Normalized derivative RMSE', 'FontSize', 9);
xline(ax2, templateTrust.trustInfo.selectedHalfWidth, 'k:', 'LineWidth', 1.0);
title(ax2, '(b) Leave-one-gap validation errors', 'FontWeight', 'normal', 'FontSize', 9);
xlabel(ax2, 'Candidate half-width L (mm)', 'FontSize', 9);

ax3 = nexttile;
plot(ax3, metrics.half_width_mm, metrics.score, '^-', ...
    'Color', [0.00 0.45 0.74], 'LineWidth', 1.3, 'MarkerSize', 4, ...
    'DisplayName', 'validation score');
hold(ax3, 'on');
xline(ax3, templateTrust.trustInfo.selectedHalfWidth, 'k:', 'LineWidth', 1.0, ...
    'DisplayName', 'selected L');
title(ax3, '(c) Selection score', 'FontWeight', 'normal', 'FontSize', 9);
xlabel(ax3, 'Candidate half-width L (mm)', 'FontSize', 9);
ylabel(ax3, 'Score', 'FontSize', 9);
legend(ax3, 'Location', 'best');

ax4 = nexttile;
xFull = templateFull.xGrid(:);
yFull = eval_gap_template(templateFull, gTrue, xFull);
xTrust = templateTrust.xGrid(:);
yTrust = eval_gap_template(templateTrust, gTrue, xTrust);
plot(ax4, xFull, yFull, '-', 'Color', [0.70 0.70 0.70], ...
    'LineWidth', 1.0, 'DisplayName', 'full static template');
hold(ax4, 'on');
plot(ax4, xTrust, yTrust, '-', 'Color', [0.85 0.33 0.10], ...
    'LineWidth', 1.5, 'DisplayName', 'trusted fitting template');
xline(ax4, templateTrust.domain(1), ':', 'Color', [0.10 0.10 0.10], ...
    'LineWidth', 1.0, 'HandleVisibility', 'off');
xline(ax4, templateTrust.domain(2), ':', 'Color', [0.10 0.10 0.10], ...
    'LineWidth', 1.0, 'HandleVisibility', 'off');
title(ax4, '(d) Full pulse versus fitting segment', 'FontWeight', 'normal', 'FontSize', 9);
xlabel(ax4, 'Spatial coordinate x (mm)', 'FontSize', 9);
ylabel(ax4, 'Voltage (V)', 'FontSize', 9);
legend(ax4, 'Location', 'best');

apply_inv_log_2_figure_style(ax1, 180, 125);
