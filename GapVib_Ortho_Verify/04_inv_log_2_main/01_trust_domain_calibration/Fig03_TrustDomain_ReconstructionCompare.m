%% Plot 03: reconstruction comparison under different fitting windows.
close all;
clc;

mainDir = fileparts(fileparts(mfilename('fullpath')));
addpath(mainDir, '-begin');
ctx = load_inv_log_2_project_context(mainDir);

if ~exist('gHoldout', 'var') || isempty(gHoldout)
    gHoldout = 0.2;
end

idxHoldout = find(abs(ctx.gapList - gHoldout) < 1e-12, 1);
modelInv = get_inv_log_2_model_def();
trustOptions = make_gap_trust_options(ctx.cfgAna);
templateTrust = make_response_template_library(ctx.gapList, ctx.xCell, ctx.yCell, ...
    gHoldout, ctx.cfgAna.xGridN, modelInv, trustOptions);
templateFull = make_response_template_library(ctx.gapList, ctx.xCell, ctx.yCell, ...
    gHoldout, ctx.cfgAna.xGridN, modelInv, struct('enable', false));

LSelected = templateTrust.trustInfo.selectedHalfWidth;
LList = unique([trustOptions.minHalfWidth, LSelected, min(abs(templateFull.domain))], 'stable');
labels = ["narrow", "selected", "full"];
colors = [0.55 0.55 0.55; 0.85 0.33 0.10; 0.00 0.45 0.74];

fig = figure('Name', 'trust domain reconstruction comparison', 'Color', 'w');
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile;
plot(ax1, ctx.xCell{idxHoldout}(:), ctx.yCell{idxHoldout}(:), 'k--', ...
    'LineWidth', 1.2, 'DisplayName', 'held-out truth');
hold(ax1, 'on');

ax2 = nexttile;
hold(ax2, 'on');

for ii = 1:numel(LList)
    domainNow = [max(templateFull.domain(1), -LList(ii)), min(templateFull.domain(2), LList(ii))];
    templateNow = restrict_gap_template_domain(templateFull, domainNow);
    xEval = templateNow.xGrid(:);
    yPred = eval_gap_template(templateNow, gHoldout, xEval);
    yTruth = interp1(ctx.xCell{idxHoldout}(:), ctx.yCell{idxHoldout}(:), xEval, 'pchip');
    plot(ax1, xEval, yPred, '-', 'Color', colors(ii, :), 'LineWidth', 1.3, ...
        'DisplayName', sprintf('%s L = %.2f mm', labels(ii), LList(ii)));
    plot(ax2, xEval, yPred - yTruth, '-', 'Color', colors(ii, :), 'LineWidth', 1.2, ...
        'DisplayName', sprintf('%s', labels(ii)));
end

yline(ax2, 0, 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
xlabel(ax1, 'Spatial coordinate x (mm)', 'FontSize', 9);
ylabel(ax1, 'Voltage (V)', 'FontSize', 9);
title(ax1, '(a) Held-out reconstruction', 'FontWeight', 'normal', 'FontSize', 9);
legend(ax1, 'Location', 'best');
xlabel(ax2, 'Spatial coordinate x (mm)', 'FontSize', 9);
ylabel(ax2, 'Residual (V)', 'FontSize', 9);
title(ax2, '(b) Reconstruction residual', 'FontWeight', 'normal', 'FontSize', 9);
legend(ax2, 'Location', 'best');
apply_inv_log_2_figure_style(ax1, 170, 120);
