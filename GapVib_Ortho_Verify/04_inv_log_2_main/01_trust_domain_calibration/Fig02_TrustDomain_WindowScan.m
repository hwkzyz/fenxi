%% Plot 02: candidate-window scan and leave-one-gap validation.
close all;
clc;

mainDir = fileparts(fileparts(mfilename('fullpath')));
addpath(mainDir, '-begin');
ctx = load_inv_log_2_project_context(mainDir);

if ~exist('gHoldout', 'var') || isempty(gHoldout)
    gHoldout = 0.2;
end

modelInv = get_inv_log_2_model_def();
trustOptions = make_gap_trust_options(ctx.cfgAna);
templateInv = make_response_template_library(ctx.gapList, ctx.xCell, ctx.yCell, ...
    gHoldout, ctx.cfgAna.xGridN, modelInv, trustOptions);
metrics = templateInv.trustInfo.metrics;

fig = figure('Name', 'trust domain window scan', 'Color', 'w');
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile;
plot(ax1, metrics.half_width_mm, metrics.template_rmse, 'o-', ...
    'Color', [0.00 0.45 0.74], 'LineWidth', 1.3, 'MarkerSize', 4, ...
    'DisplayName', 'waveform RMSE');
hold(ax1, 'on');
plot(ax1, metrics.half_width_mm, metrics.derivative_rmse, 's-', ...
    'Color', [0.85 0.33 0.10], 'LineWidth', 1.3, 'MarkerSize', 4, ...
    'DisplayName', 'derivative RMSE');
xline(ax1, templateInv.trustInfo.selectedHalfWidth, 'k:', 'LineWidth', 1.1, ...
    'DisplayName', 'selected L');
ylabel(ax1, 'Normalized error', 'FontSize', 9);
title(ax1, '(a) Leave-one-gap validation', 'FontWeight', 'normal', 'FontSize', 9);
legend(ax1, 'Location', 'northeast');

ax2 = nexttile;
yyaxis(ax2, 'left');
plot(ax2, metrics.half_width_mm, metrics.gap_sensitivity, '^-', ...
    'LineWidth', 1.3, 'MarkerSize', 4, 'DisplayName', 'gap sensitivity');
ylabel(ax2, 'Sensitivity', 'FontSize', 9);
yyaxis(ax2, 'right');
plot(ax2, metrics.half_width_mm, metrics.score, 'd-', ...
    'LineWidth', 1.3, 'MarkerSize', 4, 'DisplayName', 'score');
ylabel(ax2, 'Score', 'FontSize', 9);
xline(ax2, templateInv.trustInfo.selectedHalfWidth, 'k:', 'LineWidth', 1.1);
xlabel(ax2, 'Candidate half-width L (mm)', 'FontSize', 9);
title(ax2, '(b) Sensitivity and selection score', 'FontWeight', 'normal', 'FontSize', 9);
apply_inv_log_2_figure_style(ax1, 170, 120);
