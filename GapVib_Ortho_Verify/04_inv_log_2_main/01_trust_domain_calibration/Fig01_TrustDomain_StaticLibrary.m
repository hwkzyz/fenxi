%% Plot 01: static gap library and calibrated trusted window.
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
templateInv = make_response_template_library(ctx.gapList, ctx.xCell, ctx.yCell, ...
    gHoldout, ctx.cfgAna.xGridN, modelInv, trustOptions);

fig = figure('Name', 'trust domain static library', 'Color', 'w');
ax = axes(fig);
hold(ax, 'on');
trainGaps = templateInv.gapTrain(:);
colors = turbo(numel(trainGaps));
for ig = 1:numel(trainGaps)
    idxTrain = find(abs(ctx.gapList(:) - trainGaps(ig)) < 1e-12, 1);
    plot(ax, ctx.xCell{idxTrain}(:), ctx.yCell{idxTrain}(:), '-', ...
        'Color', colors(ig, :), 'LineWidth', 0.9, ...
        'DisplayName', sprintf('library g = %.1f mm', trainGaps(ig)));
end
plot(ax, ctx.xCell{idxHoldout}(:), ctx.yCell{idxHoldout}(:), 'k--', ...
    'LineWidth', 1.2, 'DisplayName', sprintf('held-out g = %.1f mm', gHoldout));
xline(ax, templateInv.trustInfo.effectiveDomain(1), '--', 'Color', [0.55 0.55 0.55], ...
    'LineWidth', 0.9, 'DisplayName', 'effective response domain');
xline(ax, templateInv.trustInfo.effectiveDomain(2), '--', 'Color', [0.55 0.55 0.55], ...
    'LineWidth', 0.9, 'HandleVisibility', 'off');
xline(ax, templateInv.domain(1), ':', 'Color', [0.10 0.10 0.10], ...
    'LineWidth', 1.1, 'DisplayName', 'trusted fitting window');
xline(ax, templateInv.domain(2), ':', 'Color', [0.10 0.10 0.10], ...
    'LineWidth', 1.1, 'HandleVisibility', 'off');
xlabel(ax, 'Spatial coordinate x (mm)', 'FontSize', 9);
ylabel(ax, 'Voltage (V)', 'FontSize', 9);
title(ax, 'Static library and trusted window', 'FontWeight', 'normal', 'FontSize', 9);
legend(ax, 'Location', 'eastoutside');
apply_inv_log_2_figure_style(ax, 170, 88);
