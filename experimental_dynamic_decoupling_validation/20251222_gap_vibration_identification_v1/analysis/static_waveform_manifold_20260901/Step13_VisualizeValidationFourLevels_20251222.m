%% Visualize the four-level validation framework
clear; clc; close all;
thisDir = fileparts(mfilename('fullpath'));
outDir = fullfile(thisDir, 'results');

v1 = readtable(fullfile(outDir, 'V1_SameBladeTransition_20251222.csv'));
v2 = readtable(fullfile(outDir, 'V2_LeaveOneGap_20251222.csv'));
v3 = readtable(fullfile(outDir, 'V3_LeaveOneBlade_IncrementTransfer_20251222.csv'));
v4 = readtable(fullfile(outDir, 'V4_ExperimentalClosure_PerWindow_20251222.csv'));

% Gate limits match Step10_ValidationGate_20251222.m.
v1p90 = prctile(v1.relativeRMSEpct, 90);
v2p90 = prctile(v2.relativeRMSEpct, 90);
v3p90 = prctile(v3.relativeToBaseline, 90);
v4df = mean(abs(v4.frequencyDiffHz_freeMinusNoB));

set(groot, 'defaultAxesFontName', 'Times New Roman', ...
    'defaultTextFontName', 'Times New Roman', 'defaultAxesFontSize', 8.5);
blue = [0.12 0.35 0.62]; orange = [0.85 0.38 0.10];
green = [0.16 0.52 0.28]; red = [0.70 0.14 0.12];

fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 17 12.5]);
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile; plot_sorted_metric(v1.relativeRMSEpct, 10, blue, 'V1', 'Relative RMSE (%)', v1p90, 'P90');
nexttile; plot_sorted_metric(v2.relativeRMSEpct, 10, blue, 'V2', 'Relative RMSE (%)', v2p90, 'P90');
nexttile; plot_sorted_metric(v3.relativeToBaseline, 0.50, orange, 'V3', 'Transfer / baseline', v3p90, 'P90');
nexttile; plot_sorted_metric(abs(v4.frequencyDiffHz_freeMinusNoB), 0.05, green, 'V4', '|Frequency difference| (Hz)', v4df, 'Mean');

export_figure(fig, fullfile(outDir, 'ValidationFourLevels_Summary_20251222'));

% Per-blade distributions for the three static validations.
fig2 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 17 10.5]);
tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile; grouped_box(v1.blade, v1.relativeRMSEpct, 10, 'V1', 'Relative RMSE (%)');
nexttile; grouped_box(v2.blade, v2.relativeRMSEpct, 10, 'V2', 'Relative RMSE (%)');
nexttile; grouped_box(v3.heldOutBlade, v3.relativeToBaseline, 0.50, 'V3', 'Transfer / baseline');
export_figure(fig2, fullfile(outDir, 'ValidationFourLevels_ByBlade_20251222'));

fprintf('Saved four-level validation figures under %s\n', outDir);

function plot_sorted_metric(values, limit, color, tag, yLabel, summary, summaryName)
red = [0.70 0.14 0.12];
values = values(:); values = values(isfinite(values));
values = sort(values, 'ascend');
plot(1:numel(values), values, 'o', 'Color', color, 'MarkerFaceColor', color, 'MarkerSize', 3.2);
hold on; yline(limit, '--', 'Color', red, 'LineWidth', 0.9);
yline(summary, '-', 'Color', [0.15 0.15 0.15], 'LineWidth', 0.9);
hold off; box on; grid off; set(gca, 'TickDir', 'in');
xlim([1 max([1 numel(values)])]);
xlabel('Sorted case index'); ylabel(yLabel);
title(sprintf('%s  |  %s = %.3g', tag, summaryName, summary), 'FontWeight', 'normal');
legend({'Cases', sprintf('Gate %.3g', limit), summaryName}, 'Location', 'northwest', 'Box', 'off', 'FontSize', 7.5);
end

function grouped_box(groups, values, limit, tag, yLabel)
red = [0.70 0.14 0.12];
groups = double(groups(:)); values = values(:);
valid = isfinite(groups) & isfinite(values);
boxchart(categorical(groups(valid)), values(valid), 'BoxFaceColor', blue_local(tag));
hold on; yline(limit, '--', 'Color', red, 'LineWidth', 0.9); hold off;
box on; grid off; set(gca, 'TickDir', 'in');
xlabel('Blade'); ylabel(yLabel); title(tag, 'FontWeight', 'normal');
end

function c = blue_local(tag)
if strcmp(tag, 'V3'), c = [0.85 0.38 0.10]; else, c = [0.12 0.35 0.62]; end
end

function export_figure(fig, baseName)
if exist('exportgraphics', 'file')
    exportgraphics(fig, [baseName '.png'], 'Resolution', 300);
    exportgraphics(fig, [baseName '.pdf'], 'ContentType', 'vector');
else
    print(fig, [baseName '.png'], '-dpng', '-r300');
    print(fig, [baseName '.pdf'], '-dpdf', '-painters');
end
end
