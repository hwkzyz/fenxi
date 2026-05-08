%% Plot 06: waveform fit under the fixed trusted window.
close all;
clc;

mainDir = fileparts(fileparts(mfilename('fullpath')));
addpath(mainDir, '-begin');
ctx = load_inv_log_2_project_context(mainDir);

if ~exist('gTrue', 'var') || isempty(gTrue)
    gTrue = 0.2;
end
D = compute_inv_log_2_fixed_demo(ctx, gTrue);

trainGaps = D.templateInv.gapTrain(:);
[~, nearestIdx] = min(abs(trainGaps - D.gapState.gHat));
nearestGap = trainGaps(nearestIdx);
xAligned = D.highMap.x_v(:) - D.staticState.dx0;
[xSorted, sortIdx] = sort(xAligned, 'ascend');
VNearest = eval_gap_template(D.templateInv, nearestGap, xAligned);
VStatic = eval_gap_template(D.templateInv, D.staticState.gHat, xAligned);
VProjected = eval_gap_template(D.templateInv, D.gapState.gHat, xAligned);
VJoint = D.result.VFit(:);

fig = figure('Name', 'fixed trust waveform fit', 'Color', 'w');
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile;
plot(ax1, xAligned, D.highMap.V_a(:), 'k.', 'MarkerSize', 5, ...
    'DisplayName', 'high-speed samples');
hold(ax1, 'on');
plot(ax1, xSorted, VNearest(sortIdx), '-', 'Color', [0.65 0.65 0.65], ...
    'LineWidth', 1.0, 'DisplayName', 'nearest library waveform');
plot(ax1, xSorted, VStatic(sortIdx), '-', 'Color', [0.85 0.33 0.10], ...
    'LineWidth', 1.3, 'DisplayName', 'static gap waveform');
plot(ax1, xSorted, VProjected(sortIdx), '-', 'Color', [0.00 0.45 0.74], ...
    'LineWidth', 1.4, 'DisplayName', 'projected gap waveform');
xlabel(ax1, 'Aligned spatial coordinate x - dx (mm)', 'FontSize', 9);
ylabel(ax1, 'Voltage (V)', 'FontSize', 9);
title(ax1, '(a) Static gap waveforms in trusted x-window', ...
    'FontWeight', 'normal', 'FontSize', 9);
legend(ax1, 'Location', 'best');

ax2 = nexttile;
[tSorted, tIdx] = sort(D.highMap.t_v(:), 'ascend');
plot(ax2, tSorted, D.highMap.V_a(tIdx), 'k.', 'MarkerSize', 4, ...
    'DisplayName', 'high-speed samples');
hold(ax2, 'on');
plot(ax2, tSorted, VJoint(tIdx), '-', 'Color', [0.47 0.67 0.19], ...
    'LineWidth', 1.2, 'DisplayName', 'joint vibration fit');
yline(ax2, 0, 'k--', 'LineWidth', 0.7, 'HandleVisibility', 'off');
xlabel(ax2, 'Time t (s)', 'FontSize', 9);
ylabel(ax2, 'Voltage (V)', 'FontSize', 9);
title(ax2, '(b) Final joint fit in time domain', ...
    'FontWeight', 'normal', 'FontSize', 9);
legend(ax2, 'Location', 'best');
apply_inv_log_2_figure_style(ax1, 170, 120);
