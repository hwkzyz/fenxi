%% Plot 04: high-speed samples selected by the fixed trusted window.
close all;
clc;

mainDir = fileparts(fileparts(mfilename('fullpath')));
addpath(mainDir, '-begin');
ctx = load_inv_log_2_project_context(mainDir);

if ~exist('gTrue', 'var') || isempty(gTrue)
    gTrue = 0.2;
end
D = compute_inv_log_2_fixed_demo(ctx, gTrue);

fig = figure('Name', 'fixed trust high-speed windowing', 'Color', 'w');
ax = axes(fig);
hold(ax, 'on');

revPick = D.highMap.rev_v(1);
sensorPick = D.highMap.S_v(1);
localMask = D.highMap.rev_v(:) == revPick & D.highMap.S_v(:) == sensorPick;
tCenter = mean(D.highMap.t_v(localMask));
timeHalf = 1.35 * max(abs(D.highMap.t_v(localMask) - tCenter));
idxAll = abs(D.dataHigh.t(:) - tCenter) <= timeHalf;
tAll = D.dataHigh.t(idxAll);
VAll = D.dataHigh.V_cap(idxAll);
plot(ax, tAll, VAll, '.', 'Color', [0.75 0.75 0.75], 'MarkerSize', 4, ...
    'DisplayName', 'wide pulse window');
plot(ax, D.highMap.t_v(localMask), D.highMap.V_a(localMask), 'k.', 'MarkerSize', 7, ...
    'DisplayName', 'samples used for fitting');
yline(ax, D.highMap.activeThreshold, ':', 'Color', [0.65 0.65 0.65], ...
    'LineWidth', 0.9, 'DisplayName', 'diagnostic voltage threshold');
xlabel(ax, 'Time t (s)', 'FontSize', 9);
ylabel(ax, 'Voltage (V)', 'FontSize', 9);
title(ax, 'Fixed-window high-speed sample selection', 'FontWeight', 'normal', 'FontSize', 9);
legend(ax, 'Location', 'best');
apply_inv_log_2_figure_style(ax, 170, 86);
